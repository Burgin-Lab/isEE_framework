"""
Utility functions implemented here are clearly defined unit operations. They may only be called once in the code, but
are defined separately for cleanliness and legibility.
"""

import os
import re
import sys
import copy
import math
import time
import numpy
import pytraj
import mdtraj
import parmed
import shutil
import argparse
import fileinput
import dill as pickle   # I think this is kosher!
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from isee_framework.initialize_charges import set_charges
# from main import Thread

# Two different ways to import tleap depending on, I think, paprika version
try:
    from paprika.build.system import TLeap as tleap
except ModuleNotFoundError:
    from paprika import tleap

def update_progress(progress, message='Progress', eta=0, quiet=False):
    """
    Print a dynamic progress bar to stdout.

    Credit to Brian Khuu from stackoverflow, https://stackoverflow.com/questions/3160699/python-progress-bar

    Parameters
    ----------
    progress : float
        A number between 0 and 1 indicating the fractional completeness of the bar. A value under 0 represents a 'halt'.
        A value at 1 or bigger represents 100%.
    message : str
        The string to precede the progress bar (so as to indicate what is progressing)
    eta : int
        Number of seconds to display as estimated completion time (converted into HH:MM:SS)
    quiet : bool
        If True, suppresses output entirely

    Returns
    -------
    None

    """

    if quiet:
        return None

    barLength = 10  # Modify this to change the length of the progress bar
    status = ""
    if isinstance(progress, int):
        progress = float(progress)
    if not isinstance(progress, float):
        progress = 0
        status = "error: progress var must be float\r\n"
    if progress < 0:
        progress = 0
        status = "Halt...\r\n"
    if progress >= 1:
        progress = 1
        status = "Done!          \r\n"
    block = int(round(barLength * progress))
    if eta:
        # eta is in seconds; convert into HH:MM:SS
        eta_h = str(math.floor(eta/3600))
        eta_m = str(math.floor((eta % 3600) / 60))
        eta_s = str(math.floor((eta % 3600) % 60)) + ' '
        if len(eta_m) == 1:
            eta_m = '0' + eta_m
        if len(eta_s) == 2:
            eta_s = '0' + eta_s
        eta_str = eta_h + ':' + eta_m + ':' + eta_s
        text = "\r" + message + ": [{0}] {1}% {2}".format("#" * block + "-" * (barLength - block), round(progress * 100, 2), status) + " ETA: " + eta_str
    else:
        text = "\r" + message + ": [{0}] {1}% {2}".format("#" * block + "-" * (barLength - block), round(progress * 100, 2), status)
    sys.stdout.write(text)
    sys.stdout.flush()


def mutate(coords, topology, mutation, name, settings, titrations=[]):
    """
    Apply the specified mutation to the structure given by coords and topology and return the names of the new coord
    and topology files.

    The strategy for this method is to export the solvent (including any ions) into a separate object, cast the
    remaining coordinates to the .pdb format, remove the sidechain atoms of the specified residue(s) and rename it, then
    rebuild a new coordinate file and topology using AmberTools' tleap program, which will automatically build the
    missing sidechain for the appropriate residue and resolvate it in the exported solvent (with water molecules deleted
    where there is conflict with the new model). The resulting structure will then be minimized using OpenMM directly
    and the results outputted to a new .rst7 formatted coordinate file and corresponding topology file.

    Parameters
    ----------
    coords : str
        Path to the coordinate file to mutate
    topology : str
        Path to the topology file corresponding to coords
    mutation : list
        List of mutations to apply, each given as "<resid><three-letter code>". For example, "70ASP" mutates residue 70
        to aspartate.
    name : str
        String to prepend to all filenames produced by this method
    settings : argparse.Namespace
        Settings namespace object
    titrations : list
        List of entries [resname, resid, protonation_state] to use to protonate histidine, aspartate, and glutamate
        residues. The protonation_state is the resname for the appropriate protonation state (e.g., HIP, GLH)

    Returns
    -------
    new_coords : str
        Path to the newly created, mutated coordinate file, named as name + '_min.rst7'
    new_topology : str
        Path to the newly created, mutated topology file corresponding to new_coords, named as name + '.prmtop'

    """
    # if len(name) > 200:
    #     print('WARNING: variant name is too long (>200 characters). Truncating. This may cause collisions if another '
    #           'name has the same first 200 characters. Offending name: ' + name)
    #     name = name[:200]

    # So this is dumb but this function sometimes fails at or before tleap and merely needs to be rerun, most recently
    # due to an error in writing the .mol2 files below (one of them just stopped writing mid-stream for some reason.)
    # As a failsafe I'm building in a repeat attempt.
    attempts = 0
    max_attempts = 2
    complete = False
    while not complete and attempts < max_attempts:
        attempts += 1

        if settings.SPOOF:
            return name + '_min.rst7', name + '_tleap.prmtop'

        # todo: implement format checking on 'mutation' (and coords and topology while we're at it)

        from contextlib import contextmanager

        # Define helper function to suppress unwanted output from tleap
        @contextmanager
        def suppress_stderr():
            with open(os.devnull, "w") as devnull:
                old_stderr = sys.stderr
                sys.stderr = devnull
                try:
                    yield
                finally:
                    sys.stderr = old_stderr

        # Force runtime to working directory; should not be necessary, but might be anyway...
        os.chdir(settings.working_directory)

        ### FOR DEBUGGING ###
        mutate_debug = False    # prints information for use in debugging mutate if True

        # Get all non-protein, store separately as mol2 to preserve explicit atom types and topology
        protein_resnames = ':ARG,HIS,HID,HIE,HIP,LYS,ASP,ASH,GLU,GLH,SER,THR,ASN,GLN,CYS,GLY,PRO,ALA,VAL,ILE,LEU,MET,PHE,TYR,TRP,CYX,CYM,HYP'
        if settings.treat_as_protein:   # let the user tell us something else should be considered protein too
            protein_resnames += ',' + ','.join(settings.treat_as_protein)

        traj = pytraj.load(coords, topology)

        # Take a peek at the transition state definition and identify residues that are non-protein (i.e., residues that
        # will show up in the .mol2 file in the next step) and save their coordinates
        ts_atoms = list(set(settings.ts_bonds[0] + settings.ts_bonds[1]))
        ts_atoms_full = settings.ts_bonds[0] + settings.ts_bonds[1]
        ts_xyzs = []
        if not ts_atoms == ['']:
            for atom in ts_atoms:
                atm = traj.top.select(atom)     # get atom index
                # resname = str(traj.top.residue((traj.top.atom(atm).resid)))[1:4]
                xyz = traj.xyz[0][atm]
                try:
                    assert xyz.size == 3
                except AssertionError:
                    raise RuntimeError('One or more atoms selection strings in the transition state definition did not '
                                       'match exactly one atom in the topology: ' + topology + '. The first offending '
                                       'selection found was: ' + atom + ', which matched ' + str(int(xyz.size/3)) +
                                       ' atoms.')
                ts_xyzs.append(xyz)  # coordinates of transition state atoms as a list

        if mutate_debug:
            print('DEBUG: ts_atoms = ' + str(ts_atoms))
            print('DEBUG: ts_xyzs = ' + str(ts_xyzs))

        traj.strip(protein_resnames)
        no_mol2 = False
        try:
            traj.strip(':WAT,Na+,Cl-')
        except IndexError:
            print('System is only protein and solvent; skipping generating .mol2 file.')
            no_mol2 = True    # caused if stripping the above would leave nothing left
        pytraj.write_traj(name + '_nonprot.mol2', traj, overwrite=True)

        def _all(arg_list):
            # Version of all() that returns False when passed an empty list
            if arg_list == []:
                return False
            else:
                return all(arg_list)

        def is_ts_bond(atoms, ts_atoms_mol2, ts_atoms_mol2_indices):  # todo: this is insanely ugly
            # Determine if the atoms in atoms constitute a ts bond
            #  atoms : indices of atoms from mol2 file to check
            #  ts_atoms_mol2 : indices of atoms from mol2 file corresponding to ts atoms
            #  ts_atoms_mol2_indices : indices mapping entries in ts_atoms_mol2 to entries in ts_atoms_full
            #  ts_atoms_full : strings matching ts atoms (not indices)
            if not all([atom in ts_atoms_mol2 for atom in atoms]):
                return False

            assert len(ts_atoms_full) % 2 == 0   # ts_atoms_full must be of even length
            ts_atoms_firsthalf = [ts_atoms_full[i] for i in range(0,int(len(ts_atoms_full) / 2))]
            ts_atoms_secondhalf = [ts_atoms_full[i] for i in range(int(len(ts_atoms_full) / 2), len(ts_atoms_full))]

            relevant_indices = [ts_atoms_mol2_indices[ts_atoms_mol2.index(atoms[0])], ts_atoms_mol2_indices[ts_atoms_mol2.index(atoms[1])]]
            relevant_atoms = [ts_atoms[relevant_index] for relevant_index in relevant_indices]
            for ii in range(len(ts_atoms_firsthalf)):
                if (relevant_atoms[0] == ts_atoms_firsthalf[ii] and relevant_atoms[1] == ts_atoms_secondhalf[ii]) or \
                   (relevant_atoms[1] == ts_atoms_firsthalf[ii] and relevant_atoms[0] == ts_atoms_secondhalf[ii]):
                    return True

            if mutate_debug:
                print('DEBUG: is_ts_bond between ts_atoms False: ')
                print(' DEBUG: atoms = ' + str(atoms))
                print(' DEBUG: ts_atoms_mol2 = ' + str(ts_atoms_mol2))
                print(' DEBUG: ts_atoms_mol2_indices = ' + str(ts_atoms_mol2_indices))
                print(' DEBUG: ts_atoms_full = ' + str(ts_atoms_full))

            return False

        ### Remove all bond terms between atoms with bond definitions in ts_bonds
        if mutate_debug:
            print('DEBUG: no_mol2 = ' + str(no_mol2))
        if not no_mol2:
            open(name + '_nonprot_mod.mol2', 'w').close()
            with open(name + '_nonprot_mod.mol2', 'a') as f:
                atoms_yet = False
                bonds_yet = False
                substructure_yet = False
                index_name_list = []
                ts_atoms_mol2 = []
                ts_atoms_mol2_indices = []
                removed_lines = 0
                bond_count = 0
                for line in open(name + '_nonprot.mol2', 'r').readlines():
                    if '@<TRIPOS>ATOM' in line:
                        atoms_yet = True
                        newline = line
                        f.write(newline)
                        continue
                    if '@<TRIPOS>BOND' in line:
                        bonds_yet = True
                        atoms_yet = False
                        index_name_list = list(map(list, zip(*index_name_list)))    # transpose index_name_list
                        newline = line
                        f.write(newline)
                        continue
                    if '@<TRIPOS>SUBSTRUCTURE' in line:
                        if not bonds_yet:   # happens if nonprot has no bonds, and breaks tleap
                            f.write('@<TRIPOS>BOND\n')
                        substructure_yet = True
                        bonds_yet = False
                        newline = line
                        f.write(newline)
                        continue
                    if substructure_yet:    # a substructure line
                        newline = line
                        f.write(newline)
                        continue
                    if atoms_yet and not 'WAT' in line: # an atom line not water
                        split = line.split()
                        if mutate_debug:
                            print('DEBUG: coord_compare: ' + str(split))
                        # if the coordinates for this atom match the coordinates of any of the ts atoms...
                        coord_compare = [_all([math.isclose(float(split[2]), ts_xyzs[i][0][0], abs_tol=1e-3),
                                         math.isclose(float(split[3]), ts_xyzs[i][0][1], abs_tol=1e-3),
                                         math.isclose(float(split[4]), ts_xyzs[i][0][2], abs_tol=1e-3)]) for i in range(len(ts_xyzs))]
                        if any(coord_compare):
                            if mutate_debug:
                                print('DEBUG: match!')
                            ts_atoms_mol2.append(split[0])  # add this atom index to a list of ts_atoms in the mol2
                            try:
                                assert(coord_compare.count(True) == 1)
                            except AssertionError:
                                raise RuntimeError('somehow more than one transition state atom has roughly the same '
                                                   'coordinates (to within 1e-3 nanometers in x, y, and z). This disrupts '
                                                   'isEE\'s ability to process structures. Address this and try again.')
                            ts_atoms_mol2_indices.append(coord_compare.index(True))    # to keep track of which atom is which
                        index_name_list.append([split[0], split[1]])   # index, name
                    if bonds_yet == False:  # a water atom line, or preamble
                        newline = line
                        f.write(newline)
                        continue
                    else:   # a bond line
                        atoms = line.split()[1:3]
                        try:
                            [index_name_list[1][index_name_list[0].index(atom)] for atom in atoms]
                        except ValueError:  # atom index not in list, so it's a WAT
                            bond_count += 1
                            newline = line.replace(line.split()[0], str(bond_count), 1)
                            f.write(newline)
                            continue
                        if is_ts_bond(atoms, ts_atoms_mol2, ts_atoms_mol2_indices):
                            if mutate_debug:
                                print('DEBUG: is_ts_bond True: ' + str(atoms))
                            removed_lines += 1
                            continue
                        else:
                            bond_count += 1
                            newline = line.replace(line.split()[0], str(bond_count), 1)
                            f.write(newline)
                            continue

            # Adjust number of bonds
            open(name + '_nonprot_mod_2.mol2', 'w').close()
            with open(name + '_nonprot_mod_2.mol2', 'a') as f2:
                count = 0
                for line in open(name + '_nonprot_mod.mol2', 'r').readlines():
                    if count == 2:
                        numbers = line.split()
                        newline = line[::-1].replace(numbers[1][::-1], str(int(numbers[1]) - removed_lines)[::-1], 1)[::-1] # replace once from right
                        f2.write(newline)
                    else:
                        f2.write(line)
                    count += 1

            os.remove(name + '_nonprot_mod.mol2')
            os.remove(name + '_nonprot.mol2')
            os.rename(name + '_nonprot_mod_2.mol2', name + '_nonprot.mol2')

        ### If keep_waters, save the waters to a separate pdb
        if settings.keep_waters:
            traj = pytraj.load(coords, topology)
            traj.strip('!(:WAT,HOH)')
            pytraj.write_traj(name + '_water.pdb', traj, overwrite=True)

        ### Cast remainder to separate .pdb
        traj = pytraj.load(coords, topology)
        traj.strip('!(' + protein_resnames + ')')
        pytraj.write_traj(name + '_prot.pdb', traj, overwrite=True)

        ### Mutate
        # First, reset all titrations
        # Rename all ASH -> ASP, GLH -> GLU, and HIP, HID, and HIE -> HIS
        temp_titrations = []
        for line in fileinput.input(name + '_prot.pdb', inplace=True):
            if not titrations:  # if we don't have new, explicit titrations, we want to save the old ones
                if line.split()[0] == 'ATOM':
                    temp_titrations.append([line.split()[3].replace(
                        'ASH', 'ASP').replace(
                        'GLH', 'GLU').replace(
                        'HIP', 'HIS').replace(
                        'HID', 'HIS').replace(
                        'HIE', 'HIS'), line.split()[4], line.split()[3]])

            print(line.replace(
                ' ASH ', ' ASP ').replace(
                ' GLH ', ' GLU ').replace(
                ' HIP ', ' HIS ').replace(
                ' HID ', ' HIS ').replace(
                ' HIE ', ' HIS '), end='')
        pdb_to_modify = name + '_prot.pdb'
        if temp_titrations:
            titrations = temp_titrations

        # Implement Rosetta mutator, if desired
        if settings.rosetta_mutate:
            import pyrosetta
            import rosetta
            from pyrosetta import standard_packer_task
            from pyrosetta import pose_from_file
            from pyrosetta import Pose
            from pyrosetta import create_score_function
            from rosetta.utility import vector1_bool
            from rosetta.core.chemical import aa_from_oneletter_code
            from rosetta.protocols.minimization_packing import PackRotamersMover
            from rosetta.core.pose import PDBInfo

            def mutate_residue(pose, mutant_position, mutant_aa,
                               pack_radius=0.0, pack_scorefxn='', settings=None):
                """
                Replaces the residue at  <mutant_position>  in  <pose>  with  <mutant_aa>
                    and repack any residues within  <pack_radius>  Angstroms of the mutating
                    residue's center (nbr_atom) using  <pack_scorefxn>
                note: <mutant_aa>  is the single letter name for the desired ResidueType

                example:
                    mutate_residue(pose,30,A)
                See also:
                    Pose
                    PackRotamersMover
                    MutateResidue
                    pose_from_sequence
                """
                #### a MutateResidue Mover exists similar to this except it does not pack
                ####    the area around the mutant residue (no pack_radius feature)
                # mutator = MutateResidue( mutant_position , mutant_aa )
                # mutator.apply( test_pose )
                #
                # Code adapted by Tucker E. Burgin from Evan H. Baugh, in turn adapted from Sid Chaudhury.

                if pose.is_fullatom() == False:
                    raise IOError('mutate_residue only works with fullatom poses')

                test_pose = Pose()
                test_pose.assign(pose)

                # create a beta_nov16 scorefxn by default (changes to this may require change to pyrosetta.init call)
                if not pack_scorefxn:
                    pack_scorefxn = create_score_function('beta_nov16')

                task = standard_packer_task(test_pose)

                aa_bool = rosetta.utility.vector1_bool()
                mutant_aa = aa_from_oneletter_code(mutant_aa)

                for i in range(1, 21):
                    aa_bool.append(i == mutant_aa)

                task.nonconst_residue_task(mutant_position).restrict_absent_canonical_aas(aa_bool)

                # prevent residues from packing by setting the per-residue "options" of the PackerTask
                center = pose.residue(mutant_position).nbr_atom_xyz()
                for i in range(1, pose.total_residue() + 1):
                    # only pack the mutating residue and any within the pack_radius
                    if (not i == mutant_position or center.distance_squared(
                            test_pose.residue(i).nbr_atom_xyz()) > pack_radius ** 2) \
                            or i in settings.rosetta_prevent_repacking:
                        task.nonconst_residue_task(i).prevent_repacking()

                # apply the mutation and pack nearby residues
                packer = PackRotamersMover(pack_scorefxn, task)
                packer.apply(test_pose)

                return test_pose

            def convert_3_1(resname):
                # Helper function to convert three-letter residue names to one-letter code
                all_resnames = [['ARG', 'HIS', 'LYS', 'ASP', 'GLU', 'SER', 'THR', 'ASN', 'GLN', 'CYS', 'GLY', 'PRO', 'ALA',
                                 'VAL', 'ILE', 'LEU', 'MET', 'PHE', 'TYR', 'TRP', 'GLH', 'ASH', 'HIP', 'HIE', 'HID', 'CYX',
                                 'CYM', 'HYP'],
                                ['R', 'H', 'K', 'D', 'E', 'S', 'T', 'N', 'Q', 'C', 'G', 'P', 'A', 'V', 'I', 'L', 'M', 'F',
                                 'Y', 'W', 'E', 'D', 'H', 'H', 'H', 'C', 'C', 'P']]

                try:
                    result = all_resnames[1][all_resnames[0].index(resname.upper())]
                except ValueError:
                    raise RuntimeError('got unknown residue name: ' + resname)

                return result

            init_string = '-corrections::beta_nov16'
            if not settings.rosetta_override == ['']:   # bool(['']) == True, surprisingly
                init_string += ' -PDB_components_overrides ' + ' '.join(settings.rosetta_override)
            pyrosetta.init(init_string)  # initialize with corrections for beta_nov16 weights
            opt = pyrosetta.rosetta.core.import_pose.ImportPoseOptions()    # initialize pose options object
            opt.set_keep_input_protonation_state(True)                      # don't try to protonate
            opt.set_ignore_zero_occupancy(False)                            # don't know what this does, honestly
            mypose = pose_from_file(name + '_prot.pdb', opt, False, pyrosetta.rosetta.core.import_pose.PDB_file)    # load pdb file
            for mut in mutation:        # apply each mutation
                if mut:     # specifically fixes issue with wild type (mut == '')
                    resid = int(mut[:-3])
                    target = convert_3_1(mut[-3:])
                    mypose = mutate_residue(mypose, resid, target, pack_radius=8, settings=settings)  # here's where the actual mutation is performed
            pyrosetta.rosetta.core.io.pdb.dump_pdb(mypose, name + '_rosetta.pdb')   # write output to a new .pdb file

            # Rosetta does weird things to hydrogen atoms, so we're gonna remove them all and let tleap put them back
            for line in fileinput.input(name + '_rosetta.pdb', inplace=True):
                if not line.split() or not line.split()[0] in ['ATOM', 'HETATM']:   # non-atom
                    print(line, end='')
                elif 'H' in line.split()[-1]:       # hydrogen atom
                    pass
                else:                               # non-hydogen atom
                    print(line, end='')

            pdb_to_modify = name + '_rosetta.pdb'   # set the output from this block as the input for the next one

        with open(name + '_mutated.pdb', 'w') as f:
            patterns = [re.compile('\s+[A-Z0-9]+\s+[A-Z]{3}\s+' + mutant[:-3] + '\s+') for mutant in mutation if mutant]
            for line in open(pdb_to_modify, 'r').readlines():
                if not settings.rosetta_mutate and patterns and not all(pattern.findall(line) == [] for pattern in patterns):
                    pat_index = 0
                    for pattern in patterns:
                        try:
                            if pattern.findall(line)[0].split()[0] in ['C','N','O','CA']:
                                newline = line.replace(pattern.findall(line)[0].split()[1], mutation[pat_index][-3:].upper())
                                break
                            else:
                                newline = ''
                        except IndexError:  # happens when line doesn't match pattern
                            pass
                        pat_index += 1
                elif titrations:  # if we were passed titration results to use and this line wasn't already mutated
                    try:
                        titration_index = -1
                        try:
                            resnameandindex = [line.split()[3], line.split()[4]]
                            titration_index = [[titration[0], titration[1]] for titration in titrations].index(
                                resnameandindex)
                        except IndexError:  # line doesn't have a residue listed
                            newline = line
                            f.write(newline)
                            continue
                        if titration_index >= 0 and titrations[titration_index][0] in ['HIS', 'ASP', 'GLU'] and \
                                not int(titrations[titration_index][1]) in settings.immutable and \
                                not titrations[titration_index][0] == titrations[titration_index][2] and \
                                not line.split()[2][0] == 'H':  # last statement: if this is not a hydrogen
                            newline = line.replace(titrations[titration_index][0], titrations[titration_index][2])
                        elif titration_index >= 0 and titrations[titration_index][0] in ['HIS', 'ASP', 'GLU'] and \
                                not int(titrations[titration_index][1]) in settings.immutable and \
                                not titrations[titration_index][0] == titrations[titration_index][2] and \
                                line.split()[2][0] == 'H':  # same as above but it IS a hydrogen
                            newline = ''    # remove it; tleap will add it back in if appropriate
                        else:
                            newline = line
                    except ValueError:  # this line doesn't correspond to an entry in titrations
                        newline = line
                else:
                    newline = line
                f.write(newline)

        ### Rebuild with tleap into .rst7 and .prmtop files
        # todo: there's no way to use tleap to do this without being forced into using Amber (or CHARMM?) force fields...
        # todo: I need to come up with a different strategy if I want to move away from Amber-only. Really this whole thing
        # todo: needs to be user-customizable in some way, fundamental as it is to the process as a whole.
        try:
            system = tleap()
        except TypeError:
            system = tleap.System()
        system.pbc_type = None  # turn off automatic solvation
        system.neutralize = False
        system.output_path = settings.working_directory
        system.output_prefix = name + '_tleap'
        system.template_lines = ['source ' + item + '\n' for item in settings.paths_to_forcefields if item] + \
            ['source leaprc.protein.ff19SB',
            'source leaprc.GLYCAM_06j-1',
            'source leaprc.water.opc',  # essential to load OPC last to avoid solvent model getting overwritten
            'WAT = OP3',
            'HOH = OP3',
            'mut = loadpdb ' + name + '_mutated.pdb',
            'nonprot = loadmol2 '  + name + '_nonprot.mol2',
            'model = combine { mut nonprot }',
             settings.tleap_extra,
            'solvateoct model OPC3BOX 8.0',
            'addIons model Na+ 0',
            'addIons model Cl- 0'
            # 'set model box {' + box_dimensions + '}'
        ]
        if settings.keep_waters and not no_mol2:
            system.template_lines = system.template_lines[:system.template_lines.index('model = combine { mut nonprot }')] + \
                                     ['water = loadpdb ' + name + '_water.pdb',
                                      'model = combine { mut nonprot water }'] + \
                                     system.template_lines[system.template_lines.index('model = combine { mut nonprot }') + 1:]
            system.template_lines.remove('solvateoct model OPC3BOX 8.0')
        elif not settings.keep_waters and no_mol2:
            system.template_lines.remove('nonprot = loadmol2 ' + name + '_nonprot.mol2')
            system.template_lines.remove('model = combine { mut nonprot }')
            system.template_lines[system.template_lines.index('mut = loadpdb ' + name + '_mutated.pdb')] = 'model = loadpdb ' + name + '_mutated.pdb'
        elif settings.keep_waters and no_mol2:
            system.template_lines = system.template_lines[:system.template_lines.index('model = combine { mut nonprot }')] + \
                                     ['water = loadpdb ' + name + '_water.pdb',
                                      'model = combine { mut water }'] + \
                                     system.template_lines[system.template_lines.index('model = combine { mut nonprot }') + 1:]
            system.template_lines.remove('nonprot = loadmol2 ' + name + '_nonprot.mol2')
            system.template_lines.remove('solvateoct model OPC3BOX 8.0')
        with suppress_stderr():
            try:
                system.build(clean_files=False)  # produces a ton of unwanted "WARNING" messages in stderr even when successful
            except TypeError:   # older versions don't support clean_files argument
                system.build()
        try:
            shutil.copy('leap.log', name + '_leap.log')
        except FileNotFoundError:   # encountered this once where I think it was caused by tleap being just a bit slow
            time.sleep(10)
            shutil.copy('leap.log', name + '_leap.log')

        mutated_rst = name + '_tleap.rst7'
        mutated_top = name + '_tleap.prmtop'

        if os.path.exists(mutated_top) and os.path.exists(mutated_rst):
            complete = True

    if not complete:
        raise RuntimeError('Unable to build ' + mutated_top + ' and ' + mutated_rst + ' for some reason.\n'
                           'Check ' + name + '_leap.log\n'
                           'Attempted ' + str(attempts) + ' time(s).')

    # Use ParmED to add ts_bonds to mutated_top; also handles HMR
    do_parmed(mutated_top, mutated_rst, settings)

    # If keep_waters, we won't have box info, so use parmed to copy it over from the input topology
    if settings.keep_waters:
        try:
            parmed_top = parmed.load_file(mutated_top)
            parmed_top.load_rst7(mutated_rst)
        except parmed.exceptions.FormatNotFound:
            raise RuntimeError('problem with topology file: ' + top + '\nDid something go wrong with tleap?')

        ref_top = parmed.load_file(settings.init_topology)
        parmed_top.box = ref_top.box
        parmed_top.write_parm(mutated_top)
        parmed_top.write_rst7(mutated_rst)

    # If appropriate, apply calculated charges
    if settings.initialize_charges:
        mutated_top = set_charges(mutated_top)

    if settings.min_steps > 0:
        ### Minimize with OpenMM
        # First, cast .prmtop to OpenMM topology todo: replace Amber-specific stuff with call to method of MDEngine that returns an OpenMM Simulation object
        openmm_top = AmberPrmtopFile(mutated_top)
        openmm_sys = openmm_top.createSystem(constraints=HBonds, nonbondedMethod=PME, nonbondedCutoff=0.8*nanometer)
        openmm_rst = AmberInpcrdFile(mutated_rst)
        integrator = LangevinIntegrator(300 * kelvin, 1 / picosecond, 0.002 * picoseconds)
        simulation = Simulation(openmm_top.topology, openmm_sys, integrator)
        simulation.context.setPositions(openmm_rst.positions)

        if openmm_rst.boxVectors is not None:
            simulation.context.setPeriodicBoxVectors(*openmm_rst.boxVectors)

        simulation.minimizeEnergy()
        simulation.reporters.append(PDBReporter(name + '_min.pdb', int(settings.min_steps/10), enforcePeriodicBox=False))
        simulation.reporters.append(StateDataReporter(sys.stdout, int(settings.min_steps/10), step=True, potentialEnergy=True, temperature=True))
        simulation.step(settings.min_steps)

        ### Return results
        # First, cast minimization output .pdb back to .rst7
        traj = mdtraj.load_frame(name + '_min.pdb', -1)
        traj.save_amberrst7(name + '_min.rst7')
        traj.save_pdb()

        # Some file cleanup
        os.remove(name + '_min.pdb')
        os.remove(name + '_tleap.rst7')

        to_return = name + '_min.rst7'
    else:
        to_return = mutated_rst

    ## Deprecated code to do casting with pytraj; fails to carry over box information
    # traj = pytraj.load('min.pdb', mutated_top, frame_indices=[-1])
    # pytraj.write_traj(new_name + '_min.rst7', traj)
    # if os.path.exists(new_name + '_min.rst7.1'):
    #     os.rename(new_name + '_min.rst7.1', new_name + '_min.rst7')

    # Clean up unnecessary files and return the coordinate and topology files!
    os.remove(name + '_mutated.pdb')
    os.remove(name + '_prot.pdb')
    os.remove(name + '_nonprot.mol2')
    if settings.keep_waters:
        os.remove(name + '_water.pdb')

    return to_return, mutated_top


def do_parmed(top, rst, settings):
    """
    Add settings.ts_bonds to the given topology file using parmed. The file will be modified in place.

    Also applied hydrogen mass repartitioning if settings.hmr = True

    Parameters
    ----------
    top : str
        Path to the topology file to modify
    rst : str
        Path to the restart file to modify
    settings : argparse.Namespace
        Settings namespace object

    Returns
    -------
    None

    """

    # Load topology into parmed
    try:
        parmed_top = parmed.load_file(top)
    except parmed.exceptions.FormatNotFound:
        raise RuntimeError('problem with topology file: ' + top + '\nDid something go wrong with tleap?')

    parmed_top.load_rst7(rst)

    ts_bonds = list(map(list, zip(*settings.ts_bonds)))
    if not ts_bonds == [['', '', -1, -1]]:  # default setting for no ts_bonds
        for bond in ts_bonds:
            arg = [str(item) for item in bond]
            try:
                setbond = parmed.tools.actions.setBond(parmed_top, arg[0], arg[1], arg[2], arg[3])
                setbond.execute()
            except parmed.tools.exceptions.SetParamError as e:
                raise RuntimeError('encountered parmed.tools.exceptions.SetParamError: ' + str(e) + '\n'
                                   'The offending bond and topology are: ' + str(arg) + ' and ' + top)

    if settings.hmr:
        action = parmed.tools.actions.HMassRepartition(parmed_top)
        action.execute()

    # Save the topology and coordinate files with the new bonds
    parmed_top.write_parm(top)
    parmed_top.write_rst7(rst)


def strip_and_store(traj, top, settings):
    """
    Strip water (':WAT,HOH') out of the given trajectory and topology files and store the "dry" versions in
    settings.storage_directory. Water molecules within settings.dry_distance of non-solvent atoms other than Na+ or Cl-
    are retained (none if dry_distance is 0). These distances are measured relative to the coordinates in the last frame
    of traj.

    The files are named the same as the input files, except for '_dry' inserted just before the file extension.

    Parameters
    ----------
    traj : str
        Path to trajectory file
    top : str
        Path to topology file
    settings : argparse.Namespace
        Settings namespace

    Returns
    -------
    dry_traj :
        Path to dry, stored trajectory file
    dry_top :
        Path to dry, stored topology file

    """

    if not os.path.exists(settings.storage_directory):
        os.mkdir(settings.storage_directory)

    # Trajectory
    ptraj = pytraj.iterload(traj, top)
    ptraj_ref_frame = ptraj[-1]     # save this now because the data type of ptraj is changed by strip
    ptraj.top.set_reference(ptraj_ref_frame)
    ptraj = pytraj.strip(ptraj, ':WAT,HOH & (!:WAT,HOH,Na+,Cl-)>:' + str(settings.dry_distance))
    dry_traj_name = traj[:traj.rindex('.')] + '_dry' + traj[traj.rindex('.'):]  # insert '_dry'
    if '/' in dry_traj_name:
        dry_traj_name = dry_traj_name[traj.rindex('/') + 1:]                    # remove path, leaving only filename
    pytraj.write_traj(settings.storage_directory + '/' + dry_traj_name, ptraj, overwrite=True)  # save it to storage

    # Topology
    ptop = pytraj.load_topology(top)
    ptop.set_reference(ptraj_ref_frame)
    ptop = pytraj.strip(ptop, ':WAT,HOH & (!:WAT,HOH,Na+,Cl-)>:' + str(settings.dry_distance))
    dry_top_name = top[:top.rindex('.')] + '_dry' + top[top.rindex('.'):]       # insert '_dry'
    if '/' in dry_top_name:
        dry_top_name = dry_top_name[traj.rindex('/') + 1:]                      # remove path, leaving only filename
    pytraj.write_parm(settings.storage_directory + '/' + dry_top_name, ptop, overwrite=True)    # save it to storage

    return settings.storage_directory + '/' + dry_traj_name, settings.storage_directory + '/' + dry_top_name


def muts_to_current_name(thread, settings):
    """
    From the last entry in a Thread's history.muts object, construct a string for use as its next current_name.

    Necessary to handle arbitrary-length lists of mutations without potentially running into file name length limits.

    Parameters
    ----------
    thread : Thread
        Thread object
    settings : argparse.Namespace
        Settings namespace

    Returns
    -------
    current_name : str
        Unique name derived from mutation list
    """

    def get_timestamp():
        if not os.path.exists('timestamp_name_list.txt'):
            open('timestamp_name_list.txt', 'w').close()

        for line in open('timestamp_name_list.txt', 'r').readlines():
            if thread.name + '_' + '_'.join(thread.history.muts[-1]) in line:
                return line.split()[1]

        timestamp = str(time.time())
        open('timestamp_name_list.txt', 'a').write(
            thread.name + '_' + '_'.join(thread.history.muts[-1]) + '\t' + timestamp)

        return timestamp

    if settings.name_as_timestamp:
        return get_timestamp()

    # Evaluate max length of string based on other thread attributes
    # Assume total length = thread.name + 30
    reserved_length = len(thread.name) + 30
    max_length = 256 - reserved_length
    if len('_'.join(thread.history.muts[-1])) < max_length:
        return '_'.join(thread.history.muts[-1])

    def three2one(code):
        # Convert three-letter amino acid codes to one-letter codes
        all_3 = ['ARG', 'HIS', 'LYS', 'ASP', 'GLU', 'SER', 'THR', 'ASN', 'GLN', 'CYS', 'GLY', 'PRO', 'ALA',
                 'VAL', 'ILE', 'LEU', 'MET', 'PHE', 'TYR', 'TRP']
        all_1 = ['R', 'H', 'K', 'D', 'E', 'S', 'T', 'N', 'Q', 'C', 'G', 'P', 'A', 'V', 'I', 'L', 'M', 'F',
                 'Y', 'W']
        return all_1[all_3.index(code)]

    ones = '_'.join([item[:-3] + three2one(item[-3:]) for item in thread.history.muts[-1]])
    if len(ones) < max_length:
        return ones

    return get_timestamp()


if __name__ == '__main__':
    ### This stuff is all for testing, shouldn't ever be called during an isEE run ###
    settings = argparse.Namespace()
    settings.SPOOF = False
    settings.working_directory = './'
    settings.topology = 'TmAfc_D224G_t200.prmtop'
    settings.lie_alpha = 0.18
    settings.lie_beta = 0.33
    settings.hmr = False
    settings.ts_mask = ':442,443'
    settings.paths_to_forcefields = ['171116_FPA_4NP-Xyl_ff.leaprc']
    settings.min_steps = 0

    # Residue1:AtomName1; Residue2:AtomName2; weight in kcal/mol-Å**2; equilibrium bond length in Å
    settings.ts_bonds = ([':260@OE2', ':***@O4',  ':***@O4', ':***@N1'],
                         [':***@H4O', ':***@H4O', ':***@C1', ':***@C1'],
                         [200,        200,        200,       200],
                         [1.27,       1.23,       1.9,       2.4])
    coords = 'data/one_frame.rst7'
    topology = 'data/TmAfc_D224G_t200.prmtop'

    mutate(coords, topology, ['64ASP','60ALA'], 'test', settings)

    thread = main.Thread()
    thread.history.trajs = ['']
