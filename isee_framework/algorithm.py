"""
Interface for Algorithm objects. New Algorithms can be implemented by constructing a new class that inherits from
Algorithm and implements its abstract methods.
"""

import os
import re
import sys
import abc
import copy
import time
import numpy
import shutil
import pickle
import mdtraj
import argparse
import itertools
import subprocess
import pandas as pd
from math import factorial
from filelock import Timeout, FileLock
from isee_framework import utilities
from isee_framework.infrastructure import factory


class Algorithm(abc.ABC):
    """
    Abstract base class for isEE algorithms.

    Implements methods for all of the algorithm-specific tasks that isEE might need.

    """

    def __init__(self):
        # define backlog list for algorithms that determine more than one step at a time
        if os.path.exists('backlog.pkl'):
            self.backlog = pickle.load(open('backlog.pkl', 'rb'))
        else:
            self.backlog = []

        ## Deprecated
        # ### Load algorithm_history from algorithm_history.pkl ###
        # if os.path.exists('algorithm_history.pkl'):
        #     self.algorithm_history = pickle.load(open('algorithm_history.pkl', 'rb'))
        # else:
        #     raise FileNotFoundError('algorithm_history.pkl not found; it should have been created when the first '
        #                             'thread was initialized. Did you do something unusual?')

    @staticmethod
    def build_algorithm_history(allthreads, settings):
        # First, establish lock so only one instance of build_algorithm_history will run at a time
        # If we have a shared history file, this needs to go in the same directory as that
        if settings.shared_history_file and '/' in settings.shared_history_file:
            lockfile_directory = settings.shared_history_file[:settings.shared_history_file.rindex('/') + 1]
        else:
            lockfile_directory = ''
        lock = FileLock(lockfile_directory + 'algorithm_history.lock')
        with lock:
            open(lockfile_directory + 'arbitrary_lockfile.lock', 'w').close()
        lock.acquire()
        open(lockfile_directory + 'arbitrary_lockfile.lock', 'w').close()    # this step will block until the lock is released

        # Build algorithm_history afresh from all thread history attributes
        if settings.shared_history_file:
            history_file = settings.shared_history_file
        else:
            history_file = 'algorithm_history.pkl'

        if not os.path.exists(history_file):
            algorithm_history = argparse.Namespace()
            for key in allthreads[0].history.__dict__.keys():
                algorithm_history.__dict__[key] = []    # initialize dictionary
        else:
            algorithm_history = pickle.load(open(history_file, 'rb'))

        def merge_namespace(current, new):
            # Helper function to merge contents of 'new' namespace into 'current' namespace, assuming both have the same
            # keys, that each key corresponds to a list, and that there are the same number of entries in each list
            # (i.e., all lists in 'new' have one number of entries, and all lists in 'current' have another number of
            # entries).
            # If a particular "column" of entries in 'new' is already present in 'current' then it is not merged.]))
            # # First, get lengths of current.__dict__ and new.__dict__ lists, to later restore
            # current_lengths = [len(current.__dict__[key]) for key in list(current.__dict__.keys())]
            # new_lengths = [len(new.__dict__[key]) for key in list(new.__dict__.keys())]
            # lengths = list(numpy.sum([current_lengths, new_lengths], 0))
            # This implementation accomplishes the task by converting to pandas data frames
            current_df = pd.DataFrame.from_dict(current.__dict__, orient='index').T
            # load sideways then transpose to accomodate empty entries by adding a NaN in place, to be removed later
            new_df = pd.DataFrame.from_dict(new.__dict__, orient='index').T
            merged_df = pd.concat([current_df, new_df])
            merged_df = merged_df[~merged_df.astype(str).duplicated()]
            merged_df = merged_df.sort_values('timestamps') # sort every column chronologically by timestamp attribute
            merged = merged_df.to_dict('list')
            return merged

        # Perform merge
        for thread in allthreads:
            algorithm_history.__dict__ = merge_namespace(algorithm_history, thread.history)

        # Remove extraneous NaNs
        for key in list(algorithm_history.__dict__.keys()):
            algorithm_history.__dict__[key] = [x for x in algorithm_history.__dict__[key] if not (str(x).lower() == 'nan' or str(x).lower() == 'none')]
            # algorithm_history.__dict__[key] = [x for _,x in sorted(zip(algorithm_history.__dict__['timestamps'], algorithm_history.__dict__[key]))]

        # Dump pickle file to 'algorithm_history.pkl' regardless of settings.shared_history_file, because every working
        # directory should have its own algorithm history object in addition to the shared one.
        pickle.dump(algorithm_history, open('algorithm_history.pkl.bak', 'wb'))
        if not os.path.getsize('algorithm_history.pkl.bak') == 0:
            shutil.copy('algorithm_history.pkl.bak', 'algorithm_history.pkl')  # copying after is safer
        else:
            raise RuntimeError('failed to dump algorithm history pickle file')

        # Copy to shared file if necessary
        if settings.shared_history_file:
            shutil.copy('algorithm_history.pkl', settings.shared_history_file)

        lock.release()
        return algorithm_history

    @abc.abstractmethod
    def get_first_step(self, thread, allthreads, settings):
        """
        Determine the first step for a thread. This method should be called before the first 'process' step after a new
        thread is created (as defined by having an empty thread.history.trajs attribute) to allow the algorithm to
        decide what its first step should be.

        Specifically, implementations of this method should either pass on all threads after the first one directly to
        get_next_step (for algorithms where the next step can be determined without information from the first one), or
        else idle each thread after the first one.

        Parameters
        ----------
        thread : Thread()
            Thread object to consider
        allthreads : list
            List of all thread object to consider
        settings : argparse.Namespace
            Settings namespace object

        Returns
        -------
        first_step : list or str
            List of first mutation(s) to apply to the system, formatted as a list of strings of format
            <resid><three-letter-code>; alternatively, 'TER' if terminating or 'IDLE' if doing nothing. Finally, can
            also return 'WT' to indicate that the simulation should roceed with the structure passed directly to isEE,
            without mutation.

        """

        pass

    @abc.abstractmethod
    def get_next_step(self, thread, allthreads, settings):
        """
        Determine the next step for the thread, or return 'TER' if there is none or 'IDLE' if information from other
        threads is required before the next step can be determined.

        Implementations of this method should always begin by loading a global history object from
        algorithm_history.pkl, which is a namespace object with the same attributes as a thread.history object, which
        should have been created by JobType.update_history the first time it was called. They should then end by dumping
        their copy of the algorithm history to that file. The purpose of this file is to keep track of the global
        history of the algorithm so far in order to facilitate parallelization among multiple threads.

        In general, whenever a new item is added to any of the lists in the algorithm history object, something should
        also be added to the corresponding index of every other list, even if only a placeholder. If the new item is
        being appended to the end of a list, then this requirement can be satisfied by calling buffer_history() on the
        algorithm history object after appending the desired information.

        In the case where necessary information for get_next_step is still being awaited from an assigned thread, this
        function should return the string 'IDLE', to indicate that the thread should do nothing. The functions of 'IDLE'
        and 'TER' are defined by process.process().

        Parameters
        ----------
        thread : Thread()
            Thread object to consider
        allthreads : list
            Full list of Thread() objects for this job (including thread)
        settings : argparse.Namespace
            Settings namespace object

        Returns
        -------
        next_step : list or str
            Next mutation(s) to apply to the system, formatted as a list of strings of format
            <resid><three-letter-code>; alternatively, 'TER' if terminating or 'IDLE' if doing nothing.

        """

        pass

    @abc.abstractmethod
    def reevaluate_idle(self, thread, allthreads):
        """
        Determine if the given idle thread is ready to be passed along to the next step or should remain idle. The
        result from this method can be used as a stand-in for the results of a jobtype.gatekeeper implementation when
        called on an idle thread.

        Parameters
        ----------
        thread : Thread()
            Thread object to consider
        allthreads : list
            List of all thread objects to consider

        Returns
        -------
        gatekeeper : bool
            If True, then the thread is ready for the next "interpret" step. If False, then it should remain idle.

        """

        pass


class Script(Algorithm):
    """
    Adapter class for the "algorithm" that is simply following the "script" of mutations prescribed by the user in the
    settings.mutation_script list.

    This particular implementation of Algorithm has no conditions that return 'IDLE'.

    """

    def get_first_step(self, thread, allthreads, settings):
        # if this is the first thread in allthreads and there are no history.trajs objects in any threads yet
        if thread == allthreads[0] and not any([bool(item.history.trajs) for item in allthreads]) and not settings.skip_wt:
            return 'WT'
        else:
            return Script.get_next_step(self, thread, allthreads, settings)

    def get_next_step(self, thread, allthreads, settings):
        algorithm_history = self.build_algorithm_history(allthreads, settings)

        untried = [item for item in settings.mutation_script if not item in algorithm_history.muts]
        try:
            return untried[0]     # first untried mutation
        except IndexError:  # no untried mutation remains
            return 'TER'

    def reevaluate_idle(self, thread, allthreads):
        return True    # this algorithm never returns 'IDLE', so it's always ready for the next step

class Random(Algorithm):
    """
    Adapter class that chooses new mutants to explore completely randomly, excluding already attempted mutants.
    """

    def get_first_step(self, thread, allthreads, settings):
        # if this is the first thread in allthreads and there are no history.trajs objects in any threads yet
        if thread == allthreads[0] and not any([bool(item.history.trajs) for item in allthreads]) and not settings.skip_wt:
            return 'WT'
        else:
            return Random.get_next_step(self, thread, allthreads, settings)

    def get_next_step(self, thread, allthreads, settings):
        if not settings.DEBUG_TERTIME == 0 and settings.DEBUG_TERTIME < time.time():    # just for debugging something
            return 'TER'
        # todo: implement any sort of termination criterion?

        algorithm_history = self.build_algorithm_history(allthreads, settings)

        all_resnames = ['ARG', 'HIS', 'LYS', 'ASP', 'GLU', 'SER', 'THR', 'ASN', 'GLN', 'CYS', 'GLY', 'PRO', 'ALA',
                        'VAL', 'ILE', 'LEU', 'MET', 'PHE', 'TYR', 'TRP']

        WT_seq = [str(int(str(atom).replace('-CA', '')[3:]) + 1) + str(atom)[0:3] for atom in
                  mdtraj.load_prmtop(settings.init_topology).atoms if (atom.residue.is_protein and (atom.name == 'CA'))]

        # First, build list of all possible single mutants:
        single_muts = [[str(int(resid + 1)) + res for res in all_resnames] for resid in range(len(WT_seq))
                       if not resid + 1 in settings.immutable]
        single_muts_temp = []
        for item in single_muts:  # combine lists
            single_muts_temp += item
        single_muts = single_muts_temp

        for item in WT_seq:  # remove wild type sequence from single_muts
            if item in single_muts:
                single_muts.remove(item)

        len_all_combs = int(sum([factorial(len(single_muts)) / (factorial(r) * factorial(len(single_muts) - r)) for r in
                                 range(settings.max_plurality + 1)]))
        if len(algorithm_history.muts) == len_all_combs:  # no more permissible mutants
            return 'TER'

        # Implement Monte Carlo algorithm
        comb = []
        while not comb:  # repeat until we find a valid combination
            for i in range(max([int(len(single_muts) / 10), 10])):
                comb = []
                while not comb or set(comb) in [set(item) for item in algorithm_history.muts]:  # todo: this algorithm will become slow as the sample size becomes a significant fraction of the possibility space, unsure if I should care
                    comb = []  # reset at the start of each loop
                    for null in range(numpy.random.randint(settings.min_plurality, settings.max_plurality + 1)):  # todo: distribute this non-equally?
                        random_index = -1
                        while random_index < 0 or single_muts[random_index][:-3] in [item[:-3] for item in comb]:
                            random_index = numpy.random.randint(0, len(single_muts))
                        comb.append(single_muts[random_index])

            if comb and settings.stability_model:  # check for stability of proposed mutant
                stability_model = factory.stability_model_factory(settings.stability_model)
                if settings.destabilization_cutoff >= stability_model.predict(settings.initial_coordinates[0], settings.init_topology, list(comb)):
                    comb = []

        return list(comb)

    def reevaluate_idle(self, thread, allthreads):
        return True  # this algorithm never returns 'IDLE', so it's always ready for the next step


if __name__ == '__main__':
    # Just testing stuff, this should never be used in an actual run
    jobtype = factory.jobtype_factory('isee')
    thread = argparse.Namespace()
    settings = argparse.Namespace()
    settings.working_directory = './'
    settings.initial_coordinates = ['data/one_frame.rst7']
    settings.init_topology = 'data/TmAfc_D224G_t200.prmtop'
    settings.covariance_reference_resid = 260
    settings.immutable = [218, 260] + [ii for ii in range(350, 442)]
    settings.max_plurality = 12
    settings.min_plurality = 12
    settings.plural_penalty = 1 #0.85
    settings.rmsd_covar = [2.409151127122221, 1.9845629910244407, 1.6409770071516554, 1.2069303829003601, 1.1774104996742552, 1.1414646722509711, 1.6277495210688269, 1.7178594550668773, 1.5573849264923196, 1.6222616936921075, 2.139257339920909, 2.0355790819694772, 2.0670293620816573, 2.0509967163136467, 2.3718218856798265, 2.725089368986255, 2.6368147154575774, 2.2260307731037887, 2.5657250982226993, 2.8139534330304286, 2.4156272283232267, 2.1185344711518828, 1.7966304212974096, 1.495615071178279, 1.0665349714934624, 0.62622389357105, 0.41661834767956046, 0.33238540271523664, 0.5266911884945009, 0.8603863335688027, 1.2961691493210963, 1.3149073997096503, 1.087848067605599, 1.3980104852851503, 1.7593319546319754, 1.436783423351133, 1.6735766023866119, 1.659270934828769, 1.6463508914996159, 1.3603778641613482, 1.546642862021431, 1.6936155445917012, 1.7574027783441544, 1.7342633556126716, 1.4121915763550301, 1.6524978464246118, 1.6129345105308022, 1.7510921326857187, 1.6307382119143303, 2.154613585911121, 2.156907141198287, 1.725686467726665, 1.8871688474992456, 2.288093265532053, 2.0877648683775574, 1.6154540392514545, 1.2448718430424812, 1.1084758370905747, 1.2047520664808857, 1.1741113899674331, 1.4872555630783437, 1.8148940365888924, 1.9624111544303016, 2.1034988697779524, 2.360337367399307, 2.692259080940934, 2.8470775307056306, 2.950603784872049, 3.4151379698075583, 3.474764677558553, 3.0042367509224386, 2.794936213963638, 2.5714411950890304, 3.0499310290409363, 3.2612488945150475, 2.863811282105302, 2.925291191195256, 3.4715106088872867, 3.383187628094896, 2.9925172066994925, 3.078550824373092, 3.5276035223482007, 3.6829112885280217, 3.4131480501006095, 2.975245985352162, 2.571204510725302, 2.026526280899236, 2.1705092136262616, 2.326261646306838, 1.895536239749708, 1.612370269197273, 1.7443935702592746, 1.6920798364775174, 1.2961090424218624, 1.1321792266686648, 0.7915806960860717, 0.8425653659092112, 0.6276519979344761, 0.5508744129335861, 0.768327956507145, 0.9515841437190209, 1.3683201822203757, 1.1712898660583293, 1.0250037039145583, 1.462105078021419, 1.686335245815784, 1.4334984152413268, 1.5708144546831355, 2.052314708700172, 2.038735418403611, 1.8131322137605297, 2.2361306675822448, 1.9378961600424507, 2.190625937133571, 1.787426435393106, 1.2976323362080442, 0.8843715780240702, 0.5228643898438631, 0.2753030693786261, 0.2931480175456677, 0.5604517800951593, 0.7368742584734022, 1.146118756018045, 1.3938059917073062, 1.243670202048456, 1.0929851633070546, 0.6902545909682317, 0.4509524708374205, 0.4099084244249464, 0.48631771273118757, 0.4715336762950448, 0.5811367041156211, 0.7975226590410984, 1.1169981318533977, 1.299278924163868, 1.353778763719865, 0.9631797347498909, 0.7778446523137179, 0.7216003876497649, 1.0897635260493828, 1.2557407427022462, 1.1505684470150626, 0.974161238371388, 0.765544297085648, 0.6202937686638598, 0.6299491474448252, 0.533219470284319, 0.6634134980282882, 0.9140578548659407, 1.0212212614103604, 1.0497564563853488, 1.3142009062245186, 1.5813282681733314, 1.5830934657785871, 1.6888525831972165, 2.041600833871994, 2.202213598464107, 2.2170919890886727, 2.4894435986060954, 2.1440275518530654, 2.011302037148799, 1.5824784649567594, 1.231793242869025, 0.7597221084756685, 0.36574467108503445, 0.2426180788866027, 0.5627032962183899, 0.5795986303204472, 0.8688177925717158, 1.2011175714700462, 1.4146732517566771, 1.9096741487027076, 1.9433687887027842, 1.9567698304844803, 2.1513283474511753, 2.6395192768727735, 2.7210238754183904, 2.3910420851463057, 1.9554328922013766, 2.2069190451734633, 1.8735633762533965, 1.3824968668748268, 1.484110887336474, 1.7264446496422705, 1.263819735591131, 1.137969578171706, 1.6015427254654742, 1.855641517560868, 1.6201557567389, 1.7824029538859167, 2.233427836014575, 1.972627246902554, 1.5540540530661726, 1.454894772726324, 1.2679617701951245, 0.851167820885501, 0.7074277161313157, 0.8157599686830077, 0.5664890398789291, 0.2701987764217301, 0.30107691606115244, 0.3163238299766273, 0.3685276068911975, 0.5869647407044045, 0.5324965531897539, 0.6922229012169531, 1.075642431945148, 1.1037202279417786, 1.043566513472472, 1.285622065768631, 1.5418846846427452, 1.441335928529488, 1.6879800474718196, 1.3087187529675477, 0.9074400812850182, 0.606917811710234, 0.21780464531946842, 0.2931568214704138, 0.630420891795513, 0.47699966727896803, 0.5547007426080225, 0.9079568930591633, 0.9576678785506575, 0.8452547790687734, 0.5146568706079566, 0.25019025232020536, 0.12654961819877517, 0.19407457914767376, 0.39592379771816294, 0.7609881678155445, 0.8226159250285602, 0.7910535454423352, 1.0979486806926377, 1.4144657484975933, 1.4041643681220197, 1.4528880499993324, 1.8154599187476193, 2.057022967204504, 2.0486456185722006, 2.128184414690382, 2.5047800212121523, 2.3410527214924435, 1.8597505389873465, 1.5375925603209826, 1.1361808589684976, 0.8608179507169844, 0.566523921777011, 0.13140513506918064, 0.07643223807901155, 0.10645400076481107, 0.27127477466451755, 0.6354730069052482, 0.9679355477112124, 1.0532284826810165, 1.446919662451723, 1.2985703066987582, 0.9327201827220035, 0.6363473270869764, 0.37707143679148636, 0.0, 0.11117381419196652, 0.5411263237797045, 0.5286961430415533, 0.22935841839544527, 0.275447433171813, 0.7361389457070859, 1.177071043945897, 1.0809764298283535, 1.7106017113521783, 1.4211480726764867, 1.745983440371482, 1.6773188079734158, 1.767584681909782, 1.3388175680809646, 1.1410965551353385, 0.8290482834641163, 0.4444106899569443, 0.17953802145559583, 0.21121193561551252, 0.22396349255208844, 0.47577322321209203, 0.7254588799332433, 0.7779900420863387, 0.5091329394865541, 0.704930093038628, 0.7941977888666573, 1.0859653157623244, 1.3038849488772195, 1.3662859941822716, 1.0967829645658635, 0.9715060851202962, 0.6464536128311025, 0.5987589389687789, 0.6077785965036145, 0.3611553664445629, 0.28945654938493703, 0.6287773515956873, 1.0913176823927886, 1.2927752008730398, 0.9581180586300749, 1.090772797990753, 1.6000260470249799, 1.4915620210429987, 1.2186495598768083, 1.6242335959151546, 1.9598572280225335, 1.6580632496121226, 1.5876816490809949, 2.0526459274293436, 2.205011437726407, 1.9124488341758386, 2.061931045203091, 1.7368555755010284, 1.6297440578280893, 1.3103347862305699, 1.0151235722146335, 0.8059030226851512, 0.42043910934251677, 0.3848455687139601, 0.3590899804844901, 0.4141687821613241, 0.6995170831116444, 0.9826975673263751, 0.8194318159601455, 0.5092828602216573, 0.44155123531463325, 0.40661681691716456, 0.4052022342291145, 0.42814104434780625, 0.3833280005477048, 0.48796059768077726, 0.7969440863216863, 0.9240711606558762, 0.9388317747586525, 1.1852224388384098, 1.4941374463147943, 1.5350917933650963, 1.6461872916744764, 1.9756075334098921, 2.2140541337950506, 2.2657139881707336, 2.4478347857196536, 2.7996263736786267, 3.0125933434981746, 2.9837737832734907, 2.901736796335611, 3.328834576772897, 3.1183658185321623, 2.73079868631526, 2.9248641891616955, 3.202118688821212, 3.1508518551592273, 3.0757607896187644, 2.9124350891931727, 3.1606175487987334, 3.0644915574986578, 2.68605399609867, 2.3542348162145617, 2.1849715516617287, 2.275185568753259, 2.096041955897397, 2.3356478796927447, 2.2685770315015517, 1.989714490289354, 1.655795013026901, 1.8393096927883619, 1.8222212723938744, 2.2501436955325262, 2.34359750439921, 2.7869629794529933, 2.994581454486375, 3.4884471102193877, 3.8101881861875384, 4.262079394553044, 4.548642999489229, 4.238979005642821, 3.7733912293831864, 3.393830151211838, 3.1279857232921846, 2.674127307105069, 2.4175490279309866, 1.8687369311056565, 1.8869166551159962, 2.2636483829618825, 2.4023400423920345, 2.359951525860288, 2.8270711837342333, 3.1686124840254863, 3.4165226314876707, 3.17498211022827, 2.972637284797088, 2.8439145950543123, 2.5982878703543366, 2.6852835930288164, 3.138381512010268, 3.6813803788278388, 4.023955267760979, 4.448302859331009, 4.606420642764985, 5.105594724983221, 4.858767352349876, 4.353596739482308, 4.284761409624414, 3.897795345829884, 3.7587067732613675, 3.712495038883498, 4.157025759828298, 4.5783732940487525, 4.581539121601806, 4.712942631238383, 4.418594026092707, 4.566824163913762, 4.32212482190961, 4.29595069676135, 4.136138750950603, 3.9951211532035513, 3.7980376683223147, 3.3070358891562206, 3.2361462063231747, 3.5892744274502157, 3.6740648671037346, 3.739058741077101, 3.8511284283921143, 3.5419133044670237, 3.594028848231803, 3.196600532786006, 3.356354056977835, 3.402522717011724, 2.9375519598739452, 2.7466527424402343, 2.963536795599139, 2.773755925825151, 2.2682060292216804, 1.9490444166797825, 2.3281984080601004, 2.7692411688918215, 3.0911847946978086, 3.546368506384973, 3.9041936075483354, 4.165817729092229, 4.584574723822292]
    settings.algorithm = 'predictor_guided' #'monte_carlo'
    settings.stability_model = '' # 'ddgunmean'
    settings.destabilization_cutoff = -3
    settings.shared_history_file = 'shared.pkl'
    settings.imported_history_file = 'data/resampled_history.pkl'
    jobtype.update_history(thread, settings, **{'initialize': True})
    # thread.history.trajs = ['extant']
    for i in range(1):
        test = PredictorGuided()
        next = test.get_next_step(thread, [thread], settings)
        print(next)
