.. _configuration_file:

The Configuration File
======================

The configuration file is the primary means of controlling the behavior of isEE_framework. In order to support the wide array of functionality that any given user may need, the configuration file supports many options and can be quite long; however, in most cases a relatively short configuration file will be sufficient. This page provides some recommendations for building the configuration file for a handful of common use cases, and provides detailed documentation for each setting.

The contents of the configuration file are read line-by-line into isEE_framework as literal python code, which enables invocation of python built-in functions as well as methods of pytraj and numpy (and anything else you may wish to import). This means comments can be included in-line or on their own lines preceded by a '#' character, and blank lines are simply ignored. **Warning**: This input is not sanitized in any way. For this reason among others, "shutil.rmtree('/')" makes for a poor working directory!!

.. toctree::
   :maxdepth: 2
   :caption: Contents of this page:

   configuration_file.rst


Core Settings
-------------

Certain settings should be given for every job. The following settings should be set in every configuration file:

``batch_system``

    Indicates the type of batch system on which ATESA is running. Supported options are 'slurm' and 'pbs' (the latter is also known as TORQUE).

``restart``

    Indicates whether this is a new job (False) or a continuation of an old one *in the same working directory* (True). Restarting is not supported for job_type = 'find_ts'.

``overwrite``

    Indicates whether to delete the existing working directory (if one exists) and create a new one *if and only if* restart = False (has no effect otherwise). If restart = False, overwrite = True will *always* delete the working directory if it exists; conversely, overwrite = False will *never* delete it (regardless of the "restart" setting).

``init_topology``

    An absolute or relative path given as a string and pointing to the simulation topology file.

``initial_coordinates``

    Path to initial coordinate file, corresponding to the topology in ``init_topology``. This model is treated as the template from which all mutants are formed. Formatted as a list of strings, though isEE_framework does not currently support multiple independent coordinate files.

``working_directory``

    An absolute or relative path given as a string and pointing to the desired working directory (this can be omitted if the working directory is set in the command line). This is the directory in which all of the simulations will be performed. It will be created if it does not exist.

``md_engine``

    The MD engine to use for running simulations. This also controls which template and input files ATESA looks for (see :ref:`SettingUpSimulationFiles` Supported options are 'amber' or 'cp2k'. Note that support for CP2K is experimental. Default = 'amber'

``algorithm``

    The algorithm used to select mutations to apply and simulate. Options are: ``random`` and ``script``. See below for additional options for controlling each.


.. _FilePathSettings:

File Path Settings
~~~~~~~~~~~~~~~~~~

These settings define the paths where isEE_framework will search for user-defined input files and template files.

``path_to_input_files``

    Absolute path (as a string enclosed in quotes) to the directory containing the input files. The default is the directory 'data/input_files' located inside the ATESA installation directory.

``path_to_templates``

    Absolute path (as a string enclosed in quotes) to the directory containing the template files. The default is the directory 'data/templates' located inside the ATESA installation directory.


Other Settings
--------------

The following options may be important for your application, but they have sensible defaults and will not be necessary in every configuration file.


Basic Optional Settings
~~~~~~~~~~~~~~~~~~~~~~~

``degeneracy``

    An integer indicating how many independent threads to spawn from each set of initial coordinates. For example, if one set of initial coordinates is given and degeneracy is set to 3, then three independent threads of simulations will begin from the same initial coordinates. Use this option to parallelize isEE_framework across available resources. Default = 0 (just one thread per coordinate file)

``skip_wt``

    Indicates whether a simulation should be run with no mutations applied (False means do run it, True means don't.) This setting functions independent of the choice of ``algorithm``. Note that despite the name of this option, the simulation run in this way may or not be the "wild type" of any given protein; it will just be whatever is provided in the ``init_topology`` and ``init_coordinates`` files. Default = False

``pH``

    A float indicating the pH at which to protonate titratable residues using PROPKA. Default = 7.0

``keep_waters``

    A boolean. If True, all waters in the initial structure are kept during mutation and re-solvation of mutants is skipped. Default = False

``hmr``

    A boolean. Governs whether to apply hydrogen mass repartioning (with the default factor of 3) to each model after mutating. If set to True, MD simulation timesteps up to 4 fs are likely to be stable, but this can break certain models. Default = False

``treat_as_protein``

    A list of strings. Used to force isEE_framework to treat non-standard residue names as part of the protein for the purposes of applying mutations and building models. Useful for models containing modified residues or prosthetic groups, almost always in conjunction with ``paths_to_forcefields``. Default = ['']

``paths_to_forcefields``

    List of paths to Amber force field files to ``source`` in tleap during each model build step. This option enables isEE_framework to handle non-standard or custom force field files, which can be important for models containing non-standard residues, including small molecules. By default, ff19SB, OPC, and GLYCAM_06j-1 are loaded without needing to be specified.

``tleap_extra``

    A string containing any extra lines to execute after building the tleap ``model`` object when building each model. Default = ''

``rosetta_prevent_repacking``

    A list of integers corresponding to residue indices that should not be repacked by Rosetta during mutation. This can be important if you want to ensure that the initial positions of certain residues are not moved from one model to another for whatever reason. Default = []

``shared_history_file``

    Path to an ``algorithm_history.pkl`` file created by an isEE_framework job. If this option is specified, then rather than creating a new ``algorithm_history.pkl`` for itself, this job will treat the specified file as its own (and write its own history to that file, as well). This option enables multiple independent runs of isEE_framework to work cooperatively on the same model without repeating any simulations between them. Default = ''

``name_as_timestamp``

    A boolean. If True, forces each file for each simulation to be named using a timestamp set when the model is built, rather than using the list of mutations applied. This can be useful for very long lists of mutations that might clash with OS limits on file names. Note that even if this option is set to False, timestamps will be used whenever a filename might be longer than 256 characters. Default = False

``restart_terminated_threads``

    A boolean. If this is True, and if ``restart = True``, threads that had individually terminated already (e.g., due to ``max_steps_per_thread``) will be restarted as well as if they were fresh threads. Default = True.

``resubmit_on_failure``

    An integer indicating the number of times that a simulation job should be reattempted if it fails to run. This can be useful if your batch system is occasionally unstable or rejects valid jobs, though excessively high values can cause broken jobs to be repeatedly resubmitted. Default = 1


Algorithm Settings
~~~~~~~~~~~~~~~~~~

These settings are specific to individual selections for the ``algorithm`` option.

For the ``random`` algorithm, which applies randomly selected mutations. A random number of mutations are selected with equal probability across every position from the canonical amino acids, excluding the existing identity of the selected residue (e.g., no ALA to ALA mutations allowed).

``min_plurality``

    An integer setting the minimum number of mutations to apply per simulation, inclusive. Default = 1

``max_plurality``

    An integer setting the maximum number of mutations to apply per simulation, inclusive. Default = 3

``immutable``

    A list of integer positions where mutations will never be applied. The integers should correspond to the positions as they are indexed in the topology file. Use this setting to prevent mutations in positions that you aren't interested in exploring. Default = []

``max_steps_per_thread``

    The maximum number of simulations that each thread may run before terminating. Use this to put an upper limit on the number of simulations performed; otherwise, isEE_framework runs with ``algorithm = random`` will continue indefinitely until the main job terminates.


For the ``script`` algorithm, which applies user-specified mutations in the order that they are specified and then terminates.

``mutation_script``

    A list of lists of strings indicating the sequence of mutations to apply. Each sublist in the top-level list corresponds to an individual simulation, whereas each string in those sublists corresponds to a single mutation to apply, formatted as the integer position of the residue to be mutated followed by the three-letter code of the target amino acid to mutate to. For example,

    .. code-block:: python

        mutation_script = [['12ALA'], ['12ASP', '9ARG']]


    would run two simulations: one where residue 12 is mutated to alanine, and another where both residue 12 is mutated to aspartate and residue 9 is mutated to arginine.


.. _BatchTemplateSettings:

Batch Template Settings
~~~~~~~~~~~~~~~~~~~~~~~

These settings are used to fill in the template slots in the user-provided template files. If you do not wish for isEE_framework to use an option, you can simply omit its template slot from the appropriate file and leave it unset in the configuration file.

Technically, if you modify your batch templates appropriately, you can use each of these template slots however you wish (so long as data types are respected, e.g., ``nodes`` must always be an integer). This documentation only refers to how each option is intended to be used.

``nodes``

    The number of compute nodes to request for simulations, given as an integer. Default = 1

``ppn``

    The number of cores or processes to request per node (ppn: "processes per node") for simulations, given as an integer. Default = 1

``mem``

    The amount of RAM to request for simulations, given as a string of appropriate format for the batch system. Depending on the batch system, this may be interpreted as total memory, or as memory per core. Default = '4000mb'

``walltime``

    The amount of walltime (real time limit for the batch job) to request for simulations, given as a string of appropriate format for the batch system. Default = '02:00:00'

``solver``

    The name of the executable to use to perform simulations, given as a string. Default = 'pmemd.cuda' (which is specific to Amber)

``extra``

    An additional template slot for simulations to be used however the user sees fit. This option is provided in case a user has an unforseen need to template something other than the above options. Default = '' (an empty string)


.. _NVIDIA_MPS:

NVIDIA MPS Settings
~~~~~~~~~~~~~~~~~~~

isEE_framework is configured to maximize the overall throughput of MD simulations (i.e., total ns of simulation across all simulations per real hour) on NVIDIA GPUs by taking advantage of NVIDIA's MPS (Multi-Process Service) framework, available on cards with compute capability version 3.5 or higher. MPS interleaves multiple independent CUDA processes on a single GPU so as to avoid inefficient downtime while waiting for other hardware (such as the CPU, RAM, or writes to disk). This slows down each individual simulation but speeds up the throughput of simulations overall. Note that if you want to use this option, you also need to engage the MPS daemon appropriate in each batch job; the example batch template file in ``isee_framework/data/example_amber_slurm_mps.tpl`` shows one way to do this.

``nvidia_mps``

    An integer controlling the number of jobs that should be co-localized on a single GPU. If this is set to 1, MPS will not be used. Although the relationship between this number and the total simulation throughput  depends on the available hardware and on the speed of each individual simulation, it generally increases asymptotically from the throughput without MPS to an upper limit representing 100% GPU usage. For most protein simulations, 3 is a reasonable choice. Default = 1 (no MPS)

``mps_patient``

    A boolean indicating whether each thread of simulations should execute all of its jobs in each step right away, or if sets of simulations fewer than ``nvidia_mps`` should hold off until another thread produces additional steps that can be used to fill the batch job up to ``nvidia_mps`` independent simulations. In other words, setting this to ``True`` prioritizes maximum efficiency in resource usage at the cost of potentially leaving some resources idle while other jobs finish, whereas setting it to ``False`` prioritizes not making any threads wait on one another at the cost of potentially using less than 100% of a given GPU at a time. In making your selection for this option, consider whether ``degeneracy`` is divisible by ``nvidia_mps``. Default = True


Charge Initialization Settings
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

isEE_framework supports dynamic charge assignments for specified atoms using a single step of QM simulation. This option can be useful if you expect that mutations applied by isEE_framework might polarize small molecules bound to your protein to an extent significant enough to be worth the additional overhead.

``initialize_charges``

    A boolean. If True, a single step of QM simulation is used to initialize charges in ``ic_qm_mask`` before each simulation. Default = False

``ic_qm_mask``

    A string formatted as an Amber-style mask (e.g., ":1@CA" matches the atom named CA in residue 1). Corresponds to the ``qm_mask`` option in Amber. This is the QM mask used during the single step of QM and corresponds to the atoms whose charges will be set dynamically. Default = ''

``ic_qm_theory``

    A string indicating Amber's ``qm_theory`` option. Default = 'DFTB3'


``ic_qm_cut``

    A float indicating Amber's ``qm_cut`` option. The non-bonded cutoff in angstroms. Default = 12.0

``ic_qm_charge``

    An integer indicating Amber's ``qm_charge`` option. The charge of the ``qm_mask``. Default = 0

``ic_dftb_telec``

    A float indicating Amber's ``dftb_telec`` option. The electronic temperature used during DFTB. Default = 0
