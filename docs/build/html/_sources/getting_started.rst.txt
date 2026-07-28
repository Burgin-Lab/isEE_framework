Getting Started
===============

Installation
------------

isEE_framework has a couple of unclean dependencies that have not been resolved yet. I have streamlined it as much as possible without refactoring large portions of the code for now. The following instructions are the best option for now until pAPRika and pytraj can be excised from the dependency chain.

First, ensure that Amber22 is sourced (so that pytraj is available):

``source /path/to/amber22/amber.sh``

Then, create a new conda environment that uses exactly Python version 3.11:

``conda create --prefix=/path/to/new-environment python=3.11``

isEE_framework itself can be built from the pre-packaged wheel file located in ``isEE_framework/dist``:

``pip install isee_framework-[version_number]-py3-none-any.whl``

After, you must also install the pAPRika package from conda-forge. Note that this package has a name clash with an unrelated package on PyPI, so you cannot use pip to install it.

``conda install -c conda-forge paprika``

You may need to follow this up with updating the version of numpy:

``conda update numpy``

Finally, install pyrosetta from its dedicated conda channel:

``conda install pyrosetta -c https://levinthal:paradox@conda.graylab.jhu.edu``


Invoking isEE_framework
-----------------------

isEE_framework is invoked on the command line as:

``python /path/to/isee_framework.main.py [config] [[working_directory]]``

The first argument, which points to a configuration file (see :ref:`configuration_file`) is required. The second, which defines the directory in which isEE_framework will run its simulations and analyses, is optional (it can also be defined in the configuration file), but setting it on the command line at runtime may be necessary on certain compute cluster configurations where working/scratch space is defined only when a job begins.

On batch-managed systems it is typical to run isee_framework as a batch job, but there is no need to allocate extensive resources to that job. isEE_framework will submit additional jobs with their own resource requests. You need only ensure that the main isee_framework job has sufficient walltime for all the simulations that you want to run, since it will persist during the entire runtime of the program.

Settings up input and template files
------------------------------------

In addition to the configuration file, the most important files to setup in advance are MD input files and batch job template files. These files should be located in two separate directories that will be indicated in the configuration file.

The input files directory should contain all of the MD input files for your simulations. For simulations using Amber, typically this consists of three files: ``min_amber.in`` (an energy minimization), ``heat_amber.in`` (an NVT heating simulation), and ``isee_amber.in`` (an equilibration/production simulation). These should be configured however you desire for your simulations.

The template files directory is used by isee_framework to set up batch jobs. Template slots denoted by double curly braces (``{{ example }}``) are parsed and replaced as appropriate for each simulation by isEE_framework. Template files are named as ``[md_engine]_[batch_system].tpl``; so for example, for Amber simulations on a Slurm batch system, the template file should be named: ``amber_slurm.tpl``. There is also a separate, optional file for use with the charge initialization option (see options starting with ``ic_`` in the configuration file documentation), named as ``[md_engine]_initialize_charges_[batch_system].tpl`` (so, e.g., ``amber_init_charges_slurm.tpl``)

The recommended way to put together these files is to modify the examples provided in ``isEE_framework/isee_framework/data``. Two Slurm examples are provided: one configured for use with NVIDIA MPS (see :ref:`NVIDIA_MPS`) and one without. In either case, you should modify the preamble (lines starting with #SLURM) and add anything else necessary for compatibility with your batch system, but leave the rest of the content untouched unless you are sure that you know what you're doing.