#####################
CLI Options in Depth
#####################

While the :doc:`Quickstart Guide <../quickstart>` shows a few examples of the common usage
of Cmai, there are many more options available in our package. Namely, you can preciesly
control what's being run and how the various steps are steps at your finger tip. Here,
we details of each argument and provide some examples on common usages.

----------------------------

**********************
Arguments Overview
**********************

Cmai Arguments
-----------------

================================== ================== ========== ===========================================================================================================================
Argument                            Purpose            Input       Functionality                                                         
---------------------------------- ------------------ ---------- ---------------------------------------------------------------------------------------------------------------------------
``-h`` or ``--help``                Documentation      None       Show the help message and all CLI options.
``--code``                          Config             ``str``    The absolute path to the Cmai directory.
``--input``                         Config             ``str``    The absolute/relative path to the ``input.csv`` file.
``--out``                           Config             ``str``    The absolute/relative path to the output directory.
``--env_path``                      Config             ``str``    The file saving the directory of the Conda environments: ``runEmbed``, ``runBind``, and RoseTTA enviroments in order.
``--rf_data``                       RoseTTAFold        ``str``    The directory that store the RosettaFold structure and sequence databases. 
``--fasta``                         RoseTTAFold        ``str``    The path to the FASTA files for antigen sequences. 
``--pre_dir``                       RoseTTAFold        ``str``    The directory to save the preprocessed data, defaults to the output directory.
``--npy_dir``                       RoseTTAFold        ``str``    The direcrtory for ``npy`` files if it is different from the ``pre_dir``. 
``--cpu``                           RoseTTAFold        ``int``    The maximum number of CPU cores for embedding.
``--mem``                           RoseTTAFold        ``int``    The amount of system RAM available.
``--use_cpu``                       RoseTTAFold        None       Whether to use CPU instead of GPU.
``--seed``                          Config             ``int``    The seed used in random number generator to ensure reproducibility.
``--min_size_background_bcr``       Config             ``int``    The initial and minimum sample size of background BCRs, defaults to 100.
``--max_size_background_bcr``       Config             ``int``    The maximum size for subsample of background BCRs, which should no more than 1000000, defaults to 10000.
``--rf_para``                       RoseTTAFold        None       Use the parameters from ``paras/rf_para.txt`` for antigen embedding, defaults to False.
``--gen_msa``                       Mode               None       Generate MSA only and exit.
``--run_rf``                        Mode               None       Skip generating msa and running embedding prediction, defaults to False.
``--skip_preprocess``               Mode               None       Skip preprocessing of antigen embedding, defaults to False.
``--skip_extract``                  Mode               None       Skip extracting ``NPY`` for antigen embedding, defaults to False.
``--runEmbed``                      Mode               None       Run the embedding step only, defaults to False.
``--runBind``                       Mode               None       Run the binding step only, defaults to False.
``--skip_check``                    Mode               None       Skip the proprocessing step and data-checking step (use only if already done), defaults to False.
``--species``                       Mode               None       Match trhe species of background BCR to the target BCR.
``--suffix``                        Config             None       Add suffix to antigen id, which is used to distinguish antigens with the same name, defaults to False.
``--no_rank``                       Config             None       Only export the predicted score but without ranking with background BCRs, defaults to False.
``--verbose``                       Config             None       Enable verbose output, defaults to False.
``--merge``                         Config             None       Merge the output to the ``input.csv`` dataset, defaults to False.
================================== ================== ========== ===========================================================================================================================

As you may have noticed, there are in general four types of arguments (or called flags),
and here are their overall functions:

1. Documentation: This renders the CLI help message for Cmai.
2. Config: The overall configurations of Cmai regarding specific values and IOs.
3. RoseTTAFold: The parameters used for RoseTTAFold's embedding steps.
4. Mode: These parameters are toggles to enable or disable different modes of Cmai.

For those comments that do not have an input type, they are known as flags, which
will toggle the underlying switch if included. Otherwise, ``str`` or ``int`` inputs
are to directly follow the flag. All parameters are briefly explained in the 
"Functionality" column with similar language to the default help message. Some of the
parameters will be further explained in the tutorials and examples below.


RoseTTAFold Parameters
-----------------------

Many of the parameters of Cmai are passed to RoseTTAFold for computing embeddings.
There are two main ways of passing such parameters: the first is to use CLI arguments
as shown above; the second is the usage of a ``rf_para.txt`` configuration file that
stores the necessary parameters. Here, we focus a bit on the latter.

First of all, to enable Cmai to read the customized ``rf_para.txt``, the ``--rf_para``
flag should be specified. The file is usually stored in the ``paras/`` directory,
and it contains two self-explanatory columns: ``para`` for the name of parameters
and ``arg`` for the argument passed into the parameters. Here is a table of all the
parameters:

================== =========== ============================= ================================================================================
Argument             Cmai CLI     Default                       Functionality                                                         
------------------ ----------- ----------------------------- --------------------------------------------------------------------------------
runtime.cpu          Yes         ``32``                        The maximum number of CPU cores for embedding.
runtime.mem          Yes         ``128``                       The amount of system RAM available.
exe.hhsuite-path     No          Automatically generated       The path to ``hhsuite``.
exe.psipred-path     No          Automatically generated       The path to ``psipred``.
exe.csblast-path     No          Automatically generated       The path to ``csblast``.
exe.blast-path       No          Automatically generated       The path to ``blast``.
exe.CSBLAST-DATA     No          Automatically generated       The path to ``csblast`` data.
exe.PSIPRED-DATA     No          Automatically generated       The path to ``psipred`` data.
================== =========== ============================= ================================================================================

Notice that two of the parameters are exposed in the Cmai CLI with their own arguments.
Cmai CLI arguments take priority in each run, but they will not overwrite the configurations
stored here.

All the "automatically generated" parameters are generated upon installation via the
``./get_env_path.sh`` command. In general, there is no need to change the paths to RoseTTAFold
binaries unless the installation on the system has changed. In the latter case, re-running
the aforementioned command is recommended.