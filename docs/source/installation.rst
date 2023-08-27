######################
Installation Guide
######################

To install Cmai and successfully run the prediction pipeline, you will need the following components:

- Conda for package management
- Cmai installation
- RoseTTAFold from RosettaCommons

In this tutorial, we will walk you through each step of the way with an emphasis on the software installation
procedures. For the most efficient prediction, we will also briefly mention hardware requirements.

**********************
Hardware Requirements
**********************

Since Cmai is a deep-learning-based tool and we use RoseTTAFold in our pipeline, we **recommend** a GPU for
computation. Here is a list of GPUs that we tested to have worked well:

- Nvidia A100
- Nvidia V100
- Nvidia P40
- Nvidia P4 (Needs workaround with CPU)

While these are all enterprise accelerator cards, the main thing you should look out for is VRAM.
Comprable consumer-grade Nvidia GPU with sufficient VRAM should work as well. While cards like P40 with 24GB
of VRAM has no issue, Nvidia P4 with 8GB VRAM has some problems with RoseTTAFold (See :doc:`Working with Limited VRAM <tutorial/vram>`
for a workaround using CPU). Thus, we recommend at least 8GB of VRAM for Cmai to work at all, and more VRAM
for the best efficiency.

As for system RAM requirement, we recommend at least as much system RAM as your GPU VRAM. For most use cases,
32GB or 64GB will be safe.

--------------------------------

********************
Conda
********************

Conda is a popular package management system that is especially good at virtual environments for workflows such
as data science and deep learning. If you already have conda installed (miniconda or anaconda), you can skip
this step. Otherwise, you will need to install conda on your system. We recommend
`this guide <https://conda.io/projects/conda/en/latest/user-guide/install/index.html>`_ for detailed instructions
on installing conda on your system, and major operating systems are all supported. From our experiences,
miniconda is sufficient for all our needs.

.. note::

    We will need quite a few conda environments along the way, which can take up more than 10GB of disk space.
    If you are running on a server with restricted space in your home directory, please consider specifying a
    work directory for installing conda.


********************
Cmai
********************

Once you have conda installed, you can start setting up Cmai. First, you will need to pull the Git repository:

.. code-block:: shell

    git clone https://github.com/ice4prince/Cmai.git
    cd Cmai

Here, you will find yourself in our Cmai project. We will then create a conda environment called `runEmbed`,
which is used to run embeddings, and another called `runBind` to run the binding step:

.. code-block:: shell

    conda env create -f models/runEmbed.yml
    conda env create -f models/runBind.yml

The good news is that we have the files all set up for you to ensure reproducibility and ease of installation.
With these environments set up, we've officially installed all our components of Cmai. Next, we will move onto
RoseTTAFold, which is a bit more involved.

*******************
RoseTTAFold
*******************

Since RoseTTAFold is not our software, how it is installed depends on their implementation. To install it,
we will follow their `GitHub README <https://github.com/RosettaCommons/RoseTTAFold>`_. However, not all steps
are needed: namely, we don't PyRosetta for Cmai to work (thus, you can skip Step 6 on their tutorial). The following
steps are all you need:

1. Clone the package
2. Make conda environment
3. Download network weights
4. Install dependencies
5. Install structure databases

While you can get the most up-to-date instructions there, we will still provide the command needed here as well.

.. note::

    For researchers and others who are using a server or cloud instance, it is always best to check whether RoseTTAFold
    is available from the institution. It is in general not very practical to install the databases on home machines given
    the intensive need for storage.

First, you will need to clone the GitHub repository:

.. code-block:: shell

    cd scripts/rfold
    git clone https://github.com/RosettaCommons/RoseTTAFold.git
    cd RoseTTAFold

Then, you can go on to install conda environments:

.. code-block:: shell

    conda env create -f RoseTTAFold-linux.yml
    conda env create -f folding-linux.yml

.. note:: If your GPU is not compatible with CUDA11, please see RoseTTAFold's documentation for different environmnets.

Next, network weights are needed:

.. code-block::

    wget https://files.ipd.uw.edu/pub/RoseTTAFold/weights.tar.gz
    tar xfz weights.tar.gz

The third-party dependencies are need as well:

.. code-block::

    ./install_dependencies.sh

Finally, you will need to install sequence and structure databases. These files are massive:
more than **2TB** of storage space is needed to fully install all the databases after untar.
Again, we encourage the use of existing databases if available. If not, then you can download
the following:

.. code-block::

    # uniref30 [46G]
    wget http://wwwuser.gwdg.de/~compbiol/uniclust/2020_06/UniRef30_2020_06_hhsuite.tar.gz
    mkdir -p UniRef30_2020_06
    tar xfz UniRef30_2020_06_hhsuite.tar.gz -C ./UniRef30_2020_06

    # BFD [272G]
    wget https://bfd.mmseqs.com/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt.tar.gz
    mkdir -p bfd
    tar xfz bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt.tar.gz -C ./bfd

    # structure templates [over 100G]
    wget https://files.ipd.uw.edu/pub/RoseTTAFold/pdb100_2021Mar03.tar.gz
    tar xfz pdb100_2021Mar03.tar.gz

**Importantly**, you will need to take note of the directory in which these databases
are installed since we will need them later for Cmai.


******************
Environment Setup
******************


After you have both Cmai and RosettaFold installed, we need to save the path of your
python interpreters from the environment into a file that Cmai can use:

.. code-block:: shell

    ./get_env_path.sh

You will see a `paras/env_path.txt` file. Now, you're all done! You can start by reading the
:doc:`Quickstart Guide <quickstart>` and trying the included example we have!
