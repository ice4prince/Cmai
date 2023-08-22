####################
Quickstart Guide
####################

In this guide, we will walk you through an included example while offering some insights
on how you can use this interface to use with your own antibodies and antigens. Again, to
get things started, you will need to have Cmai and RoseTTAFold installed as detailed in
our :doc:`Installation Guide <installation>`.

-------------------------

****************
Input Data
****************

For predicting BCR-antigne binding (Antibody-antigen binding, and here we use BCR and antibody
interchangeably), you will need two main ingredients: BCRs and antigens. We require all the
information to be stored in a ``csv`` file named ``input.csv`` so that you can conveniently
run Cmai without the hassle of configuring multiple input files. However, given that we have
a one-file-run-them-all approach, it is critical to have the correct format. Below, we list
all the necessary columns in order:

============ ======================================================== ===================
Column Name     Contents                                                 Example         
------------ -------------------------------------------------------- -------------------
Antigen_id     An identifier (name) for the antigen in this entry.      ``ova``
Antigen_seq    The sequence for the antigen.                            See below.
BCR_id         An identifier for the BCR given the antigen.             ``BCR_11``
BCR_Vh         The **variable (V) region** of the given BCR.            See below.
BCR_CDR3h      The **CDR3 region** of the given BCR.                    See below.
BCR_species    The species that corresponds to the BCR.                 ``Human``
============ ======================================================== ===================

Now that the general format is given, let's unpack this a bit. As we can see, the first two
columns of ``input.csv`` (Yes, I know: we are using rows in the above table) correspond to the
antigen whereas the last four are for BCRs. The id is more or less arbitary and defined by the
user where it makes sense. The antigen sequence is a raw sequence of protein, and since we used
``ova`` antigen as an example, here we will show its protein sequence as an example:

.. code-block:: text

    mgsigaasmefcfdvfkelkvhhanenifycpiaimsalamvylgakdstrtqinkvvrfdklpgfgdsieaqcgtsvnvhsslrdilnqitkpndvysfslasrlyaeerypilpeylqcvkelyrgglepinfqtaadqarelinswvesqtngiirnvlqpssvdsqtamvlvnaivfkglwekafkdedtqampfrvteqeskpvqmmyqiglfrvasmasekmkilelpfasgtmsmlvllpdevsgleqlesiinfekltewtssnvmeerkikvylprmkmeekynltsvlmamgitdvfsssanlsgissaeslkisqavhaahaeineagrevvgsaeagvdaasvseefradhpflfcikhiatnavlffgrcvsp

Of course, your antigen will likely be different, but you get the idea of how to make this
work (See :doc:`CLI Options in Depth <tutorial/cli>` for other options of entering antigen sequences).
Then, you will need to specify your input BCR. Again, the id is arbitary. For example,
``BCR_11`` can stand for the 11th BCR of interest. The most important things you need to
specify are the Vh and CDR3h of the BCR. If you already have the information, you can directly
use them. The Vh region for ``BCR_11`` is:

.. code-block:: text

    QVQLVQSGAEVKKPGASVKVSCKASGYTFTDYYMHWVRQAPGQGLEWMGKVNPNKRGTTYNQKFEGRVTMTTDTSTSTAYMELRSLRSDDTAVYY	

whereas the CDR3h region of the BCR is

.. code-block::

    CARSNALDAW

Cmai makes the following assumptions (requirements) about Vh and CDR3h inputs:

1. The Vh gene sequence must be complete at the N terminal. If the full N terminal is not covered by sequencing,
   corresponding parts of the germline sequence of the Vh gene can be extracted to pad the N terminal. For one
   example on how to do this, please refer to imputation of uncovered V-D-J regions in the
   `mixcr software <https://github.com/milaboratory/mixcr>`_.
2. The C Terminal does not include the "C" amino acid, which is included in the CDR3h sequence.

.. note::

    Cmai allows both uppser- and lower-case letters to represent amino acids in BCRs and antigens.

In the case that you don't have the information or you don't know how to get it from raw sequence,
don't worry: you can use the IgBlast tool which will find the regions for you. Below, we provide a
small detour on introducing the IgBlast pipeline.


IgBlast 
----------

`IgBLAST <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3692102/>`_ is a tool for immunoglobulin
variable domain seuquence analysis, which can be used to find the Vh and CDR3h regions of the BCRs.
Given that you already have the sequence for your BCR, you have two options:

1. Use the `webtool <https://www.ncbi.nlm.nih.gov/igblast/index.cgi>`_: this is a good option for a small scale analysis with a few BCRs.
2. Use a local installation: We recommend the `PyIR <https://github.com/crowelab/PyIR>`_ package, which is not only easy to install but also well implemented with parallel execution.

Given the large-scale nature of antibody-antigen interaction analyses, we recommend the second option
for better efficiency and less stress on the web tool.

.. note::

    The IgBLAST web tool is orders of magnitude slower than a local version with a multicore CPU.
    For large scale computations, we do not recommend the web tool.

----------------------------

**********************
The ``Cmai`` Pipeline
**********************


Once you have your ``input.csv`` ready and all the necessary installations, you are well poised to run the pipeline.
Here are the command-line arguments you will need:


============== ==================================================================================
Argument         Functionality                                                         
-------------- ----------------------------------------------------------------------------------
``--code``       The absolute path to the Cmai directory.
``--input``      The absolute/relative path to the ``input.csv`` file.
``--out``        The absolute/relative path to the output directory.
``--rf_data``    The directory that store the RosettaFold structure and sequence databases. 
============== ==================================================================================

Note from above that the ``--code`` option must be absolute path, where as ``--input`` and ``--out`` options can
support relative paths if they are stored within the Cmai directory. Here is a full example using the included
dataset from Cmai:

.. code-block:: shell

    conda activate runBind
    python Cmai.py --code '<path_to_Cmai_dir>' --input '<path_to_input.csv>' --out '<path_to_output_dir>' --rf_data '<path_to_RoseTTAFold_database_dir>'


As long as all the configurations are correct, Cmai will run and provide you with detailed console outputs
on its progress. Once finished, you will be able to access the results as shown in the section below.

----------------------------

**********************
Output Format
**********************


All the outputs are stored in the directory you specified in the Cmai call with the ``-out`` argument.
If the inference process is successful, there are two things in the directory:

1. A ``binding_results.csv`` that stores the antibody-antigen binding results.
2. An ``RFoutputs/`` folder with intermediate outputs from RoseTTAFold. 

For this quickstart guide and most workflows, the main interest is in the ``binding_results.csv``.
The format of the CSV is the following:

============ ======================================================================================================================================================================================================
Column Name     Contents                                                   
------------ ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
record_id     An identifier (name) for this entry.
Antigen       The antigen name that corresponds to the ``antigen_id`` from the input.
BCR_id        An identifier for the BCR.
Score         The binding score between the given antigen and antibody: smaller means stronger binding.
Rank          The percentile rank of the predicted binding affinity of the query antibody towards the query antigen, as compared to those of a large population of background antibodies against the same antigen.
============ ======================================================================================================================================================================================================

The ``record_id``, ``antigen``, and ``BCR_id`` are used to identify the antigen and antibody
used in the model. What is of interest for most users are the ``Score`` and ``Rank``. The
score is the direct output from the model, and a smaller score indicates better binding. To
provide an easy interpretation of the predicted binding score, we compare this binding score
to the predicted binding scores of a large population of naturally occurring BCRs as background,
and derived a metric called the percentile rank. The exact steps of this procedure is as follows: 

1. First, the antibody is compared against ``min_size_background_bcr`` or 100 background BCRs.
2. If results from Step 1 have met a predetermined threshold of percentile rank (*i.e.* desired performance),
   the antibody is then compared against 1000 background BCRs. (Note that the number always increases by a
   factor of 10.)
3. This process is repeated until the number of background BCRs reach the ``max_size_background_bcr`` or 1,000,000 BCRs, whichever is reached first.

The details of setting the ``max_size_background_bcr`` and threshold are further discussed in the :doc:`CLI Options in Depth <tutorial/cli>`
section.
