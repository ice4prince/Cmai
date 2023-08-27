###########################
Working with Limited VRAM
###########################

Since Cmai both uses RoseTTAFold and has a deep learning model under the hood, it is necessary
to have a GPU for inference. In the :doc:`Installation Guide <../installation>`, we briefly
described the requirements: a Nvidia GPU with *sufficient* VRAM. The definition of "sufficient"
here is a bit vague: while state-of-the-art GPUs such as A100 with 40G will definitely be
an ideal chocie, consumer GPUs with anything less than 8GB will be a poor choice.
While this is not a problem for servers or cloud nodes with modern enterprise GPUs, average
users at home may not have access to a high-end consumer-grade GPU. Thus, in this tutorial,
we go more in depth on a work-around to enable the inference of Cmai using a GPU with
8GB of VRAM.

.. note::

    While we were not able to test all GPUs on our end, you are certainly welcomed to try your
    configuration using Cmai. If your GPU has more than 8GB of VRAM (*e.g.* 12GB or 16GB), you
    can start with the standard pipeline. If your GPU has less than 8GB of VRAM, follow this
    guide, but we can't guarantee that it will work.


-------------------------

************************************
Standard Pipeline: Out of Memory
************************************

If you follow the example in the :doc:`Quickstart Guide <../quickstart>`, you will be running
both the embedding step and binding step together. Under this setting, GPU is used throughout,
which can cause problems for the embedding step if your VRAM is 8GB or less. While we won't
go through the example again, we will instead show you the error you will receive:

.. code-block:: text

    RuntimeError: CUDA out of memory. Tried to allocate 3.38 GiB (GPU 0; 7.43 GiB total capacity; 3.22 GiB already allocated; 2.75 GiB free; 3.78 GiB reserved in total by PyTorch)


As seen, it is very clear that your GPU is running out of memory if you receive this error.
In fact, this error stems from RoseTTAFold. If you encounter something different, then it's
still worth a shot to try the pipeline here, but you may want to investigate further or
head over to GitHub to talk about the issue.


***************************************
Workaround: Using CPU and Limiting RAM
***************************************

The core idea of enabling Cmai on limited-VRAM GPU is through a two-step process:

1. Run the embedding step using CPU and system RAM only
2. Run the binding step using GPU

This works because the latter is less VRAM intensive, and we have shown that it can work
for only 8GB of VRAM. To get started, you will also need to know the amount of system RAM
so that RoseTTAFold knows what it is working with.

The first `runEmbed` step looks like this:

.. code-block:: shell

    python Cmai.py --use_cpu --mem 32 --runEmbed --code '<path_to_Cmai_dir>' --input '<path_to_input.csv>' --out '<path_to_output_dir>' --rf_data '<path_to_RoseTTAFold_database_dir>'


In the above code snippet, you will need to replace the desired paths. Also, you will want to change the ``--mem``
value to slightly less than your acutal system RAM avaiable so that RoseTTAFold can know the amount to use. Once
this is successful, you can run the `runBind` step:

.. code-block:: shell

    python Cmai.py --code '<path_to_Cmai_dir>' --out '<path_to_output_dir>' --skip_check --runBind

Note: it is important to use the same ``--out`` for both steps so that the results from the first step can be
used in the second.

If all is good, you have successfully run Cmai with limited VRAM!