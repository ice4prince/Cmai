###########################
The Step-by-Step Pipeline
###########################

In the :doc:`Quickstart Guide <../quickstart>`, we showcased the all-in-one approach by
running both ``runEmbed`` and ``runBind`` at the same time. However, we also support
running each of the step separately. In this tutorial, we will go through an example
of this approach.

--------------------------------

************************************
The Embedding Step
************************************

The embedding step uses RoseTTAFold to embed the XXX. This can be done in two ways:
either in one step or two steps. The one-step process is as follows:

.. code-block:: shell

    python Cmai.py --code '<path_to_Cmai_dir>' --input '<path_to_input.csv>' --out '<path_to_output_dir>' --rf_data '<path_to_RoseTTAFold_database_dir>'  --runEmbed


Note that the difference between this procedure and the all-in-one approach is that
the ``--runEmbed`` flag is included.

If you wish to have a more detailed breakdown of the pipeline, we can do it in two steps:

.. code-block:: shell

    python Cmai.py --code '<path_to_Cmai_dir>' --input '<path_to_input.csv>' --out '<path_to_output_dir>' --rf_data '<path_to_RoseTTAFold_database_dir>'  --runEmbed --gen_msa --use_cpu
    python Cmai.py --code '<path_to_Cmai_dir>' --input '<path_to_input.csv>' --out '<path_to_output_dir>' --rf_data '<path_to_RoseTTAFold_database_dir>'  --runEmbed --run_rf

Here, we first use the ``--gen_msa`` flag and then the ``--run_rf`` step. Running
both is functionally equivalent to the one-step embedding step above.

************************************
The Binding Step
************************************

After you have successfully run the embedding step, you can proceed to run the binding step:

.. code-block:: shell

    python Cmai.py --code '<path_to_Cmai_dir>' --out '<path_to_output_dir>' --skip_check --runBind


Here, we use ``--skip_check`` because we already did all the checking in the previous step.
Also, it is important to have the same output directory as the previous step as supplied
in the ``--out`` parameter.

After the binding step, you will see the ``binding_results.csv`` file in your directory,
and our pipeline is officially finished!
