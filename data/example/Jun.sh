#!/bin/bash
#SBATCH --job-name=Jun1
#SBATCH -p 256GB
#SBATCH --nodes=1                                                    # number of nodes requested by user
#SBATCH --ntasks=1                                                  # number of total tasks
#SBATCH --ntasks-per-node=1
#SBATCH --time=5-00:00:00                                            # run time, format: D-H:M:S (max wallclock time)
#SBATCH --output=./sbatch_output_%j                                 # redirect both standard output and erro output to the same file
#SBATCH --error=./sbatch_error_%j                                 # standard error output file name
#SBATCH --mail-type=ALL                             # send email when job status change (start, end, abortion and etc.)
#SBATCH --mail-user bing.song@utsouthwestern.edu          # specify an email address #SBATCH --gres=gpu:1
###########JOBSTART########################
source ~/bashrc_condaset.txt
module load python/3.8.x-anaconda
conda activate /project/DPDS/Wang_lab/shared/BCR_antigen/code/DeLAnO
#export CUDA_VISIBLE_DEVICES=0
cd /project/DPDS/Wang_lab/shared/BCR_antigen/code/DeLAnO
python -u runEmbed.py --use_cpu cpu --cpu 48 --mem 256 --fasta /project/DPDS/Wang_lab/shared/BCR_antigen/code/CLAnO/data/example/cleaned.antigan.Jun.fasta --path.out  /project/DPDS/Wang_lab/shared/BCR_antigen/code/CLAnO/data/example/output --verbose True #&
#cp ~/output/results/pred/*pair.npy intermediates/NPY
#srun --ntasks=1 --nodes=1 python runEmbed.py --use_cpu cpu --fasta /project/DPDS/Wang_lab/shared/BCR_antigen/data/rfscript/Antigen_input/update_0321/combined_93.fasta --path.out /project/DPDS/Wang_lab/shared/BCR_antigen/data/rfscript/Antigen_embed/update0321_combined_9_3 --verbose True &
#wait
