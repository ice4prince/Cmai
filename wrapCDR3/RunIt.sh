#!/bin/bash
#SBATCH --job-name=CDR3_wrap
#SBATCH -p GPUv100s
#SBATCH --nodes=1                                                    # number of nodes requested by user
#SBATCH --ntasks=1                                                  # number of total tasks
#SBATCH --time=10-00:00:00                                            # run time, format: D-H:M:S (max wallclock time)
#SBATCH --output=./sbatch_output_%j                                 # redirect both standard output and erro output to the same file
#SBATCH --error=./sbatch_error_%j                                 # standard error output file name
#SBATCH --mail-user= yourUsername@utsouthwestern.edu                     # specify an email address
#SBATCH --mail-type=ALL                             # send email when job status change (start, end, abortion and etc.)
#SBATCH --gres=gpu:1
###########JOBSTART########################
source ~/.bashrc_conda
module load python/3.8.x-anaconda
# conda create --name myenv --file spec-file.txt # for the FRIST TIME, create your environment using the spec-file.txt
conda activate myenv #Replace myenv with your specifed environment name
# conda env create -f torchEnv.yml
# conda activate torchEnv.yml #OR create your environment from the yml file and activate it.
export CUDA_VISIBLE_DEVICES=0
python3.8 ~/CDR3_wrapup.py \
	-i ~/inputfile \
	# input the sequence or gene id file. Usually in txt.
	-o ~/outputFolder \
	# set the directory of outputs
	> ~/log.txt
