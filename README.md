README

## Installation and Prerequisite
Download Cmai and install the environments.
```sh
cd /path/to/Cmai
git clone git@github.com:ice4prince/Cmai.git
conda env create -f models/runEmbed.yml
conda env create -f models/runBind.yml
```
In the folder scripts/rfold, Install RoseTTAFold from [RoseTTAFold's git](https://github.com/RosettaCommons/RoseTTAFold).
```
cd scripts/rfold
git clone git@github.com:RosettaCommons/RoseTTAFold.git
```
Follow the README in RoseTTAFold OR go through  

	1. create the environment
	
	# If your NVIDIA driver compatible with cuda11
	conda env create -f RoseTTAFold-linux.yml
	# If not (but compatible with cuda10)
	conda env create -f RoseTTAFold-linux-cu101.yml
	# create conda environment for pyRosetta folding & running DeepAccNet
	conda env create -f folding-linux.yml
	
	2. download network weights
	
	wget https://files.ipd.uw.edu/pub/RoseTTAFold/weights.tar.gz
	tar xfz weights.tar.gz
	
	3. download and install third-party software
	
	./install_dependencies.sh
	 
	4. download sequence and structure databases
	
	# uniref30 [46G]
	wget http://wwwuser.gwdg.de/~compbiol/uniclust/2020_06/UniRef30_2020_06_hhsuite.tar.gz
	mkdir -p UniRef30_2020_06
	tar xfz UniRef30_2020_06_hhsuite.tar.gz -C ./UniRef30_2020_06
	
	# BFD [272G]
	wget https://bfd.mmseqs.com/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt.tar.gz
	mkdir -p bfd
	tar xfz bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt.tar.gz -C ./bfd
	
	# structure templates (including *_a3m.ffdata, *_a3m.ffindex) [over 100G]
	wget https://files.ipd.uw.edu/pub/RoseTTAFold/pdb100_2021Mar03.tar.gz
	tar xfz pdb100_2021Mar03.tar.gz
	# for CASP14 benchmarks, we used this one: https://files.ipd.uw.edu/pub/RoseTTAFold/pdb100_2020Mar11.tar.gz


Save the environment path to paras/env_path IN ORDER
```
./get_env_path.sh
```
## Data Preparation

Required columns:  

	| Antigen_id | Antigen_seq | BCR_Vh | BCR_CDR3h |  
	**To be noticed:**  
	If there is no 'Antigen_seq' column, a fasta file MUST be provided.  
Optional columns:  
| BCR_id | Score | test | ... |  

## Pipeline

```sh
python Cmai.py --code '/path/to/Cmai' --input 'data/example/input.csv' --out '/path/to/Cmai/example/data/example/output' --rf_data 'path/to/RoseTTAFold_database'
```
## Step-by-step Pipeline

```sh
# Antigen Embedding:
# In 1-step:
python Cmai.py --code '/path/to/Cmai' --input 'data/example/binary_example.csv' --out 'data/example/output' --rf_data 'path/to/RoseTTAFold_database'  --runEmbed
# In 2 steps:
	python Cmai.py --code '/path/to/Cmai' --input 'data/example/binary_example.csv' --out 'data/example/output' --rf_data 'path/to/RoseTTAFold_database'  --runEmbed --gen_msa --use_cpu
	python Cmai.py --code '/path/to/Cmai' --input 'data/example/binary_example.csv' --out 'data/example/output' --rf_data 'path/to/RoseTTAFold_database'  --runEmbed --run_rf

# Binding Predict:
python Cmai.py --code '/path/to/Cmai' --out 'data/example/output' --skip_check --runBind
```

## Usage

```sh
usage: Cmai.py [-h] [--code CODE] [--input INPUT] [--out OUT] [--env_path ENV_PATH]
               [--rf_data RF_DATA] [--fasta FASTA] [--pre_dir PRE_DIR]
               [--npy_dir NPY_DIR] [--cpu CPU] [--mem MEM] [--use_cpu] [--seed SEED]
               [--min_size_background_bcr MIN_SIZE_BACKGROUND_BCR]
               [--max_size_background_bcr MAX_SIZE_BACKGROUND_BCR] [--rf_para]
               [--gen_msa] [--run_rf] [--skip_preprocess] [--skip_extract]
               [--runEmbed] [--runBind] [--skip_check] [--suffix] [--no_rank]
               [--verbose] [--merge] [--move_npy] [--embedBCR] [--bcr_heatmap]
               [--debug]
```

Parameters for the interface script.

```sh
optional arguments:
  -h, --help            show this help message and exit
  --code CODE           the Cmai directory
  --input INPUT         the input files in csv which should include
                        Antigen_id,BCR_Vh,BCR_CDR3h
  --out OUT             the directory for output files. An absolute path is required.

  --env_path ENV_PATH   the file saving the directory of the Conda environments-
                        python of runEmbed, python of runBind, and RoseTTAFold in
                        order.
  --rf_data RF_DATA     the database folder for RoseTTAFold
  --fasta FASTA         The fasta file entering runEbed. When no sequence included in
                        the input, the separate fasta file of antigens is required
  --pre_dir PRE_DIR     the directory to save the preprocessed data. If not defiend,
                        same with output directory.
  --npy_dir NPY_DIR     the npy folder if different with preprocess folder
  --cpu CPU             the maximum of cpus for antigen embedding. If not defined, use
                        the value of paras/rf_para.txt
  --mem MEM             the maximum of memory in GB for antigen embedding. If not
                        defined, use the value of paras/rf_para.txt
  --use_cpu             the option to use cpu or gpu.
  --seed SEED           the seed for the first 100 background BCRs. To use the
                        prepared embeded 100 BCRs, keep the seed to default 1
  --min_size_background_bcr MIN_SIZE_BACKGROUND_BCR
                        the initial and minimum sample size of background BCRs. The
                        default is 100

  --max_size_background_bcr MAX_SIZE_BACKGROUND_BCR
                        the maximum size for subsample of background BCRs, which
                        should no more than 1000000. The default is 10000
  --rf_para             use the parameters from paras/rf_para.txt for antigen
                        embedding. Default is False
  --gen_msa             only run generating msa and exit. Default is False
  --run_rf              skip generating msa and running embedding prediction. Default
                        is False
  --skip_preprocess     skip preprocess of antigen_embedding. Default is False
  --skip_extract        skip extracting NPY for antigen embedding. Default is False
  --runEmbed            only run antigen embedding. Default is False
  --runBind             only run binding. Default is False
  --skip_check          skip check and preprocess of input data, only use when it has
                        been done before. Default is False
  --suffix              Adding suffix to antigen id. Only use to distinguish same-name
                        antigens. The default is False.
  --no_rank             Only export the predicted score but no rank in background
                        BCRs, default is False.
  --verbose             Enable verbose output, default is False.
  --merge               Enable merging output to input, default is False.
  --move_npy            only move npy files to the desired directory. Default is False
  --embedBCR            extract the bcr sequences and embeddings to the folder of
                        preprocessed data. Default is False
  --bcr_heatmap         export full embedding results including the heatmap
                        comparison. Default is False
  --debug               Switch to the debug mode and print output step by step.
                        Default is False

```
