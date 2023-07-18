README

#  Installation and Prerequisite
```sh
cd /path/to/Cmai
git clone git@github.com:ice4prince/Cmai.git
conda env create -f models/runEmbed.yml
conda activate runEmbed
which python >paras/runEmbed_path #save your path of runEmbed env
conda deactivate
conda env create -f models/runBC.yml
cd scripts/rfscripts
rm -r RoseTTAFold
#  Install RoseTTAFold from RoseTTAFold's git.
git clone git@github.com:RosettaCommons/RoseTTAFold.git
# Follow the README in RoseTTAFold OR go through
	# create the environment
	# download network weights 
	# download and install third-party softwareand
	# download sequence and structure databases
```


## Pipeline
```sh
conda activate runBC
python Cmai.py --code '/path/to/Cmai' --input 'data/example/binary_example.csv' --out 'data/example/output' --rf_data 'path/to/RoseTTAFold_database'
```
## Step-by-step Pipeline

```sh
conda activate runBC

python Cmai.py --code '/path/to/Cmai' --input 'data/example/binary_example.csv' --out 'data/example/output' --rf_data 'path/to/RoseTTAFold_database'  --runEmbed --gen_msa --use_cpu 'cpu'

python Cmai.py --code '/path/to/Cmai' --input 'data/example/binary_example.csv' --out 'data/example/output' --rf_data 'path/to/RoseTTAFold_database'  --runEmbed --run_rf --use_cpu 'gpu'

python Cmai.py --code '/path/to/Cmai' --out 'data/example/output' --skip_check --runBC
```

## Usage

```sh
usage: Cmai.py [-h] [--code CODE] [--input INPUT] [--out OUT] [--env_path ENV_PATH] [--rf_data RF_DATA]
               [--fasta FASTA] [--pre_dir PRE_DIR] [--seed SEED] [--subsample SUBSAMPLE] [--bottomline BOTTOMLINE]
               [--continuous] [--rf_para] [--gen_msa] [--run_rf] [--skip_preprocess] [--skip_extract] [--runEmbed]
               [--runBC] [--species] [--verbose]
```

Parameters for the interface script.

```sh
optional arguments:
  -h, --help            show this help message and exit
  --code CODE           the Cmai directory
  --input INPUT         the input files in csv which should include Antigen_id,BCR_Vh,BCR_CDR3h
  --out OUT             the directory for output files
  --env_path ENV_PATH   the file saving the directory of the runEmbed Conda environment
  --rf_data RF_DATA     the database folder for RoseTTAFold
  --fasta FASTA         The fasta file entering runEbed. When no sequence included in the input, the separate fasta
                        file of antigens is required
  --pre_dir PRE_DIR     the directory to save the preprocessed data.
  --seed SEED           the seed for the first 100 background BCRs. To use the prepared embeded 100 BCRs, keep the
                        seed to default 1
  --subsample SUBSAMPLE
                        the initial sample size of background BCRs. The default is 100
  --bottomline BOTTOMLINE
                        the maximum size for subsample of background BCRs, which should no more than 1000000. The
                        default is 10000
  --continuous          swtich the mode from binary to continuous, default mode is binary.
  --rf_para             use the parameters from paras/rf_para.txt for antigen embedding. Default is False
  --gen_msa             only run generating msa and exit. Default is False
  --run_rf              skip generating msa and running RoseTTAFold. Default is False
  --skip_preprocess     skip preprocess of antigen_embedding. Default is False
  --skip_extract        skip extracting NPY for antigen embedding. Default is False
  --runEmbed            only run antigen embedding. Default is False
  --runBC               only run binding or comparing. Default is False
  --species             match the species of background BCR to the target BCR. NOTE: the species MUST BE specified and
                        unique in the target BCR input.
  --verbose             Enable verbose output, default is False.
```
