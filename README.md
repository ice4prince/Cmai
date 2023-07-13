## README
#  install multiple environment
	cd /path/to/Cmai
	git clone git@github.com:ice4prince/Cmai.git
	conda env create -f models/runEmbed.yml
	conda env create -f models/runBC.yml
	cd scripts/rfscripts
	rm -r RoseTTAFold
#  Install RoseTTAFold
Clone the package  

	git clone git@github.com:RosettaCommons/RoseTTAFold.git
	# Follow the README in RoseTTAFold OR go through
	cd RoseTTAFold
Create conda environment using RoseTTAFold-linux.yml file and folding-linux.yml file. The latter is required to run a pyrosetta version only (run_pyrosetta_ver.sh)
	
	# create conda environment for RoseTTAFold
	conda env create -f RoseTTAFold-linux.yml  #If your NVIDIA driver compatible with cuda11
	conda env create -f RoseTTAFold-linux-cu101.yml #If not (but compatible with cuda10)
 Download network weights
	
	wget https://files.ipd.uw.edu/pub/RoseTTAFold/weights.tar.gz
	tar xfz weights.tar.gz
Download and install third-party software.
	
	./install_dependencies.sh
Download sequence and structure databases  
	
	# uniref30 [46G]
	wget http://wwwuser.gwdg.de/~compbiol/uniclust/2020_06/UniRef30_2020_06_hhsuite.tar.gz
	mkdir -p UniRef30_2020_06
	tar xfz UniRef30_2020_06_hhsuite.tar.gz -C ./UniRef30_2020_06
	# BFD [272G]
	wget https://bfd.mmseqs.com/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt.tar.gz
	mkdir -p bfd
	tar xfz
	bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt.tar.gz -C ./bfd
	# structure templates (including *_a3m.ffdata, *_a3m.ffindex) [over 100G]
	wget https://files.ipd.uw.edu/pub/RoseTTAFold/pdb100_2021Mar03.tar.gz
	tar xfz pdb100_2021Mar03.tar.gz
	# for CASP14 benchmarks, we used this one:https://files.ipd.uw.edu/pub/RoseTTAFold/pdb100_2020Mar11.tar.gz`


## Pipeline
	conda activate runBC
	python Cmai.py --code '/path/to/Cmai' --input 'data/example/binary_example.csv' --out 'data/example/output' --rf_data 'path/to/RoseTTAFold_database' 
## Usage
	usage: Cmai.py [-h] [--code CODE] [--input INPUT] [--out OUT] [--rf_data RF_DATA]
               [--fasta FASTA] [--pre_dir PRE_DIR] [--seed SEED]
               [--subsample SUBSAMPLE] [--bottomline BOTTOMLINE] [--continuous]
               [--rf_para] [--gen_msa] [--run_rf] [--skip_preprocess]
               [--skip_extract] [--runEmbed] [--runBC] [--species] [--verbose]

Parameters for the interface script.

optional arguments:
  -h, --help            show this help message and exit
  --code CODE           the CLAnO directory
  --input INPUT         the input files in csv which should include
                        Antigen_id,BCR_Vh,BCR_CDR3h
  --out OUT             the directory for output files
  --rf_data RF_DATA     the database folder for RoseTTAFold
  --fasta FASTA         The fasta file entering runEbed. When no sequence included
                        in the input, the seperate fasta file of antigens is
                        required
  --pre_dir PRE_DIR     the directory to save the preprocessed data.
  --seed SEED           the seed for the first 100 background BCRs. To use the
                        prepared embeded 100 BCRs, keep the seed to default 1
  --subsample SUBSAMPLE
                        the initial sample size of background BCRs. The default is
                        100
  --bottomline BOTTOMLINE
                        the maximum size for subsample of background BCRs, which
                        should no more than 1000000. The deafult is 10000
  --continuous          swtich the mode from binary to continuous, default mode is
                        binary.
  --rf_para             use the parameters from paras/rf_para.txt for antigen
                        embedding. Default is False
  --gen_msa             only run generating msa and exit. Default is False
  --run_rf              skip generating msa and running RoseTTAFold. Default is
                        False
  --skip_preprocess     skip preprocess of antigen_embedding. Default is False
  --skip_extract        skip extracting NPY for antigen embedding. Default is
                        False
  --runEmbed            only run antigen embedding. Default is False
  --runBC               only run binding or comparing. Default is False
  --species             match the species of background BCR to the target BCR.
                        NOTE: the species MUST BE specified and unique in the
                        target BCR input.
  --verbose             Enable verbose output, default is False.
