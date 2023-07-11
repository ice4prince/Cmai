## README
#  install multiple environment
	cd /path/to/Cmai
	git clone git@github.com:ice4prince/Cmai.git
	conda env create -f models/runEmbed.yml
	conda env create -f models/torch_a100.yml
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
	python Cmai.py ....
## Usage
	usage: runEmbed.py [-h] [--fasta FASTA] [--verbose VERBOSE] [--cpu CPU] [--mem MEM] [--use_cpu USE_CPU]
			   [--skip_gen_msa SKIP_GEN_MSA] [--skip_run_rf SKIP_RUN_RF] [--skip_preprocess SKIP_PREPROCESS]
			   [--skip_extract SKIP_EXTRACT] [--path.data PATH.DATA] [--path.scripts PATH.SCRIPTS]
			   [--path.rf PATH.RF] [--path.out PATH.OUT] [--path.logs PATH.LOGS] [--exe.hhsearch EXE.HHSEARCH]
			   [--exe.hhblits EXE.HHBLITS] [--exe.hhfilter EXE.HHFILTER] [--exe.psipred EXE.PSIPRED]
			   [--exe.psipass2 EXE.PSIPASS2]

	Config for runEmbed script.

	options:
	  -h, --help            show this help message and exit
	  --fasta FASTA         Input fasta files. (default:
				/project/DPDS/Wang_lab/shared/BCR_antigen/code/DeLAnO/example/input.fa)
	  --verbose VERBOSE     Print verbose messages or not. (default: False)
	  --cpu CPU             Max CPUs. (default: 32)
	  --mem MEM             Max Memory (in GB). (default: 128)
	  --use_cpu USE_CPU     cpu|gpu (default: cpu)
	  --skip_gen_msa SKIP_GEN_MSA
				Skip generate_msa (default: False)
	  --skip_run_rf SKIP_RUN_RF
				Skip run_rf (default: False)
	  --skip_preprocess SKIP_PREPROCESS
				Skip preprocess (default: False)
	  --skip_extract SKIP_EXTRACT
				Skip extraction (default: False)

	path:
	  --path.data PATH.DATA
				Data base folder. (default: /project/DPDS/Wang_lab/shared/BCR_antigen/data)
	  --path.scripts PATH.SCRIPTS
				Scripts folder. (default: /project/DPDS/Wang_lab/shared/BCR_antigen/code/DeLAnO)
	  --path.rf PATH.RF     RosettaFold folder. (default: rfscripts/RoseTTAFold)
	  --path.out PATH.OUT   Output folder. (default:
				/project/DPDS/Wang_lab/shared/BCR_antigen/code/DeLAnO/example/outputs)
	  --path.logs PATH.LOGS
				Logs folder. (default: logs)

	exe:
	  --exe.hhsearch EXE.HHSEARCH
				HHsearch executable. (default: rfscripts/bin/hhsearch)
	  --exe.hhblits EXE.HHBLITS
				HHblits executable. (default: rfscripts/bin/hhblits)
	  --exe.hhfilter EXE.HHFILTER
				HHfilter executable. (default: rfscripts/bin/hhfilter)
	  --exe.psipred EXE.PSIPRED
				PSIPRED executable. (default: psipred)
	  --exe.psipass2 EXE.PSIPASS2
				PSIPASS2 executable. (default: psipass2)
