README

## Installation and Prerequisite
Download Cmai.
```sh
cd /Parent/Dir
git clone git@github.com:ice4prince/Cmai.git
```
In the folder scripts/rfold, Install RoseTTAFold from [RoseTTAFold's git](https://github.com/RosettaCommons/RoseTTAFold).

```
cd Cmai/scripts/rfold
git clone git@github.com:RosettaCommons/RoseTTAFold.git
cd RoseTTAFold
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

**NOTE**: the paths for the sequence and structure databases (**UniRef30_2020_06, bfd, and pdb100_2021Mar03**) should be recored and saved for downstream analysis.

Go back to Cmai and Install the environments:

```
cd /Parent/Dir/Cmai
conda env create -f models/runEmbed.yml
conda env create -f models/runBind.yml
```

If a path is specified to install the conda environments:

```
conda env create -f models/runEmbed.yml -p /path/to/cmai_envs/runEmbed
conda env create -f models/runBind.yml -p /path/to/cmai_envs/runBind
```

please remember to add the path to the conda env_dirs

``` 
conda config --append envs_dirs /path/to/cmai_envs
```

and make sure that all three env names are listed in your conda environment list:

```
conda env list
# conda environments:
#
base                  *  /path/to/python/3.8.x-anaconda
runBind                  /path/to/cmai_envs/runBind
runEmbed                 /path/to/cmai_envs/runEmbed
RoseTTAFold              /path/to/.conda/envs/RoseTTAFold
```

*(RoseTTAFold can also been installed in the /path/to/cmai_envs if you like)*

Then run the script provided to save the environment path to paras/env_path.

```bash
./get_env_path.sh
```

## Test with the example data

```
python Cmai.py --input 'data/example/input.csv' --out '/path/to/Cmai/example/data/example/output' --rf_data 'path/to/RoseTTAFold_database'
```

## Data Preparation and Examples

**Required** columns:  

| Antigen_id | Antigen_seq | BCR_Vh | BCR_CDR3h |
| ---------- | ----------- | ------ | --------- |

**Note:**  

If there is NO 'Antigen_seq' column, a fasta file **MUST** be provided.  

**Optional** columns:  

| BCR_id | Score | test | ...  |
| ------ | ----- | ---- | ---- |

NOTE: Rows with **NA** will be **REMOVED**!!!!

**Test data**

| Antigen_id | Antigen_seq                                                  | BCR_id | BCR_Vh                                                       | BCR_CDR3h            | BCR_species |
| ---------- | ------------------------------------------------------------ | ------ | ------------------------------------------------------------ | -------------------- | ----------- |
| H4Hubei    | mlsivilfllvaenssqnytgnpvicmghhavangtmvktltddqvevvaaqelvesqnlpelcpsplrlvdgqtcdiingalgspgcdhlngaewdvfierpnamdtcypfdvpdyqslrsilasngkfefiaeefqwttvkqdgksgackranvndffnrlnwlvksdgnayplqnltkvnngdyarlyiwgvhhpstdteqtnlyknnpggvtvstktsqtsvvpniggrpwvrgqsgrisfywtivepgdlivfntignliaprghyklnnqkkstilntaipigscvskchtdkgslsttkpfqnisriaigncpkyvkqgslklatgmrnipekasrglfgaiagfiengwqglidgwygfrhqnaegtgtaadlkstqaaidqingklnrliektnekyhqiekefeqvegriqdlekyvedtkidlwsynaellvalenqhtidvtdsemnklfervrrqlrenaedkgngcfeifhkcdnnciesirngtydhdiyrdeainnrfqiqgvkltqgykdtilwisfsiscfllvalllafvlwacqngnircqici | BCR_6  | QVQLQQSGPGLVKPSQTLSLTCAISGDSVSSYNAVWNWIRQSPSRGLEWLGRTYYRSGWYNDYAESVKSRITINPDTSKNQFSLQLNSVTPEDTAVYY | CARSGHITVFGVNVDAFDMW |             |
| H4Hubei    | mlsivilfllvaenssqnytgnpvicmghhavangtmvktltddqvevvaaqelvesqnlpelcpsplrlvdgqtcdiingalgspgcdhlngaewdvfierpnamdtcypfdvpdyqslrsilasngkfefiaeefqwttvkqdgksgackranvndffnrlnwlvksdgnayplqnltkvnngdyarlyiwgvhhpstdteqtnlyknnpggvtvstktsqtsvvpniggrpwvrgqsgrisfywtivepgdlivfntignliaprghyklnnqkkstilntaipigscvskchtdkgslsttkpfqnisriaigncpkyvkqgslklatgmrnipekasrglfgaiagfiengwqglidgwygfrhqnaegtgtaadlkstqaaidqingklnrliektnekyhqiekefeqvegriqdlekyvedtkidlwsynaellvalenqhtidvtdsemnklfervrrqlrenaedkgngcfeifhkcdnnciesirngtydhdiyrdeainnrfqiqgvkltqgykdtilwisfsiscfllvalllafvlwacqngnircqici | BCR_7  | QVQLQQSGPGLVKPSQTLSLTCAISGFSVSSYNAVWNWIRQSPSRGLEWLGRTYYRSGWYNDYAESVKSRITINPDTSKNQFSLQLNSVTPEDTAVYY | CARSGHITVFGVNVDAFDMW |             |
| ova        | mgsigaasmefcfdvfkelkvhhanenifycpiaimsalamvylgakdstrtqinkvvrfdklpgfgdsieaqcgtsvnvhsslrdilnqitkpndvysfslasrlyaeerypilpeylqcvkelyrgglepinfqtaadqarelinswvesqtngiirnvlqpssvdsqtamvlvnaivfkglwekafkdedtqampfrvteqeskpvqmmyqiglfrvasmasekmkilelpfasgtmsmlvllpdevsgleqlesiinfekltewtssnvmeerkikvylprmkmeekynltsvlmamgitdvfsssanlsgissaeslkisqavhaahaeineagrevvgsaeagvdaasvseefradhpflfcikhiatnavlffgrcvsp | BCR_11 | QVQLVQSGAEVKKPGASVKVSCKASGYTFTDYYMHWVRQAPGQGLEWMGKVNPNKRGTTYNQKFEGRVTMTTDTSTSTAYMELRSLRSDDTAVYY | CARSNALDAW           | human       |
| ova        | mgsigaasmefcfdvfkelkvhhanenifycpiaimsalamvylgakdstrtqinkvvrfdklpgfgdsieaqcgtsvnvhsslrdilnqitkpndvysfslasrlyaeerypilpeylqcvkelyrgglepinfqtaadqarelinswvesqtngiirnvlqpssvdsqtamvlvnaivfkglwekafkdedtqampfrvteqeskpvqmmyqiglfrvasmasekmkilelpfasgtmsmlvllpdevsgleqlesiinfekltewtssnvmeerkikvylprmkmeekynltsvlmamgitdvfsssanlsgissaeslkisqavhaahaeineagrevvgsaeagvdaasvseefradhpflfcikhiatnavlffgrcvsp | BCR_12 | QVQLVQSGAEVKKPGASVKVSCKASGYTFTDYYMHWVRQAPGQGLEWMGRVNPNGRGTTYNQKFEGRVTMTTDTSTSTAYMELRSLRSDDTAVYY | CARSNALDAW           | human       |

**Preprocessing and Intermediate Data**

Intermediate data will be generated and saved in the **preprocessed directory** defined by '**--pre_dir**' which will be the same as **output directory** if not specified.

**Preprocessed  Data** 

| Antigen_id | Antigen_seq                                                  |      |      |      |      | BCR_id | BCR_Vh                                                       | BCR_CDR3h            | BCR_species | record_id |
| ---------- | ------------------------------------------------------------ | ---- | ---- | ---- | ---- | ------ | ------------------------------------------------------------ | :------------------- | ----------- | --------- |
| H4Hubei    | MLSIVILFLLVAENSSQNYTGNPVICMGHHAVANGTMVKTLTDDQVEVVAAQELVESQNLPELCPSPLRLVDGQTCDIINGALGSPGCDHLNGAEWDVFIERPNAMDTCYPFDVPDYQSLRSILASNGKFEFIAEEFQWTTVKQDGKSGACKRANVNDFFNRLNWLVKSDGNAYPLQNLTKVNNGDYARLYIWGVHHPSTDTEQTNLYKNNPGGVTVSTKTSQTSVVPNIGGRPWVRGQSGRISFYWTIVEPGDLIVFNTIGNLIAPRGHYKLNNQKKSTILNTAIPIGSCVSKCHTDKGSLSTTKPFQNISRIAIGNCPKYVKQGSLKLATGMRNIPEKASRGLFGAIAGFIENGWQGLIDGWYGFRHQNAEGTGTAADLKSTQAAIDQINGKLNRLIEKTNEKYHQIEKEFEQVEGRIQDLEKYVEDTKIDLWSYNAELLVALENQHTIDVTDSEMNKLFERVRRQLRENAEDKGNGCFEIFHKCDNNCIESIRNGTYDHDIYRDEAINNRFQIQGVKLTQGYKDTILWISFSISCFLLVALLLAFVLWACQNGNIRCQICI |      |      |      |      | BCR_6  | QVQLQQSGPGLVKPSQTLSLTCAISGDSVSSYNAVWNWIRQSPSRGLEWLGRTYYRSGWYNDYAESVKSRITINPDTSKNQFSLQLNSVTPEDTAVYY | CARSGHITVFGVNVDAFDMW |             | record_0  |
| H4Hubei    | MLSIVILFLLVAENSSQNYTGNPVICMGHHAVANGTMVKTLTDDQVEVVAAQELVESQNLPELCPSPLRLVDGQTCDIINGALGSPGCDHLNGAEWDVFIERPNAMDTCYPFDVPDYQSLRSILASNGKFEFIAEEFQWTTVKQDGKSGACKRANVNDFFNRLNWLVKSDGNAYPLQNLTKVNNGDYARLYIWGVHHPSTDTEQTNLYKNNPGGVTVSTKTSQTSVVPNIGGRPWVRGQSGRISFYWTIVEPGDLIVFNTIGNLIAPRGHYKLNNQKKSTILNTAIPIGSCVSKCHTDKGSLSTTKPFQNISRIAIGNCPKYVKQGSLKLATGMRNIPEKASRGLFGAIAGFIENGWQGLIDGWYGFRHQNAEGTGTAADLKSTQAAIDQINGKLNRLIEKTNEKYHQIEKEFEQVEGRIQDLEKYVEDTKIDLWSYNAELLVALENQHTIDVTDSEMNKLFERVRRQLRENAEDKGNGCFEIFHKCDNNCIESIRNGTYDHDIYRDEAINNRFQIQGVKLTQGYKDTILWISFSISCFLLVALLLAFVLWACQNGNIRCQICI |      |      |      |      | BCR_7  | QVQLQQSGPGLVKPSQTLSLTCAISGFSVSSYNAVWNWIRQSPSRGLEWLGRTYYRSGWYNDYAESVKSRITINPDTSKNQFSLQLNSVTPEDTAVYY | CARSGHITVFGVNVDAFDMW |             | record_1  |
| ova        | MGSIGAASMEFCFDVFKELKVHHANENIFYCPIAIMSALAMVYLGAKDSTRTQINKVVRFDKLPGFGDSIEAQCGTSVNVHSSLRDILNQITKPNDVYSFSLASRLYAEERYPILPEYLQCVKELYRGGLEPINFQTAADQARELINSWVESQTNGIIRNVLQPSSVDSQTAMVLVNAIVFKGLWEKAFKDEDTQAMPFRVTEQESKPVQMMYQIGLFRVASMASEKMKILELPFASGTMSMLVLLPDEVSGLEQLESIINFEKLTEWTSSNVMEERKIKVYLPRMKMEEKYNLTSVLMAMGITDVFSSSANLSGISSAESLKISQAVHAAHAEINEAGREVVGSAEAGVDAASVSEEFRADHPFLFCIKHIATNAVLFFGRCVSP |      |      |      |      | BCR_11 | QVQLVQSGAEVKKPGASVKVSCKASGYTFTDYYMHWVRQAPGQGLEWMGKVNPNKRGTTYNQKFEGRVTMTTDTSTSTAYMELRSLRSDDTAVYY | CARSNALDAW           | human       | record_2  |
| ova        | MGSIGAASMEFCFDVFKELKVHHANENIFYCPIAIMSALAMVYLGAKDSTRTQINKVVRFDKLPGFGDSIEAQCGTSVNVHSSLRDILNQITKPNDVYSFSLASRLYAEERYPILPEYLQCVKELYRGGLEPINFQTAADQARELINSWVESQTNGIIRNVLQPSSVDSQTAMVLVNAIVFKGLWEKAFKDEDTQAMPFRVTEQESKPVQMMYQIGLFRVASMASEKMKILELPFASGTMSMLVLLPDEVSGLEQLESIINFEKLTEWTSSNVMEERKIKVYLPRMKMEEKYNLTSVLMAMGITDVFSSSANLSGISSAESLKISQAVHAAHAEINEAGREVVGSAEAGVDAASVSEEFRADHPFLFCIKHIATNAVLFFGRCVSP |      |      |      |      | BCR_12 | QVQLVQSGAEVKKPGASVKVSCKASGYTFTDYYMHWVRQAPGQGLEWMGRVNPNGRGTTYNQKFEGRVTMTTDTSTSTAYMELRSLRSDDTAVYY | CARSNALDAW           | human       | record_3  |

(processed_input.csv)

The output directory for the test data is '**data/example/output**', where the intermediate file for antigen embeddings should be saved in the folder '**RFoutputs**', which can be found at **[RFoutputsLink](https://qbrc.swmed.edu/labs/wanglab/download/cmai/RFoutputs.zip)**.

If there is no 'antigen embedding directory' specified using **'--npy_dir'**, the antigen embeddings will be extracted and saved in the folder '**NPY**' in the 'preprocessed directory' as **'pred_dir/NPY',**  which is not included in the github but can be found at **[NPYLink](https://qbrc.swmed.edu/labs/wanglab/download/cmai/NPY.zip)**.

**Expected Output**

The final output **WITHOUT** merged with the input data (using **--no_merge**):

| record_id | Antigen | BCR_id_y | Score    | Rank     |
| --------- | ------- | -------- | -------- | -------- |
| record_0  | H4Hubei | BCR_6    | -3.97083 | 0.581942 |
| record_1  | H4Hubei | BCR_7    | -3.97073 | 0.584042 |
| record_2  | ova     | BCR_11   | -4.16114 | 0.025797 |
| record_3  | ova     | BCR_12   | -4.16291 | 0.022598 |

If **MERGED** with the input data (Default):

| Antigen_id | Antigen_seq                                                  | BCR_id_x | BCR_Vh                                                       | BCR_CDR3h            | BCR_species | record_id | Antigen | BCR_id_y | Score    | Rank     |
| ---------- | ------------------------------------------------------------ | -------- | ------------------------------------------------------------ | -------------------- | ----------- | --------- | ------- | -------- | -------- | -------- |
| H4Hubei    | mlsivilfllvaenssqnytgnpvicmghhavangtmvktltddqvevvaaqelvesqnlpelcpsplrlvdgqtcdiingalgspgcdhlngaewdvfierpnamdtcypfdvpdyqslrsilasngkfefiaeefqwttvkqdgksgackranvndffnrlnwlvksdgnayplqnltkvnngdyarlyiwgvhhpstdteqtnlyknnpggvtvstktsqtsvvpniggrpwvrgqsgrisfywtivepgdlivfntignliaprghyklnnqkkstilntaipigscvskchtdkgslsttkpfqnisriaigncpkyvkqgslklatgmrnipekasrglfgaiagfiengwqglidgwygfrhqnaegtgtaadlkstqaaidqingklnrliektnekyhqiekefeqvegriqdlekyvedtkidlwsynaellvalenqhtidvtdsemnklfervrrqlrenaedkgngcfeifhkcdnnciesirngtydhdiyrdeainnrfqiqgvkltqgykdtilwisfsiscfllvalllafvlwacqngnircqici | BCR_6    | QVQLQQSGPGLVKPSQTLSLTCAISGDSVSSYNAVWNWIRQSPSRGLEWLGRTYYRSGWYNDYAESVKSRITINPDTSKNQFSLQLNSVTPEDTAVYY | CARSGHITVFGVNVDAFDMW |             | record_0  | H4Hubei | BCR_6    | -3.97083 | 0.581942 |
| H4Hubei    | mlsivilfllvaenssqnytgnpvicmghhavangtmvktltddqvevvaaqelvesqnlpelcpsplrlvdgqtcdiingalgspgcdhlngaewdvfierpnamdtcypfdvpdyqslrsilasngkfefiaeefqwttvkqdgksgackranvndffnrlnwlvksdgnayplqnltkvnngdyarlyiwgvhhpstdteqtnlyknnpggvtvstktsqtsvvpniggrpwvrgqsgrisfywtivepgdlivfntignliaprghyklnnqkkstilntaipigscvskchtdkgslsttkpfqnisriaigncpkyvkqgslklatgmrnipekasrglfgaiagfiengwqglidgwygfrhqnaegtgtaadlkstqaaidqingklnrliektnekyhqiekefeqvegriqdlekyvedtkidlwsynaellvalenqhtidvtdsemnklfervrrqlrenaedkgngcfeifhkcdnnciesirngtydhdiyrdeainnrfqiqgvkltqgykdtilwisfsiscfllvalllafvlwacqngnircqici | BCR_7    | QVQLQQSGPGLVKPSQTLSLTCAISGFSVSSYNAVWNWIRQSPSRGLEWLGRTYYRSGWYNDYAESVKSRITINPDTSKNQFSLQLNSVTPEDTAVYY | CARSGHITVFGVNVDAFDMW |             | record_1  | H4Hubei | BCR_7    | -3.97073 | 0.584042 |
| ova        | mgsigaasmefcfdvfkelkvhhanenifycpiaimsalamvylgakdstrtqinkvvrfdklpgfgdsieaqcgtsvnvhsslrdilnqitkpndvysfslasrlyaeerypilpeylqcvkelyrgglepinfqtaadqarelinswvesqtngiirnvlqpssvdsqtamvlvnaivfkglwekafkdedtqampfrvteqeskpvqmmyqiglfrvasmasekmkilelpfasgtmsmlvllpdevsgleqlesiinfekltewtssnvmeerkikvylprmkmeekynltsvlmamgitdvfsssanlsgissaeslkisqavhaahaeineagrevvgsaeagvdaasvseefradhpflfcikhiatnavlffgrcvsp | BCR_11   | QVQLVQSGAEVKKPGASVKVSCKASGYTFTDYYMHWVRQAPGQGLEWMGKVNPNKRGTTYNQKFEGRVTMTTDTSTSTAYMELRSLRSDDTAVYY | CARSNALDAW           | human       | record_2  | ova     | BCR_11   | -4.16114 | 0.025797 |
| ova        | mgsigaasmefcfdvfkelkvhhanenifycpiaimsalamvylgakdstrtqinkvvrfdklpgfgdsieaqcgtsvnvhsslrdilnqitkpndvysfslasrlyaeerypilpeylqcvkelyrgglepinfqtaadqarelinswvesqtngiirnvlqpssvdsqtamvlvnaivfkglwekafkdedtqampfrvteqeskpvqmmyqiglfrvasmasekmkilelpfasgtmsmlvllpdevsgleqlesiinfekltewtssnvmeerkikvylprmkmeekynltsvlmamgitdvfsssanlsgissaeslkisqavhaahaeineagrevvgsaeagvdaasvseefradhpflfcikhiatnavlffgrcvsp | BCR_12   | QVQLVQSGAEVKKPGASVKVSCKASGYTFTDYYMHWVRQAPGQGLEWMGRVNPNGRGTTYNQKFEGRVTMTTDTSTSTAYMELRSLRSDDTAVYY | CARSNALDAW           | human       | record_3  | ova     | BCR_12   | -4.16291 | 0.022598 |

(merged_results.csv)

**Expected Time**

The expected runime for processing the test data is around 55m55.012s.



## Pipeline

```sh
python Cmai.py --input 'data/example/input.csv' --out '/path/to/Cmai/example/data/example/output' --rf_data 'path/to/RoseTTAFold_database'
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

Cmai.py is the main interface for users to execute Cmai after installation is successful. 

```shi
usage: Cmai.py [-h] [--code CODE] [--input INPUT] [--out OUT] [--env_path ENV_PATH] [--rf_data RF_DATA]
               [--fasta FASTA] [--pre_dir PRE_DIR] [--npy_dir NPY_DIR] [--cpu CPU] [--mem MEM] [--use_cpu]
               [--seed SEED] [--min_size_background_bcr MIN_SIZE_BACKGROUND_BCR]
               [--max_size_background_bcr MAX_SIZE_BACKGROUND_BCR]
               [--min_size_background_antigen MIN_SIZE_BACKGROUND_ANTIGEN]
               [--max_size_background_antigen MAX_SIZE_BACKGROUND_ANTIGEN] [--export_background] [--add_rank]
               [--background_score BACKGROUND_SCORE] [--e_values E_VALUES] [--background_npy BACKGROUND_NPY]
               [--rf_para] [--gen_msa] [--run_rf] [--skip_preprocess] [--skip_extract] [--runEmbed] [--runBind]
               [--skip_check] [--suffix] [--no_rank] [--verbose] [--no_merge] [--move_npy] [--gen_npy] [--embedBCR]
               [--bcr_heatmap] [--debug] [--backBCR_only] [--backAntigen_only]

Parameters for the interface script.

optional arguments:
  -h, --help            show this help message and exit
  --code CODE           the Cmai directory
  --input INPUT         the input files in csv which should include Antigen_id,BCR_Vh,BCR_CDR3h
  --out OUT             the directory for output files. An absolute path is required.
  --env_path ENV_PATH   the file saving the directory of the Conda environments- python of runEmbed, python of
                        runBind, and RoseTTAFold in order.
  --rf_data RF_DATA     the database folder for RoseTTAFold
  --fasta FASTA         The fasta file entering runEbed. When no sequence included in the input, the separate fasta
                        file of antigens is required
  --pre_dir PRE_DIR     the directory to save the preprocessed data. If not defiend, same with output directory.
  --npy_dir NPY_DIR     the npy folder if different with "preprocess folder/NPY".
  --cpu CPU             the maximum of cpus for antigen embedding. If not defined, use the value of paras/rf_para.txt
  --mem MEM             the maximum of memory in GB for antigen embedding. If not defined, use the value of
                        paras/rf_para.txt
  --use_cpu             the option to use cpu or gpu.
  --seed SEED           the seed for the first 100 background BCRs. To use the prepared embeded 100 BCRs, keep the
                        seed to default 1
  --min_size_background_bcr MIN_SIZE_BACKGROUND_BCR
                        the initial and minimum sample size of background BCRs. The default is 100
  --max_size_background_bcr MAX_SIZE_BACKGROUND_BCR
                        the maximum size for subsample of background BCRs, which should no more than 1000000. The
                        default is 10000
  --min_size_background_antigen MIN_SIZE_BACKGROUND_ANTIGEN
                        the initial sample size of background Antigens. The default is 10
  --max_size_background_antigen MAX_SIZE_BACKGROUND_ANTIGEN
                        the maximum size of subsample of background antigens, which should no more than 100. The
                        deafult is 100
  --export_background   Only export the score dict for background BCRs of quantity of the max_size_background_bcr
                        number, default is False.
  --add_rank            Only add ranks from background BCR scores to no_ranked results, default is False.
  --background_score BACKGROUND_SCORE
                        the pkl file of the score dictionary of background BCRs
  --e_values E_VALUES   E-value cutoff for inclusion in result alignment. Default is '1e-30 1e-10 1e-6 1e-3'
  --background_npy BACKGROUND_NPY
                        the npy folder of background antigens if different with “preprocess folder/NPY”.
  --rf_para             use the parameters from paras/rf_para.txt for antigen embedding. Default is False
  --gen_msa             only run generating msa and exit. Default is False
  --run_rf              skip generating msa and running embedding prediction. Default is False
  --skip_preprocess     skip preprocess of antigen_embedding. Default is False
  --skip_extract        skip extracting NPY for antigen embedding. Default is False
  --runEmbed            only run antigen embedding. Default is False
  --runBind             only run binding. Default is False
  --skip_check          skip check and preprocess of input data, only use when it has been done before. Default is
                        False
  --suffix              Adding suffix to antigen id. Only use to distinguish same-name antigens. The default is False.
  --no_rank             Only export the predicted score but no rank in background BCRs, default is False.
  --verbose             Enable verbose output, default is False.
  --no_merge            Unable merging output to input, default is False.
  --move_npy            only move npy files to the desired directory. Default is False
  --gen_npy             extract npy from npz files. Default is False
  --embedBCR            extract the bcr sequences and embeddings to the folder of preprocessed data. Default is False
  --bcr_heatmap         export full embedding results including the heatmap comparison. Default is False
  --debug               Switch to the debug mode and print output step by step. Default is False
  --backBCR_only        Only get the rank in background BCRs. Default is False
  --backAntigen_only    Only get the rank in background antigens. Default is False

```

## Runtime and Memory Usage

Cmai is designed for large-scale inference on the binding properties between antigens and antibodies. Speed is thus an important factor in its scalability. We randomly sampled 50 BCR-antigen binding pairs (50 unique antigens) from SabDab and ran our model. In the Antigen-anchor mode, the model took an average of ~12 minutes to finish computing for each pair, which included both the antigen embedding computation step by RosettaFold (average time=*11.84 minutes*), and also the binding prediction/rank percentile computation step against 1 million BCRs (average time=*32.67 seconds*). In practice, the most common scenario is to test the binding between many BCRs against a specific antigen of interest, and in this case, the embedding step needs to be run **only once**, enabling efficient inference on massive datasets. As a contrast, the model took an average of 42.06 seconds on the the binding prediction/rank percentile computation step against 10 to 100 background antigens in the BCR-anchor mode. Both mode share the same step of antigen embedding, so it should take less than 13 minutes on average to run both the embedding and binding prediction among background antigens. Since the embedding step runs only once, two types of ranks can be both predicted for one BCR-antigen pair within 14 minutes on average.

In terms of computational hardware requirements, GPUs are required for the binding prediction phase (executed with `--runBind`), and are strongly recommended for the antigen embedding phase (`--run_rf`). For optimal performance, RoseTTAFold generally requires a GPU with at least 40 GB of memory to prevent out-of-memory errors. On the other hand, the MSA generation phase (`--gen_msa`) runs on CPUs. If GPU resources are limited, it is advisable to run the MSA generation phase separately on CPUs before proceeding with the other phases.

Our training and prediction computations were executed on the A100 GPU nodes of our UT Southwestern BioHPC server (https://portal.biohpc.swmed.edu/content/). 

## Round-2 Training ##

In version 1.0, we only trained Cmai with the antigen as the anchor and contrasted between positive and negative binding BCRs. When Cmai 1.0 generates predictions, it also always uses the query antigen as the anchor, and tests different BCRs against the antigen (this is how we calculated rank%). In version 1.1.0, we updated Cmai to include with the other direction of contrast, where the BCR serves as the anchor, and the antigens will be compared for their binding towards the query BCR. Again, the training, validation and application of Cmai in this direction all follow this pattern. This updated version incorporates the original trained model (antigen as anchor) alongside a newly trained model (BCR as anchor). By default, predictions are performed for both directions, and the merged outputs include rank% (antigen anchor), rank% (BCR anchor), maximum/average rank% of both directions. Additionally, options are available for just generating predictions for each direction: **--backBCR_only** for antigen-anchored data and **--backAntigen_only** for BCR-anchored data.

Validation of Cmai version 1.1 was conducted, with the results visualized in the ROC curves presented below.

![Antigen-anchored data](Z:\shared\BCR_antigen\code\Cmai\archive\R2_Antigen_ROC_Curve.png)

 

![BCR-anchored data](Z:\shared\BCR_antigen\code\Cmai\archive\R2_BCR_ROC_Curve.png)

Strikingly and very interestingly, we found that the branch of Cmai that is trained with the antigen as the anchor would also perform better for the validation data where the antigens serve as the anchor (query an antigen’ s binding towards a set of BCRs), and that the branch of Cmai that is trained with the BCR as the anchor would also perform better for the validation data where the BCRs serve as the anchor (query a BCR’ s binding towards a set of antigens). This result makes sense as the behavior of Cmai in application would be influenced by how it is trained. And it also gives us a hint of how Cmai (and its future improved versions) should be applied. Namely, the antigen-anchor Cmai will be more helpful for antibody design and optimization, where we have one query antigen and need to find a good/better binding BCR. On the other hand, the BCR-anchor Cmai would be more helpful for antigen discovery and vaccine development, where the goal is to identify the binding target of an antibody in question. Additionally, the average or maximum of rank%s of both directions should probably be employed in situations where the objective is to identify the binding antigen-BCR pairs from a set of candidate BCRs and a set of candidate antigens.











