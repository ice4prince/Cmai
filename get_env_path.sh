CONDA_BASE=$(conda info | grep -i 'base environment' | awk '{print $4}')

# Source the conda.sh script using the value of CONDA_BASE
source "$CONDA_BASE/etc/profile.d/conda.sh"
# source /cm/shared/apps/python/3.8.x-anaconda/etc/profile.d/conda.sh
conda activate runEmbed
which python >paras/env_path
conda activate runBind
which python >>paras/env_path
conda activate RoseTTAFold
echo $CONDA_PREFIX >>paras/env_path
sed -i "s|/home2/s205236/.conda/envs/RoseTTAFold|$(echo $CONDA_PREFIX)|g" paras/rf_para.txt
conda deactivate
conda deactivate
conda deactivate
