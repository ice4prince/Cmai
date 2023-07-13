#!/bin/bash

# make the script stop when error (non-true exit code) is occured
set -e

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('conda' 'shell.bash' 'hook' 2> /dev/null)"
eval "$__conda_setup"
unset __conda_setup
# <<< conda initialize <<<

conda create -n py310 python=3.10
conda activate py310
pip install rich smile_config biopython -U
conda deactivate


cd RoseTTAFold
conda env create -f RoseTTAFold-linux.yml

./install_dependencies.sh
wget https://files.ipd.uw.edu/pub/RoseTTAFold/weights.tar.gz
tar xfz weights.tar.gz

cd ..

conda activate py310
