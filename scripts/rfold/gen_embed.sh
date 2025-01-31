#!/bin/bash

# make the script stop when error (non-true exit code) is occured
set -e

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('conda' 'shell.bash' 'hook' 2> /dev/null)"
eval "$__conda_setup"
unset __conda_setup
# <<< conda initialize <<<

conda activate RoseTTAFold

RFD=$1

DB=$2/pdb100_2021Mar03/pdb100_2021Mar03

IN="$3/data"
OUTD="$3/pred"

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

mkdir -p $IN
mkdir -p $OUTD

if [ ! -f $IN/IDX.uniq ]; then
  echo "IDX.uniq does not exist, sort and uniq IDX..."
  sort $IN/IDX | uniq | sort > $IN/IDX.uniq
else
  echo "IDX.uniq exists, skipping..."
fi


#python "$SCRIPT_DIR"/$RFD/network/predict_smile.py $RFD/weights $DB $IN $OUTD $4
cp $SCRIPT_DIR/model.py $RFD/network/RoseTTAFoldModel.py
cp $SCRIPT_DIR/predict.py $RFD/network/predict_smile.py
python $RFD/network/predict_smile.py $RFD/weights $DB $IN $OUTD $4
