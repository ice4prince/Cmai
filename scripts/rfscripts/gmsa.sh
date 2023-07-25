#!/bin/bash

# make the script stop when error (non-true exit code) is occured
set -e

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('conda' 'shell.bash' 'hook' 2> /dev/null)"
eval "$__conda_setup"
unset __conda_setup
# <<< conda initialize <<<

export PIPEDIR=$1

CPU=$2
MEM=$3

WDIR=$4/work
DDIR=$4/data

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

LOGD=$WDIR/$5

DB=$6

IN=$7
OUT=$8

HHb=$9
HHs=${10}
HHf=${11}

mkdir -p $4
mkdir -p $WDIR
mkdir -p $DDIR
mkdir -p $LOGD

conda activate RoseTTAFold
############################################################
# 1. generate MSAs
############################################################
if [ ! -s $WDIR/"$OUT"_.msa0.a3m ]
then
    echo "Running HHblits"
    "$SCRIPT_DIR"/make_msa.sh $IN $WDIR $CPU $MEM $DB "$OUT" $HHb $HHf > $LOGD/make_msa.stdout 2> $LOGD/make_msa.stderr
fi


############################################################
# 2. predict secondary structure for HHsearch run
############################################################
if [ ! -s $WDIR/"$OUT"_.ss2 ]
then
    echo "Running PSIPRED"
    "$SCRIPT_DIR"/make_ss.sh $WDIR/"$OUT"_.msa0.a3m $WDIR/"$OUT"_.ss2 > $LOGD/make_ss.stdout 2> $LOGD/make_ss.stderr
fi


############################################################
# 3. search for templates
############################################################
DB="$DB/pdb100_2021Mar03/pdb100_2021Mar03"
if [ ! -s $WDIR/"$OUT"_.hhr ]
then
    echo "Running hhsearch"
    HH="$HHs -b 50 -B 500 -z 50 -Z 500 -mact 0.05 -cpu $CPU -maxmem $MEM -aliw 100000 -e 100 -p 5.0 -d $DB"
    cat $WDIR/"$OUT"_.ss2 $WDIR/"$OUT"_.msa0.a3m > $WDIR/"$OUT"_.msa0.ss2.a3m
    $HH -i $WDIR/"$OUT"_.msa0.ss2.a3m -o $WDIR/"$OUT"_.hhr -atab $WDIR/"$OUT"_.atab -v 0 > $LOGD/hhsearch.stdout 2> $LOGD/hhsearch.stderr
fi

cp $WDIR/"$OUT"_.{msa0.a3m,hhr,atab} $DDIR
echo $OUT >> $DDIR/IDX

cp "$SCRIPT_DIR"/predict.py $PIPEDIR/network/predict_smile.py
cp "$SCRIPT_DIR"/model.py $PIPEDIR/network/RoseTTAFoldModel.py
