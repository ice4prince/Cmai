#!/bin/bash

# make the script stop when error (non-true exit code) is occured
set -e

# inputs
in_fasta="$1"
out_dir="$2"

# resources
CPU="$3"
MEM="$4"

DB=$5
OUT=$6

HHb=$7
HHf=$8

# sequence databases
declare -a DATABASES=( \
    "$DB/UniRef30_2020_06/UniRef30_2020_06" \
    "$DB/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt")

# setup hhblits command
HHBLITS="$HHb -o /dev/null -mact 0.35 -maxfilt 100000000 -neffmax 20 -cov 25 -cpu $CPU -nodiff -realign_max 100000000 -maxseq 1000000 -maxmem $MEM -n 4"
echo $HHBLITS

tmp_dir="$out_dir/hhblits"
mkdir -p $tmp_dir
out_prefix="$out_dir/""$OUT""_"

# perform iterative searches
prev_a3m="$in_fasta"
for (( i=0; i<${#DATABASES[@]}; i++ ))
do
    for e in 1e-30 1e-10 1e-6 1e-3
    do
        echo db:$i e:$e
        $HHBLITS -i $prev_a3m -oa3m $tmp_dir/"$OUT"_.db$i.$e.a3m -e $e -v 0 -d ${DATABASES[$i]}

        $HHf -id 90 -cov 75 -i $tmp_dir/"$OUT"_.db$i.$e.a3m -o $tmp_dir/"$OUT"_.db$i.$e.id90cov75.a3m
        prev_a3m="$tmp_dir/"$OUT"_.db$i.$e.id90cov50.a3m"
        n75=`grep -c "^>" $tmp_dir/"$OUT"_.db$i.$e.id90cov75.a3m`
        echo " -- cov75 results: $n75"
        if (($n75>2000))
        then
            cp $tmp_dir/"$OUT"_.db$i.$e.id90cov75.a3m ${out_prefix}.msa0.a3m
            break 2
        fi

        $HHf -id 90 -cov 50 -i $tmp_dir/"$OUT"_.db$i.$e.a3m -o $tmp_dir/"$OUT"_.db$i.$e.id90cov50.a3m
        n50=`grep -c "^>" $tmp_dir/"$OUT"_.db$i.$e.id90cov50.a3m`
        echo " -- cov50 results: $n50"
        if (($n50>4000))
        then
            cp $tmp_dir/"$OUT"_.db$i.$e.id90cov50.a3m ${out_prefix}.msa0.a3m
            break 2
        fi
    done
done

if [ ! -s ${out_prefix}.msa0.a3m ]
then
    cp $prev_a3m ${out_prefix}.msa0.a3m
fi
