#!/usr/bin/env python
# -*- coding: utf-8 -*-

r"""
Python ♡ Nasy.

    |             *         *
    |                  .                .
    |           .                              登
    |     *                      ,
    |                   .                      至
    |
    |                               *          恖
    |          |\___/|
    |          )    -(             .           聖 ·
    |         =\ -   /=
    |           )===(       *
    |          /   - \
    |          |-    |
    |         /   -   \     0.|.0
    |  NASY___\__( (__/_____(\=/)__+1s____________
    |  ______|____) )______|______|______|______|_
    |  ___|______( (____|______|______|______|____
    |  ______|____\_|______|______|______|______|_
    |  ___|______|______|______|______|______|____
    |  ______|______|______|______|______|______|_
    |  ___|______|______|______|______|______|____

author   : Nasy https://nasy.moe
date     : Nov 13, 2022
email    : Nasy <nasyxx+python@gmail.com>
filename : tofastas.py
project  : rosettafold
license  : GPL-3.0+

Convert single fasta to mulitple fastas.
"""
# Others
from Bio import SeqIO
from fire import Fire
from tqdm import tqdm


def main(inf: str, out: str):
    """From input fasta file to out dir fasta files."""
    for key, seq in tqdm(SeqIO.index(inf, "fasta").items()):
        SeqIO.write(seq, f"{out}/{key}.fasta", "fasta")


if __name__ == "__main__":
    Fire(main)
