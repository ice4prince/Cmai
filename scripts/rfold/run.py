#!/usr/bin/env python
# -*- coding: utf-8 -*-

r"""Python ♡ Nasy.

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
date     : Aug  3, 2023
email    : Nasy <nasyxx+python@gmail.com>
filename : run.sh
project  : Cmai
license  : GPL-3.0+

Run rfpd.py
"""
from config import Conf, config
from rfpd import preprocess, make_msa, make_ss, make_hhr_atab, to_data_dir
from glob import iglob
from Bio import SeqIO, SeqRecord
import os


def make_temp_fasta(key: str, seq: SeqRecord.SeqRecord, conf: Conf) -> str:
    """Make temp fasta file."""
    fasta = os.path.join(conf.env.RF_RUNTIME_BASE, "temp.fasta", f"{key}.fasta")
    seq.seq = seq.seq.upper()  # type: ignore[attr-defined]
    SeqIO.write(seq, fasta, "fasta")
    return fasta


def main() -> None:
    """Run main function."""
    conf = config()
    preprocess(conf.env)
    for fasta in iglob(conf.in_fasta):
        for key, seq in SeqIO.index(fasta, "fasta").items():
            fpath = make_temp_fasta(key, seq, conf)
            a3m = make_msa(fpath, key, conf)
            ss2 = make_ss(a3m, conf)
            hhr, atab = make_hhr_atab(a3m, ss2, conf)
            to_data_dir(key, a3m, ss2, hhr, atab, conf)


if __name__ == "__main__":
    main()
