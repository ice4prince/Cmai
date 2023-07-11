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
date     : Jul  5, 2022
email    : Nasy <nasyxx+python@gmail.com>
filename : utils.py
project  : rosettafold
license  : GPL-3.0+

Utils
"""
# Types
from typing import Iterator


def read_fasta(path: str) -> Iterator[tuple[str, str]]:
    """Read fasta file to dictionary."""
    with open(path) as f:
        name, seq = "", ""
        for line in f:
            if line.startswith(">"):
                yield name, seq
                name = line.strip()[1:]
                seq = ""
            else:
                seq += line.strip()
        yield name, seq
