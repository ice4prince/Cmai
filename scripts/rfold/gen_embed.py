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
date     : Aug  7, 2023
email    : Nasy <nasyxx+python@gmail.com>
filename : gen_embed.py
project  : Cmai
license  : GPL-3.0+

Generate embeddings.
"""

# Standard Library
import subprocess
import os

# Others
from .rconfig import Conf, config
import tyro


def gen_embed(conf: Conf) -> None:
    """Generate embeddings from rosettafold."""
    script_dir = os.path.dirname(os.path.realpath(__file__))
    subprocess.run([
        "bash",
        f"{script_dir}/gen_embed.sh",
        conf.env.RF_BASE,
        conf.env.RF_DATA_BASE,
        conf.env.RF_RUNTIME_BASE,
        f"{conf.runtime.use_cpu}"
    ])

if __name__ == '__main__':
    gen_embed(tyro.cli(Conf))
