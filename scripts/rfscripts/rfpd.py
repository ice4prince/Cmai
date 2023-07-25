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
date     : Jun 29, 2022
email    : Nasy <nasyxx+python@gmail.com>
filename : rfpd.py
project  : rosettafold
license  : GPL-3.0+

Rosettafold scripts
"""
# Standard Library
import os, shutil
from glob import iglob

# Config
from config import Conf, config

# Others
from Bio import SeqIO
from rich import print
from rich.prompt import Confirm
from rich.traceback import install as __install

__install()

def generate_msa(conf: Conf) -> None:
    """Generate msas."""
    InFasta = conf.path.out + '/preprocess/'+ conf.fasta.rsplit('/', 1)[1]
    PathOut = conf.path.out + '/results'
    RF_dir = os.path.join(conf.path.scripts,conf.path.rf)
    for fasta in iglob(InFasta):
        for key, seq in SeqIO.index(fasta, "fasta").items():
            print(f"Processing:\n{seq}")
            os.makedirs(f"{conf.path.data}/temp/", exist_ok=True)
            path = f"{conf.path.data}/temp/{key}"
           if ps := len(list(iglob(f"{path}*"))):
               path = f"{path}_{ps:03}"
               key = f"{key}_{ps:03}"
            SeqIO.write(seq, path, "fasta")
            script_dir = os.path.dirname(os.path.realpath(__file__))
            print(
                f"{script_dir}/gmsa.sh "
                f"{RF_dir} {conf.cpu} {conf.mem} "
                f"{PathOut} {conf.path.logs} {conf.path.data} "
                f"{path} {key} "
                f"{conf.exe.hhblits} {conf.exe.hhsearch} {conf.exe.hhfilter}"
            )
            os.system(
                f"{script_dir}/gmsa.sh "
                f"{RF_dir} {conf.cpu} {conf.mem} "
                f"{PathOut} {conf.path.logs} {conf.path.data} "
                f"{path} {key} "
                f"{conf.exe.hhblits} {conf.exe.hhsearch} {conf.exe.hhfilter}"
            )


def check_exe(exe: str) -> None:
    """Check exe in the PATH."""
    print(f"Checking {exe}...", end="")
    if shutil.which(exe):
        print("[green]OK!")
        return
    raise OSError(f"{exe} not found in PATH.")


def check_exes(conf: Conf) -> None:
    """Check exes in the PATH."""
    for exe in conf.exe.__dataclass_fields__:
        check_exe(exe)


def run_rf(conf: Conf):
    """Run rosetta fold."""
    script_dir = os.path.dirname(os.path.realpath(__file__))
    PathOut = conf.path.out + '/results'
    RF_dir = os.path.join(conf.path.scripts,conf.path.rf)
    cmd = f"{script_dir}/run.sh {RF_dir} {conf.path.data} {PathOut} {conf.use_cpu}"
    print(cmd)
    os.system(cmd)

def run_msa(conf: Conf):
#        if Confirm.ask("Is the config correct?", default=True):
            # check_exes(config)
    generate_msa(conf)
def run_rfpd(conf: Conf):
    run_rf(conf)



if __name__ == "__main__":
    print(config)
    if not config.skip_gen_msa:
        run_msa(config)
#        if Confirm.ask("Is the config correct?", default=True):
            # check_exes(config)
#       generate_msa(config)
    if not config.skip_run_rf:
        run_rfpd(config)
#    if config.run_rf:
#        if Confirm.ask("Run rosetta?", default=True):
#            run_rf(config)
