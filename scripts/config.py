#!/usr/bin/env python
# -*- coding: utf-8 -*-

r"""

author   : Bing Song https://nasy.moe
date     : Dec 9, 2022
email    : Sunny Song <bing.song@utsouthwestern.edu>
filename : settings.py
project  : DeLAnO
license  : GPL-3.0+

DeLAnO scripts config.
"""
# Standard Library
from dataclasses import asdict, dataclass, field

# Config
from smile_config import from_dataclass

# Types
from typing import Annotated, cast

# Others
from rich import print


@dataclass
class Path:

    data: Annotated[str, "Data base folder."] = "/project/DPDS/Wang_lab/shared/BCR_antigen/data"

    # bfd: str = "bfd"
    # pdb: str = "pdb100_2021Mar03"
    # msa: str = "UniRef30_2020_06"

    scripts: Annotated[str, "Scripts folder."] = "/project/DPDS/Wang_lab/shared/BCR_antigen/code/CLAnO/scripts"

    rf: Annotated[str, "RosettaFold folder."] = "rfscripts/RoseTTAFold"

    out: Annotated[str, "Output folder."] = "/project/DPDS/Wang_lab/shared/BCR_antigen/code/CLAnO/data/example/output/RFoutputs"

    logs: Annotated[str, "Logs folder."] = "logs"

@dataclass
class Exe:
    """Executables."""

    hhsearch: Annotated[str, "HHsearch executable."] = "rfscripts/bin/hhsearch"
    hhblits: Annotated[str, "HHblits executable."] = "rfscripts/bin/hhblits"
    hhfilter: Annotated[str, "HHfilter executable."] = "rfscripts/bin/hhfilter"
    psipred: Annotated[str, "PSIPRED executable."] = "psipred"
    psipass2: Annotated[str, "PSIPASS2 executable."] = "psipass2"

# @dataclass
# class Wrap:
#     """Previously wrapped"""
#
#     Vinput: Annotated[str, "BCR_Vh input"] = "/project/DPDS/Wang_lab/shared/BCR_antigen/code/DeLAnO/example/outputs/preprocess/V.txt"
#     Voutput: Annotated[str, "BCR_Vh output"] = "/project/DPDS/Wang_lab/shared/BCR_antigen/code/DeLAnO/example/outputs/V"
#     Vtype: Annotated[str,"type of V input"] = "sequence"
#     Vspecies: Annotated[str,"species"] = "human"
#     CDR3input: Annotated[str, "BCR_CDR3h input"] = "/project/DPDS/Wang_lab/shared/BCR_antigen/code/DeLAnO/example/outputs/preprocess/CDR3.txt"
#     CDR3output: Annotated[str, "BCR_CDR3h output"] = "/project/DPDS/Wang_lab/shared/BCR_antigen/code/DeLAnO/example/outputs/CDR3"


@dataclass
class Conf:
    """Config for DeLAnO script."""

    fasta: Annotated[str, "Input fasta files."] = "/project/DPDS/Wang_lab/shared/BCR_antigen/code/CLAnO/data/intermediates/antigens.fasta"

    verbose: Annotated[bool,"Print verbose messages or not."] = False
    cpu: Annotated[str, "Max CPUs."] = "32"
    mem: Annotated[str, "Max Memory (in GB)."] = "64"
    use_cpu: Annotated[str, "cpu|gpu"] = "gpu"

    gen_msa: Annotated[bool, "run generate_msa and exit"] = False
    run_rf: Annotated[bool, "Skip generate_msa and run RoseTTAFold"] = False
    skip_preprocess: Annotated[bool,"Skip preprocess"] = False
    skip_extract: Annotated[bool,"Skip extraction"] = False

    path: Path = Path()
    exe: Exe = Exe()
    # wrap: Wrap = Wrap()

config = cast(Conf, from_dataclass(Conf()).config)

if __name__ == "__main__":
    print("Show help info by adding `--help.`")
    if config.verbose:
        print(config)
