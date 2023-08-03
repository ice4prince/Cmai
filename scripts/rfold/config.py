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
date     : Jan 12, 2023
email    : Nasy <nasyxx+python@gmail.com>
filename : config.py
project  : rfold
license  : GPL-3.0+

RoseTTAFold in SMILE.
"""

# Standard Library
import os, shutil
from collections.abc import Callable
from dataclasses import asdict, dataclass, field, fields

# Types
from typing import Annotated, cast

# Utils
from rich import print
from rich.console import Console

# Config
import tyro


@dataclass
class Env:
    """Envs."""

    RF_BASE: Annotated[str, "RoseTTAFold repo path"] = os.environ.get(
        "RF_BASE", "./RoseTTAFold"
    )
    RF_DATA_BASE: Annotated[str, "RoseTTAFold data base path"] = os.environ.get(
        "RF_DATA_BASE", "./data"
    )
    RF_RUNTIME_BASE: Annotated[str, "RoseTTAFold runtime base path"] = os.environ.get(
        "RF_RUNTIME_BASE", "./runtime"
    )
    CSBLAST_DATA: Annotated[str, "CSBLAST data path"] = os.environ.get(
        "CSBLAST_DATA", "./deps/csblast-2.2.3/data"
    )
    PSIPRED_DATA: Annotated[str, "PSIPRED data path"] = os.environ.get(
        "PSIPRED_DATA", "./deps/psipred/data"
    )

    auto_gen_env: bool = True

    def gen_env(self) -> None:
        """Generate envs."""
        with open("env.sh", "w") as f:
            f.write("#!/bin/bash\n\n")
            for k, v in asdict(self).items():
                if k.isupper():
                    f.write(f"export {k}='{v}'\n")

    def __post_init__(self) -> None:
        """Post init."""
        os.environ.update(
            dict(filter(lambda kv: kv[0].isupper(), asdict(self).items()))
        )

        if self.auto_gen_env:
            self.gen_env()


@dataclass
class RuntimeEnv:
    """Runtime envs."""

    cpu: Annotated[int, "Max CPUs per process"] = 8
    processes: Annotated[int, "Max processes"] = 8
    mem: Annotated[int, "Max memory"] = 32
    use_cpu: Annotated[bool, "Use CPU or GPU"] = False


@dataclass
class Exe:
    """Executables."""

    hhsearch: Annotated[
        str, "HHsearch executable."
    ] = "deps/hh-suite/build/bin/hhsearch"
    hhblits: Annotated[str, "HHblits executable."] = "deps/hh-suite/build/bin/hhblits"
    hhfilter: Annotated[
        str, "HHfilter executable."
    ] = "deps/hh-suite/build/bin/hhfilter"

    psipred: Annotated[str, "PSIPRED executable."] = "deps/psipred/bin/psipred"
    psipass2: Annotated[str, "PSIPASS2 executable."] = "deps/psipred/bin/psipass2"

    csbuild: Annotated[str, "csbuild executable."] = "deps/csblast-2.2.3/bin/csbuild"

    makemat: Annotated[str, "makemat executable."] = "deps/blast-legacy/bin/makemat"


@dataclass
class DB:
    """DB prefix."""

    bfd: str = "bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt"
    uniref90: str = "UniRef30_2020_06/UniRef30_2020_06"
    pdb: str = "pdb100_2021Mar03/pdb100_2021Mar03"


@dataclass
class Conf:
    """Config."""

    in_fasta: str = "./*.fasta"

    env: Env = field(default_factory=Env)
    runtime: RuntimeEnv = field(default_factory=RuntimeEnv)

    exe: Exe = field(default_factory=Exe)

    db: DB = field(default_factory=DB)

    _conf_path: Annotated[str, "Don't change this!"] = field(init=False)

    def __post_init__(self) -> None:
        """Post init."""
        self._conf_path = __file__
        self._console = Console()

    @property
    def log(self) -> Callable[[str], None]:
        """Log function"""
        return self._console.log


def config() -> Conf:
    """Get the Config."""
    return tyro.cli(Conf)


def check_exe(exe: str) -> None:
    """Check exe in the PATH."""
    print(f"Checking {exe}...", end="")
    if shutil.which(exe):
        print("[green]OK![/green]")
        return
    raise OSError(f"{exe} not found in PATH.")


def check_exes(conf: Conf) -> None:
    """Check exes in the PATH."""
    for exe in fields(conf.exe):
        check_exe(getattr(conf.exe, exe.name))


if __name__ == "__main__":
    print(config())
