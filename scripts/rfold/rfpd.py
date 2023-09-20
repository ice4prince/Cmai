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
filename : rfpd.py
project  : rfold
license  : GPL-3.0+

Rfold Preprocessor for Data.
"""
from .rconfig import Conf, Env, config

# Standard Library
import shlex, shutil, subprocess
from os import listdir, makedirs, remove, environ
from os.path import abspath, basename, dirname, exists, join

# Types
from typing import cast


def preprocess(env: Env) -> None:
    """Preprocess runtime env."""
    runtime = env.RF_RUNTIME_BASE
    makedirs(runtime, exist_ok=True)
    makedirs(join(runtime, "work"), exist_ok=True)
    makedirs(join(runtime, "data"), exist_ok=True)
    makedirs(join(runtime, "pred"), exist_ok=True)
    makedirs(join(runtime, "temp.fasta"), exist_ok=True)
    makedirs(join(runtime, "work/logs"), exist_ok=True)
    makedirs(join(runtime, "work/temp"), exist_ok=True)
    makedirs(join(runtime, "work/hhblits"), exist_ok=True)


def hhblits(in_fasta: str, a3m_prefix: str, db: str, e: str, conf: Conf) -> list[str]:
    """Build hhblits cmd."""
    cmd = shlex.split(
        f"{conf.exe.hhblits} -o /dev/null -mact 0.35 -maxfilt 100000000 -neffmax 20 "
        f"-cov 25 -cpu {conf.runtime.cpu} -nodiff "
        f"-realign_max 100000000 -maxseq 1000000 -maxmem {conf.runtime.mem} -n 4 "
        f"-i {in_fasta} -oa3m {a3m_prefix}.a3m -e {e} -v 0 -d {db}"
    )
    conf.log(f"HHBLITS cmds: {cmd}")
    conf.log(" ".join(cmd))
    return cmd


def hhfilter(a3m_prefix: str, cov: int, conf: Conf) -> list[str]:
    """Build hhfilter cmd."""
    cmd = shlex.split(
        f"{conf.exe.hhfilter} -id 90 -cov {cov} "
        f"-i {a3m_prefix}.a3m -o {a3m_prefix}.id90cov{cov}.a3m"
    )
    conf.log(f"HHFILTER cmds: {cmd}")
    conf.log(" ".join(cmd))
    return cmd


def one_search(
    in_fasta: str, out: str, dbi: int, db: str, e: str, conf: Conf
) -> tuple[bool, str]:
    """One search.

    Return (success, out_a3m_path)
    """
    conf.log(f"HHBLITS: {in_fasta}, {out}, db:{dbi}:{db}, e:{e}")

    out_dir = join(conf.env.RF_RUNTIME_BASE, "work")
    a3m_prefix = join(out_dir, "hhblits", f"{out}_.db{dbi}.{e}")

    logs = join(out_dir, "logs")

    with (
        open(join(logs, "make_msa.stdout"), "ab") as stdout,
        open(join(logs, "make_msa.stderr"), "ab") as stderr,
    ):
        # hhblits
        res = subprocess.run(
            hhblits(in_fasta, a3m_prefix, db, e, conf),
            stdout=stdout,
            stderr=stderr,
        )
        if res.returncode:
            raise RuntimeError(
                f"HHBLITS failed with {in_fasta=},{out=},{dbi=},{db=},{e=}."
            )

        # hhfilter
        # cov 75
        res = subprocess.run(
            hhfilter(a3m_prefix, 75, conf),
            stdout=stdout,
            stderr=stderr,
        )
        if res.returncode:
            raise RuntimeError(
                f"HHFILTER 75 failed with {in_fasta=},{out=},{dbi=},{db=},{e=}."
            )
        with open(f"{a3m_prefix}.id90cov75.a3m") as f75:
            n75 = len(list(filter(lambda _l: _l.startswith(">"), f75)))
            conf.log(f"Cov 75 result: {n75}")
        if n75 > 2000:
            return True, f"{a3m_prefix}.id90cov75.a3m"

        # cov 50
        res = subprocess.run(
            hhfilter(a3m_prefix, 50, conf),
            stdout=stdout,
            stderr=stderr,
        )
        if res.returncode:
            raise RuntimeError(
                f"HHFILTER 50 failed with {in_fasta=},{out=},{dbi=},{db=},{e=}."
            )
        with open(f"{a3m_prefix}.id90cov50.a3m") as f50:
            n50 = len(list(filter(lambda _l: _l.startswith(">"), f50)))
            conf.log(f"Cov 75 result: {n75}")
        if n50 > 4000:
            return True, f"{a3m_prefix}.id90cov50.a3m"
    return False, f"{a3m_prefix}.id90cov75.a3m"


def make_msa(path: str, key: str, conf: Conf) -> str:
    """Generate MSAs.

    Make MSA for in_fasta PATH and output fasta KEY.
    """
    dbs = (
        join(conf.env.RF_DATA_BASE, conf.db.uniref90),
        join(conf.env.RF_DATA_BASE, conf.db.bfd),
    )

    if exists(a3m_target := join(conf.env.RF_RUNTIME_BASE, "work", f"{key}_.msa0.a3m")):
        conf.log(f"MSA exists: {a3m_target}")
        return a3m_target
    for dbi, db in enumerate(dbs):
        for e in "1e-30 1e-10 1e-6 1e-3".split():
            p, a3m = one_search(path, key, dbi, db, e, conf)
            if p:
                return cast(str, shutil.copy2(a3m, f"{a3m_target}"))
            shutil.copy2(a3m, f"{a3m_target}_bk")
    return cast(str, shutil.move(f"{a3m_target}_bk", a3m_target))


def ss_cmd(
    name: str,
    cmd: str,
    conf: Conf,
    *,
    stdout: None | str = None,
    env: None | dict = None,
) -> None:
    """Run ss CMD."""
    logs = join(conf.env.RF_RUNTIME_BASE, "work", "logs")
    cmds = shlex.split(cmd)
    if env is None:
        env = {}
    with (
        open(stdout or join(logs, "make_ss.stdout"), "ab") as stdoutf,
        open(join(logs, "make_ss.stderr"), "ab") as stderr,
    ):
        conf.log(f"{name}: {cmds}")
        conf.log(" ".join(cmds))
        res = subprocess.run(
            cmds, stdout=stdoutf, stderr=stderr, env={**environ, **env}
        )
        if res.returncode:
            raise RuntimeError(f"make_ss failed with {name}, {cmd=}")


def make_ss(a3m: str, conf: Conf) -> str:
    """Predict secondary structure for HHsearch run."""
    ss2 = a3m.replace(".msa0.a3m", ".ss2")
    if exists(ss2):
        conf.log(f"SS exists: {ss2}")
        return ss2

    ID = basename(a3m).replace(".a3m", ".tmp")

    ss_cmd(
        "csbuild",
        f"{conf.exe.csbuild} -i {a3m} -I a3m -D {conf.env.CSBLAST_DATA}/K4000.crf "
        f"-o {ID}.chk -O chk",
        conf,
    )

    with open(a3m) as fr, open(f"{ID}.fasta", "w") as gr:
        gr.write(fr.readline())
        gr.write(fr.readline())
    with open(f"{ID}.pn", "w") as fw, open(f"{ID}.sn", "w") as gw:
        fw.write(f"{ID}.chk\n")
        gw.write(f"{ID}.fasta\n")

    ss_cmd(
        "makemat",
        f"{conf.exe.makemat} -P {ID}",
        conf,
        env={"BLASTMAT": conf.env.BLASTMAT},
    )

    ddir = conf.env.PSIPRED_DATA
    ss_cmd(
        "psipred",
        f"{conf.exe.psipred} {ID}.mtx "
        f"{ddir}/weights.dat {ddir}/weights.dat2 {ddir}/weights.dat3",
        conf,
        stdout=f"{ID}.ss",
    )
    ss_cmd(
        "psipass",
        f"{conf.exe.psipass2} {ddir}/weights_p2.dat 1 1.0 1.0 "
        f"{a3m}.csb.hhblits.ss2 {ID}.ss",
        conf,
        stdout=f"{ID}.horiz",
    )

    with open(f"{ID}.horiz") as fr2, open(ss2, "w") as gr2:
        gr2.write(">ss_pred\n")
        for lp in fr2:

            if lp.startswith("Pred") and len(lp.split()) == 2:

                gr2.write(lp.split()[1])
                gr2.write("\n")

        gr2.write(">ss_conf\n")
        fr2.seek(0)
        for lc in fr2:

            if lc.startswith("Conf") and len(lc.split()) == 2:

                gr2.write(lc.split()[1])
                gr2.write("\n")

    for fs in filter(lambda _f: _f.startswith(f"{ID}"), listdir(".")):
        remove(fs)
    remove(f"{a3m}.csb.hhblits.ss2")

    return ss2


def make_hhr_atab(a3m: str, ss2: str, conf: Conf) -> tuple[str, str]:
    """Search for templates."""
    hhr = a3m.replace(".msa0.a3m", ".hhr")
    atab = ss2.replace(".ss2", ".atab")
    if exists(hhr) and exists(atab):
        conf.log(f"HHR exists:  {hhr}")
        conf.log(f"atab exists: {atab}")
        return hhr, atab

    ss2a3m = a3m.replace(".a3m", ".ss2.a3m")
    db = join(conf.env.RF_DATA_BASE, conf.db.pdb)

    with open(a3m) as f, open(ss2) as g, open(ss2a3m, "w") as h:
        h.write(g.read())
        h.write(f.read())

    logs = join(conf.env.RF_RUNTIME_BASE, "work", "logs")
    with (
        open(join(logs, "hhsearch.stdout"), "ab") as stdout,
        open(join(logs, "hhsearch.stderr"), "ab") as stderr,
    ):
        hhcmd = shlex.split(
            f"{conf.exe.hhsearch} -b 50 -B 500 -z 50 -Z 500 -mact 0.05 "
            f"-cpu {conf.runtime.cpu} -maxmem {conf.runtime.mem} "
            f"-aliw 100000 -e 100 -p 5.0 -d {db} "
            f"-i {ss2a3m} -o {hhr} -atab {atab} -v 0"
        )
        conf.log(f"HHsearch: {hhcmd}")
        conf.log(" ".join(hhcmd))
        res = subprocess.run(hhcmd, stdout=stdout, stderr=stderr)
        if res.returncode:
            raise RuntimeError(f"HHsearch failed with {a3m=},{ss2=}")

    return hhr, atab


def to_data_dir(key: str, a3m: str, ss2: str, hhr: str, atab: str, conf: Conf) -> None:
    """Move files to data dir."""
    data_dir = join(conf.env.RF_RUNTIME_BASE, "data")
    for ori in (a3m, ss2, hhr, atab):
        shutil.copy2(ori, data_dir)
    with open(join(data_dir, "IDX"), "a") as f:
        f.write(f"{key}\n")


def test() -> None:
    """Run test."""
    conf = config()
    preprocess(conf.env)
    a3m = make_msa(
        f"{dirname(abspath(__file__))}/../../xx.fasta", "Eta", conf
    )  # passed
    ss2 = make_ss(a3m, conf)
    hhr, atab = make_hhr_atab(a3m, ss2, conf)
    to_data_dir("Eta", a3m, ss2, hhr, atab, conf)


if __name__ == "__main__":
    test()
