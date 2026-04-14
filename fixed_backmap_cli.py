#!/usr/bin/env python3
import argparse
import pickle
import shlex
import subprocess
import re
from pathlib import Path

try:
    from backmapping_funcs_with_cg2all_v2 import backmap_ca_chain_multimer, constraints_to_residue_pairs
except ImportError:
    from backmapping_funcs import backmap_ca_chain_multimer, constraints_to_residue_pairs

_CHI_RE = re.compile(r"Chi\^2\s*=\s*([0-9.eE+-]+)")


def run_foxs(pyfoxs_cmd: str, aa_pdb: Path, saxs_dat: Path, max_q: float,
             outdir: Path, summary_file: Path | None = None):
    outdir.mkdir(parents=True, exist_ok=True)
    foxs_log = outdir / (aa_pdb.stem + "_foxs.log")

    cmd = shlex.split(pyfoxs_cmd) + [
        str(aa_pdb),
        str(saxs_dat),
        "--max_q", str(max_q),
    ]

    p = subprocess.run(cmd, capture_output=True, text=True)

    foxs_log.write_text(
        f"$ {' '.join(cmd)}\n\n=== STDOUT ===\n{p.stdout}\n\n=== STDERR ===\n{p.stderr}\n"
    )

    chisq = None
    m = _CHI_RE.search(p.stdout)
    if m:
        chisq = m.group(1)

    if summary_file is not None:
        summary_file.parent.mkdir(parents=True, exist_ok=True)
        with open(summary_file, "a") as f:
            if p.returncode == 0 and chisq is not None:
                f.write(f"{aa_pdb} {chisq}\n")
            else:
                f.write(f"{aa_pdb} ERROR\n")

    return p.returncode, foxs_log, chisq

def main():
    ap = argparse.ArgumentParser(description="Carbonara coords -> all-atom backmapping + optional FoXS")
    ap.add_argument("--coords", required=True)
    ap.add_argument("--fingerprint", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--name", required=True)
    ap.add_argument("--backend", choices=["modeller", "cg2all"], default="modeller")
    ap.add_argument("--scenario-root", required=True)
    ap.add_argument("--disulfide-file", default=None)
    ap.add_argument("--cg2all-exec", default=None)
    ap.add_argument("--do-foxs", action="store_true")
    ap.add_argument("--foxs-py", default="pyFoXS/pyFoXS/foxs.py")
    ap.add_argument("--saxs", default=None)
    ap.add_argument("--max-q", type=float, default=None)
    ap.add_argument("--foxs-out", default=None)
    args = ap.parse_args()

    if args.do_foxs and (args.saxs is None or args.max_q is None):
        ap.error("--do-foxs requires --saxs and --max-q")

    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    scenario_root = Path(args.scenario_root).resolve()

    with open(scenario_root / "chainLengths.dat", "rb") as f:
        chain_lengths = pickle.load(f)
    lengths = list(chain_lengths.values())

    disulfides = None
    if args.disulfide_file:
        disulfides = constraints_to_residue_pairs(args.fingerprint, args.disulfide_file)

    cg2all_exec = None
    if args.backend == "cg2all":
        if args.cg2all_exec is None:
            ap.error("--backend cg2all requires --cg2all-exec")
        cg2all_exec = shlex.split(args.cg2all_exec)

    backmap_ca_chain_multimer(
        args.coords,
        args.fingerprint,
        str(outdir),
        args.name,
        lengths,
        disulfides=disulfides,
        method=args.backend,
        cg2all_exec=cg2all_exec,
    )

    aa_pdb_path = outdir / f"{args.name}_AA.pdb"

    if args.do_foxs:
        summary = Path(args.foxs_out) if args.foxs_out else None
        rc, log_path, chisq = run_foxs(
            pyfoxs_cmd=args.foxs_py,
            aa_pdb=aa_pdb_path,
            saxs_dat=Path(args.saxs).resolve(),
            max_q=args.max_q,
            outdir=outdir,
            summary_file=summary,
        )
        if rc != 0:
            print(f"FoXS failed (exit {rc}). See: {log_path}")
        else:
            print(f"FoXS complete (Chi^2={chisq}). Log: {log_path}")

if __name__ == "__main__":
    main()
