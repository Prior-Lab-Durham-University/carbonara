#!/usr/bin/env python3
import argparse
from pathlib import Path
import sys
import subprocess

# IMPORTANT: this import must resolve exactly the same way it does
# in your notebooks / scripts
from backmapping_funcs import backmap_ca_chain
import subprocess

import re

_CHI_RE = re.compile(r"Chi\^2\s*=\s*([0-9.eE+-]+)")

def run_foxs(pyfoxs_script: Path, aa_pdb: Path, saxs_dat: Path, max_q: float,
             outdir: Path, summary_file: Path | None = None):
    outdir.mkdir(parents=True, exist_ok=True)
    foxs_log = outdir / (aa_pdb.stem + "_foxs.log")

    cmd = [
        "python3",
        str(pyfoxs_script),
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

    # Append one-line result (real-time)
    if summary_file is not None:
        summary_file.parent.mkdir(parents=True, exist_ok=True)
        with open(summary_file, "a") as f:
            if p.returncode == 0 and chisq is not None:
                f.write(f"{aa_pdb} {chisq}\n")
            else:
                f.write(f"{aa_pdb} ERROR\n")

    return p.returncode, foxs_log, chisq



def main():
    ap = argparse.ArgumentParser(
        description="Carbonara CA → all-atom backmapping (monomer)"
    )

    ap.add_argument("--coords", required=True,
                    help="Carbonara coords (.dat)")
    ap.add_argument("--fingerprint", required=True,
                    help="Fingerprint / secondary structure file")
    ap.add_argument("--outdir", required=True,
                    help="Directory to write outputs into")
    ap.add_argument("--name", required=True,
                    help="Base name for outputs (no extension)")
    ap.add_argument("--rate", default="fast",
                    choices=["fast", "slow"],
                    help="Backmapping rate (default: fast)")
    ap.add_argument("--no-ss", action="store_true",
                    help="Disable secondary structure constraints")
    ap.add_argument("--do-foxs", action="store_true")
    ap.add_argument("--foxs-py", default="pyFoXS/pyFoXS/foxs.py")
    ap.add_argument("--saxs", default=None)
    ap.add_argument("--max-q", type=float, default=None)
    ap.add_argument("--foxs-out", default=None,
                help="Append FoXS chi^2 results to this file (one line per AA PDB)")



    args = ap.parse_args()

    Path(args.outdir).mkdir(parents=True, exist_ok=True)

    backmap_ca_chain(
        coords_file=args.coords,
        fingerprint_file=args.fingerprint,
        write_directory=args.outdir,
        name=args.name,
        ss_constraint=not args.no_ss,
        rate=args.rate
    )
    
    aa_pdb_path = Path(args.outdir) / f"{args.name}_AA__{args.rate}.pdb"
    
    aa_pdb_path = Path(args.outdir) / f"{args.name}_AA__{args.rate}.pdb"

    if args.do_foxs:
        summary = Path(args.foxs_out) if args.foxs_out else None

        rc, log_path, chisq = run_foxs(
            pyfoxs_script=Path(args.foxs_py).resolve(),
            aa_pdb=aa_pdb_path,
            saxs_dat=Path(args.saxs).resolve(),
            max_q=args.max_q,
            outdir=Path(args.outdir),
            summary_file=summary,
        )
        if rc != 0:
            print(f"FoXS failed (exit {rc}). See: {log_path}")
        else:
            print(f"FoXS complete (Chi^2={chisq}). Log: {log_path}")


if __name__ == "__main__":
    main()
