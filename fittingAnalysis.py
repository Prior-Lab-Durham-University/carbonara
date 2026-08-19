from matplotlib import pyplot as plt
import numpy as np
import mdtraj as md
import re
from pathlib import Path
import json
import pandas as pd
import math
import os
import glob
import CarbonaraDataTools as CDT



import base64
from pathlib import Path
import subprocess
import contextlib

import matplotlib.pyplot as plt
from IPython.display import HTML, display

# handling the possibility py3D didn't automatically install

try:
    import py3Dmol
    HAS_PY3DMOL = True
    _PY3DMOL_IMPORT_ERROR = None
except ImportError as e:
    py3Dmol = None
    HAS_PY3DMOL = False
    _PY3DMOL_IMPORT_ERROR = e
import string
import io
from Bio.PDB import PDBParser, PDBIO, Superimposer
from tqdm import tqdm
from collections import defaultdict
from typing import List, Tuple, Optional
from tempfile import NamedTemporaryFile
from pdbfixer import PDBFixer
from openmm.app import PDBFile

#####################################################
## Warning function for the visulisation routines which require pymol3d
#####################################################
_PY3DMOL_WARNED = False

def _warn_missing_py3dmol():
    global _PY3DMOL_WARNED
    if _PY3DMOL_WARNED:
        return
    print(
        "Warning: py3Dmol is not installed, so inline notebook visualisation is unavailable. "
        "The rest of fittingAnalysis has imported correctly."
    )
    _PY3DMOL_WARNED = True

#####################################################
## count the number of fitting runs available to a user
#####################################################


def count_fitting_runs(folder_path):
    """
    Counts fitLog{i}.dat files in the given folder and returns a message.
    """
    pattern = os.path.join(folder_path, "fitLog*.dat")
    files = glob.glob(pattern)
    n_runs = len(files)
    return f"There are {n_runs} independent fitting runs for this molecule."

#####################################################
## To Read a fitlog file
#####################################################

def read_fitlog_runs(path):
    """
    Reads the NDJSON (JSON-lines) log file and returns a list of DataFrames,
    one per 'run'.

    A "run" is defined as a block of ImprovementIndex records following a
    'Run' metadata line.

    Each returned DataFrame has:
        ImprovementIndex      - 0..N-1 (run-local index, useful for merging)
        ImprovementIndexOrig  - original index from the log (for reference)
        FitStep, ScatterFitFirst, OverlapPenalty, ElapsedTimeMin
    """
    runs = []
    current_run = []

    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            d = json.loads(line)

            # Starting a new run: push the previous one (if any)
            if "Run" in d and "ImprovementIndex" not in d:
                if current_run:
                    runs.append(current_run)
                    current_run = []
                continue

            # Skip any other non-improvement lines
            if "ImprovementIndex" not in d:
                continue

            # Find the elapsed-time key robustly (handles µ/μ etc.)
            elapsed_key = None
            for k in d.keys():
                if k.lower().startswith("elapsedtime"):
                    elapsed_key = k
                    break
            if elapsed_key is None:
                raise KeyError("Could not find ElapsedTime key in a log record.")

            elapsed_us = d[elapsed_key]
            elapsed_min = float(elapsed_us) / 1e6 / 60.0  # µs -> minutes

            rec = {
                "ImprovementIndexOrig": int(d["ImprovementIndex"]),
                "FitStep": int(d.get("FitStep", -1)),
                "ScatterFitFirst": float(d.get("ScatterFitFirst", float("nan"))),
                "OverlapPenalty": float(d.get("OverlapPenalty", float("nan"))),
                "ElapsedTimeMin": elapsed_min,
            }
            current_run.append(rec)

    # Add the final run if there was one
    if current_run:
        runs.append(current_run)

    if not runs:
        raise ValueError("No improvement records found in the fit log.")

    # Convert each run to a DataFrame and give it a run-local ImprovementIndex
    run_dfs = []
    for r in runs:
        df = pd.DataFrame(r)
        # Sort by the original ImprovementIndex, then assign 0..N-1 local index
        df = df.sort_values("ImprovementIndexOrig").reset_index(drop=True)
        df["ImprovementIndex"] = range(len(df))  # run-local
        run_dfs.append(df)

    return run_dfs

#####################################################
## Ensure the ordering of a run
#####################################################

def _scatter_sort_key(fp: str) -> int:
    """
    Sort key for scattering filenames so that:
      initial -> step_0 -> step_1 -> step_2 -> ... -> end

    Assumes names like:
      ...initial...,
      ...step_0..., step_1..., etc,
      ...end...
    """
    name = os.path.basename(fp).lower()

    # initial first
    if "initial" in name:
        return 0

    # step_<n> in numeric order
    m = re.search(r"step[_-]?(\d+)", name)
    if m:
        # Put steps after initial; step_0, step_1, ...
        return int(m.group(1)) + 1

    # end last
    if "end" in name:
        return 10**9

    # fallback (odd names) near the end but before 'end'
    return 10**8

_num = re.compile(r"^[+-]?(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?$")

#####################################################
## To read a scattering file
#####################################################

def read_scattering_list(path):
    """
    Reads 'path value' list. If a line has no numeric value (e.g. ends with ERROR),
    it uses the previous valid value (carry-forward).
    """
    rows = []
    last_value = None

    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for i, raw in enumerate(f, 1):
            line = raw.strip()
            if not line:
                continue
            parts = line.split()

            # find the last numeric token anywhere on the line
            j = None
            for k in range(len(parts) - 1, -1, -1):
                if _num.match(parts[k]):
                    j = k
                    break

            if j is not None:
                value = float(parts[j])
                file_path = " ".join(parts[:j])
                last_value = value
            else:
                # no number on this line (e.g., "... ERROR"): carry forward
                file_path = " ".join(parts[:-1]) if parts and parts[-1].upper() == "ERROR" else " ".join(parts)
                if last_value is None:
                    raise ValueError(f"Line {i}: no numeric value and no previous value to carry forward: {raw!r}")
                value = last_value

            rows.append({"FilePath": file_path, "ScatteringValue": value})

    if not rows:
        raise ValueError("No entries found in scattering values list.")

    df = pd.DataFrame(rows)
    df["SortKey"] = df["FilePath"].apply(_scatter_sort_key)
    df = df.sort_values("SortKey").reset_index(drop=True)
    df["ImprovementIndex"] = range(len(df))

    return df[["ImprovementIndex", "FilePath", "ScatteringValue"]]

#####################################################
## Read a fit lof and (all atomistis) scattering list file
#####################################################

def getRunData(FITLOG_FILE, SCATTER_LIST_FILE):
    """
    - Read ALL runs from FITLOG_FILE.
    - Read scattering list from SCATTER_LIST_FILE.
    - Select the *last* run whose number of entries matches the scattering file.
    - Merge by run-local ImprovementIndex (0..N-1).
    """
    run_dfs = read_fitlog_runs(FITLOG_FILE)
    scat_df = read_scattering_list(SCATTER_LIST_FILE)

    n_scat = len(scat_df)

    matching = [df for df in run_dfs if len(df) == n_scat]
    if not matching:
        lengths = [len(df) for df in run_dfs]
        raise ValueError(
            f"No run in {FITLOG_FILE} has {n_scat} entries. "
            f"Run lengths found: {lengths}"
        )

    # your policy: if multiple runs match, take the last one
    fit_df = matching[-1]

    merged = pd.merge(
        fit_df,
        scat_df,
        on="ImprovementIndex",
        how="inner",
        validate="one_to_one",
    ).sort_values("ImprovementIndex").reset_index(drop=True)

    merged = merged[[
        "ImprovementIndex",        # 0..N-1, same logical time index for both
        "ImprovementIndexOrig",    # as logged by FOXS/fit code
        "FitStep",
        "ScatterFitFirst",
        "OverlapPenalty",
        "ElapsedTimeMin",
        "ScatteringValue",
        "FilePath",
    ]]

    return merged

def load_runs_DF(baseFolder):
    pattern = os.path.join(baseFolder, "fitLog*.dat")
    files = glob.glob(pattern)
    n_runs = len(files)
    FITLOG_FILES = [baseFolder+"/fitLog"+str(i)+".dat" for i in range(1,n_runs+1)]
    SCATTER_LIST_FILES = [baseFolder+"/allAtomRun"+str(i)+"/foxsFull.dat" for i in range(1,n_runs+1)] # the file containing "path value" per line
    return [getRunData(FITLOG_FILES[i],SCATTER_LIST_FILES[i]) for i in range(100)]

#####################################################
## Plot a vriable
#####################################################

def plot_vars(
    df,
    x,
    y,
    color=None,
    title=None,
    kind="scatter",      
    fit_line=False,      
    annotate=False,      
    figsize=(7, 5),
    marker="o",
    # NEW:
    xlim=None,           
    ylim=None,           
    x_label=None,        
    y_label=None,        
    color_label=None,    
    save_path=None,      
    dpi=300,
    # --- NEW style controls ---
    font_size=16,        # base font size
    title_size=18,       # title font size
    legend_size=14,      # legend font size
    marker_size=120,     # scatter marker size
    line_width=2.5,       # line thickness
    tick_size=14, 
    legend_labels=None  
):
    """
    Flexible plot from one or multiple DataFrames (publication-ready).
    """
    import numpy as np
    import matplotlib.pyplot as plt

    dfs = df if isinstance(df, (list, tuple)) else [df]
    fig, ax = plt.subplots(figsize=figsize)

    if kind == "scatter":
        sc = None
        for d in dfs:
            sc = ax.scatter(
                d[x], d[y],
                c=d[color] if color else "C0",
                cmap="viridis", s=marker_size,
                edgecolor="k", alpha=0.8, marker=marker
            )

            if fit_line:
                coeffs = np.polyfit(d[x], d[y], 1)
                poly = np.poly1d(coeffs)
                xx = np.linspace(d[x].min(), d[x].max(), 200)
                ax.plot(xx, poly(xx), "--", lw=line_width,
                        label=f"fit: y={coeffs[0]:.2f}x+{coeffs[1]:.2f}",)
                ax.legend(fontsize=legend_size)

            if annotate and "ImprovementIndex" in d.columns:
                for _, row in d.iterrows():
                    ax.text(row[x], row[y], str(row["ImprovementIndex"]),
                            fontsize=font_size-2, ha="right", va="bottom")

        if color and sc is not None:
            cbar = plt.colorbar(sc, ax=ax)
            cbar.set_label(color_label if color_label else color, fontsize=font_size)

    elif kind == "line":
        if color:
            for d in dfs:
                for key, grp in d.groupby(color):
                    ax.plot(grp[x], grp[y], marker=marker,
                            lw=line_width, label=f"{color}={key}")
            ax.legend(title=(color_label if color_label else color),
                      fontsize=legend_size, title_fontsize=legend_size)
        else:
            for i, d in enumerate(dfs):
                label = legend_labels[i] if legend_labels and i < len(legend_labels) else f"df{i+1}"
                ax.plot(d[x], d[y], marker=marker, lw=line_width, label=label)
            if len(dfs) > 1:
                ax.legend(fontsize=legend_size)

    else:
        raise ValueError("kind must be 'scatter' or 'line'")

    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)

    # Labels & title with bigger fonts
    ax.set_xlabel(x_label if x_label else x, fontsize=font_size)
    ax.set_ylabel(y_label if y_label else y, fontsize=font_size)
    ax.tick_params(axis="both", which="major", labelsize=tick_size)
    if title is None:
        title = f"{y} vs {x}"
    ax.set_title(title, fontsize=title_size, weight="bold")

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        print(f"Plot saved to {save_path}")
    else:
        plt.show()
        plt.close(fig)


def _chains_present_in_pdb(pdb_text: str):
    """
    Return sorted unique chain IDs present in ATOM/HETATM records of a PDB string.
    """
    chains = set()
    for line in pdb_text.splitlines():
        if line.startswith(("ATOM  ", "HETATM")) and len(line) >= 22:
            ch = line[21].strip()
            if ch:  # ignore blank chain IDs
                chains.add(ch)
    # Sort with A,B,C... first if present
    return sorted(chains, key=lambda c: (c not in string.ascii_uppercase, c))

def _read_structure_for_viewer(structure_path):
    """
    Read a structure file for py3Dmol and return:
        (structure_text, format_string)

    Supports .pdb, .cif, .mmcif
    """
    structure_path = str(structure_path)
    suffix = Path(structure_path).suffix.lower()

    if suffix == ".pdb":
        fmt = "pdb"
    elif suffix in [".cif", ".mmcif"]:
        fmt = "cif"
    else:
        raise ValueError(f"Unsupported structure format: {suffix}")

    with open(structure_path, "r") as f:
        structure_data = f.read()

    return structure_data, fmt


def _chains_present_in_structure(structure_text: str, fmt: str):
    """
    Return chain IDs when easily available from PDB text.
    For CIF/mmCIF, return None and let callers fall back to whole-model styling.
    """
    if fmt != "pdb":
        return None
    return _chains_present_in_pdb(structure_text)


def _load_mdtraj_any(structure_path):
    """
    Load PDB or CIF/mmCIF into MDTraj.
    For CIF/mmCIF, route through PDBFixer -> temporary PDB if needed.
    Returns (traj, tmp_path_to_cleanup_or_None).
    """
    structure_path = str(structure_path)
    suffix = Path(structure_path).suffix.lower()

    if suffix == ".pdb":
        return md.load(structure_path), None

    if suffix in [".cif", ".mmcif"]:
        fixer = PDBFixer(filename=structure_path)
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens()

        with NamedTemporaryFile(mode="w", suffix=".pdb", delete=False) as tmp:
            PDBFile.writeFile(fixer.topology, fixer.positions, tmp)
            tmp_path = tmp.name

        traj = md.load(tmp_path)
        return traj, tmp_path

    raise ValueError(f"Unsupported structure format: {suffix}")


def _chain_id_for_mdtraj_chain(chain, chain_key="id"):
    if chain_key == "index":
        return chain.index
    elif chain_key == "id":
        cid = getattr(chain, "chain_id", None)
        if cid is None or str(cid).strip() == "":
            letters = string.ascii_uppercase
            if chain.index < len(letters):
                return letters[chain.index]
            return str(chain.index)
        return str(cid).strip()
    else:
        raise ValueError("chain_key must be 'index' or 'id'")


def _ca_map_mdtraj(traj, chain_key="id"):
    """
    Return:
        ca_map: dict[(chain_id, resSeq, icode)] -> xyz(3,)
    Coordinates are returned in Å.
    """
    top = traj.topology
    xyz = traj.xyz[0] * 10.0  # MDTraj nm -> Å

    ca_map = {}

    for atom in top.atoms:
        if atom.name != "CA":
            continue

        res = atom.residue
        ch = _chain_id_for_mdtraj_chain(res.chain, chain_key=chain_key)
        resseq = int(res.resSeq)
        icode = getattr(res, "insertion_code", "") or ""
        key = (ch, resseq, icode)

        ca_map[key] = xyz[atom.index]

    return ca_map


def _kabsch_transform(P, Q):
    """
    Find rotation/translation that maps Q onto P.
    P, Q : (N,3)
    Returns R, t such that Q @ R + t matches P
    """
    Pc = P.mean(axis=0)
    Qc = Q.mean(axis=0)

    P0 = P - Pc
    Q0 = Q - Qc

    C = Q0.T @ P0
    V, S, Wt = np.linalg.svd(C)
    d = np.sign(np.linalg.det(V @ Wt))
    R = V @ np.diag([1.0, 1.0, d]) @ Wt
    t = Pc - Qc @ R
    return R, t


def _traj_to_pdb_string_with_transform(traj, R=None, t=None):
    """
    Apply optional rigid transform to an MDTraj trajectory and return PDB text.
    """
    xyz = traj.xyz.copy()  # nm
    if R is not None and t is not None:
        xyzA = xyz[0] * 10.0
        xyzA = xyzA @ R + t
        xyz[0] = xyzA / 10.0

    tmp = NamedTemporaryFile(mode="w", suffix=".pdb", delete=False)
    tmp_path = tmp.name
    tmp.close()

    try:
        traj2 = traj[:]
        traj2.xyz = xyz
        traj2.save_pdb(tmp_path)
        with open(tmp_path, "r") as f:
            pdb_text = f.read()
    finally:
        try:
            os.remove(tmp_path)
        except OSError:
            pass

    return pdb_text


def superimpose_structure_files_by_ca(ref_path, mob_path, chain_key="id"):
    """
    Superimpose mobile structure onto reference using matched Cα atoms.

    Supports PDB and CIF/mmCIF inputs.

    Returns
    -------
    aligned_mob_pdb_str : str
        Mobile structure transformed and written as PDB text
    rmsd : float
        RMSD over matched Cα atoms in Å
    n_matched : int
        Number of matched Cα atoms
    """
    ref_traj, ref_tmp = _load_mdtraj_any(ref_path)
    mob_traj, mob_tmp = _load_mdtraj_any(mob_path)

    try:
        ref_ca = _ca_map_mdtraj(ref_traj, chain_key=chain_key)
        mob_ca = _ca_map_mdtraj(mob_traj, chain_key=chain_key)

        common_keys = sorted(set(ref_ca.keys()) & set(mob_ca.keys()))
        if len(common_keys) < 3:
            raise ValueError(f"Not enough matched Cα atoms for superposition (matched={len(common_keys)}).")

        P = np.array([ref_ca[k] for k in common_keys], float)
        Q = np.array([mob_ca[k] for k in common_keys], float)

        R, t = _kabsch_transform(P, Q)
        Q_aln = Q @ R + t
        diff = P - Q_aln
        rmsd = float(np.sqrt(np.mean(np.sum(diff**2, axis=1))))

        aligned_mob_pdb_str = _traj_to_pdb_string_with_transform(mob_traj, R=R, t=t)
        return aligned_mob_pdb_str, rmsd, len(common_keys)

    finally:
        for tmp_path in [ref_tmp, mob_tmp]:
            if tmp_path is not None:
                try:
                    os.remove(tmp_path)
                except OSError:
                    pass

def _resolve_latest_prediction(directory, runNo, subNo=0,subRun=False):
    """
    Resolve the latest prediction for a given runNo and subNo.

    Priority:
    1. use *_end__AA.pdb / *_end__CA.pdb if present
    2. otherwise use the largest numeric step_<predNo>

    Returns
    -------
    pred_tag : str or int
        Either "end" or an integer predNo
    aa_path : str
    ca_path : str
    """
    if subRun:
        run_dir = os.path.join(directory, f"allAtomRun{runNo}")
    else:
        run_dir =directory

    # First check for explicit "end" files
    aa_end = os.path.join(run_dir, f"mol{runNo}_sub_{subNo}_end__AA.pdb")
    ca_end = os.path.join(run_dir, f"mol{runNo}_sub_{subNo}_end__CA.pdb")
    if os.path.exists(aa_end) and os.path.exists(ca_end):
        return "end", aa_end, ca_end

    # Otherwise scan numeric step files
    aa_pattern = os.path.join(run_dir, f"mol{runNo}_sub_{subNo}_step_*__AA.pdb")
    aa_files = glob.glob(aa_pattern)

    step_re = re.compile(
        rf"mol{re.escape(str(runNo))}_sub_{re.escape(str(subNo))}_step_(\d+)__AA\.pdb$"
    )

    step_nums = []
    for aa_path in aa_files:
        fname = os.path.basename(aa_path)
        m = step_re.match(fname)
        if not m:
            continue

        predNo = int(m.group(1))
        ca_path = os.path.join(run_dir, f"mol{runNo}_sub_{subNo}_step_{predNo}__CA.pdb")
        if os.path.exists(ca_path):
            step_nums.append(predNo)

    if not step_nums:
        raise FileNotFoundError(
            f"No matching prediction files found for runNo={runNo}, subNo={subNo} in {run_dir}"
        )

    predNo = max(step_nums)
    aa_path = os.path.join(run_dir, f"mol{runNo}_sub_{subNo}_step_{predNo}__AA.pdb")
    ca_path = os.path.join(run_dir, f"mol{runNo}_sub_{subNo}_step_{predNo}__CA.pdb")
    return predNo, aa_path, ca_path


def visualisePrediction(directory, runNo, predNo=None, subNo=0, subRun=False):
    if not HAS_PY3DMOL:
        _warn_missing_py3dmol()
        return None

    view = py3Dmol.view(width=800, height=600)

    if subRun:
        run_dir = os.path.join(directory, f"allAtomRun{runNo}")
    else:
        run_dir = directory

    # Resolve file names
    if predNo is None:
        pred_tag, aa_path, ca_path = _resolve_latest_prediction(directory, runNo, subNo=subNo)
        print(f"Using latest prediction: {pred_tag}")
    else:
        if predNo == "end":
            aa_fname = f"mol{runNo}_sub_{subNo}_end__AA.pdb"
            ca_fname = f"mol{runNo}_sub_{subNo}_end__CA.pdb"
        else:
            aa_fname = f"mol{runNo}_sub_{subNo}_step_{predNo}__AA.pdb"
            ca_fname = f"mol{runNo}_sub_{subNo}_step_{predNo}__CA.pdb"

        aa_path = os.path.join(run_dir, aa_fname)
        ca_path = os.path.join(run_dir, ca_fname)

        if not os.path.exists(aa_path):
            raise FileNotFoundError(f"AA file not found: {aa_path}")
        if not os.path.exists(ca_path):
            raise FileNotFoundError(f"CA file not found: {ca_path}")

    # ---- load AA model ----
    aa_data, aa_fmt = _read_structure_for_viewer(aa_path)
    view.addModel(aa_data, aa_fmt)   # model 0
    aa_chains = _chains_present_in_structure(aa_data, aa_fmt)

    # ---- load CA model ----
    ca_data, ca_fmt = _read_structure_for_viewer(ca_path)
    view.addModel(ca_data, ca_fmt)   # model 1
    ca_chains = _chains_present_in_structure(ca_data, ca_fmt)

    palette = ["blue", "green", "red", "yellow", "cyan", "magenta",
               "orange", "purple", "lime", "gray"]

    # Style AA model
    if aa_chains:
        for i, ch in enumerate(aa_chains):
            color = palette[i % len(palette)]
            view.setStyle({"model": 0, "chain": ch}, {"cartoon": {"color": color}})
    else:
        view.setStyle({"model": 0}, {"cartoon": {"color": "lightgray"}})

    # Style CA model
    if ca_chains:
        for i, ch in enumerate(ca_chains):
            color = palette[i % len(palette)]
            view.setStyle({"model": 1, "chain": ch}, {"sphere": {"color": color, "opacity": 0.5}})
    else:
        view.setStyle({"model": 1}, {"sphere": {"color": "red", "opacity": 0.5}})

    view.zoomTo()
    view.show()
    return view
    
def _structure_from_pdb_string(pdb_str, struct_id="X"):
    parser = PDBParser(QUIET=True)
    return parser.get_structure(struct_id, io.StringIO(pdb_str))

def _pdb_string_from_structure(structure):
    out = io.StringIO()
    io_writer = PDBIO()
    io_writer.set_structure(structure)
    io_writer.save(out)
    return out.getvalue()

def _ca_map(structure):
    """
    Map (chain_id, resseq, icode) -> CA atom
    """
    model = next(structure.get_models())
    ca = {}
    for chain in model:
        ch = chain.id.strip()
        for res in chain:
            # skip hetero/water etc if you want: res.id[0] != ' ' means hetero
            if "CA" not in res:
                continue
            hetflag, resseq, icode = res.id
            icode = (icode or "").strip()
            ca[(ch, int(resseq), icode)] = res["CA"]
    return ca

def superimpose_pdb_strings_by_ca(ref_pdb_str, mob_pdb_str):
    """
    Returns (aligned_mob_pdb_str, rmsd, n_matched)
    mob is superimposed onto ref using matched CA atoms.
    """
    ref = _structure_from_pdb_string(ref_pdb_str, "REF")
    mob = _structure_from_pdb_string(mob_pdb_str, "MOB")

    ref_ca = _ca_map(ref)
    mob_ca = _ca_map(mob)

    common_keys = sorted(set(ref_ca.keys()) & set(mob_ca.keys()))
    if len(common_keys) < 3:
        raise ValueError(f"Not enough matched Cα atoms for superposition (matched={len(common_keys)}).")

    ref_atoms = [ref_ca[k] for k in common_keys]
    mob_atoms = [mob_ca[k] for k in common_keys]

    sup = Superimposer()
    sup.set_atoms(ref_atoms, mob_atoms)

    # Apply transform to *all atoms* in the mobile structure
    mob_all_atoms = list(mob.get_atoms())
    sup.apply(mob_all_atoms)

    aligned_mob_str = _pdb_string_from_structure(mob)
    return aligned_mob_str, float(sup.rms), len(common_keys)

def visualisePredictionComparison(directory, runNo1, runNo2, predNo1, predNo2, subNo1, subNo2, do_superpose=True):
    if not HAS_PY3DMOL:
        _warn_missing_py3dmol()
        return None

    view = py3Dmol.view(width=800, height=600)

    aa_fname = f"mol{runNo1}_sub_{subNo1}_step_{predNo1}__AA.pdb"
    ca_fname = f"mol{runNo2}_sub_{subNo2}_step_{predNo2}__AA.pdb"

    aa_path = os.path.join(directory, "allAtomRun" + str(runNo1), aa_fname)
    ca_path = os.path.join(directory, "allAtomRun" + str(runNo2), ca_fname)

    aa_data, aa_fmt = _read_structure_for_viewer(aa_path)

    if do_superpose:
        ca_data_to_show, rmsd, nmatch = superimpose_structure_files_by_ca(aa_path, ca_path)
        ca_fmt = "pdb"
        print(f"Superposed model 1 onto model 0 using {nmatch} matched Cα atoms. RMSD = {rmsd:.3f} Å")
    else:
        ca_data_to_show, ca_fmt = _read_structure_for_viewer(ca_path)

    view.addModel(aa_data, aa_fmt)
    view.addModel(ca_data_to_show, ca_fmt)

    aa_chains = _chains_present_in_pdb(aa_data) if aa_fmt == "pdb" else None
    ca_chains = _chains_present_in_pdb(ca_data_to_show) if ca_fmt == "pdb" else None

    palette = ["blue", "green", "red", "yellow", "cyan", "magenta", "orange", "purple", "lime", "gray"]

    if aa_chains:
        for i, ch in enumerate(aa_chains):
            color = palette[i % len(palette)]
            view.setStyle({"model": 0, "chain": ch}, {"cartoon": {"color": color}})
    else:
        view.setStyle({"model": 0}, {"cartoon": {"color": "lightgray"}})

    if ca_chains:
        for i, ch in enumerate(ca_chains):
            color = palette[i % len(palette)]
            view.setStyle({"model": 1, "chain": ch}, {"sphere": {"color": color, "opacity": 0.5}})
    else:
        view.setStyle({"model": 1}, {"sphere": {"color": "red", "opacity": 0.5}})

    view.zoomTo()
    view.show()
    return view

def _infer_length_unit_from_ca(coords):
    """
    coords: (N,3)
    Returns "nm" or "A" based on typical CA-CA spacing.
    """
    if coords.shape[0] < 5:
        # fallback: use scale of coords
        span = np.linalg.norm(coords.max(axis=0) - coords.min(axis=0))
        return "nm" if span < 50 else "A"  # rough fallback
    d = np.linalg.norm(coords[1:] - coords[:-1], axis=1)
    med = float(np.median(d[np.isfinite(d)]))
    # CA-CA ~0.38 nm or ~3.8 Å
    return "nm" if med < 1.0 else "A"

def read_ca_coords(pdb_or_cif, chain_key="index", force_angstrom=True, unit="auto"):
    """
    Read Cα coordinates from PDB or mmCIF using MDTraj/PDBFixer.

    MDTraj coordinates are handled in nm internally and converted to Å if requested.

    unit: "auto" | "nm" | "A"
      - "auto": infer from CA-CA spacing
      - "nm": treat coords as nm
      - "A":  treat coords as Å
    force_angstrom: if True, returns coords in Å.
    """
    traj, tmp_path = _load_mdtraj_any(pdb_or_cif)
    try:
        top = traj.topology
        xyz = traj.xyz[0]  # nm

        coords = []
        keys = []

        for atom in top.atoms:
            if atom.name != "CA":
                continue

            res = atom.residue
            chain = res.chain

            if chain_key == "index":
                ch = chain.index
            elif chain_key == "id":
                ch = _chain_id_for_mdtraj_chain(chain, chain_key="id")
            else:
                raise ValueError("chain_key must be 'index' or 'id'")

            icode = getattr(res, "insertion_code", "") or ""
            coords.append(xyz[atom.index])
            keys.append((ch, res.resSeq, icode))

        coords = np.array(coords, float)

        if unit == "auto":
            unit = _infer_length_unit_from_ca(coords)

        if force_angstrom:
            if unit == "nm":
                coords = coords * 10.0
            elif unit == "A":
                pass
            else:
                raise ValueError("unit must be 'auto', 'nm', or 'A'")

        return coords, keys
    finally:
        if tmp_path is not None:
            try:
                os.remove(tmp_path)
            except OSError:
                pass


def kabsch_align_Q_to_P(P, Q):
    # Returns Q_aligned in P frame
    Pc = P.mean(axis=0)
    Qc = Q.mean(axis=0)
    P0 = P - Pc
    Q0 = Q - Qc

    C = Q0.T @ P0
    V, S, Wt = np.linalg.svd(C)
    d = np.sign(np.linalg.det(V @ Wt))
    R = V @ np.diag([1.0, 1.0, d]) @ Wt

    return (Q0 @ R) + Pc
    
def best_resno_offset(keysP, keysQ, max_abs_offset=10):
    """
    Find integer offset o such that P resno r matches Q resno (r + o),
    maximizing overlap count, per chain+icode.
    """
    # Work per (chain, icode) so we don't cross chains or insertion-codes
    Pset = set(keysP)
    Qset = set(keysQ)

    best_o, best_n = 0, -1
    for o in range(-max_abs_offset, max_abs_offset + 1):
        n = 0
        for (ch, r, ic) in Pset:
            if (ch, r + o, ic) in Qset:
                n += 1
        if n > best_n:
            best_n, best_o = n, o
    return best_o, best_n

def compare_structures(pdb1, pdb2, max_abs_offset=10):
    P, keysP = read_ca_coords(pdb1)
    Q, keysQ = read_ca_coords(pdb2)

    # Detect best residue-number offset
    off, n_overlap = best_resno_offset(keysP, keysQ, max_abs_offset=max_abs_offset)
    if n_overlap <= 0:
        raise ValueError("No overlapping residues found (even after offset search).")

    # Build index maps for fast lookup (avoid list.index O(N^2))
    idxP = {k: i for i, k in enumerate(keysP)}
    idxQ = {k: i for i, k in enumerate(keysQ)}

    commonP = []
    for (ch, r, ic) in idxP.keys():
        kq = (ch, r + off, ic)
        if kq in idxQ:
            commonP.append((ch, r, ic))

    Pm = P[[idxP[k] for k in commonP]]
    Qm = Q[[idxQ[(k[0], k[1] + off, k[2])] for k in commonP]]
    N = len(Pm)

    Q_aln = kabsch_align_Q_to_P(Pm, Qm)
    diff = Pm - Q_aln
    per = np.linalg.norm(diff, axis=1)
    rmsd = np.sqrt((diff**2).sum() / N)

    d0 = 1.24 * (N - 15)**(1/3) - 1.8
    tm_score = (1.0/N) * np.sum(1.0 / (1.0 + (per/d0)**2))

    cutoffs = [1.0, 2.0, 4.0, 8.0]
    gdt = [np.mean(per <= c) for c in cutoffs]
    gdt_ts = 100 * np.mean(gdt)

    print(f"Detected residue-number offset (Q = P + off): off = {off}  (overlap={N})")
    print(f"Matched residues: {N}")
    print(f"Global Cα RMSD : {rmsd:.3f} Å")
    print(f"TM-score       : {tm_score:.3f}")
    print(f"GDT-TS         : {gdt_ts:.1f}% "
          f"(within 1/2/4/8 Å: " + ", ".join(f"{100*x:.1f}%" for x in gdt) + ")")

    plt.figure(figsize=(8,3))
    plt.plot(range(1, N+1), per, lw=1)
    plt.xlabel("Residue index (matched set)")
    plt.ylabel("Per-residue deviation (Å)")
    plt.title("Per-residue Cα deviation (after offset-corrected match)")
    plt.tight_layout()
    plt.show()

def compare_structures_vals(pdb1, pdb2, max_abs_offset=10):
    P, keysP = read_ca_coords(pdb1)
    Q, keysQ = read_ca_coords(pdb2)

    # Detect best residue-number offset
    off, n_overlap = best_resno_offset(keysP, keysQ, max_abs_offset=max_abs_offset)
    if n_overlap <= 0:
        raise ValueError("No overlapping residues found (even after offset search).")

    # Build index maps for fast lookup (avoid list.index O(N^2))
    idxP = {k: i for i, k in enumerate(keysP)}
    idxQ = {k: i for i, k in enumerate(keysQ)}

    commonP = []
    for (ch, r, ic) in idxP.keys():
        kq = (ch, r + off, ic)
        if kq in idxQ:
            commonP.append((ch, r, ic))

    Pm = P[[idxP[k] for k in commonP]]
    Qm = Q[[idxQ[(k[0], k[1] + off, k[2])] for k in commonP]]
    N = len(Pm)

    Q_aln = kabsch_align_Q_to_P(Pm, Qm)
    diff = Pm - Q_aln
    per = np.linalg.norm(diff, axis=1)
    rmsd = np.sqrt((diff**2).sum() / N)

    d0 = 1.24 * (N - 15)**(1/3) - 1.8
    tm_score = (1.0/N) * np.sum(1.0 / (1.0 + (per/d0)**2))

    cutoffs = [1.0, 2.0, 4.0, 8.0]
    gdt = [np.mean(per <= c) for c in cutoffs]
    gdt_ts = 100 * np.mean(gdt)
    return np.array([rmsd, tm_score, gdt_ts])



def read_ca_carbonara(dat_file):
    """
    Read Carbonara backbone coordinates from a plain text file
    with one xyz triplet per line.

    Returns
    -------
    coords : (N, 3) np.ndarray
    """
    coords = np.loadtxt(dat_file, dtype=float)

    if coords.ndim != 2 or coords.shape[1] != 3:
        raise ValueError(
            f"Expected an Nx3 coordinate file, got shape {coords.shape} from {dat_file}"
        )

    return coords


def compare_structures_vals_carbonara(aa_pdb, carbonara_dat):
    """
    Compare an AA prediction PDB against a Carbonara backbone coordinate file.

    Assumptions
    -----------
    - The AA PDB residue order is the reference order
    - The Carbonara coordinate file contains one CA-like backbone point per residue
      in the same order as the AA model
    - No residue-number offset search is needed

    Parameters
    ----------
    aa_pdb : str
        Path to all-atom prediction PDB.
    carbonara_dat : str
        Path to Carbonara backbone coordinates file.

    Returns
    -------
    np.ndarray
        [rmsd, tm_score, gdt_ts]
    """
    P, keysP = read_ca_coords(aa_pdb)      # AA model CA coordinates + keys
    Q = read_ca_carbonara(carbonara_dat)   # Carbonara backbone coordinates only

    if len(P) != len(Q):
        raise ValueError(
            f"Length mismatch: AA model has {len(P)} CA atoms but "
            f"Carbonara file has {len(Q)} coordinates."
        )

    N = len(P)

    # Align Carbonara coords onto AA coords
    Q_aln = kabsch_align_Q_to_P(P, Q)

    diff = P - Q_aln
    per = np.linalg.norm(diff, axis=1)
    rmsd = np.sqrt((diff**2).sum() / N)

    # safer d0 for short structures
    d0 = 1.24 * max(N - 15, 1)**(1/3) - 1.8
    d0 = max(d0, 0.5)

    tm_score = (1.0 / N) * np.sum(1.0 / (1.0 + (per / d0)**2))

    cutoffs = [1.0, 2.0, 4.0, 8.0]
    gdt = [np.mean(per <= c) for c in cutoffs]
    gdt_ts = 100.0 * np.mean(gdt)

    return np.array([rmsd, tm_score, gdt_ts])

def collect_allatom_end_pdbs(directory):
    """
    Find all mol<runNo>_sub_<i>_end__AA.pdb files inside allAtomRun<runNo> folders.

    Returns
    -------
    paths : list of str
        Full paths to matching PDB files.
    """
    paths = []

    # Find all directories matching allAtomRun*
    for run_dir in glob.glob(os.path.join(directory, "allAtomRun*")):
        if not os.path.isdir(run_dir):
            continue

        run_name = os.path.basename(run_dir)
        m = re.match(r"allAtomRun(\d+)$", run_name)
        if not m:
            continue

        runNo = m.group(1)

        pattern = os.path.join(
            run_dir,
            f"mol{runNo}_sub_*_end__AA.pdb"
        )

        paths.extend(sorted(glob.glob(pattern)))

    return paths

def collect_final_predictions(directory):
    """
    Find the final AA prediction for each (runNo, subNo) found directly in `directory`.

    Priority
    --------
    1. *_end__AA.pdb if present
    2. otherwise largest step_<N>

    Returns
    -------
    paths : list of str
        Full paths to selected AA PDB files, one per (runNo, subNo).
    """

    best = {}

    for fname in os.listdir(directory):

        if not fname.endswith("__AA.pdb"):
            continue

        parts = fname.split("_")

        # expected:
        # molX_sub_Y_step_N__AA.pdb
        # molX_sub_Y_end__AA.pdb

        if len(parts) < 5:
            continue
        if not parts[0].startswith("mol"):
            continue
        if parts[1] != "sub":
            continue

        try:
            runNo = int(parts[0][3:])   # from "molX"
            subNo = int(parts[2])       # from "..._sub_Y_..."
        except ValueError:
            continue

        key = (runNo, subNo)
        fullpath = os.path.join(directory, fname)

        if parts[3] == "end":
            best[key] = ("end", fullpath)

        elif parts[3] == "step":
            if len(parts) < 6:
                continue
            try:
                step = int(parts[4])
            except ValueError:
                continue

            if key in best and best[key][0] == "end":
                continue

            if key not in best or step > best[key][0]:
                best[key] = (step, fullpath)

    return [best[k][1] for k in sorted(best)]


def collect_allatom_end_pdb_sets(directory):
    """
    Return per-run sets of mol<runNo>_sub_<i>_end__AA.pdb files.

    Returns
    -------
    sets_per_run : list[list[str]]
        Each inner list is the files for one runNo, sorted by sub index.
        The outer list is sorted by runNo.
    """
    run_to_files = defaultdict(list)

    for run_dir in glob.glob(os.path.join(directory, "allAtomRun*")):
        if not os.path.isdir(run_dir):
            continue

        run_name = os.path.basename(run_dir)
        m = re.match(r"allAtomRun(\d+)$", run_name)
        if not m:
            continue

        runNo = m.group(1)
        pattern = os.path.join(run_dir, f"mol{runNo}_sub_*_end__AA.pdb")

        for path in glob.glob(pattern):
            base = os.path.basename(path)
            msub = re.match(rf"mol{runNo}_sub_(\d+)_end__AA\.pdb$", base)
            if not msub:
                continue
            sub_i = int(msub.group(1))
            run_to_files[int(runNo)].append((sub_i, path))

    # Sort runs and sort files within each run by sub_i
    sets_per_run = []
    for runNo in sorted(run_to_files.keys()):
        files = [p for sub_i, p in sorted(run_to_files[runNo], key=lambda t: t[0])]
        sets_per_run.append(files)

    return sets_per_run


def collect_pdb_sets_ordered(directory, min_run=None, max_run=None):
    from collections import defaultdict
    found = defaultdict(list)

    for run_dir in glob.glob(os.path.join(directory, "allAtomRun*")):
        if not os.path.isdir(run_dir):
            continue
        m = re.match(r"allAtomRun(\d+)$", os.path.basename(run_dir))
        if not m:
            continue
        runNo = int(m.group(1))

        pattern = os.path.join(run_dir, f"mol{runNo}_sub_*_end__AA.pdb")
        for path in glob.glob(pattern):
            msub = re.search(r"_sub_(\d+)_end__AA\.pdb$", path)
            if msub:
                found[runNo].append((int(msub.group(1)), path))

    if not found:
        return [], []

    if min_run is None:
        min_run = min(found)
    if max_run is None:
        max_run = max(found)

    runNos = list(range(min_run, max_run + 1))
    pdb_sets = []

    for r in runNos:
        if r in found:
            pdb_sets.append([p for _, p in sorted(found[r])])
        else:
            pdb_sets.append(None)

    return pdb_sets, runNos


def collect_end_scatter_mixtures_ordered(directory, min_run=None, max_run=None):
    """
    Return mixtures ordered by run number, with None for missing runs.

    Returns
    -------
    mixtures : list[list[float] | None]
        Index i corresponds to runNo = start + i
    runNos : list[int]
        The run numbers represented.
    """
    found = {}

    pattern = os.path.join(directory, "mol*_end_scatter.dat")
    for path in glob.glob(pattern):
        base = os.path.basename(path)
        m = re.match(r"mol(\d+)_end_scatter\.dat$", base)
        if not m:
            continue

        runNo = int(m.group(1))

        with open(path, "r") as f:
            lines = f.read().splitlines()

        for line in reversed(lines):
            s = line.strip()
            if s:
                found[runNo] = [float(x) for x in s.split()]
                break

    if not found:
        return [], []

    if min_run is None:
        min_run = min(found)
    if max_run is None:
        max_run = max(found)

    runNos = list(range(min_run, max_run + 1))
    mixtures = [found.get(r) for r in runNos]

    return mixtures, runNos


def plot_overlaid_histograms(
    data_list,
    metric,
    labels,
    xlabel,
    title,
    bins=20,
    figsize=(4.5, 3.5),
    colors=None,
    alpha=0.5,
    density=False,
    savepath=None,
    dpi=300,
):
    """
    Plot multiple histograms overlaid on the same axes.

    Parameters
    ----------
    data_list : list of lists
        Each element is a list of [ref_pdb, model_pdb, array([v0,v1,v2])].
    value_index : int
        Index into the metric array (0, 1, or 2).
    labels : list of str
        Legend labels for each dataset.
    xlabel : str
        X-axis label.
    title : str
        Figure title.
    bins : int or array
        Number of bins or explicit bin edges.
    figsize : tuple
        Figure size in inches.
    colors : list of str
        Matplotlib color specs (defaults to muted palette).
    alpha : float
        Transparency for fills.
    density : bool
        If True, plot probability densities instead of counts.
    savepath : str or None
        If given, save figure to this path (e.g. "figS3.pdf").
    dpi : int
        DPI for raster formats (ignored for PDF/SVG).
    """

    assert len(data_list) == len(labels), "labels must match data_list length"

    # Extract values
    values = [
        np.array([entry[metric] for entry in data], dtype=float)
        for data in data_list
    ]

    # Shared bins across all datasets
    all_vals = np.concatenate(values)
    bin_edges = (
        np.histogram_bin_edges(all_vals, bins=bins)
        if not hasattr(bins, "__len__")
        else bins
    )

    if colors is None:
        colors = plt.cm.Greys(np.linspace(0.3, 0.8, len(values)))

    fig, ax = plt.subplots(figsize=figsize)

    for vals, lab, col in zip(values, labels, colors):
        ax.hist(
            vals,
            bins=bin_edges,
            histtype="stepfilled",
            alpha=alpha,
            color=col,
            edgecolor="black",
            linewidth=0.8,
            density=density,
            label=lab,
        )

    ax.set_xlabel(xlabel, fontsize=11)
    ax.set_ylabel("Density" if density else "Count", fontsize=11)
    ax.set_title(title, fontsize=12)

    ax.tick_params(axis="both", labelsize=10)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.legend(frameon=False, fontsize=9)

    plt.tight_layout()

    if savepath is not None:
        plt.savefig(savepath, dpi=dpi, bbox_inches="tight")

    plt.show()


def radius_of_gyration(pdb_or_cif, atom_selection="all", mass_weighted=False):
    """
    Compute radius of gyration (Rg) from a structure file.

    Parameters
    ----------
    pdb_or_cif : str
        Path to PDB or mmCIF.
    atom_selection : {"all","CA","backbone"} or MDTraj selection string
        Which atoms to use.
    mass_weighted : bool
        If True, use atomic masses (approx) for weighting; else uniform weights.

    Returns
    -------
    rg_A : float
        Radius of gyration in Å.
    """
    traj = md.load(pdb_or_cif)           # traj.xyz in nm
    top = traj.topology
    xyz = traj.xyz[0]                   # (N,3) in nm

    if atom_selection == "all":
        idx = np.arange(top.n_atoms)
    elif atom_selection == "CA":
        idx = top.select("name CA")
    elif atom_selection == "backbone":
        idx = top.select("backbone")
    else:
        # allow arbitrary mdtraj selection language
        idx = top.select(atom_selection)

    X = xyz[idx]                        # (M,3) in nm

    if X.shape[0] == 0:
        raise ValueError("Selection contained zero atoms.")

    if mass_weighted:
        # Simple mass lookup by element symbol (MDTraj usually has atom.element)
        mass_map = {"H": 1.008, "C": 12.011, "N": 14.007, "O": 15.999, "S": 32.06, "P": 30.974}
        masses = []
        atoms = list(top.atoms)
        for ai in idx:
            el = atoms[ai].element.symbol if atoms[ai].element is not None else None
            masses.append(mass_map.get(el, 12.011))  # default ~carbon if unknown
        w = np.asarray(masses, float)
    else:
        w = np.ones(X.shape[0], float)

    wsum = w.sum()
    r_cm = (X * w[:, None]).sum(axis=0) / wsum
    rg2 = (w * np.sum((X - r_cm) ** 2, axis=1)).sum() / wsum

    rg_nm = np.sqrt(rg2)
    rg_A = rg_nm * 10.0                 # nm -> Å
    return float(rg_A)



# -----------------------------------------------------------------------------
# Prediction-record helpers for single-structure and mixture-aware analysis
# -----------------------------------------------------------------------------

def _is_pathlike_object(x):
    """True for strings/path objects, but not for lists/tuples/dicts."""
    return isinstance(x, (str, os.PathLike, Path))


def _is_path_sequence(x):
    """True for a non-empty list/tuple whose entries are all path-like."""
    return isinstance(x, (list, tuple)) and len(x) > 0 and all(_is_pathlike_object(v) for v in x)


_PRED_STEP_PDB_RE = re.compile(
    r"^mol(?P<run>\d+)_sub_(?P<sub>\d+)_step_(?P<step>\d+)(?:_xyz)?_+AA\.pdb$"
)
_PRED_END_PDB_RE = re.compile(
    r"^mol(?P<run>\d+)_sub_(?P<sub>\d+)_end(?:_xyz)?_+AA\.pdb$"
)


def _parse_prediction_pdb_metadata(path):
    """
    Parse both historical and current Carbonara AA PDB names.

    Supported examples
    ------------------
    mol7_sub_0_step_12__AA.pdb
    mol7_sub_0_step_12_xyz_AA.pdb
    mol7_sub_0_end__AA.pdb
    mol7_sub_0_end_xyz_AA.pdb
    """
    name = Path(path).name
    m = _PRED_STEP_PDB_RE.match(name)
    if m:
        run_no = int(m.group("run"))
        sub_no = int(m.group("sub"))
        step_no = int(m.group("step"))
        return {
            "runNo": run_no,
            "run_no": run_no,
            "subNo": sub_no,
            "sub_no": sub_no,
            "predTag": f"step_{step_no}",
            "pred_tag": f"step_{step_no}",
            "label": f"mol{run_no}_step_{step_no}",
        }

    m = _PRED_END_PDB_RE.match(name)
    if m:
        run_no = int(m.group("run"))
        sub_no = int(m.group("sub"))
        return {
            "runNo": run_no,
            "run_no": run_no,
            "subNo": sub_no,
            "sub_no": sub_no,
            "predTag": "end",
            "pred_tag": "end",
            "label": f"mol{run_no}_end",
        }

    return {
        "runNo": None,
        "run_no": None,
        "subNo": None,
        "sub_no": None,
        "predTag": None,
        "pred_tag": None,
        "label": Path(path).stem,
    }


def _infer_fitdata_dir_from_pdb_paths(pdb_paths):
    """
    Infer the fitdata directory from component paths such as
    .../fitdata/allAtomRun7/mol7_sub_0_step_12_xyz_AA.pdb.
    """
    for p in pdb_paths:
        p = Path(p)
        for parent in [p.parent] + list(p.parents):
            if re.fullmatch(r"allAtomRun\d+", parent.name):
                return parent.parent
    return None


def _match_record_by_pdb_names(candidates, pdb_paths, chi2=None):
    """Find a FoXS record whose component PDB basenames match ``pdb_paths``."""
    wanted_names = {Path(p).name for p in pdb_paths}
    wanted_resolved = set()
    for p in pdb_paths:
        try:
            wanted_resolved.add(Path(p).resolve())
        except Exception:
            pass

    for cand in candidates:
        cand_paths = [Path(p) for p in cand.get("pdb_paths", [])]
        have_names = {p.name for p in cand_paths}
        have_resolved = set()
        for p in cand_paths:
            try:
                have_resolved.add(p.resolve())
            except Exception:
                pass

        same_paths = bool(wanted_names) and wanted_names == have_names
        if wanted_resolved and have_resolved:
            same_paths = same_paths or (wanted_resolved == have_resolved)

        if not same_paths:
            continue

        if chi2 is not None:
            try:
                if abs(float(cand.get("chi2")) - float(chi2)) > 1e-8:
                    continue
            except Exception:
                continue
        return cand

    return None


def _try_recover_prediction_record_from_files(pdb_paths, fitdata_dir=None, chi2=None):
    """
    Recover a full record, including mixture weights and fit curve, from
    foxs_mixture_results.txt/foxs_results.txt when only PDB paths were supplied.
    """
    search_dirs = []
    if fitdata_dir is not None:
        search_dirs.append(Path(fitdata_dir))

    inferred = _infer_fitdata_dir_from_pdb_paths(pdb_paths)
    if inferred is not None:
        search_dirs.append(Path(inferred))

    seen = set()
    for fd in search_dirs:
        fd = Path(fd)
        key = str(fd.resolve()) if fd.exists() else str(fd)
        if key in seen:
            continue
        seen.add(key)
        try:
            candidates = read_foxs_prediction_records(fd, mode="auto", require_exists=False)
            hit = _match_record_by_pdb_names(candidates, pdb_paths, chi2=chi2)
            if hit is not None:
                return hit
        except Exception:
            pass

    return None


def _normalise_prediction_collection(predictions, fitdata_dir: str | Path | None = None):
    """
    Convert a mixed collection of records, legacy tuples, paths and mixture path
    lists into record dictionaries. ``None`` entries are skipped.

    Important compatibility rule: a top-level flat list of paths is treated as
    the historical input style, i.e. a list of single structures. A nested list
    of paths is treated as one mixture prediction.
    """
    if predictions is None:
        return []

    if isinstance(predictions, dict) or _is_pathlike_object(predictions):
        return [_normalise_prediction_record(predictions, fitdata_dir=fitdata_dir)]

    if isinstance(predictions, tuple):
        # Legacy form: (path_or_paths, chi2)
        if len(predictions) >= 2 and (_is_pathlike_object(predictions[0]) or _is_path_sequence(predictions[0])):
            return [_normalise_prediction_record(predictions, fitdata_dir=fitdata_dir)]

    if isinstance(predictions, (list, tuple)):
        if len(predictions) == 0:
            return []

        # Historical flat list of PDB files: [pdb1, pdb2, ...]
        if all(_is_pathlike_object(x) for x in predictions):
            return [_normalise_prediction_record(p, fitdata_dir=fitdata_dir) for p in predictions]

        records = []
        for item in predictions:
            if item is None:
                continue
            records.append(_normalise_prediction_record(item, fitdata_dir=fitdata_dir))
        return records

    return [_normalise_prediction_record(predictions, fitdata_dir=fitdata_dir)]


def flatten_prediction_pdbs(predictions, fitdata_dir: str | Path | None = None, keep_metadata: bool = False):
    """
    Flatten prediction records/mixtures into individual component PDBs.

    This is the appropriate representation for RMSD, TM-score and GDT analyses:
    mixtures are not averaged structurally; each component structure is compared
    one-by-one.
    """
    records = _normalise_prediction_collection(predictions, fitdata_dir=fitdata_dir)
    rows = []
    for pred_i, rec in enumerate(records):
        pdb_paths = [Path(p) for p in rec.get("pdb_paths", [])]
        weights = rec.get("weights")
        if weights is None or len(weights) != len(pdb_paths):
            weights = [None] * len(pdb_paths)

        for comp_i, (pdb, weight) in enumerate(zip(pdb_paths, weights)):
            meta = _parse_prediction_pdb_metadata(pdb)
            rows.append({
                "pdb": pdb,
                "prediction_i": pred_i,
                "component_i": comp_i,
                "weight": None if weight is None else float(weight),
                "chi2": rec.get("chi2"),
                "type": rec.get("type", "single"),
                "label": rec.get("label") or meta.get("label"),
                "runNo": rec.get("run_no", meta.get("runNo")),
                "run_no": rec.get("run_no", meta.get("run_no")),
                "subNo": meta.get("subNo"),
                "sub_no": meta.get("sub_no"),
                "predTag": meta.get("predTag"),
                "pred_tag": meta.get("pred_tag"),
                "record": rec,
            })

    if keep_metadata:
        return rows
    return [r["pdb"] for r in rows]


def _weights_for_prediction_record(rec, pred_i, n_components, mixtures=None, default_mixture_weights="error"):
    """
    Resolve weights for an Rg weighted average.

    Priority:
    1. explicit ``mixtures`` argument, if supplied;
    2. weights stored in the FoXS mixture record;
    3. [1] for a single structure;
    4. equal weights only if ``default_mixture_weights='equal'``.
    """
    weights = None

    if mixtures is not None:
        run_no = rec.get("run_no")
        pred_tag = rec.get("pred_tag") or rec.get("predTag")
        if pred_tag is None and rec.get("pdb_paths"):
            pred_tag = _parse_prediction_pdb_metadata(rec["pdb_paths"][0]).get("predTag")
        if run_no is None and rec.get("pdb_paths"):
            run_no = _parse_prediction_pdb_metadata(rec["pdb_paths"][0]).get("runNo")

        if isinstance(mixtures, dict):
            if (run_no, pred_tag) in mixtures:
                weights = mixtures[(run_no, pred_tag)]
            elif run_no in mixtures:
                weights = mixtures[run_no]
        else:
            try:
                weights = mixtures[pred_i]
            except Exception:
                weights = None

    if weights is None:
        weights = rec.get("weights")

    if weights is None or len(weights) != n_components:
        if n_components == 1:
            weights = [1.0]
        elif default_mixture_weights == "equal":
            weights = np.ones(n_components, dtype=float) / float(n_components)
        else:
            raise ValueError(
                "No valid mixture weights were available for an Rg weighted average. "
                "Use collect_good_prediction_files(..., return_records=True) or "
                "collect_best_prediction_per_run_closest_to_one(..., return_records=True), "
                "or pass fitdata_dir so weights can be recovered from foxs_mixture_results.txt. "
                "For a fallback only, set default_mixture_weights='equal'."
            )

    weights = np.asarray(weights, dtype=float)
    if len(weights) != n_components:
        raise ValueError(f"Expected {n_components} weights, got {len(weights)}")
    if np.any(~np.isfinite(weights)):
        raise ValueError("Mixture weights contain NaN or infinite values")
    if np.any(weights < 0):
        raise ValueError("Mixture weights contain negative values")
    s = float(weights.sum())
    if s <= 0:
        raise ValueError("Mixture weights sum to zero")
    return weights / s

def pairwise_structure_metrics(pdb_files, compare_func, fitdata_dir: str | Path | None = None):
    """
    Compute pairwise RMSD / TM / GDT for all unique structure pairs.

    Mixture-aware behaviour
    -----------------------
    ``pdb_files`` may now be a flat list of PDB paths, records returned by
    ``collect_*`` with ``return_records=True``, legacy ``(path(s), chi2)``
    tuples, or a list containing mixture component lists. Mixture entries are
    flattened and compared component-by-component.
    """
    components = flatten_prediction_pdbs(pdb_files, fitdata_dir=fitdata_dir, keep_metadata=True)
    n = len(components)
    n_pairs = n * (n - 1) // 2

    results = []

    with tqdm(total=n_pairs, desc="Pairwise comparisons", unit="pair") as pbar:
        for i in range(n):
            for j in range(i + 1, n):
                c1 = components[i]
                c2 = components[j]
                p1 = c1["pdb"]
                p2 = c2["pdb"]

                rmsd, tm, gdt = compare_func(p1, p2)

                results.append({
                    "i": i,
                    "j": j,
                    "pdb1": p1,
                    "pdb2": p2,
                    "prediction_i1": c1.get("prediction_i"),
                    "prediction_i2": c2.get("prediction_i"),
                    "component_i1": c1.get("component_i"),
                    "component_i2": c2.get("component_i"),
                    "weight1": c1.get("weight"),
                    "weight2": c2.get("weight"),
                    "chi2_1": c1.get("chi2"),
                    "chi2_2": c2.get("chi2"),
                    "label1": c1.get("label"),
                    "label2": c2.get("label"),
                    "rmsd": float(rmsd),
                    "tm": float(tm),
                    "gdt_ts": float(gdt),
                })

                pbar.update(1)

    return results


def structure_metrics_vs_reference(pdb_files, ref_pdb, compare_func, fitdata_dir: str | Path | None = None):
    """
    Compute RMSD / TM / GDT for each structure against ``ref_pdb``.

    Mixture entries are flattened and compared component-by-component; no
    structural weighted average is attempted.
    """
    components = flatten_prediction_pdbs(pdb_files, fitdata_dir=fitdata_dir, keep_metadata=True)
    results = []

    with tqdm(total=len(components), desc="Comparisons vs reference", unit="pdb") as pbar:
        for i, comp in enumerate(components):
            pdb = comp["pdb"]
            rmsd, tm, gdt = compare_func(pdb, ref_pdb)

            results.append({
                "i": i,
                "pdb": pdb,
                "ref_pdb": ref_pdb,
                "prediction_i": comp.get("prediction_i"),
                "component_i": comp.get("component_i"),
                "weight": comp.get("weight"),
                "chi2": comp.get("chi2"),
                "label": comp.get("label"),
                "rmsd": float(rmsd),
                "tm": float(tm),
                "gdt_ts": float(gdt),
            })

            pbar.update(1)

    return results


def structure_metrics_vs_carbonara(pdb_files, carbonara_dir, compare_func, fitdata_dir: str | Path | None = None):
    """
    Compare each component PDB against Carbonara coordinate files.

    Mixture entries are flattened and compared component-by-component.
    """

    coord_files = sorted(glob.glob(os.path.join(carbonara_dir, "coordinates*.dat")))

    if not coord_files:
        raise ValueError(f"No coordinates*.dat files found in {carbonara_dir}")

    components = flatten_prediction_pdbs(pdb_files, fitdata_dir=fitdata_dir, keep_metadata=True)
    total = len(components) * len(coord_files)
    results = []

    with tqdm(total=total, desc="Comparisons vs Carbonara", unit="pair") as pbar:

        for dat in coord_files:
            for i, comp in enumerate(components):
                pdb = comp["pdb"]
                rmsd, tm, gdt = compare_func(pdb, dat)

                results.append({
                    "i": i,
                    "pdb": pdb,
                    "carbonara_dat": dat,
                    "prediction_i": comp.get("prediction_i"),
                    "component_i": comp.get("component_i"),
                    "weight": comp.get("weight"),
                    "chi2": comp.get("chi2"),
                    "label": comp.get("label"),
                    "rmsd": float(rmsd),
                    "tm": float(tm),
                    "gdt_ts": float(gdt),
                })

                pbar.update(1)

    return results


    
import os
import multiprocessing as mp
from tqdm import tqdm

def _pair_indices(n):
    for i in range(n):
        for j in range(i + 1, n):
            yield i, j

def _worker_compare_pair(args):
    """
    Worker must be top-level for multiprocessing pickling.
    args = (i, j, pdb_files)
    """
    i, j, pdb_files = args
    p1 = pdb_files[i]
    p2 = pdb_files[j]

    # Call the module-level function directly (most robust for pickling)
    vals = compare_structures_vals(p1, p2)  # returns array-like [rmsd, tm, gdt_ts]
    return i, j, float(vals[0]), float(vals[1]), float(vals[2])


def pairwise_structure_metrics_mp(pdb_files, nprocs=None, chunksize=20, start_method=None,
                                  fitdata_dir: str | Path | None = None):
    """
    Parallel pairwise RMSD/TM/GDT over all unique component pairs.

    Mixture entries are flattened first. The returned rows include prediction
    and component indices so the component comparisons can be traced back to
    their mixture records.
    """
    components = flatten_prediction_pdbs(pdb_files, fitdata_dir=fitdata_dir, keep_metadata=True)
    pdb_flat = [str(c["pdb"]) for c in components]

    n = len(pdb_flat)
    n_pairs = n * (n - 1) // 2
    if n_pairs == 0:
        return []

    if nprocs is None:
        nprocs = os.cpu_count() or 1

    if start_method is None:
        # Best default on Linux for this kind of numeric work
        start_method = "fork" if os.name == "posix" else "spawn"

    ctx = mp.get_context(start_method)
    tasks = ((i, j, pdb_flat) for i, j in _pair_indices(n))

    results = []
    with ctx.Pool(processes=nprocs) as pool:
        it = pool.imap_unordered(_worker_compare_pair, tasks, chunksize=chunksize)
        for i, j, rmsd, tm, gdt in tqdm(it, total=n_pairs, desc="Pairwise comparisons", unit="pair"):
            c1 = components[i]
            c2 = components[j]
            results.append({
                "i": i,
                "j": j,
                "pdb1": c1["pdb"],
                "pdb2": c2["pdb"],
                "prediction_i1": c1.get("prediction_i"),
                "prediction_i2": c2.get("prediction_i"),
                "component_i1": c1.get("component_i"),
                "component_i2": c2.get("component_i"),
                "weight1": c1.get("weight"),
                "weight2": c2.get("weight"),
                "chi2_1": c1.get("chi2"),
                "chi2_2": c2.get("chi2"),
                "label1": c1.get("label"),
                "label2": c2.get("label"),
                "rmsd": rmsd,
                "tm": tm,
                "gdt_ts": gdt,
            })

    results.sort(key=lambda d: (d["i"], d["j"]))
    return results


def canonical_mixture_key(mix, decimals=2):
    """
    Convert a mixture list into an order-invariant, rounded tuple key.

    Examples:
      [0.8, 0.2]     -> (0.2, 0.8)
      [0.79, 0.21]   -> (0.21, 0.79)  (with decimals=2)
      [0.33,0.33,0.34] -> (0.33,0.33,0.34)
    """
    return tuple(sorted(round(float(x), decimals) for x in mix))


def bucket_rg_differences_by_mixture(
    mixtures,
    rg_sets,
    key_decimals=2,
    diff_mode="pair_abs",
    skip_none=True,
):
    """
    Bucket Rg differences by unordered mixture composition.

    mixtures : list[list[float] | None]
    rg_sets  : list[list[float] | None]
    """
    if len(mixtures) != len(rg_sets):
        raise ValueError("mixtures and rg_sets must have same length")

    buckets = defaultdict(list)

    for mix, rgs in zip(mixtures, rg_sets):
        if mix is None or rgs is None:
            if skip_none:
                continue
            else:
                raise ValueError("Found None entry")

        if len(mix) != len(rgs):
            raise ValueError(f"Mixture length {len(mix)} != Rg-set length {len(rgs)}")

        # --- canonical unordered mixture key ---
        key = canonical_mixture_key(mix, decimals=key_decimals)

        # --- reduce Rg set ---
        if diff_mode == "pair_abs":
            if len(rgs) != 2:
                raise ValueError("pair_abs requires exactly 2 Rg values")
            value = abs(float(rgs[0]) - float(rgs[1]))
        elif diff_mode == "range":
            rr = list(map(float, rgs))
            value = max(rr) - min(rr)
        elif diff_mode == "all":
            value = list(map(float, rgs))
        else:
            raise ValueError("diff_mode must be 'pair_abs', 'range', or 'all'")

        buckets[key].append(value)

    return dict(buckets)


def calc_rg_distribution(
    pdb_files,
    rg_func,
    weighted=False,
    mixtures=None,
    keep_components=True,
    fitdata_dir: str | Path | None = None,
    default_mixture_weights="error",
):
    """
    Calculate radius-of-gyration values from single predictions or mixtures.

    Parameters
    ----------
    pdb_files : sequence
        May be a flat list of PDB paths, records returned by
        ``collect_good_prediction_files(..., return_records=True)`` or
        ``collect_best_prediction_per_run_closest_to_one(..., return_records=True)``,
        legacy ``(path(s), chi2)`` tuples, or nested lists of component PDBs.
    rg_func : callable
        Function like ``radius_of_gyration(pdb_path) -> float``.
    weighted : bool
        If False, return one Rg per component structure. This is the direct
        analogue of RMSD/TM one-by-one analysis.
        If True, return one Rg per prediction: for mixtures this is
        ``sum_i weight_i * Rg_i`` using the approximate MultiFoXS weights.
    mixtures : dict or list, optional
        Optional explicit weights. This overrides weights stored in records.
        Supported forms are dict keyed by ``(runNo, predTag)`` or ``runNo``, or a
        list aligned with the sorted/normalised prediction list.
    keep_components : bool
        In weighted mode, include component PDBs, component Rg values and weights
        in each returned record.
    fitdata_dir : str or Path, optional
        Used to recover weights/fit metadata when only raw component PDB paths
        were supplied.
    default_mixture_weights : {"error", "equal"}
        What to do if a multi-component prediction has no recoverable weights.
        The default is to raise an error rather than silently report an
        unphysical average.

    Returns
    -------
    list of dict
        Unweighted mode: one row per component, with key ``"rg"``.
        Weighted mode: one row per prediction, with key ``"rg"`` equal to the
        weighted average; component values are in ``"rg_components"`` when
        ``keep_components=True``.
    """
    records = _normalise_prediction_collection(pdb_files, fitdata_dir=fitdata_dir)
    results = []

    # --------------------------------------------------
    # Unweighted mode: one output per component PDB
    # --------------------------------------------------
    if not weighted:
        components = flatten_prediction_pdbs(records, fitdata_dir=fitdata_dir, keep_metadata=True)
        with tqdm(total=len(components), desc="Calculating Rg", unit="pdb") as pbar:
            for i, comp in enumerate(components):
                pdb = comp["pdb"]
                rg = float(rg_func(str(pdb)))

                results.append({
                    "i": i,
                    "prediction_i": comp.get("prediction_i"),
                    "component_i": comp.get("component_i"),
                    "pdb": pdb,
                    "rg": rg,
                    "weight": comp.get("weight"),
                    "chi2": comp.get("chi2"),
                    "label": comp.get("label"),
                    "type": comp.get("type"),
                    "runNo": comp.get("runNo"),
                    "subNo": comp.get("subNo"),
                    "predTag": comp.get("predTag"),
                })

                pbar.update(1)

        return results

    # --------------------------------------------------
    # Weighted mode: one output per prediction/mixture
    # --------------------------------------------------
    with tqdm(total=len(records), desc="Calculating weighted Rg", unit="pred") as pbar:
        for i, rec in enumerate(records):
            pdb_group = [Path(p) for p in rec.get("pdb_paths", [])]
            if not pdb_group:
                pbar.update(1)
                continue

            rg_components = [float(rg_func(str(pdb))) for pdb in pdb_group]
            weights = _weights_for_prediction_record(
                rec,
                pred_i=i,
                n_components=len(rg_components),
                mixtures=mixtures,
                default_mixture_weights=default_mixture_weights,
            )

            weighted_rg = float(np.sum(weights * np.asarray(rg_components, dtype=float)))
            meta = _parse_prediction_pdb_metadata(pdb_group[0])

            rec_out = {
                "i": i,
                "runNo": rec.get("run_no", meta.get("runNo")),
                "run_no": rec.get("run_no", meta.get("run_no")),
                "predTag": meta.get("predTag"),
                "pred_tag": meta.get("pred_tag"),
                "label": rec.get("label") or meta.get("label"),
                "type": rec.get("type", "single"),
                "chi2": rec.get("chi2"),
                "rg": weighted_rg,
            }

            if keep_components:
                rec_out.update({
                    "pdbs": pdb_group,
                    "pdb_paths": pdb_group,
                    "rg_components": rg_components,
                    "weights": list(map(float, weights)),
                })

            results.append(rec_out)
            pbar.update(1)

    return results


#########################################################
#
#   For the live update fornt end this will grab all the "good" predictions (user specified chi squared)
#
#########################################################
    

def _sort_allatom_run_dir(run_dir: Path) -> int:
    """Sort allAtomRun<N> directories by N, with unknown names at the end."""
    m = re.fullmatch(r"allAtomRun(\d+)", Path(run_dir).name)
    return int(m.group(1)) if m else 10**12


_SINGLE_FOXS_LINE_NUM_RE = re.compile(r"^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$")
_MIX_CHI_RE = re.compile(r"\bchi2=([0-9.eE+-]+)")
_MIX_SCALE_RE = re.compile(r"\bscale=([0-9.eE+-]+)")
_MIX_WEIGHTS_RE = re.compile(r"\bweights=([^\s]+)")
_MIX_PDBS_RE = re.compile(r"\bpdbs=([^\s]+)")
_MIX_PROFILES_RE = re.compile(r"\bprofiles=([^\s]+)")
_MIX_C1_RE = re.compile(r"\bc1=([0-9.eE+-]+)")
_MIX_C2_RE = re.compile(r"\bc2=([0-9.eE+-]+)")
_MIX_BEST_COMPONENT_CHI2_RE = re.compile(r"\bbest_component_chi2=([0-9.eE+-]+)")
_MIX_LABEL_STEP_RE = re.compile(r"^mol(\d+)_step_(\d+)$")
_MIX_LABEL_END_RE = re.compile(r"^mol(\d+)_end$")
_MIX_LABEL_INITIAL_RE = re.compile(r"^mol(\d+)_initial$")


def _resolve_recorded_path(path_text: str, cwd: Path | None = None, require_exists: bool = True) -> Path | None:
    """
    Resolve a path recorded in a FoXS summary file.

    Older summary files may contain absolute paths from another checkout/session.
    If the path contains a ``carbonara_runs`` component, rebuild it relative to
    the current working directory as a fallback.
    """
    cwd = Path.cwd() if cwd is None else Path(cwd)
    p = Path(path_text)

    candidates = []
    candidates.append(p)
    if not p.is_absolute():
        candidates.append(cwd / p)

    try:
        idx = p.parts.index("carbonara_runs")
        candidates.append(cwd / Path(*p.parts[idx:]))
    except ValueError:
        pass

    seen = set()
    for cand in candidates:
        cand = cand.resolve() if cand.exists() else cand
        key = str(cand)
        if key in seen:
            continue
        seen.add(key)
        if cand.exists():
            return cand

    # If existence is not required, return the most portable candidate when possible.
    if not require_exists:
        try:
            idx = p.parts.index("carbonara_runs")
            return cwd / Path(*p.parts[idx:])
        except ValueError:
            return p

    return None


def _parse_run_no_from_allatom_dir(run_dir: Path) -> int | None:
    m = re.fullmatch(r"allAtomRun(\d+)", Path(run_dir).name)
    return int(m.group(1)) if m else None


def _sub_sort_key_from_path(path: Path):
    m = re.search(r"_sub_(\d+)_", Path(path).name)
    return int(m.group(1)) if m else 10**9


def _infer_mixture_pdbs_from_label(run_dir: Path, label: str) -> list[Path]:
    """
    Infer component AA PDBs for a mixture label such as ``mol7_step_12``.

    Supports both the current watcher naming convention
    ``mol7_sub_0_step_12_xyz_AA.pdb`` and older double-underscore names such as
    ``mol7_sub_0_step_12__AA.pdb``.
    """
    label = str(label).strip()
    patterns = []

    m = _MIX_LABEL_STEP_RE.match(label)
    if m:
        run_no, step = int(m.group(1)), int(m.group(2))
        patterns.extend([
            f"mol{run_no}_sub_*_step_{step}_xyz_AA.pdb",
            f"mol{run_no}_sub_*_step_{step}__AA.pdb",
            f"mol{run_no}_sub_*_step_{step}_*_AA.pdb",
        ])

    m = _MIX_LABEL_END_RE.match(label)
    if m:
        run_no = int(m.group(1))
        patterns.extend([
            f"mol{run_no}_sub_*_end_xyz_AA.pdb",
            f"mol{run_no}_sub_*_end__AA.pdb",
            f"mol{run_no}_sub_*_end*_AA.pdb",
        ])

    m = _MIX_LABEL_INITIAL_RE.match(label)
    if m:
        run_no = int(m.group(1))
        patterns.extend([
            f"mol{run_no}_sub_*_initial_xyz_AA.pdb",
            f"mol{run_no}_sub_*_initial__AA.pdb",
            f"mol{run_no}_sub_*_initial*_AA.pdb",
        ])

    found = []
    seen = set()
    for pat in patterns:
        for p in run_dir.glob(pat):
            if p not in seen:
                found.append(p)
                seen.add(p)

    return sorted(found, key=_sub_sort_key_from_path)


def _parse_float_list_csv(text: str | None) -> list[float] | None:
    if text is None:
        return None
    text = text.strip()
    if not text:
        return []
    vals = []
    for item in text.split(','):
        item = item.strip()
        if not item:
            continue
        vals.append(float(item))
    return vals


def _parse_pdb_list_csv(text: str | None, cwd: Path, require_exists: bool) -> list[Path]:
    if text is None:
        return []
    paths = []
    for item in text.split(','):
        item = item.strip()
        if not item:
            continue
        p = _resolve_recorded_path(item, cwd=cwd, require_exists=require_exists)
        if p is not None:
            paths.append(p)
    return sorted(paths, key=_sub_sort_key_from_path)


def _read_single_foxs_records(run_dir: Path, require_exists: bool = True, cwd: Path | None = None) -> list[dict]:
    """Read ordinary single-structure FoXS records from ``foxs_results.txt``."""
    cwd = Path.cwd() if cwd is None else Path(cwd)
    foxs_file = run_dir / "foxs_results.txt"
    if not foxs_file.exists():
        return []

    run_no = _parse_run_no_from_allatom_dir(run_dir)
    records = []

    for line in foxs_file.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        parts = line.split()
        if len(parts) < 2:
            continue

        pdb_path_str, chi_str = parts[0], parts[1]
        if chi_str.upper() == "ERROR":
            continue
        if not _SINGLE_FOXS_LINE_NUM_RE.match(chi_str):
            continue

        try:
            chi2 = float(chi_str)
        except ValueError:
            continue

        pdb_path = _resolve_recorded_path(pdb_path_str, cwd=cwd, require_exists=require_exists)
        if pdb_path is None:
            continue

        records.append({
            "type": "single",
            "run_no": run_no,
            "label": Path(pdb_path).stem,
            "chi2": chi2,
            "pdb_path": pdb_path,
            "pdb_paths": [pdb_path],
            "weights": [1.0],
            "scale": None,
            "fit_file": None,
            "summary_file": foxs_file,
            "line": line,
        })

    return records


def _read_mixture_foxs_records(run_dir: Path, require_exists: bool = True, cwd: Path | None = None) -> list[dict]:
    """
    Read approximate MultiFoXS-style records from ``foxs_mixture_results.txt``.

    Expected watcher line format:
        mol7_step_12 chi2=<...> scale=<...> weights=w0,w1,... pdbs=p0,p1,...
    """
    cwd = Path.cwd() if cwd is None else Path(cwd)
    mix_file = run_dir / "foxs_mixture_results.txt"
    if not mix_file.exists():
        return []

    records = []
    fallback_run_no = _parse_run_no_from_allatom_dir(run_dir)

    for line in mix_file.read_text().splitlines():
        line = line.strip()
        if not line:
            continue

        parts = line.split(maxsplit=1)
        if not parts:
            continue
        label = parts[0]

        m_chi = _MIX_CHI_RE.search(line)
        if not m_chi:
            continue
        try:
            chi2 = float(m_chi.group(1))
        except ValueError:
            continue

        m_scale = _MIX_SCALE_RE.search(line)
        scale = None
        if m_scale:
            try:
                scale = float(m_scale.group(1))
            except ValueError:
                scale = None

        m_c1 = _MIX_C1_RE.search(line)
        c1 = None
        if m_c1:
            try:
                c1 = float(m_c1.group(1))
            except ValueError:
                c1 = None

        m_c2 = _MIX_C2_RE.search(line)
        c2 = None
        if m_c2:
            try:
                c2 = float(m_c2.group(1))
            except ValueError:
                c2 = None

        m_best = _MIX_BEST_COMPONENT_CHI2_RE.search(line)
        best_component_chi2 = None
        if m_best:
            try:
                best_component_chi2 = float(m_best.group(1))
            except ValueError:
                best_component_chi2 = None

        m_weights = _MIX_WEIGHTS_RE.search(line)
        weights = _parse_float_list_csv(m_weights.group(1) if m_weights else None)

        # Important: in current partial-profile mixture summaries the line has
        #     ... pdbs=p0,p1,p2 profiles=p0.dat,p1.dat,p2.dat
        # so pdbs must stop at the next whitespace, not run to end-of-line.
        m_pdbs = _MIX_PDBS_RE.search(line)
        pdbs = _parse_pdb_list_csv(m_pdbs.group(1) if m_pdbs else None, cwd=cwd, require_exists=require_exists)
        # Guard against old parser artefacts or malformed lines: component
        # structures are AA PDBs, while FoXS partial profiles end in .pdb.dat.
        pdbs = [p for p in pdbs if Path(p).suffix.lower() == ".pdb"]

        m_profiles = _MIX_PROFILES_RE.search(line)
        profiles = _parse_pdb_list_csv(m_profiles.group(1) if m_profiles else None, cwd=cwd, require_exists=False)

        if not pdbs:
            pdbs = _infer_mixture_pdbs_from_label(run_dir, label)
            if require_exists:
                pdbs = [p for p in pdbs if p.exists()]

        if require_exists and not pdbs:
            continue

        if weights is not None and pdbs and len(weights) != len(pdbs):
            # Keep the record, but make the mismatch visible rather than silently
            # assigning incorrect component weights.
            weight_mismatch = True
        else:
            weight_mismatch = False

        run_no = fallback_run_no
        m_label = re.match(r"^mol(\d+)_", label)
        if m_label:
            run_no = int(m_label.group(1))

        fit_file = run_dir / f"{label}_foxs_mixture_fit.dat"

        records.append({
            "type": "mixture",
            "run_no": run_no,
            "label": label,
            "chi2": chi2,
            "pdb_path": None,
            "pdb_paths": pdbs,
            "weights": weights,
            "scale": scale,
            "c1": c1,
            "c2": c2,
            "best_component_chi2": best_component_chi2,
            "profile_paths": profiles,
            "fit_file": fit_file if fit_file.exists() else None,
            "summary_file": mix_file,
            "line": line,
            "weight_mismatch": weight_mismatch,
        })

    return records


def read_foxs_prediction_records(
    fitdata_dir: str | Path,
    mode: str = "auto",
    require_exists: bool = True,
) -> list[dict]:
    """
    Read FoXS scoring records from a Carbonara ``fitdata`` directory.

    Parameters
    ----------
    fitdata_dir : str or Path
        Directory containing ``allAtomRun*/`` folders.
    mode : {"auto", "single", "mixture", "both"}
        ``single`` reads ``foxs_results.txt``.
        ``mixture`` reads ``foxs_mixture_results.txt``.
        ``auto`` reads mixture results for a run when present, otherwise single
        results. This is the safest default for mixed old/new analyses.
        ``both`` reads both files if both exist.
    require_exists : bool
        If True, discard records whose PDB paths cannot be found.

    Returns
    -------
    list[dict]
        Each record has at least:
        ``type`` ("single" or "mixture"), ``run_no``, ``label``, ``chi2``,
        ``pdb_paths``, ``weights``, ``fit_file`` and ``summary_file``.
    """
    fitdata_dir = Path(fitdata_dir)
    mode = str(mode).lower()
    if mode not in {"auto", "single", "mixture", "both"}:
        raise ValueError("mode must be 'auto', 'single', 'mixture', or 'both'")

    records = []
    cwd = Path.cwd()
    for run_dir in sorted(fitdata_dir.glob("allAtomRun*"), key=_sort_allatom_run_dir):
        if not run_dir.is_dir():
            continue

        has_mix = (run_dir / "foxs_mixture_results.txt").exists()
        has_single = (run_dir / "foxs_results.txt").exists()

        if mode == "single":
            records.extend(_read_single_foxs_records(run_dir, require_exists=require_exists, cwd=cwd))
        elif mode == "mixture":
            records.extend(_read_mixture_foxs_records(run_dir, require_exists=require_exists, cwd=cwd))
        elif mode == "both":
            if has_single:
                records.extend(_read_single_foxs_records(run_dir, require_exists=require_exists, cwd=cwd))
            if has_mix:
                records.extend(_read_mixture_foxs_records(run_dir, require_exists=require_exists, cwd=cwd))
        else:  # auto
            if has_mix:
                records.extend(_read_mixture_foxs_records(run_dir, require_exists=require_exists, cwd=cwd))
            elif has_single:
                records.extend(_read_single_foxs_records(run_dir, require_exists=require_exists, cwd=cwd))

    records.sort(key=lambda r: (
        10**12 if r.get("run_no") is None else int(r.get("run_no")),
        str(r.get("label", "")),
        float(r.get("chi2", np.inf)),
    ))
    return records


def _legacy_prediction_tuple(record: dict):
    """
    Convert a prediction record into the historical return style.

    Single records become ``(Path, chi2)``.
    Mixture records become ``([Path, ...], chi2)`` because a mixture prediction
    has several component structures.
    """
    if record.get("type") == "mixture":
        return (list(record.get("pdb_paths", [])), float(record["chi2"]))
    return (Path(record["pdb_paths"][0]), float(record["chi2"]))


def collect_good_prediction_files(
    fitdata_dir: str | Path,
    chi2_threshold: float,
    require_exists: bool = True,
    sort_by_chi2: bool = True,
    mode: str = "auto",
    return_records: bool = False,
):
    """
    Collect predictions whose FoXS chi^2 is <= ``chi2_threshold``.

    This now supports both Carbonara output modes:

    - ordinary single-structure scoring in ``allAtomRun*/foxs_results.txt``;
    - approximate MultiFoXS-style mixture scoring in
      ``allAtomRun*/foxs_mixture_results.txt``.

    Parameters
    ----------
    fitdata_dir : str or Path
        Path to the ``fitdata`` directory containing ``allAtomRun*`` folders.
    chi2_threshold : float
        Maximum chi^2 to accept.
    require_exists : bool
        If True, only return predictions whose PDB files can be found.
    sort_by_chi2 : bool
        If True, sort results by increasing chi^2.
    mode : {"auto", "single", "mixture", "both"}
        Which result files to read. In ``auto`` mode, a run uses mixture results
        when ``foxs_mixture_results.txt`` is present, otherwise ordinary single
        FoXS results.
    return_records : bool
        If False, preserve the old style as far as possible:
        single predictions return ``(Path, chi2)`` and mixture predictions return
        ``([Path, ...], chi2)``.
        If True, return dictionaries with metadata including weights and the
        mixture fit curve file.
    """
    records = [
        rec for rec in read_foxs_prediction_records(
            fitdata_dir,
            mode=mode,
            require_exists=require_exists,
        )
        if float(rec["chi2"]) <= float(chi2_threshold)
    ]

    if sort_by_chi2:
        records.sort(key=lambda r: float(r["chi2"]))

    if return_records:
        return records
    return [_legacy_prediction_tuple(rec) for rec in records]


def _normalise_prediction_record(prediction, fitdata_dir: str | Path | None = None) -> dict:
    """
    Convert a prediction specification into the internal record dictionary form.

    Accepted inputs
    ---------------
    - record dicts returned by ``collect_*`` with ``return_records=True``;
    - legacy ``(Path, chi2)`` or ``([Path, ...], chi2)`` tuples;
    - a single PDB path;
    - a raw list of component PDB paths for one mixture prediction.
    """
    if isinstance(prediction, dict):
        return prediction

    # Raw list of component PDB paths for one mixture prediction.
    if _is_path_sequence(prediction):
        pdb_paths = [Path(p) for p in prediction]
        recovered = _try_recover_prediction_record_from_files(
            pdb_paths,
            fitdata_dir=fitdata_dir,
            chi2=None,
        )
        if recovered is not None:
            return recovered

        meta = _parse_prediction_pdb_metadata(pdb_paths[0]) if pdb_paths else {}
        rec_type = "mixture" if len(pdb_paths) > 1 else "single"
        return {
            "type": rec_type,
            "run_no": meta.get("runNo"),
            "label": meta.get("label"),
            "chi2": None,
            "pdb_path": pdb_paths[0] if rec_type == "single" and pdb_paths else None,
            "pdb_paths": pdb_paths,
            "weights": [1.0] if rec_type == "single" else None,
            "scale": None,
            "fit_file": None,
            "summary_file": None,
            "line": "",
        }

    # Legacy collect_* style: (path_or_paths, chi2)
    if isinstance(prediction, tuple) and len(prediction) >= 2:
        paths, chi2 = prediction[0], float(prediction[1])
        if _is_path_sequence(paths):
            pdb_paths = [Path(p) for p in paths]
            rec_type = "mixture" if len(pdb_paths) > 1 else "single"
        else:
            pdb_paths = [Path(paths)]
            rec_type = "single"

        recovered = _try_recover_prediction_record_from_files(
            pdb_paths,
            fitdata_dir=fitdata_dir,
            chi2=chi2,
        )
        if recovered is not None:
            return recovered

        meta = _parse_prediction_pdb_metadata(pdb_paths[0]) if pdb_paths else {}
        return {
            "type": rec_type,
            "run_no": meta.get("runNo"),
            "label": meta.get("label"),
            "chi2": chi2,
            "pdb_path": pdb_paths[0] if rec_type == "single" and pdb_paths else None,
            "pdb_paths": pdb_paths,
            "weights": [1.0] if rec_type == "single" else None,
            "scale": None,
            "fit_file": None,
            "summary_file": None,
            "line": "",
        }

    if _is_pathlike_object(prediction):
        p = Path(prediction)
        recovered = _try_recover_prediction_record_from_files([p], fitdata_dir=fitdata_dir, chi2=None)
        if recovered is not None:
            return recovered

        meta = _parse_prediction_pdb_metadata(p)
        return {
            "type": "single",
            "run_no": meta.get("runNo"),
            "label": meta.get("label") or p.stem,
            "chi2": None,
            "pdb_path": p,
            "pdb_paths": [p],
            "weights": [1.0],
            "scale": None,
            "fit_file": None,
            "summary_file": None,
            "line": "",
        }

    raise TypeError(
        "prediction must be a record dict, a legacy (path(s), chi2) tuple, "
        "a path, or a raw list of component paths"
    )


def load_foxs_fit_curve(fit_file: str | Path, max_q: float | None = None):
    """
    Load a FoXS or mixture-fit curve file and return plotting arrays.

    Supported column layouts:
    - mixture curve: ``q I_exp sigma I_fit``;
    - FoXS-like 4-column curve: ``q I_exp sigma I_fit``;
    - FoXS-like 3-column curve: ``q I_exp I_fit``.
    """
    fit_file = Path(fit_file)
    data = np.loadtxt(fit_file)
    if data.ndim == 1:
        data = data[None, :]
    if data.shape[1] < 3:
        raise ValueError(f"Fit curve must have at least 3 columns: {fit_file}")

    q = data[:, 0].astype(float)
    if data.shape[1] >= 4:
        i_exp = data[:, 1].astype(float)
        sigma = data[:, 2].astype(float)
        i_fit = data[:, 3].astype(float)
        sigma = np.where(sigma <= 0, 1.0, sigma)
        residual = (i_exp - i_fit) / sigma
    else:
        i_exp = data[:, 1].astype(float)
        sigma = None
        i_fit = data[:, 2].astype(float)
        residual = i_exp - i_fit

    if max_q is not None:
        mask = q <= float(max_q)
        q = q[mask]
        i_exp = i_exp[mask]
        i_fit = i_fit[mask]
        residual = residual[mask]
        if sigma is not None:
            sigma = sigma[mask]

    return {
        "q": q,
        "i_exp": i_exp,
        "sigma": sigma,
        "i_fit": i_fit,
        "residual": residual,
        "fit_file": fit_file,
    }


def plot_foxs_fit_curve(
    fit_file: str | Path,
    chi2: float | None = None,
    title: str | None = None,
    max_q: float | None = None,
    figsize=(6.0, 6.0),
    show=True,
    save_path=None,
):
    """
    Plot a stored FoXS fit curve, including approximate MultiFoXS curves written
    by the watcher as ``*_foxs_mixture_fit.dat``.
    """
    curve = load_foxs_fit_curve(fit_file, max_q=max_q)
    q = curve["q"]
    i_exp = curve["i_exp"]
    i_fit = curve["i_fit"]
    residual = curve["residual"]

    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(2, 1, height_ratios=[3, 1], hspace=0.08)
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1], sharex=ax1)

    ax1.plot(q, i_exp, "o", ms=4, label="Experimental")
    ax1.plot(q, i_fit, "-", lw=2, label="Approx. MultiFoXS fit" if "mixture" in Path(fit_file).name else "FoXS fit")
    ax1.set_yscale("log")
    ax1.set_ylabel("Intensity")

    if title is None:
        title = "FoXS fit"
    if chi2 is not None:
        title += f"  (chi² = {float(chi2):.4g})"
    ax1.set_title(title)
    ax1.legend()
    ax1.tick_params(axis="x", labelbottom=False)

    ax2.axhline(0.0, lw=1)
    ax2.plot(q, residual, "o", ms=3)
    ax2.set_xlabel("q")
    ax2.set_ylabel("Residual")

    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches="tight")
    if show:
        plt.show()
    return fig


def plot_prediction_foxs_fit(prediction, fitdata_dir: str | Path | None = None, max_q=None, **kwargs):
    """
    Plot the stored SAXS/FoXS fit for a prediction record.

    For mixture records this uses the watcher's precomputed
    ``*_foxs_mixture_fit.dat`` file, i.e. the approximate MultiFoXS fit. For
    ordinary single records, a stored fit curve is only available if ``fit_file``
    is present in the record; otherwise use ``show_structure_and_foxs_side_by_side``
    to rerun FoXS for that single PDB.
    """
    rec = _normalise_prediction_record(prediction, fitdata_dir=fitdata_dir)
    fit_file = rec.get("fit_file")
    if fit_file is None:
        raise ValueError(
            "No stored fit curve is associated with this prediction. "
            "For single-PDB predictions, use show_structure_and_foxs_side_by_side(...) "
            "to rerun FoXS, or pass a record with fit_file set."
        )
    return plot_foxs_fit_curve(
        fit_file,
        chi2=rec.get("chi2"),
        title=("Approx. MultiFoXS fit: " + str(rec.get("label"))) if rec.get("type") == "mixture" else None,
        max_q=max_q,
        **kwargs,
    )


def visualisePredictionMixture(pdb_paths, weights=None, ncols=3, panel_width=350, panel_height=300, show_labels=True):
    """
    Visualise the component structures of one mixture prediction.

    ``pdb_paths`` may be a list of paths or a prediction record returned by
    ``collect_good_prediction_files(..., return_records=True)``.
    """
    if isinstance(pdb_paths, dict):
        rec = pdb_paths
        weights = rec.get("weights") if weights is None else weights
        pdb_paths = rec.get("pdb_paths", [])
    elif isinstance(pdb_paths, tuple) and len(pdb_paths) >= 1:
        # Legacy collect_good_prediction_files tuple: ([pdbs], chi2)
        first = pdb_paths[0]
        if isinstance(first, (list, tuple)):
            pdb_paths = first

    pdb_paths = [Path(p) for p in pdb_paths]
    if not pdb_paths:
        raise ValueError("No component PDB files supplied for mixture visualisation.")

    if show_labels and weights is not None and len(weights) == len(pdb_paths):
        print("Mixture weights:")
        for p, w in zip(pdb_paths, weights):
            print(f"  {Path(p).name}: {float(w):.4g}")

    return visualisePrediction_panel(
        pdb_paths,
        ncols=ncols,
        panel_width=panel_width,
        panel_height=panel_height,
        show_labels=show_labels,
        color_by_chain=True,
    )


############################################

## Visulaisation routine for a specific path, used in live notebook

############################################


def visualisePredictionIndividual(aa_path):
    """
    Visualise a single AA PDB, or the component PDBs of a mixture prediction.

    Backwards compatible behaviour:
      visualisePredictionIndividual("model_AA.pdb")

    New mixture-aware behaviour:
      visualisePredictionIndividual(record)
      visualisePredictionIndividual([pdb0, pdb1, ...])
      visualisePredictionIndividual(([pdb0, pdb1, ...], chi2))
    """
    if isinstance(aa_path, dict):
        if aa_path.get("type") == "mixture" or len(aa_path.get("pdb_paths", [])) > 1:
            return visualisePredictionMixture(aa_path)
        paths = aa_path.get("pdb_paths", [])
        if paths:
            aa_path = paths[0]

    elif isinstance(aa_path, tuple) and len(aa_path) >= 1:
        first = aa_path[0]
        if isinstance(first, (list, tuple)):
            return visualisePredictionMixture(first)
        aa_path = first

    elif isinstance(aa_path, (list, tuple)) and not isinstance(aa_path, (str, bytes)):
        return visualisePredictionMixture(aa_path)

    if not HAS_PY3DMOL:
        _warn_missing_py3dmol()
        return None

    view = py3Dmol.view(width=800, height=600)

    structure_data, fmt = _read_structure_for_viewer(aa_path)
    view.addModel(structure_data, fmt)   # model 0
    chains = _chains_present_in_structure(structure_data, fmt)

    palette = ["blue", "green", "red", "yellow", "cyan", "magenta",
               "orange", "purple", "lime", "gray"]

    if chains:
        for i, ch in enumerate(chains):
            color = palette[i % len(palette)]
            view.setStyle({"model": 0, "chain": ch}, {"cartoon": {"color": color}})
    else:
        view.setStyle({"model": 0}, {"cartoon": {"color": "lightgray"}})

    view.zoomTo()
    view.show()
    return view


def visualisePredictionComp(pdb1, pdb2, do_superpose=True):
    if not HAS_PY3DMOL:
        _warn_missing_py3dmol()
        return None

    view = py3Dmol.view(width=800, height=600)

    data1, fmt1 = _read_structure_for_viewer(pdb1)

    if do_superpose:
        data2_to_show, rmsd, nmatch = superimpose_structure_files_by_ca(pdb1, pdb2)
        fmt2 = "pdb"
        print(f"Superposed model 1 onto model 0 using {nmatch} matched Cα atoms. RMSD = {rmsd:.3f} Å")
    else:
        data2_to_show, fmt2 = _read_structure_for_viewer(pdb2)

    view.addModel(data1, fmt1)
    view.addModel(data2_to_show, fmt2)

    chains1 = _chains_present_in_pdb(data1) if fmt1 == "pdb" else None
    chains2 = _chains_present_in_pdb(data2_to_show) if fmt2 == "pdb" else None

    palette = ["blue", "green", "red", "yellow", "cyan", "magenta",
               "orange", "purple", "lime", "gray"]

    if chains1:
        for i, ch in enumerate(chains1):
            color = palette[i % len(palette)]
            view.setStyle({"model": 0, "chain": ch}, {"cartoon": {"color": color}})
    else:
        view.setStyle({"model": 0}, {"cartoon": {"color": "lightgray"}})

    if chains2:
        for i, ch in enumerate(chains2):
            color = palette[(i + 1) % len(palette)]
            view.setStyle({"model": 1, "chain": ch}, {"cartoon": {"color": color, "opacity": 0.6}})
    else:
        view.setStyle({"model": 1}, {"cartoon": {"color": "red"}})

    view.zoomTo()
    view.show()



def visualise_linker_sections(structure_path, selected_secs, chain_id_map=None, ignore_chain=False):
    structure_path = str(structure_path)
    suffix = Path(structure_path).suffix.lower()

    if suffix == ".pdb":
        fmt = "pdb"
    elif suffix in [".cif", ".mmcif"]:
        fmt = "cif"
    else:
        raise ValueError(f"Unsupported structure format: {suffix}")

    with open(structure_path, "r") as f:
        structure_data = f.read()

    view = py3Dmol.view(width=900, height=650)
    view.addModel(structure_data, fmt)

    view.setStyle(
        {"model": 0},
        {"cartoon": {"color": "lightgray", "opacity": 0.5}}
    )

    colours = ["red", "orange", "yellow", "cyan", "magenta", "lime", "blue"]

    # Try to infer PDB chain IDs if none supplied
    inferred_chain_ids = []
    if fmt == "pdb":
        seen = set()
        for line in structure_data.splitlines():
            if line.startswith(("ATOM", "HETATM")) and len(line) > 21:
                ch = line[21].strip()
                if ch and ch not in seen:
                    seen.add(ch)
                    inferred_chain_ids.append(ch)

    if chain_id_map is None and inferred_chain_ids:
        chain_id_map = {i + 1: ch for i, ch in enumerate(inferred_chain_ids)}

    print("chain_id_map =", chain_id_map)
    print("ignore_chain =", ignore_chain)

    parsed = []

    for i, sec in enumerate(selected_secs):
        m = re.search(r"Chain\s+(\d+)\s+ResID:\s*(\d+)-(\d+)", sec)
        if not m:
            print(f"Could not parse section label: {sec}")
            continue

        chain_num = int(m.group(1))
        start = int(m.group(2))
        end = int(m.group(3))
        colour = colours[i % len(colours)]

        sel = {"model": 0, "resi": f"{start}-{end}"}

        if not ignore_chain and chain_id_map is not None and chain_num in chain_id_map:
            sel["chain"] = chain_id_map[chain_num]

        print("Applying selection:", sel, "for", sec)

        view.addStyle(
            sel,
            {"cartoon": {"color": colour, "opacity": 1.0}}
        )

        parsed.append((chain_num, start, end))

    print("Parsed sections:", parsed)

    view.zoomTo()
    view.show()


def visualisePredictionComp_panel(
    file_list,
    reference_file,
    do_superpose=True,
    ncols=3,
    panel_width=350,
    panel_height=300,
    max_panels=None,
    show_labels=True,
):
    """
    Show a grid of pairwise comparisons against a fixed reference.

    model 0 = reference
    model 1 = one member of file_list
    """
    if not HAS_PY3DMOL:
        _warn_missing_py3dmol()
        return None

    if max_panels is not None:
        file_list = file_list[:max_panels]

    n = len(file_list)
    if n == 0:
        print("No files to display.")
        return None

    ncols = max(1, int(ncols))
    nrows = (n + ncols - 1) // ncols

    view = py3Dmol.view(
        viewergrid=(nrows, ncols),
        width=ncols * panel_width,
        height=nrows * panel_height,
        linked=False,
    )

    ref_data, ref_fmt = _read_structure_for_viewer(reference_file)

    for k, mobile_file in enumerate(file_list):
        r = k // ncols
        c = k % ncols
        viewer = (r, c)

        try:
            # reference
            view.addModel(ref_data, ref_fmt, viewer=viewer)

            # mobile
            if do_superpose:
                mob_data_to_show, rmsd, nmatch = superimpose_structure_files_by_ca(
                    reference_file, mobile_file
                )
                mob_fmt = "pdb"
                panel_title = f"{Path(mobile_file).name}\nRMSD={rmsd:.2f} Å"
            else:
                mob_data_to_show, mob_fmt = _read_structure_for_viewer(mobile_file)
                panel_title = Path(mobile_file).name

            view.addModel(mob_data_to_show, mob_fmt, viewer=viewer)

            # style reference
            ref_chains = _chains_present_in_structure(ref_data, ref_fmt)
            if ref_chains:
                for ch in ref_chains:
                    view.setStyle(
                        {"model": 0, "chain": ch},
                        {"cartoon": {"color": "lightgray", "opacity": 0.85}},
                        viewer=viewer,
                    )
            else:
                view.setStyle(
                    {"model": 0},
                    {"cartoon": {"color": "lightgray", "opacity": 0.85}},
                    viewer=viewer,
                )

            # style mobile
            mob_chains = _chains_present_in_structure(mob_data_to_show, mob_fmt)
            if mob_chains:
                for ch in mob_chains:
                    view.setStyle(
                        {"model": 1, "chain": ch},
                        {"cartoon": {"color": "red", "opacity": 0.85}},
                        viewer=viewer,
                    )
            else:
                view.setStyle(
                    {"model": 1},
                    {"cartoon": {"color": "red", "opacity": 0.85}},
                    viewer=viewer,
                )

            view.zoomTo(viewer=viewer)

            if show_labels:
                view.addLabel(
                    panel_title,
                    {
                        "fontSize": 10,
                        "backgroundColor": "white",
                        "backgroundOpacity": 0.7,
                        "fontColor": "black",
                        "borderThickness": 0,
                        "inFront": True,
                    },
                    viewer=viewer,
                )

        except Exception as e:
            print(f"Failed for {mobile_file}: {e}")
            view.addLabel(
                f"Failed:\n{Path(mobile_file).name}",
                {
                    "fontSize": 12,
                    "backgroundColor": "mistyrose",
                    "backgroundOpacity": 0.8,
                    "fontColor": "black",
                    "borderThickness": 0,
                    "inFront": True,
                },
                viewer=viewer,
            )

    view.show()
    return view




def show_structure_and_foxs_side_by_side(
    pdb_name,
    saxs_name,
    foxs_cmd="pyfoxs",
    max_q=None,
    structure_width=480,
    structure_height=420,
    plot_width=520,
    print_summary=False,
):
    """
    Display:
      left  = structure viewer
      right = FoXS fit + residuals

    Returns
    -------
    dict with keys:
        chi2, stdout, stderr, fit_file, view_html
    """

    pdb_path = Path(pdb_name).resolve()
    saxs_path = Path(saxs_name).resolve()

    if not pdb_path.exists():
        raise FileNotFoundError(f"Structure file not found: {pdb_path}")
    if not saxs_path.exists():
        raise FileNotFoundError(f"SAXS file not found: {saxs_path}")

    temp_pdb_to_clean = None

    try:
        pdb_for_foxs, temp_pdb_to_clean = CDT._convert_cif_to_pdb_for_foxs(pdb_path)
        pdb_for_foxs = Path(pdb_for_foxs)

        base_cmd = CDT._normalise_cmd(foxs_cmd)
        cmd = base_cmd + [str(pdb_for_foxs), str(saxs_path)]

        if max_q is not None:
            cmd += ["--max_q", str(max_q)]

        proc = subprocess.run(cmd, capture_output=True, text=True)
        stdout = proc.stdout
        stderr = proc.stderr

        if proc.returncode != 0:
            raise RuntimeError(
                f"pyFoXS failed with exit code {proc.returncode}\n\nSTDERR:\n{stderr}\n\nSTDOUT:\n{stdout}"
            )

        combined_text = stdout + "\n" + stderr
        chi2 = None
        chi_patterns = [
            r"Chi(?:\^?2| square)\s*[:=]\s*([0-9.eE+-]+)",
            r"chi(?:\^?2| square)\s*[:=]\s*([0-9.eE+-]+)",
            r"\bchi\s*=\s*([0-9.eE+-]+)",
            r"\bChi\s*=\s*([0-9.eE+-]+)",
            r"\bchi2\s*[:=]\s*([0-9.eE+-]+)",
            r"\bChi2\s*[:=]\s*([0-9.eE+-]+)",
        ]
        for pat in chi_patterns:
            m = re.search(pat, combined_text)
            if m:
                try:
                    chi2 = float(m.group(1))
                    break
                except ValueError:
                    pass

        fit_file = CDT._find_foxs_fit_file(pdb_for_foxs, saxs_path)

        if fit_file is None:
            raise FileNotFoundError("Could not identify a FoXS fit file automatically.")

        fit = CDT.load_numeric_table_loose(fit_file, min_cols=2)
        q = fit[:, 0]

        if fit.shape[1] >= 4:
            i_exp = fit[:, 1]
            sigma = fit[:, 2]
            i_fit = fit[:, 3]
            residual = (i_exp - i_fit) / sigma
        elif fit.shape[1] == 3:
            i_exp = fit[:, 1]
            i_fit = fit[:, 2]
            residual = i_exp - i_fit
        else:
            raise ValueError("Fit file must have at least 3 columns for residual plotting.")

        # -----------------------------
        # Build matplotlib figure -> HTML image
        # -----------------------------
        fig = plt.figure(figsize=(6.0, 6.0))
        gs = fig.add_gridspec(2, 1, height_ratios=[3, 1], hspace=0.08)

        ax1 = fig.add_subplot(gs[0])
        ax2 = fig.add_subplot(gs[1], sharex=ax1)

        ax1.plot(q, i_exp, "o", ms=4, label="Experimental")
        ax1.plot(q, i_fit, "-", lw=2, label="FoXS fit")
        ax1.set_yscale("log")
        ax1.set_ylabel("Intensity")

        title = "Initial FoXS check"
        if chi2 is not None:
            title += f"  (chi² = {chi2:.4g})"
        if max_q is not None:
            title += f", max_q={max_q}"
        ax1.set_title(title)
        ax1.legend()
        ax1.tick_params(axis="x", labelbottom=False)

        ax2.axhline(0.0, lw=1)
        ax2.plot(q, residual, "o", ms=3)
        ax2.set_xlabel("q")
        ax2.set_ylabel("Residual")

        plt.tight_layout()

        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=160, bbox_inches="tight")
        plt.close(fig)
        buf.seek(0)
        plot_b64 = base64.b64encode(buf.read()).decode("utf-8")
        plot_html = f'<img src="data:image/png;base64,{plot_b64}" style="width:{plot_width}px; max-width:100%;">'

        # -----------------------------
        # Build py3Dmol viewer -> HTML quietly
        # -----------------------------
        if not HAS_PY3DMOL:
            _warn_missing_py3dmol()
            viewer_html = "<div style='padding:20px;border:1px solid #ddd;border-radius:6px;'>py3Dmol is not available.</div>"
        else:
            view = py3Dmol.view(width=structure_width, height=structure_height)

            structure_data, fmt = _read_structure_for_viewer(pdb_name)
            view.addModel(structure_data, fmt)

            chains = _chains_present_in_structure(structure_data, fmt)
            palette = ["blue", "green", "red", "yellow", "cyan", "magenta",
                       "orange", "purple", "lime", "gray"]

            if chains:
                for i, ch in enumerate(chains):
                    color = palette[i % len(palette)]
                    view.setStyle({"model": 0, "chain": ch}, {"cartoon": {"color": color}})
            else:
                view.setStyle({"model": 0}, {"cartoon": {"color": "lightgray"}})

            view.zoomTo()

            silent_out = io.StringIO()
            silent_err = io.StringIO()
            with contextlib.redirect_stdout(silent_out), contextlib.redirect_stderr(silent_err):
                viewer_html = view._make_html()

        # -----------------------------
        # Display side by side
        # -----------------------------
        html = f"""
        <div style="
            display:flex;
            flex-wrap:wrap;
            gap:20px;
            align-items:flex-start;
            margin-top:10px;
            margin-bottom:10px;
        ">
            <div style="flex:0 0 auto;">
                <div style="font-weight:600; margin-bottom:8px;">Structure</div>
                {viewer_html}
            </div>
            <div style="flex:0 0 auto;">
                <div style="font-weight:600; margin-bottom:8px;">Initial FoXS fit</div>
                {plot_html}
            </div>
        </div>
        """
        display(HTML(html))

    finally:
        if temp_pdb_to_clean is not None:
            try:
                os.remove(temp_pdb_to_clean)
            except OSError:
                pass


def show_prediction_record_and_foxs_side_by_side(
    prediction,
    saxs_name=None,
    foxs_cmd="pyfoxs",
    max_q=None,
    fitdata_dir: str | Path | None = None,
    structure_width=480,
    structure_height=420,
    plot_width=520,
):
    """
    Display a prediction plus its SAXS fit.

    For ordinary single-PDB predictions this delegates to
    ``show_structure_and_foxs_side_by_side`` and reruns FoXS unless a stored
    ``fit_file`` is supplied in the record.

    For mixture predictions this displays all component structures in a panel
    and plots the stored approximate MultiFoXS curve
    ``*_foxs_mixture_fit.dat``.
    """
    rec = _normalise_prediction_record(prediction, fitdata_dir=fitdata_dir)

    if rec.get("type") != "mixture" and len(rec.get("pdb_paths", [])) <= 1:
        pdb = rec.get("pdb_path") or rec.get("pdb_paths", [None])[0]
        if rec.get("fit_file") is None:
            if saxs_name is None:
                raise ValueError("saxs_name is required to rerun FoXS for a single-PDB prediction.")
            return show_structure_and_foxs_side_by_side(
                pdb,
                saxs_name,
                foxs_cmd=foxs_cmd,
                max_q=max_q,
                structure_width=structure_width,
                structure_height=structure_height,
                plot_width=plot_width,
            )

    pdbs = [Path(p) for p in rec.get("pdb_paths", [])]
    if not pdbs:
        raise ValueError("Prediction record contains no PDB paths.")

    fit_file = rec.get("fit_file")
    if fit_file is None:
        raise ValueError(
            "No stored mixture fit curve found for this prediction. "
            "Use collect_good_prediction_files(..., return_records=True) so the "
            "record includes fit_file, or pass fitdata_dir to recover it."
        )

    # Build structure panel HTML.
    if not HAS_PY3DMOL:
        _warn_missing_py3dmol()
        viewer_html = "<div style='padding:20px;border:1px solid #ddd;border-radius:6px;'>py3Dmol is not available.</div>"
    else:
        n = len(pdbs)
        ncols = min(3, max(1, n))
        nrows = (n + ncols - 1) // ncols
        view = py3Dmol.view(
            viewergrid=(nrows, ncols),
            width=ncols * structure_width,
            height=nrows * structure_height,
            linked=False,
        )
        weights = rec.get("weights")
        palette = ["blue", "green", "red", "yellow", "cyan", "magenta", "orange", "purple", "lime", "gray"]

        for k, pdb in enumerate(pdbs):
            r = k // ncols
            c = k % ncols
            viewer = (r, c)
            structure_data, fmt = _read_structure_for_viewer(pdb)
            view.addModel(structure_data, fmt, viewer=viewer)
            chains = _chains_present_in_structure(structure_data, fmt)
            if chains:
                for chain_i, ch in enumerate(chains):
                    view.setStyle(
                        {"model": 0, "chain": ch},
                        {"cartoon": {"color": palette[chain_i % len(palette)], "opacity": 0.9}},
                        viewer=viewer,
                    )
            else:
                view.setStyle({"model": 0}, {"cartoon": {"color": palette[k % len(palette)], "opacity": 0.9}}, viewer=viewer)
            label = Path(pdb).name
            if weights is not None and k < len(weights):
                label += f"\nw={float(weights[k]):.3g}"
            view.addLabel(
                label,
                {"fontSize": 10, "backgroundColor": "white", "backgroundOpacity": 0.7, "fontColor": "black", "borderThickness": 0, "inFront": True},
                viewer=viewer,
            )
            view.zoomTo(viewer=viewer)

        silent_out = io.StringIO()
        silent_err = io.StringIO()
        with contextlib.redirect_stdout(silent_out), contextlib.redirect_stderr(silent_err):
            viewer_html = view._make_html()

    # Build SAXS plot HTML from stored mixture curve.
    curve = load_foxs_fit_curve(fit_file, max_q=max_q)
    fig = plt.figure(figsize=(6.0, 6.0))
    gs = fig.add_gridspec(2, 1, height_ratios=[3, 1], hspace=0.08)
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1], sharex=ax1)
    ax1.plot(curve["q"], curve["i_exp"], "o", ms=4, label="Experimental")
    ax1.plot(curve["q"], curve["i_fit"], "-", lw=2, label="Approx. MultiFoXS fit")
    ax1.set_yscale("log")
    ax1.set_ylabel("Intensity")
    title = f"Approx. MultiFoXS fit: {rec.get('label', '')}"
    if rec.get("chi2") is not None:
        title += f"  (chi² = {float(rec['chi2']):.4g})"
    ax1.set_title(title)
    ax1.legend()
    ax1.tick_params(axis="x", labelbottom=False)
    ax2.axhline(0.0, lw=1)
    ax2.plot(curve["q"], curve["residual"], "o", ms=3)
    ax2.set_xlabel("q")
    ax2.set_ylabel("Residual")
    plt.tight_layout()
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=160, bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    plot_b64 = base64.b64encode(buf.read()).decode("utf-8")
    plot_html = f'<img src="data:image/png;base64,{plot_b64}" style="width:{plot_width}px; max-width:100%;">'

    html = f"""
    <div style="display:flex; flex-wrap:wrap; gap:20px; align-items:flex-start; margin-top:10px; margin-bottom:10px;">
        <div style="flex:0 0 auto;">
            <div style="font-weight:600; margin-bottom:8px;">Mixture component structures</div>
            {viewer_html}
        </div>
        <div style="flex:0 0 auto;">
            <div style="font-weight:600; margin-bottom:8px;">Approx. MultiFoXS fit</div>
            {plot_html}
        </div>
    </div>
    """
    display(HTML(html))
    return {"record": rec, "fit_file": fit_file}


def visualisePrediction_panel(
    file_list,
    ncols=3,
    panel_width=350,
    panel_height=300,
    max_panels=None,
    show_labels=True,
    color_by_chain=True,
):
    """
    Display a set of structures side by side in a py3Dmol viewer grid.

    Each panel contains one independent structure. No alignment or
    comparison against a reference structure is performed.

    Parameters
    ----------
    file_list : sequence of str or Path
        Structure files to display.

    ncols : int
        Number of columns in the viewer grid.

    panel_width, panel_height : int
        Approximate dimensions of each panel in pixels.

    max_panels : int or None
        Maximum number of structures to display.

    show_labels : bool
        Show the filename in each panel.

    color_by_chain : bool
        If True, assign a different colour to each chain.
        If False, display the whole structure in one colour.

    Returns
    -------
    py3Dmol.view or None
    """
    if not HAS_PY3DMOL:
        _warn_missing_py3dmol()
        return None

    file_list = list(file_list)

    if max_panels is not None:
        file_list = file_list[:max_panels]

    n = len(file_list)

    if n == 0:
        print("No files to display.")
        return None

    ncols = max(1, int(ncols))
    nrows = (n + ncols - 1) // ncols

    view = py3Dmol.view(
        viewergrid=(nrows, ncols),
        width=ncols * panel_width,
        height=nrows * panel_height,
        linked=False,
    )

    palette = [
        "blue",
        "green",
        "red",
        "yellow",
        "cyan",
        "magenta",
        "orange",
        "purple",
        "lime",
        "gray",
    ]

    for k, structure_file in enumerate(file_list):
        r = k // ncols
        c = k % ncols
        viewer = (r, c)

        try:
            structure_data, structure_fmt = _read_structure_for_viewer(
                structure_file
            )

            view.addModel(
                structure_data,
                structure_fmt,
                viewer=viewer,
            )

            chains = _chains_present_in_structure(
                structure_data,
                structure_fmt,
            )

            if color_by_chain and chains:
                for chain_i, chain_id in enumerate(chains):
                    view.setStyle(
                        {
                            "model": 0,
                            "chain": chain_id,
                        },
                        {
                            "cartoon": {
                                "color": palette[chain_i % len(palette)],
                                "opacity": 0.9,
                            }
                        },
                        viewer=viewer,
                    )
            else:
                view.setStyle(
                    {"model": 0},
                    {
                        "cartoon": {
                            "color": palette[k % len(palette)],
                            "opacity": 0.9,
                        }
                    },
                    viewer=viewer,
                )

            view.zoomTo(viewer=viewer)

            if show_labels:
                view.addLabel(
                    Path(structure_file).name,
                    {
                        "fontSize": 10,
                        "backgroundColor": "white",
                        "backgroundOpacity": 0.7,
                        "fontColor": "black",
                        "borderThickness": 0,
                        "inFront": True,
                    },
                    viewer=viewer,
                )

        except Exception as e:
            print(f"Failed for {structure_file}: {e}")

            view.addLabel(
                f"Failed:\n{Path(structure_file).name}",
                {
                    "fontSize": 12,
                    "backgroundColor": "mistyrose",
                    "backgroundOpacity": 0.8,
                    "fontColor": "black",
                    "borderThickness": 0,
                    "inFront": True,
                },
                viewer=viewer,
            )

    view.show()
    return view

def collect_best_prediction_per_run_closest_to_one(
    fitdata_dir: str | Path,
    require_exists: bool = True,
    mode: str = "auto",
    return_records: bool = False,
    target_chi2: float = 1.0,
):
    """
    For each run number from 1 up to the maximum detected ``allAtomRun*``,
    select the prediction whose FoXS chi^2 is closest to ``target_chi2``.

    This now supports both ordinary single-structure FoXS summaries
    (``foxs_results.txt``) and mixture summaries
    (``foxs_mixture_results.txt``). In ``mode='auto'`` a run uses the mixture
    summary when it exists, otherwise the ordinary single summary.

    Missing runs or runs with no valid prediction return ``None``.

    Returns
    -------
    list
        If ``return_records=False``:
            single run entry   -> ``(Path, chi2)``
            mixture run entry  -> ``([Path, ...], chi2)``
        If ``return_records=True``:
            each non-missing entry is a metadata dictionary containing
            ``pdb_paths``, ``weights``, ``fit_file`` and ``chi2``.
    """
    fitdata_dir = Path(fitdata_dir)

    run_dirs = {}
    max_run = 0
    for p in fitdata_dir.glob("allAtomRun*"):
        if not p.is_dir():
            continue
        m = re.fullmatch(r"allAtomRun(\d+)", p.name)
        if not m:
            continue
        run_no = int(m.group(1))
        run_dirs[run_no] = p
        max_run = max(max_run, run_no)
    for p in fitdata_dir.glob("fitLog*.dat"):
        m = re.fullmatch(r"fitLog(\d+)\.dat", p.name)
        if not m:
            continue
        max_run = max(max_run, int(m.group(1)))

    if max_run == 0:
        return []

    all_records = read_foxs_prediction_records(
        fitdata_dir,
        mode=mode,
        require_exists=require_exists,
    )
    by_run = defaultdict(list)
    for rec in all_records:
        run_no = rec.get("run_no")
        if run_no is not None:
            by_run[int(run_no)].append(rec)

    results = []

    for run_no in range(1, max_run + 1):
        recs = by_run.get(run_no, [])
        if not recs:
            results.append(None)
            continue

        best = min(recs, key=lambda r: abs(float(r["chi2"]) - float(target_chi2)))
        results.append(best if return_records else _legacy_prediction_tuple(best))

    return results

#########################################################
#
# Weighted, diversity-constrained selection of structures
# for seeding a new Carbonara run
#
#########################################################


def pairwise_results_to_matrices(pdb_files, pairwise_results):
    """Convert ``pairwise_structure_metrics`` output into dense matrices.

    Parameters
    ----------
    pdb_files : sequence of path-like
        Structures in the same order used for the pairwise calculation.
    pairwise_results : sequence of dict
        Output from ``pairwise_structure_metrics`` or
        ``pairwise_structure_metrics_mp``.

    Returns
    -------
    rmsd_matrix, tm_matrix, gdt_matrix : np.ndarray
        Symmetric ``(N, N)`` matrices. Their diagonals are 0, 1 and 100,
        respectively.
    """
    n = len(pdb_files)
    rmsd = np.full((n, n), np.nan, dtype=float)
    tm = np.full((n, n), np.nan, dtype=float)
    gdt = np.full((n, n), np.nan, dtype=float)

    np.fill_diagonal(rmsd, 0.0)
    np.fill_diagonal(tm, 1.0)
    np.fill_diagonal(gdt, 100.0)

    for rec in pairwise_results:
        i = int(rec["i"])
        j = int(rec["j"])
        if not (0 <= i < n and 0 <= j < n and i != j):
            raise ValueError(f"Invalid pair indices ({i}, {j}) for N={n}.")

        rmsd[i, j] = rmsd[j, i] = float(rec["rmsd"])
        tm[i, j] = tm[j, i] = float(rec["tm"])
        gdt[i, j] = gdt[j, i] = float(rec["gdt_ts"])

    missing = np.argwhere(np.isnan(rmsd))
    if len(missing):
        i, j = map(int, missing[0])
        raise ValueError(
            "The pairwise results are incomplete; "
            f"the first missing RMSD is for pair ({i}, {j})."
        )

    return rmsd, tm, gdt


def _passes_seed_separation(
    candidate,
    selected,
    rmsd_matrix,
    tm_matrix,
    min_rmsd=None,
    max_tm=None,
    separation_rule="both",
):
    """Return True when a candidate is sufficiently distinct from all seeds."""
    if separation_rule not in {"both", "either"}:
        raise ValueError("separation_rule must be 'both' or 'either'.")

    for other in selected:
        tests = []
        if min_rmsd is not None:
            tests.append(rmsd_matrix[candidate, other] >= float(min_rmsd))
        if max_tm is not None:
            tests.append(tm_matrix[candidate, other] <= float(max_tm))

        if not tests:
            continue

        pair_ok = all(tests) if separation_rule == "both" else any(tests)
        if not pair_ok:
            return False

    return True


def _medoid_objective(distance_matrix, weights, medoids, distance_power=2.0):
    """Weighted distance-to-nearest-medoid objective."""
    medoids = np.asarray(medoids, dtype=int)
    nearest = np.min(distance_matrix[:, medoids], axis=1)
    return float(np.sum(weights * nearest**float(distance_power)))


def _constrained_weighted_kmedoids(
    distance_matrix,
    weights,
    n_select,
    rmsd_matrix,
    tm_matrix,
    min_rmsd=None,
    max_tm=None,
    separation_rule="both",
    distance_power=2.0,
    n_starts=20,
    max_swap_passes=50,
):
    """PAM-like weighted k-medoids with pairwise seed-separation constraints.

    This is deliberately dependency-free. It uses several deterministic greedy
    starts followed by medoid/non-medoid swap refinement.
    """
    distance_matrix = np.asarray(distance_matrix, dtype=float)
    weights = np.asarray(weights, dtype=float)
    n = len(weights)

    if distance_matrix.shape != (n, n):
        raise ValueError("distance_matrix has the wrong shape.")
    if not 1 <= int(n_select) <= n:
        raise ValueError("n_select must lie between 1 and the retained ensemble size.")

    n_select = int(n_select)
    n_starts = max(1, min(int(n_starts), n))

    # Prefer good one-medoid solutions as deterministic starting points.
    singleton_cost = np.sum(
        weights[:, None] * distance_matrix**float(distance_power), axis=0
    )
    start_candidates = np.argsort(singleton_cost, kind="stable")[:n_starts]

    best_medoids = None
    best_objective = np.inf

    for first in start_candidates:
        medoids = [int(first)]

        # Greedily add the candidate giving the greatest objective reduction.
        while len(medoids) < n_select:
            best_add = None
            best_add_obj = np.inf

            for candidate in range(n):
                if candidate in medoids:
                    continue
                if not _passes_seed_separation(
                    candidate,
                    medoids,
                    rmsd_matrix,
                    tm_matrix,
                    min_rmsd=min_rmsd,
                    max_tm=max_tm,
                    separation_rule=separation_rule,
                ):
                    continue

                trial = medoids + [candidate]
                obj = _medoid_objective(
                    distance_matrix, weights, trial, distance_power=distance_power
                )
                if obj < best_add_obj - 1e-12:
                    best_add_obj = obj
                    best_add = candidate

            if best_add is None:
                medoids = None
                break
            medoids.append(int(best_add))

        if medoids is None:
            continue

        # Standard PAM-style local swap refinement, respecting separation.
        current_obj = _medoid_objective(
            distance_matrix, weights, medoids, distance_power=distance_power
        )

        for _ in range(int(max_swap_passes)):
            swap_medoids = None
            swap_obj = current_obj
            medoid_set = set(medoids)

            for pos in range(n_select):
                fixed = medoids[:pos] + medoids[pos + 1 :]
                for candidate in range(n):
                    if candidate in medoid_set:
                        continue
                    if not _passes_seed_separation(
                        candidate,
                        fixed,
                        rmsd_matrix,
                        tm_matrix,
                        min_rmsd=min_rmsd,
                        max_tm=max_tm,
                        separation_rule=separation_rule,
                    ):
                        continue

                    trial = list(medoids)
                    trial[pos] = candidate
                    obj = _medoid_objective(
                        distance_matrix,
                        weights,
                        trial,
                        distance_power=distance_power,
                    )
                    if obj < swap_obj - 1e-12:
                        swap_obj = obj
                        swap_medoids = trial

            if swap_medoids is None:
                break

            medoids = swap_medoids
            current_obj = swap_obj

        if current_obj < best_objective - 1e-12:
            best_objective = current_obj
            best_medoids = list(map(int, medoids))

    if best_medoids is None:
        raise ValueError(
            f"Could not find {n_select} mutually separated structures. "
            "Reduce n_select, lower min_rmsd, raise max_tm, or use "
            "separation_rule='either'."
        )

    medoids = np.asarray(best_medoids, dtype=int)
    labels = np.argmin(distance_matrix[:, medoids], axis=1)

    return medoids, labels, best_objective


def _highest_weight_representatives_with_separation(
    labels,
    medoids,
    weights,
    rmsd_matrix,
    tm_matrix,
    min_rmsd=None,
    max_tm=None,
    separation_rule="both",
):
    """Choose one high-weight member per cluster without breaking separation.

    The globally best feasible combination is found by a small branch-and-bound
    search. For the intended M=3--4 this is inexpensive. The medoid combination
    is always a feasible fallback because the clustering itself was constrained.
    """
    labels = np.asarray(labels, dtype=int)
    medoids = np.asarray(medoids, dtype=int)
    weights = np.asarray(weights, dtype=float)
    n_clusters = len(medoids)

    candidates = []
    for cluster in range(n_clusters):
        members = np.flatnonzero(labels == cluster)
        # Ensure the medoid remains available even under pathological ties.
        if medoids[cluster] not in members:
            members = np.append(members, medoids[cluster])
        members = sorted(
            set(map(int, members)),
            key=lambda i: (-weights[i], i),
        )
        candidates.append(members)

    # Search clusters with fewer choices first, retaining original cluster IDs.
    cluster_order = sorted(range(n_clusters), key=lambda c: len(candidates[c]))
    upper_best = [max(weights[i] for i in candidates[c]) for c in cluster_order]
    remaining_upper = np.cumsum(upper_best[::-1])[::-1]

    best_score = -np.inf
    best_by_cluster = None
    chosen = []
    chosen_by_cluster = {}

    def recurse(depth, score):
        nonlocal best_score, best_by_cluster

        if depth == n_clusters:
            if score > best_score + 1e-15:
                best_score = score
                best_by_cluster = dict(chosen_by_cluster)
            return

        if score + remaining_upper[depth] <= best_score + 1e-15:
            return

        cluster = cluster_order[depth]
        for candidate in candidates[cluster]:
            if not _passes_seed_separation(
                candidate,
                chosen,
                rmsd_matrix,
                tm_matrix,
                min_rmsd=min_rmsd,
                max_tm=max_tm,
                separation_rule=separation_rule,
            ):
                continue

            chosen.append(candidate)
            chosen_by_cluster[cluster] = candidate
            recurse(depth + 1, score + weights[candidate])
            chosen.pop()
            del chosen_by_cluster[cluster]

    recurse(0, 0.0)

    if best_by_cluster is None:
        # This should not normally occur, but gives a safe deterministic fallback.
        return medoids.copy()

    return np.asarray(
        [best_by_cluster[c] for c in range(n_clusters)], dtype=int
    )


def select_weighted_carbonara_seeds(
    pdb_files,
    weights,
    n_select=4,
    cumulative_weight=0.95,
    compare_func=compare_structures_vals,
    pairwise_results=None,
    use_multiprocessing=False,
    nprocs=None,
    chunksize=20,
    distance_metric="rmsd",
    distance_power=2.0,
    min_rmsd=None,
    max_tm=None,
    separation_rule="both",
    representative="highest_weight",
    n_starts=20,
):
    """Select diverse, SAXS-weighted structures for a new Carbonara run.

    Workflow
    --------
    1. Normalize the SAXS weights.
    2. Retain the smallest high-weight subset containing ``cumulative_weight``.
    3. Compute C-alpha RMSD/TM/GDT matrices for that retained subset.
    4. Run weighted k-medoids subject to a hard minimum-separation rule.
    5. Return either the weighted medoids or the highest-weight feasible member
       of each cluster.

    Parameters
    ----------
    pdb_files : sequence of str/path-like
        The MD snapshots, ordered exactly as ``weights``.
    weights : sequence of float
        Non-negative SAXS/BME weights.
    n_select : int, default 4
        Number of Carbonara seed structures.
    cumulative_weight : float or None, default 0.95
        Retain the smallest descending-weight subset carrying this fraction of
        the full posterior weight. Use None to retain every positive-weight
        snapshot.
    pairwise_results : sequence of dict or None
        Optional precomputed pairwise results for the *full* ``pdb_files`` list.
        If omitted, pairwise metrics are calculated only for retained snapshots.
    distance_metric : {'rmsd', 'tm'}, default 'rmsd'
        Clustering distance. TM uses ``1 - TM``. RMSD is recommended here.
    distance_power : float, default 2
        Power in sum_i w_i min_k d(i,k)^p.
    min_rmsd : float or None
        Require chosen seeds to be at least this many Angstrom apart.
    max_tm : float or None
        Require chosen seeds to have pairwise TM-score no larger than this.
    separation_rule : {'both', 'either'}, default 'both'
        With both thresholds active, 'both' requires every pair to pass both
        tests. This is the conservative redundancy guard.
    representative : {'medoid', 'highest_weight'}, default 'highest_weight'
        Which real snapshot to return from each cluster.

    Returns
    -------
    dict
        Includes ``selected_pdbs``, ``selection_table``, ``assignments``,
        pairwise matrices, medoids, cluster weights and effective sample sizes.
    """
    pdb_files = [str(p) for p in pdb_files]
    weights = np.asarray(weights, dtype=float)
    n = len(pdb_files)

    if len(weights) != n:
        raise ValueError("pdb_files and weights must have the same length.")
    if n == 0:
        raise ValueError("No structures were supplied.")
    if not np.all(np.isfinite(weights)) or np.any(weights < 0):
        raise ValueError("weights must be finite and non-negative.")
    if weights.sum() <= 0:
        raise ValueError("At least one weight must be positive.")
    if representative not in {"medoid", "highest_weight"}:
        raise ValueError("representative must be 'medoid' or 'highest_weight'.")
    if distance_metric not in {"rmsd", "tm"}:
        raise ValueError("distance_metric must be 'rmsd' or 'tm'.")

    weights_full = weights / weights.sum()
    n_eff_full = float(1.0 / np.sum(weights_full**2))

    positive = np.flatnonzero(weights_full > 0)
    weight_order = positive[np.argsort(-weights_full[positive], kind="stable")]

    if cumulative_weight is None:
        retained_indices = weight_order
    else:
        cumulative_weight = float(cumulative_weight)
        if not 0 < cumulative_weight <= 1:
            raise ValueError("cumulative_weight must lie in (0, 1].")
        cumulative = np.cumsum(weights_full[weight_order])
        n_keep = int(np.searchsorted(cumulative, cumulative_weight, side="left") + 1)
        n_keep = max(int(n_select), n_keep)
        n_keep = min(n_keep, len(weight_order))
        retained_indices = weight_order[:n_keep]

    if len(retained_indices) < int(n_select):
        raise ValueError(
            f"Only {len(retained_indices)} positive-weight snapshots remain, "
            f"fewer than n_select={n_select}."
        )

    # Keep retained structures in descending posterior-weight order. This makes
    # output deterministic and keeps local index 0 as the highest-weight member.
    retained_indices = np.asarray(retained_indices, dtype=int)
    retained_pdbs = [pdb_files[i] for i in retained_indices]
    retained_mass = float(weights_full[retained_indices].sum())
    retained_weights = weights_full[retained_indices] / retained_mass
    n_eff_retained = float(1.0 / np.sum(retained_weights**2))

    if pairwise_results is None:
        if use_multiprocessing:
            retained_pairwise = pairwise_structure_metrics_mp(
                retained_pdbs,
                nprocs=nprocs,
                chunksize=chunksize,
            )
        else:
            retained_pairwise = pairwise_structure_metrics(
                retained_pdbs,
                compare_func,
            )
        rmsd, tm, gdt = pairwise_results_to_matrices(
            retained_pdbs, retained_pairwise
        )
    else:
        full_rmsd, full_tm, full_gdt = pairwise_results_to_matrices(
            pdb_files, pairwise_results
        )
        ix = np.ix_(retained_indices, retained_indices)
        rmsd = full_rmsd[ix]
        tm = full_tm[ix]
        gdt = full_gdt[ix]
        retained_pairwise = None

    distance = rmsd if distance_metric == "rmsd" else 1.0 - tm

    medoids_local, labels, objective = _constrained_weighted_kmedoids(
        distance,
        retained_weights,
        n_select=n_select,
        rmsd_matrix=rmsd,
        tm_matrix=tm,
        min_rmsd=min_rmsd,
        max_tm=max_tm,
        separation_rule=separation_rule,
        distance_power=distance_power,
        n_starts=n_starts,
    )

    if representative == "medoid":
        selected_local = medoids_local.copy()
    else:
        selected_local = _highest_weight_representatives_with_separation(
            labels,
            medoids_local,
            retained_weights,
            rmsd,
            tm,
            min_rmsd=min_rmsd,
            max_tm=max_tm,
            separation_rule=separation_rule,
        )

    selected_global = retained_indices[selected_local]
    medoids_global = retained_indices[medoids_local]

    assignment_rows = []
    selection_rows = []

    for cluster in range(int(n_select)):
        members_local = np.flatnonzero(labels == cluster)
        cluster_full_weight = float(weights_full[retained_indices[members_local]].sum())
        cluster_retained_weight = float(retained_weights[members_local].sum())
        medoid = int(medoids_local[cluster])
        selected = int(selected_local[cluster])
        member_rmsd = rmsd[members_local, medoid]

        selection_rows.append({
            "cluster": cluster,
            "selected_index": int(retained_indices[selected]),
            "selected_pdb": pdb_files[int(retained_indices[selected])],
            "selected_weight": float(weights_full[int(retained_indices[selected])]),
            "medoid_index": int(retained_indices[medoid]),
            "medoid_pdb": pdb_files[int(retained_indices[medoid])],
            "medoid_weight": float(weights_full[int(retained_indices[medoid])]),
            "cluster_weight_full": cluster_full_weight,
            "cluster_weight_retained": cluster_retained_weight,
            "n_members": int(len(members_local)),
            "mean_rmsd_to_medoid": float(np.average(
                member_rmsd, weights=retained_weights[members_local]
            )),
            "max_rmsd_to_medoid": float(np.max(member_rmsd)),
        })

        for member in members_local:
            global_i = int(retained_indices[member])
            assignment_rows.append({
                "cluster": cluster,
                "index": global_i,
                "pdb": pdb_files[global_i],
                "weight_full": float(weights_full[global_i]),
                "weight_within_retained": float(retained_weights[member]),
                "rmsd_to_medoid": float(rmsd[member, medoid]),
                "tm_to_medoid": float(tm[member, medoid]),
                "is_medoid": bool(member == medoid),
                "is_selected": bool(member == selected),
            })

    pair_rows = []
    for a in range(int(n_select)):
        for b in range(a + 1, int(n_select)):
            ia = int(selected_local[a])
            ib = int(selected_local[b])
            pair_rows.append({
                "cluster_a": a,
                "cluster_b": b,
                "index_a": int(retained_indices[ia]),
                "index_b": int(retained_indices[ib]),
                "rmsd": float(rmsd[ia, ib]),
                "tm": float(tm[ia, ib]),
                "gdt_ts": float(gdt[ia, ib]),
            })

    selection_table = pd.DataFrame(selection_rows).sort_values(
        "cluster_weight_full", ascending=False
    ).reset_index(drop=True)
    assignments = pd.DataFrame(assignment_rows).sort_values(
        ["cluster", "weight_full"], ascending=[True, False]
    ).reset_index(drop=True)
    selected_pair_metrics = pd.DataFrame(pair_rows)

    return {
        "selected_indices": list(map(int, selected_global)),
        "selected_pdbs": [pdb_files[i] for i in selected_global],
        "selected_weights": [float(weights_full[i]) for i in selected_global],
        "medoid_indices": list(map(int, medoids_global)),
        "medoid_pdbs": [pdb_files[i] for i in medoids_global],
        "retained_indices": list(map(int, retained_indices)),
        "retained_pdbs": retained_pdbs,
        "retained_weight_mass": retained_mass,
        "n_eff_full": n_eff_full,
        "n_eff_retained": n_eff_retained,
        "objective": float(objective),
        "distance_metric": distance_metric,
        "selection_table": selection_table,
        "assignments": assignments,
        "selected_pair_metrics": selected_pair_metrics,
        "rmsd_matrix_retained": rmsd,
        "tm_matrix_retained": tm,
        "gdt_matrix_retained": gdt,
        "pairwise_results_retained": retained_pairwise,
    }




def read_and_align_saxs_weights(
    weights_file,
    pdb_files,
    base_dir=None,
    normalize=True,
    require_all=True,
):
    """
    Read a two-column SAXS weights file and align the weights with an
    independently generated list of PDB files.

    Parameters
    ----------
    weights_file : str or Path
        Text file containing:
            path/to/frame1.pdb   weight
            path/to/frame2.pdb   weight

    pdb_files : sequence of str or Path
        PDB files whose order should be preserved.

    base_dir : str or Path or None
        Directory relative to which paths in the weights file are interpreted.
        If None, the current working directory is used.

    normalize : bool
        Normalize the aligned weights so they sum to one.

    require_all : bool
        If True, raise an error when a PDB has no matching weight.
        If False, unmatched PDBs are omitted.

    Returns
    -------
    aligned_pdbs : list[Path]
        PDB files in the same order as the input `pdb_files`, excluding
        unmatched files when require_all=False.

    aligned_weights : np.ndarray
        Weight corresponding to each returned PDB.

    unmatched_pdbs : list[Path]
        PDB files for which no weight was found.
    """
    weights_file = Path(weights_file)
    base_dir = Path.cwd() if base_dir is None else Path(base_dir)

    # --------------------------------------------------
    # Read the weights file
    # --------------------------------------------------
    weight_by_resolved_path = {}
    weight_by_name = {}

    with weights_file.open("r", encoding="utf-8") as f:
        for line_no, raw_line in enumerate(f, start=1):
            line = raw_line.strip()

            if not line or line.startswith("#"):
                continue

            # Split from the right, so paths containing spaces still work.
            try:
                path_text, weight_text = line.rsplit(maxsplit=1)
                weight = float(weight_text)
            except ValueError as exc:
                raise ValueError(
                    f"Could not parse line {line_no} of {weights_file}:\n"
                    f"{raw_line.rstrip()}"
                ) from exc

            path = Path(path_text).expanduser()

            if not path.is_absolute():
                path = base_dir / path

            resolved = path.resolve()

            if resolved in weight_by_resolved_path:
                raise ValueError(
                    f"Duplicate path in weights file: {resolved}"
                )

            weight_by_resolved_path[resolved] = weight

            # Basename fallback, e.g. frame57.pdb.
            weight_by_name.setdefault(path.name, []).append(
                (resolved, weight)
            )

    # --------------------------------------------------
    # Align weights to the PDB-list order
    # --------------------------------------------------
    aligned_pdbs = []
    aligned_weights = []
    unmatched_pdbs = []

    for pdb in pdb_files:
        pdb = Path(pdb).expanduser()
        resolved_pdb = pdb.resolve()

        # First choice: exact resolved-path match.
        if resolved_pdb in weight_by_resolved_path:
            weight = weight_by_resolved_path[resolved_pdb]

        else:
            # Fallback: match by filename only, provided it is unique.
            basename_matches = weight_by_name.get(pdb.name, [])

            if len(basename_matches) == 1:
                weight = basename_matches[0][1]

            elif len(basename_matches) > 1:
                raise ValueError(
                    f"Ambiguous basename match for {pdb.name}: "
                    f"{len(basename_matches)} entries occur in the weights file."
                )

            else:
                unmatched_pdbs.append(pdb)
                continue

        aligned_pdbs.append(pdb)
        aligned_weights.append(weight)

    if unmatched_pdbs and require_all:
        missing = "\n".join(f"  {p}" for p in unmatched_pdbs)
        raise ValueError(
            f"{len(unmatched_pdbs)} PDB files have no matching SAXS weight:\n"
            f"{missing}"
        )

    aligned_weights = np.asarray(aligned_weights, dtype=float)

    if np.any(~np.isfinite(aligned_weights)):
        raise ValueError("The aligned weights contain NaN or infinite values.")

    if np.any(aligned_weights < 0):
        raise ValueError("The aligned weights contain negative values.")

    if normalize:
        total = aligned_weights.sum()

        if total <= 0:
            raise ValueError("The aligned weights sum to zero.")

        aligned_weights = aligned_weights / total

    return aligned_pdbs, aligned_weights, unmatched_pdbs

# -----------------------------------------------------------------------------
# Export helpers
# -----------------------------------------------------------------------------

def zip_prediction_pdbs(
    predictions,
    zip_name="best_prediction_pdbs.zip",
    output_dir=None,
    fitdata_dir: str | Path | None = None,
    overwrite: bool = True,
    include_manifest: bool = True,
    require_exists: bool = True,
    flat: bool = False,
):
    """
    Put all PDB files referenced by prediction records into a zip archive.

    Parameters
    ----------
    predictions : object
        Usually the output of
        collect_good_prediction_files(..., return_records=True) or
        collect_best_prediction_per_run_closest_to_one(..., return_records=True).
        Legacy inputs are also accepted: a single path, a list of paths,
        (path_or_paths, chi2), or nested mixture path lists.

    zip_name : str or Path
        Name of the zip file to create. If no .zip suffix is supplied, one is
        added. Relative paths are interpreted relative to output_dir, or the
        current working directory if output_dir is None.

    output_dir : str or Path or None
        Directory in which to create the zip file. Defaults to the current
        working directory.

    fitdata_dir : str or Path or None
        Optional fitdata directory. This is useful if predictions are legacy
        path lists and the full mixture records need to be recovered from
        foxs_mixture_results.txt.

    overwrite : bool
        If False, raise FileExistsError if the target zip already exists.

    include_manifest : bool
        If True, add best_predictions_manifest.csv to the zip with component
        metadata: label, chi2, weights, c1/c2 when available, source path, and
        archive name.

    require_exists : bool
        If True, raise FileNotFoundError if any referenced PDB does not exist.
        If False, missing PDBs are skipped and listed in the returned manifest
        rows with status='missing'.

    flat : bool
        If False, keep the default organised archive layout:
            allAtomRunN/label/component_i_filename.pdb

        If True, write PDB files at the top level of the zip, with no
        subfolders. This is useful for tools such as MultiFoXS that expect a
        flat directory of input PDB files. Use include_manifest=False as well
        if the zip should contain only PDB files.

    Returns
    -------
    zip_path : pathlib.Path
        Path to the created zip file.

    rows : list[dict]
        Manifest rows describing what was added or skipped.
    """
    import csv
    import zipfile
    from io import StringIO

    if output_dir is None:
        output_dir = Path.cwd()
    else:
        output_dir = Path(output_dir)

    zip_path = Path(zip_name)
    if zip_path.suffix.lower() != ".zip":
        zip_path = zip_path.with_suffix(".zip")
    if not zip_path.is_absolute():
        zip_path = output_dir / zip_path

    zip_path.parent.mkdir(parents=True, exist_ok=True)
    if zip_path.exists() and not overwrite:
        raise FileExistsError(zip_path)

    records = _normalise_prediction_collection(predictions, fitdata_dir=fitdata_dir)

    rows = []
    used_arcnames = set()
    added_resolved = set()

    def _safe_text(x):
        if x is None:
            return "NA"
        s = str(x)
        s = re.sub(r"[^A-Za-z0-9_.+-]+", "_", s)
        s = s.strip("_")
        return s or "NA"

    def _archive_name_for_pdb(pdb_path: Path, rec, component_i: int, pred_i: int):
        label = _safe_text(rec.get("label", f"prediction_{pred_i}"))
        run_no = rec.get("run_no", rec.get("runNo"))
        if run_no is None:
            meta = _parse_prediction_pdb_metadata(pdb_path)
            run_no = meta.get("run_no") or meta.get("runNo")
        run_part = f"allAtomRun{run_no}" if run_no is not None else f"prediction_{pred_i}"

        stem = pdb_path.stem
        suffix = pdb_path.suffix

        if flat:
            # MultiFoXS-style export: only top-level PDB files, no folders.
            # Keep the original filename unless it would collide with another
            # different PDB in the same archive.
            arc = pdb_path.name
            if arc not in used_arcnames:
                used_arcnames.add(arc)
                return arc

            k = 2
            while True:
                arc2 = f"{stem}_{k}{suffix}"
                if arc2 not in used_arcnames:
                    used_arcnames.add(arc2)
                    return arc2
                k += 1

        # Default organised, collision-resistant layout.
        arc = f"{run_part}/{label}/component_{component_i}_{pdb_path.name}"
        if arc not in used_arcnames:
            used_arcnames.add(arc)
            return arc

        k = 2
        while True:
            arc2 = f"{run_part}/{label}/component_{component_i}_{stem}_{k}{suffix}"
            if arc2 not in used_arcnames:
                used_arcnames.add(arc2)
                return arc2
            k += 1

    with zipfile.ZipFile(zip_path, mode="w", compression=zipfile.ZIP_DEFLATED) as zf:
        for pred_i, rec in enumerate(records):
            pdb_paths = [Path(p) for p in rec.get("pdb_paths", [])]
            weights = rec.get("weights")
            if weights is None or len(weights) != len(pdb_paths):
                weights = [None] * len(pdb_paths)

            for component_i, (pdb_path, weight) in enumerate(zip(pdb_paths, weights)):
                row = {
                    "prediction_i": pred_i,
                    "component_i": component_i,
                    "type": rec.get("type"),
                    "run_no": rec.get("run_no", rec.get("runNo")),
                    "label": rec.get("label"),
                    "chi2": rec.get("chi2"),
                    "weight": weight,
                    "scale": rec.get("scale"),
                    "c1": rec.get("c1"),
                    "c2": rec.get("c2"),
                    "best_component_chi2": rec.get("best_component_chi2"),
                    "source_pdb": str(pdb_path),
                    "archive_name": "",
                    "status": "",
                }

                if not pdb_path.exists():
                    row["status"] = "missing"
                    rows.append(row)
                    if require_exists:
                        raise FileNotFoundError(pdb_path)
                    continue

                try:
                    resolved = pdb_path.resolve()
                except Exception:
                    resolved = pdb_path

                if resolved in added_resolved:
                    row["status"] = "duplicate_skipped"
                    rows.append(row)
                    continue

                arcname = _archive_name_for_pdb(pdb_path, rec, component_i, pred_i)
                zf.write(pdb_path, arcname=arcname)
                added_resolved.add(resolved)
                row["archive_name"] = arcname
                row["status"] = "added"
                rows.append(row)

        if include_manifest:
            fieldnames = [
                "prediction_i", "component_i", "type", "run_no", "label",
                "chi2", "weight", "scale", "c1", "c2", "best_component_chi2",
                "source_pdb", "archive_name", "status",
            ]
            sio = StringIO()
            writer = csv.DictWriter(sio, fieldnames=fieldnames)
            writer.writeheader()
            for row in rows:
                writer.writerow({k: row.get(k, "") for k in fieldnames})
            zf.writestr("best_predictions_manifest.csv", sio.getvalue())

    print(f"Wrote {sum(r['status'] == 'added' for r in rows)} PDB files to {zip_path}")
    return zip_path, rows


# =============================================================================
# Visualisation colour-mode patch
# -----------------------------------------------------------------------------
# These definitions intentionally override the earlier visualisation functions in
# this module.  They keep the old call patterns working, while adding:
#   - colour_mode="chain"      : existing per-chain colouring
#   - colour_mode="model"      : one colour for the whole model
#   - colour_mode="length"     : N-to-C / residue-position spectrum colouring
#   - colour_mode="secondary"  : secondary-structure colouring using 3Dmol's
#                                 ssPyMOL colour scheme
# and n-structure support for visualisePredictionComp.
# =============================================================================

_VIS_PALETTE = [
    "blue", "green", "red", "yellow", "cyan", "magenta",
    "orange", "purple", "lime", "gray",
]


def _normalise_visual_color_mode(color_mode="chain", color_by_chain=None):
    """Normalise visualisation colour-mode aliases."""
    if color_mode is None:
        if color_by_chain is None:
            color_mode = "chain"
        else:
            color_mode = "chain" if color_by_chain else "model"

    mode = str(color_mode).strip().lower().replace("-", "_")
    aliases = {
        "chain": "chain",
        "chains": "chain",
        "by_chain": "chain",
        "model": "model",
        "single": "model",
        "uniform": "model",
        "one": "model",
        "length": "length",
        "residue": "length",
        "residue_index": "length",
        "sequence": "length",
        "spectrum": "length",
        "n_to_c": "length",
        "ntoc": "length",
        "secondary": "secondary",
        "secondary_structure": "secondary",
        "ss": "secondary",
        "ss_pymol": "secondary",
        "sspymol": "secondary",
    }
    if mode not in aliases:
        raise ValueError(
            "Unknown colour mode {!r}. Use one of: 'chain', 'model', "
            "'length', or 'secondary'.".format(color_mode)
        )
    return aliases[mode]


def _set_style_for_viewer(view, selection, style, viewer=None):
    """py3Dmol wrapper so the same helper works inside and outside grids."""
    if viewer is None:
        view.setStyle(selection, style)
    else:
        view.setStyle(selection, style, viewer=viewer)


def _add_model_for_viewer(view, data, fmt, viewer=None):
    if viewer is None:
        view.addModel(data, fmt)
    else:
        view.addModel(data, fmt, viewer=viewer)


def _zoom_for_viewer(view, viewer=None):
    if viewer is None:
        view.zoomTo()
    else:
        view.zoomTo(viewer=viewer)


def _add_label_for_viewer(view, label, style, viewer=None):
    if viewer is None:
        view.addLabel(label, style)
    else:
        view.addLabel(label, style, viewer=viewer)


def _representation_style(rep="cartoon", *, color=None, colorscheme=None, opacity=0.9):
    rep_dict = {}
    if color is not None:
        rep_dict["color"] = color
    if colorscheme is not None:
        rep_dict["colorscheme"] = colorscheme
    if opacity is not None:
        rep_dict["opacity"] = float(opacity)
    return {rep: rep_dict}


def _style_structure_model(
    view,
    model_index,
    structure_data,
    fmt,
    *,
    viewer=None,
    color_mode="chain",
    representation="cartoon",
    opacity=0.9,
    model_color=None,
    palette=None,
):
    """
    Apply one of Carbonara's standard structure colour modes to one model.

    Parameters
    ----------
    color_mode : {'chain', 'model', 'length', 'secondary'}
        chain      - distinct colour per chain, falling back to a single colour
        model      - one colour for the whole model
        length     - 3Dmol spectrum colouring along residue/model order
        secondary  - 3Dmol ssPyMOL secondary-structure colouring
    """
    palette = list(palette or _VIS_PALETTE)
    mode = _normalise_visual_color_mode(color_mode)
    base_sel = {"model": int(model_index)}

    if mode == "secondary":
        _set_style_for_viewer(
            view,
            base_sel,
            _representation_style(representation, colorscheme="ssPyMOL", opacity=opacity),
            viewer=viewer,
        )
        return

    if mode == "length":
        _set_style_for_viewer(
            view,
            base_sel,
            _representation_style(representation, colorscheme="spectrum", opacity=opacity),
            viewer=viewer,
        )
        return

    if mode == "chain":
        chains = _chains_present_in_structure(structure_data, fmt)
        if chains:
            for i, ch in enumerate(chains):
                _set_style_for_viewer(
                    view,
                    {"model": int(model_index), "chain": ch},
                    _representation_style(
                        representation,
                        color=palette[i % len(palette)],
                        opacity=opacity,
                    ),
                    viewer=viewer,
                )
            return

    # mode == 'model', or chain mode with no chain IDs available.
    color = model_color or palette[int(model_index) % len(palette)]
    _set_style_for_viewer(
        view,
        base_sel,
        _representation_style(representation, color=color, opacity=opacity),
        viewer=viewer,
    )


def _prediction_path_list(obj):
    """Return a list of PDB paths from records/lists/legacy tuples, or None."""
    if obj is None:
        return None

    if isinstance(obj, dict):
        paths = obj.get("pdb_paths")
        if paths:
            return [Path(p) for p in paths]
        path = obj.get("pdb_path")
        if path:
            return [Path(path)]
        return []

    if isinstance(obj, tuple) and len(obj) >= 1:
        first = obj[0]
        if isinstance(first, (list, tuple)) and not isinstance(first, (str, bytes, os.PathLike)):
            return [Path(p) for p in first]
        if isinstance(first, (str, bytes, os.PathLike, Path)):
            return [Path(first)]

    if isinstance(obj, (list, tuple)) and not isinstance(obj, (str, bytes, os.PathLike)):
        return [Path(p) for p in obj]

    return None


def visualisePrediction(
    directory,
    runNo,
    predNo=None,
    subNo=0,
    subRun=False,
    color_mode="chain",
    ca_color_mode=None,
):
    """
    Visualise one Carbonara prediction, with AA cartoon plus CA spheres.

    colour modes: 'chain', 'model', 'length', 'secondary'.
    """
    if not HAS_PY3DMOL:
        _warn_missing_py3dmol()
        return None

    view = py3Dmol.view(width=800, height=600)

    if subRun:
        run_dir = os.path.join(directory, f"allAtomRun{runNo}")
    else:
        run_dir = directory

    if predNo is None:
        pred_tag, aa_path, ca_path = _resolve_latest_prediction(
            directory,
            runNo,
            subNo=subNo,
            subRun=subRun,
        )
        print(f"Using latest prediction: {pred_tag}")
    else:
        if predNo == "end":
            aa_fname = f"mol{runNo}_sub_{subNo}_end__AA.pdb"
            ca_fname = f"mol{runNo}_sub_{subNo}_end__CA.pdb"
        else:
            aa_fname = f"mol{runNo}_sub_{subNo}_step_{predNo}__AA.pdb"
            ca_fname = f"mol{runNo}_sub_{subNo}_step_{predNo}__CA.pdb"

        aa_path = os.path.join(run_dir, aa_fname)
        ca_path = os.path.join(run_dir, ca_fname)

        if not os.path.exists(aa_path):
            raise FileNotFoundError(f"AA file not found: {aa_path}")
        if not os.path.exists(ca_path):
            raise FileNotFoundError(f"CA file not found: {ca_path}")

    aa_data, aa_fmt = _read_structure_for_viewer(aa_path)
    _add_model_for_viewer(view, aa_data, aa_fmt)
    _style_structure_model(
        view,
        0,
        aa_data,
        aa_fmt,
        color_mode=color_mode,
        representation="cartoon",
        opacity=0.9,
        model_color="lightgray" if color_mode == "model" else None,
    )

    ca_data, ca_fmt = _read_structure_for_viewer(ca_path)
    _add_model_for_viewer(view, ca_data, ca_fmt)
    _style_structure_model(
        view,
        1,
        ca_data,
        ca_fmt,
        color_mode=(ca_color_mode if ca_color_mode is not None else color_mode),
        representation="sphere",
        opacity=0.45,
        model_color="red",
    )

    _zoom_for_viewer(view)
    view.show()
    return view


def visualisePredictionComparison(
    directory,
    runNo1,
    runNo2,
    predNo1,
    predNo2,
    subNo1,
    subNo2,
    do_superpose=True,
    color_mode="chain",
    reference_color_mode=None,
    mobile_color_mode=None,
):
    """Compare two Carbonara predictions in one viewer."""
    if not HAS_PY3DMOL:
        _warn_missing_py3dmol()
        return None

    view = py3Dmol.view(width=800, height=600)

    aa_fname = f"mol{runNo1}_sub_{subNo1}_step_{predNo1}__AA.pdb"
    mob_fname = f"mol{runNo2}_sub_{subNo2}_step_{predNo2}__AA.pdb"

    aa_path = os.path.join(directory, "allAtomRun" + str(runNo1), aa_fname)
    mob_path = os.path.join(directory, "allAtomRun" + str(runNo2), mob_fname)

    aa_data, aa_fmt = _read_structure_for_viewer(aa_path)

    if do_superpose:
        mob_data_to_show, rmsd, nmatch = superimpose_structure_files_by_ca(aa_path, mob_path)
        mob_fmt = "pdb"
        print(f"Superposed model 1 onto model 0 using {nmatch} matched Cα atoms. RMSD = {rmsd:.3f} Å")
    else:
        mob_data_to_show, mob_fmt = _read_structure_for_viewer(mob_path)

    _add_model_for_viewer(view, aa_data, aa_fmt)
    _add_model_for_viewer(view, mob_data_to_show, mob_fmt)

    _style_structure_model(
        view,
        0,
        aa_data,
        aa_fmt,
        color_mode=reference_color_mode or color_mode,
        representation="cartoon",
        opacity=0.9,
        model_color="lightgray",
    )
    _style_structure_model(
        view,
        1,
        mob_data_to_show,
        mob_fmt,
        color_mode=mobile_color_mode or color_mode,
        representation="cartoon",
        opacity=0.65,
        model_color="red",
    )

    _zoom_for_viewer(view)
    view.show()
    return view


def visualisePredictionMixture(
    pdb_paths,
    weights=None,
    ncols=3,
    panel_width=350,
    panel_height=300,
    show_labels=True,
    color_mode="chain",
):
    """Visualise the component structures of one mixture prediction."""
    if isinstance(pdb_paths, dict):
        rec = pdb_paths
        weights = rec.get("weights") if weights is None else weights
        pdb_paths = rec.get("pdb_paths", [])
    elif isinstance(pdb_paths, tuple) and len(pdb_paths) >= 1:
        first = pdb_paths[0]
        if isinstance(first, (list, tuple)):
            pdb_paths = first

    pdb_paths = [Path(p) for p in pdb_paths]
    if not pdb_paths:
        raise ValueError("No component PDB files supplied for mixture visualisation.")

    if show_labels and weights is not None and len(weights) == len(pdb_paths):
        print("Mixture weights:")
        for p, w in zip(pdb_paths, weights):
            print(f"  {Path(p).name}: {float(w):.4g}")

    return visualisePrediction_panel(
        pdb_paths,
        ncols=ncols,
        panel_width=panel_width,
        panel_height=panel_height,
        show_labels=show_labels,
        color_mode=color_mode,
    )


def visualisePredictionIndividual(aa_path, color_mode="chain"):
    """
    Visualise a single AA PDB, or the component PDBs of a mixture prediction.

    colour modes: 'chain', 'model', 'length', 'secondary'.
    """
    paths = _prediction_path_list(aa_path)
    if paths is not None:
        if isinstance(aa_path, dict):
            if aa_path.get("type") == "mixture" or len(paths) > 1:
                return visualisePredictionMixture(aa_path, color_mode=color_mode)
            aa_path = paths[0]
        elif len(paths) > 1:
            return visualisePredictionMixture(paths, color_mode=color_mode)
        elif len(paths) == 1:
            aa_path = paths[0]

    if not HAS_PY3DMOL:
        _warn_missing_py3dmol()
        return None

    view = py3Dmol.view(width=800, height=600)
    structure_data, fmt = _read_structure_for_viewer(aa_path)
    _add_model_for_viewer(view, structure_data, fmt)
    _style_structure_model(
        view,
        0,
        structure_data,
        fmt,
        color_mode=color_mode,
        representation="cartoon",
        opacity=0.9,
    )

    _zoom_for_viewer(view)
    view.show()
    return view


def visualisePredictionComp(
    pdb1,
    pdb2=None,
    do_superpose=True,
    color_mode="chain",
    reference_color_mode=None,
    mobile_color_mode=None,
    ncols=3,
    panel_width=350,
    panel_height=300,
    max_panels=None,
    show_labels=True,
):
    """
    Compare structures, now including n-structure inputs.

    Backwards compatible:
        visualisePredictionComp(pdb1, pdb2)

    New forms:
        visualisePredictionComp([pdb1, pdb2, ...])
            Show n structures independently.

        visualisePredictionComp(reference_pdb, [mobile1, mobile2, ...])
            Compare each mobile against one reference.

        visualisePredictionComp([mobile1, mobile2, ...], reference_pdb)
            Same as above.

        visualisePredictionComp([ref1, ref2, ...], [mob1, mob2, ...])
            Pairwise comparison panel.
    """
    paths1 = _prediction_path_list(pdb1)
    paths2 = _prediction_path_list(pdb2)

    # One list/record only: just show the n structures.
    if paths1 is not None and pdb2 is None:
        return visualisePrediction_panel(
            paths1,
            ncols=ncols,
            panel_width=panel_width,
            panel_height=panel_height,
            max_panels=max_panels,
            show_labels=show_labels,
            color_mode=color_mode,
        )

    # List + path: compare many mobiles against one reference.
    if paths1 is not None and paths2 is None and isinstance(pdb2, (str, os.PathLike, Path)):
        return visualisePredictionComp_panel(
            paths1,
            pdb2,
            do_superpose=do_superpose,
            ncols=ncols,
            panel_width=panel_width,
            panel_height=panel_height,
            max_panels=max_panels,
            show_labels=show_labels,
            color_mode=color_mode,
            reference_color_mode=reference_color_mode,
            mobile_color_mode=mobile_color_mode,
        )

    # Path + list: compare each mobile in the list against the path reference.
    if paths1 is None and paths2 is not None and isinstance(pdb1, (str, os.PathLike, Path)):
        return visualisePredictionComp_panel(
            paths2,
            pdb1,
            do_superpose=do_superpose,
            ncols=ncols,
            panel_width=panel_width,
            panel_height=panel_height,
            max_panels=max_panels,
            show_labels=show_labels,
            color_mode=color_mode,
            reference_color_mode=reference_color_mode,
            mobile_color_mode=mobile_color_mode,
        )

    if paths1 is not None and paths2 is not None:
        if len(paths1) == 1 and len(paths2) > 1:
            return visualisePredictionComp_panel(
                paths2,
                paths1[0],
                do_superpose=do_superpose,
                ncols=ncols,
                panel_width=panel_width,
                panel_height=panel_height,
                max_panels=max_panels,
                show_labels=show_labels,
                color_mode=color_mode,
                reference_color_mode=reference_color_mode,
                mobile_color_mode=mobile_color_mode,
            )
        if len(paths2) == 1 and len(paths1) > 1:
            return visualisePredictionComp_panel(
                paths1,
                paths2[0],
                do_superpose=do_superpose,
                ncols=ncols,
                panel_width=panel_width,
                panel_height=panel_height,
                max_panels=max_panels,
                show_labels=show_labels,
                color_mode=color_mode,
                reference_color_mode=reference_color_mode,
                mobile_color_mode=mobile_color_mode,
            )
        if len(paths1) == len(paths2) and len(paths1) > 1:
            return visualisePredictionComp_pairs_panel(
                paths1,
                paths2,
                do_superpose=do_superpose,
                ncols=ncols,
                panel_width=panel_width,
                panel_height=panel_height,
                max_panels=max_panels,
                show_labels=show_labels,
                color_mode=color_mode,
                reference_color_mode=reference_color_mode,
                mobile_color_mode=mobile_color_mode,
            )
        if len(paths1) == len(paths2) == 1:
            pdb1 = paths1[0]
            pdb2 = paths2[0]
        else:
            raise ValueError(
                "For n-structure comparison, pass one reference and many mobiles, "
                "or two lists with the same length."
            )

    if pdb2 is None:
        raise ValueError("pdb2 is required unless pdb1 is a list/record of structures.")

    if not HAS_PY3DMOL:
        _warn_missing_py3dmol()
        return None

    view = py3Dmol.view(width=800, height=600)
    data1, fmt1 = _read_structure_for_viewer(pdb1)

    if do_superpose:
        data2_to_show, rmsd, nmatch = superimpose_structure_files_by_ca(pdb1, pdb2)
        fmt2 = "pdb"
        print(f"Superposed model 1 onto model 0 using {nmatch} matched Cα atoms. RMSD = {rmsd:.3f} Å")
    else:
        data2_to_show, fmt2 = _read_structure_for_viewer(pdb2)

    _add_model_for_viewer(view, data1, fmt1)
    _add_model_for_viewer(view, data2_to_show, fmt2)

    _style_structure_model(
        view,
        0,
        data1,
        fmt1,
        color_mode=reference_color_mode or color_mode,
        representation="cartoon",
        opacity=0.9,
        model_color="lightgray",
    )
    _style_structure_model(
        view,
        1,
        data2_to_show,
        fmt2,
        color_mode=mobile_color_mode or color_mode,
        representation="cartoon",
        opacity=0.65,
        model_color="red",
    )

    _zoom_for_viewer(view)
    view.show()
    return view


def visualisePredictionComp_panel(
    file_list,
    reference_file,
    do_superpose=True,
    ncols=3,
    panel_width=350,
    panel_height=300,
    max_panels=None,
    show_labels=True,
    color_mode=None,
    reference_color_mode="model",
    mobile_color_mode="model",
):
    """Show a grid of pairwise comparisons against one fixed reference."""
    if not HAS_PY3DMOL:
        _warn_missing_py3dmol()
        return None

    file_list = list(file_list)
    if max_panels is not None:
        file_list = file_list[:max_panels]

    n = len(file_list)
    if n == 0:
        print("No files to display.")
        return None

    if color_mode is not None:
        reference_color_mode = color_mode
        mobile_color_mode = color_mode

    ncols = max(1, int(ncols))
    nrows = (n + ncols - 1) // ncols

    view = py3Dmol.view(
        viewergrid=(nrows, ncols),
        width=ncols * panel_width,
        height=nrows * panel_height,
        linked=False,
    )

    ref_data, ref_fmt = _read_structure_for_viewer(reference_file)

    for k, mobile_file in enumerate(file_list):
        r = k // ncols
        c = k % ncols
        viewer = (r, c)

        try:
            _add_model_for_viewer(view, ref_data, ref_fmt, viewer=viewer)

            if do_superpose:
                mob_data_to_show, rmsd, nmatch = superimpose_structure_files_by_ca(
                    reference_file,
                    mobile_file,
                )
                mob_fmt = "pdb"
                panel_title = f"{Path(mobile_file).name}\nRMSD={rmsd:.2f} Å"
            else:
                mob_data_to_show, mob_fmt = _read_structure_for_viewer(mobile_file)
                panel_title = Path(mobile_file).name

            _add_model_for_viewer(view, mob_data_to_show, mob_fmt, viewer=viewer)

            _style_structure_model(
                view,
                0,
                ref_data,
                ref_fmt,
                viewer=viewer,
                color_mode=reference_color_mode,
                representation="cartoon",
                opacity=0.75,
                model_color="lightgray",
            )
            _style_structure_model(
                view,
                1,
                mob_data_to_show,
                mob_fmt,
                viewer=viewer,
                color_mode=mobile_color_mode,
                representation="cartoon",
                opacity=0.85,
                model_color="red",
            )

            _zoom_for_viewer(view, viewer=viewer)

            if show_labels:
                _add_label_for_viewer(
                    view,
                    panel_title,
                    {
                        "fontSize": 10,
                        "backgroundColor": "white",
                        "backgroundOpacity": 0.7,
                        "fontColor": "black",
                        "borderThickness": 0,
                        "inFront": True,
                    },
                    viewer=viewer,
                )

        except Exception as e:
            print(f"Failed for {mobile_file}: {e}")
            _add_label_for_viewer(
                view,
                f"Failed:\n{Path(mobile_file).name}",
                {
                    "fontSize": 12,
                    "backgroundColor": "mistyrose",
                    "backgroundOpacity": 0.8,
                    "fontColor": "black",
                    "borderThickness": 0,
                    "inFront": True,
                },
                viewer=viewer,
            )

    view.show()
    return view


def visualisePredictionComp_pairs_panel(
    reference_files,
    mobile_files,
    do_superpose=True,
    ncols=3,
    panel_width=350,
    panel_height=300,
    max_panels=None,
    show_labels=True,
    color_mode=None,
    reference_color_mode="model",
    mobile_color_mode="model",
):
    """Show pairwise comparisons for two equal-length file lists."""
    if not HAS_PY3DMOL:
        _warn_missing_py3dmol()
        return None

    reference_files = list(reference_files)
    mobile_files = list(mobile_files)
    if len(reference_files) != len(mobile_files):
        raise ValueError("reference_files and mobile_files must have the same length.")

    if max_panels is not None:
        reference_files = reference_files[:max_panels]
        mobile_files = mobile_files[:max_panels]

    n = len(reference_files)
    if n == 0:
        print("No files to display.")
        return None

    if color_mode is not None:
        reference_color_mode = color_mode
        mobile_color_mode = color_mode

    ncols = max(1, int(ncols))
    nrows = (n + ncols - 1) // ncols
    view = py3Dmol.view(
        viewergrid=(nrows, ncols),
        width=ncols * panel_width,
        height=nrows * panel_height,
        linked=False,
    )

    for k, (ref_file, mobile_file) in enumerate(zip(reference_files, mobile_files)):
        r = k // ncols
        c = k % ncols
        viewer = (r, c)
        try:
            ref_data, ref_fmt = _read_structure_for_viewer(ref_file)
            _add_model_for_viewer(view, ref_data, ref_fmt, viewer=viewer)

            if do_superpose:
                mob_data_to_show, rmsd, nmatch = superimpose_structure_files_by_ca(
                    ref_file,
                    mobile_file,
                )
                mob_fmt = "pdb"
                panel_title = f"{Path(mobile_file).name}\nRMSD={rmsd:.2f} Å"
            else:
                mob_data_to_show, mob_fmt = _read_structure_for_viewer(mobile_file)
                panel_title = f"{Path(ref_file).name}\nvs {Path(mobile_file).name}"

            _add_model_for_viewer(view, mob_data_to_show, mob_fmt, viewer=viewer)
            _style_structure_model(
                view,
                0,
                ref_data,
                ref_fmt,
                viewer=viewer,
                color_mode=reference_color_mode,
                representation="cartoon",
                opacity=0.75,
                model_color="lightgray",
            )
            _style_structure_model(
                view,
                1,
                mob_data_to_show,
                mob_fmt,
                viewer=viewer,
                color_mode=mobile_color_mode,
                representation="cartoon",
                opacity=0.85,
                model_color="red",
            )
            _zoom_for_viewer(view, viewer=viewer)
            if show_labels:
                _add_label_for_viewer(
                    view,
                    panel_title,
                    {
                        "fontSize": 10,
                        "backgroundColor": "white",
                        "backgroundOpacity": 0.7,
                        "fontColor": "black",
                        "borderThickness": 0,
                        "inFront": True,
                    },
                    viewer=viewer,
                )
        except Exception as e:
            print(f"Failed for {ref_file} / {mobile_file}: {e}")

    view.show()
    return view


def visualisePrediction_panel(
    file_list,
    ncols=3,
    panel_width=350,
    panel_height=300,
    max_panels=None,
    show_labels=True,
    color_by_chain=True,
    color_mode=None,
):
    """
    Display a set of structures side by side in a py3Dmol viewer grid.

    colour modes: 'chain', 'model', 'length', 'secondary'.  The older
    ``color_by_chain`` argument is still accepted; it is used only when
    ``color_mode`` is not supplied.
    """
    if not HAS_PY3DMOL:
        _warn_missing_py3dmol()
        return None

    file_list = list(file_list)
    if max_panels is not None:
        file_list = file_list[:max_panels]

    n = len(file_list)
    if n == 0:
        print("No files to display.")
        return None

    mode = _normalise_visual_color_mode(color_mode, color_by_chain=color_by_chain)

    ncols = max(1, int(ncols))
    nrows = (n + ncols - 1) // ncols

    view = py3Dmol.view(
        viewergrid=(nrows, ncols),
        width=ncols * panel_width,
        height=nrows * panel_height,
        linked=False,
    )

    for k, structure_file in enumerate(file_list):
        r = k // ncols
        c = k % ncols
        viewer = (r, c)

        try:
            structure_data, structure_fmt = _read_structure_for_viewer(structure_file)
            _add_model_for_viewer(view, structure_data, structure_fmt, viewer=viewer)
            _style_structure_model(
                view,
                0,
                structure_data,
                structure_fmt,
                viewer=viewer,
                color_mode=mode,
                representation="cartoon",
                opacity=0.9,
                model_color=_VIS_PALETTE[k % len(_VIS_PALETTE)],
            )

            _zoom_for_viewer(view, viewer=viewer)

            if show_labels:
                _add_label_for_viewer(
                    view,
                    Path(structure_file).name,
                    {
                        "fontSize": 10,
                        "backgroundColor": "white",
                        "backgroundOpacity": 0.7,
                        "fontColor": "black",
                        "borderThickness": 0,
                        "inFront": True,
                    },
                    viewer=viewer,
                )

        except Exception as e:
            print(f"Failed for {structure_file}: {e}")
            _add_label_for_viewer(
                view,
                f"Failed:\n{Path(structure_file).name}",
                {
                    "fontSize": 12,
                    "backgroundColor": "mistyrose",
                    "backgroundOpacity": 0.8,
                    "fontColor": "black",
                    "borderThickness": 0,
                    "inFront": True,
                },
                viewer=viewer,
            )

    view.show()
    return view


def show_structure_and_foxs_side_by_side(
    pdb_name,
    saxs_name,
    foxs_cmd="pyfoxs",
    max_q=None,
    structure_width=480,
    structure_height=420,
    plot_width=520,
    print_summary=False,
    color_mode="chain",
):
    """
    Display a structure and its FoXS fit side by side.

    The structure viewer now accepts ``color_mode``: 'chain', 'model',
    'length', or 'secondary'.
    """
    pdb_path = Path(pdb_name).resolve()
    saxs_path = Path(saxs_name).resolve()

    if not pdb_path.exists():
        raise FileNotFoundError(f"Structure file not found: {pdb_path}")
    if not saxs_path.exists():
        raise FileNotFoundError(f"SAXS file not found: {saxs_path}")

    temp_pdb_to_clean = None

    try:
        pdb_for_foxs, temp_pdb_to_clean = CDT._convert_cif_to_pdb_for_foxs(pdb_path)
        pdb_for_foxs = Path(pdb_for_foxs)

        base_cmd = CDT._normalise_cmd(foxs_cmd)
        cmd = base_cmd + [str(pdb_for_foxs), str(saxs_path)]

        if max_q is not None:
            cmd += ["--max_q", str(max_q)]

        proc = subprocess.run(cmd, capture_output=True, text=True)
        stdout = proc.stdout
        stderr = proc.stderr

        if proc.returncode != 0:
            raise RuntimeError(
                f"pyFoXS failed with exit code {proc.returncode}\n\nSTDERR:\n{stderr}\n\nSTDOUT:\n{stdout}"
            )

        combined_text = stdout + "\n" + stderr
        chi2 = None
        chi_patterns = [
            r"Chi(?:\^?2| square)\s*[:=]\s*([0-9.eE+-]+)",
            r"chi(?:\^?2| square)\s*[:=]\s*([0-9.eE+-]+)",
            r"\bchi\s*=\s*([0-9.eE+-]+)",
            r"\bChi\s*=\s*([0-9.eE+-]+)",
            r"\bchi2\s*[:=]\s*([0-9.eE+-]+)",
            r"\bChi2\s*[:=]\s*([0-9.eE+-]+)",
        ]
        for pat in chi_patterns:
            m = re.search(pat, combined_text)
            if m:
                try:
                    chi2 = float(m.group(1))
                    break
                except ValueError:
                    pass

        if print_summary:
            print(stdout)
            if stderr:
                print(stderr)

        fit_file = CDT._find_foxs_fit_file(pdb_for_foxs, saxs_path)
        if fit_file is None:
            raise FileNotFoundError("Could not identify a FoXS fit file automatically.")

        fit = CDT.load_numeric_table_loose(fit_file, min_cols=2)
        q = fit[:, 0]

        if fit.shape[1] >= 4:
            i_exp = fit[:, 1]
            sigma = fit[:, 2]
            i_fit = fit[:, 3]
            residual = (i_exp - i_fit) / sigma
        elif fit.shape[1] == 3:
            i_exp = fit[:, 1]
            i_fit = fit[:, 2]
            residual = i_exp - i_fit
        else:
            raise ValueError("Fit file must have at least 3 columns for residual plotting.")

        fig = plt.figure(figsize=(6.0, 6.0))
        gs = fig.add_gridspec(2, 1, height_ratios=[3, 1], hspace=0.08)
        ax1 = fig.add_subplot(gs[0])
        ax2 = fig.add_subplot(gs[1], sharex=ax1)
        ax1.plot(q, i_exp, "o", ms=4, label="Experimental")
        ax1.plot(q, i_fit, "-", lw=2, label="FoXS fit")
        ax1.set_yscale("log")
        ax1.set_ylabel("Intensity")
        title = "Initial FoXS check"
        if chi2 is not None:
            title += f"  (chi² = {chi2:.4g})"
        if max_q is not None:
            title += f", max_q={max_q}"
        ax1.set_title(title)
        ax1.legend()
        ax1.tick_params(axis="x", labelbottom=False)
        ax2.axhline(0.0, lw=1)
        ax2.plot(q, residual, "o", ms=3)
        ax2.set_xlabel("q")
        ax2.set_ylabel("Residual")
        plt.tight_layout()

        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=160, bbox_inches="tight")
        plt.close(fig)
        buf.seek(0)
        plot_b64 = base64.b64encode(buf.read()).decode("utf-8")
        plot_html = f'<img src="data:image/png;base64,{plot_b64}" style="width:{plot_width}px; max-width:100%;">'

        if not HAS_PY3DMOL:
            _warn_missing_py3dmol()
            viewer_html = "<div style='padding:20px;border:1px solid #ddd;border-radius:6px;'>py3Dmol is not available.</div>"
        else:
            view = py3Dmol.view(width=structure_width, height=structure_height)
            structure_data, fmt = _read_structure_for_viewer(pdb_name)
            _add_model_for_viewer(view, structure_data, fmt)
            _style_structure_model(
                view,
                0,
                structure_data,
                fmt,
                color_mode=color_mode,
                representation="cartoon",
                opacity=0.9,
            )
            _zoom_for_viewer(view)
            silent_out = io.StringIO()
            silent_err = io.StringIO()
            with contextlib.redirect_stdout(silent_out), contextlib.redirect_stderr(silent_err):
                viewer_html = view._make_html()

        html = f"""
        <div style="display:flex; flex-wrap:wrap; gap:20px; align-items:flex-start; margin-top:10px; margin-bottom:10px;">
            <div style="flex:0 0 auto;">
                <div style="font-weight:600; margin-bottom:8px;">Structure</div>
                {viewer_html}
            </div>
            <div style="flex:0 0 auto;">
                <div style="font-weight:600; margin-bottom:8px;">Initial FoXS fit</div>
                {plot_html}
            </div>
        </div>
        """
        display(HTML(html))

        return {"chi2": chi2, "stdout": stdout, "stderr": stderr, "fit_file": fit_file, "view_html": viewer_html}

    finally:
        if temp_pdb_to_clean is not None:
            try:
                os.remove(temp_pdb_to_clean)
            except OSError:
                pass


def show_prediction_record_and_foxs_side_by_side(
    prediction,
    saxs_name=None,
    foxs_cmd="pyfoxs",
    max_q=None,
    fitdata_dir: str | Path | None = None,
    structure_width=480,
    structure_height=420,
    plot_width=520,
    color_mode="chain",
):
    """Display a prediction plus its SAXS fit, with the new colour modes."""
    rec = _normalise_prediction_record(prediction, fitdata_dir=fitdata_dir)

    if rec.get("type") != "mixture" and len(rec.get("pdb_paths", [])) <= 1:
        pdb = rec.get("pdb_path") or rec.get("pdb_paths", [None])[0]
        if rec.get("fit_file") is None:
            if saxs_name is None:
                raise ValueError("saxs_name is required to rerun FoXS for a single-PDB prediction.")
            return show_structure_and_foxs_side_by_side(
                pdb,
                saxs_name,
                foxs_cmd=foxs_cmd,
                max_q=max_q,
                structure_width=structure_width,
                structure_height=structure_height,
                plot_width=plot_width,
                color_mode=color_mode,
            )

    pdbs = [Path(p) for p in rec.get("pdb_paths", [])]
    if not pdbs:
        raise ValueError("Prediction record contains no PDB paths.")

    fit_file = rec.get("fit_file")
    if fit_file is None:
        raise ValueError(
            "No stored mixture fit curve found for this prediction. "
            "Use collect_good_prediction_files(..., return_records=True) so the "
            "record includes fit_file, or pass fitdata_dir to recover it."
        )

    if not HAS_PY3DMOL:
        _warn_missing_py3dmol()
        viewer_html = "<div style='padding:20px;border:1px solid #ddd;border-radius:6px;'>py3Dmol is not available.</div>"
    else:
        n = len(pdbs)
        ncols = min(3, max(1, n))
        nrows = (n + ncols - 1) // ncols
        view = py3Dmol.view(
            viewergrid=(nrows, ncols),
            width=ncols * structure_width,
            height=nrows * structure_height,
            linked=False,
        )
        weights = rec.get("weights")

        for k, pdb in enumerate(pdbs):
            r = k // ncols
            c = k % ncols
            viewer = (r, c)
            structure_data, fmt = _read_structure_for_viewer(pdb)
            _add_model_for_viewer(view, structure_data, fmt, viewer=viewer)
            _style_structure_model(
                view,
                0,
                structure_data,
                fmt,
                viewer=viewer,
                color_mode=color_mode,
                representation="cartoon",
                opacity=0.9,
                model_color=_VIS_PALETTE[k % len(_VIS_PALETTE)],
            )
            label = Path(pdb).name
            if weights is not None and k < len(weights):
                label += f"\nw={float(weights[k]):.3g}"
            _add_label_for_viewer(
                view,
                label,
                {"fontSize": 10, "backgroundColor": "white", "backgroundOpacity": 0.7, "fontColor": "black", "borderThickness": 0, "inFront": True},
                viewer=viewer,
            )
            _zoom_for_viewer(view, viewer=viewer)

        silent_out = io.StringIO()
        silent_err = io.StringIO()
        with contextlib.redirect_stdout(silent_out), contextlib.redirect_stderr(silent_err):
            viewer_html = view._make_html()

    curve = load_foxs_fit_curve(fit_file, max_q=max_q)
    fig = plt.figure(figsize=(6.0, 6.0))
    gs = fig.add_gridspec(2, 1, height_ratios=[3, 1], hspace=0.08)
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1], sharex=ax1)
    ax1.plot(curve["q"], curve["i_exp"], "o", ms=4, label="Experimental")
    ax1.plot(curve["q"], curve["i_fit"], "-", lw=2, label="Approx. MultiFoXS fit")
    ax1.set_yscale("log")
    ax1.set_ylabel("Intensity")
    title = f"Approx. MultiFoXS fit: {rec.get('label', '')}"
    if rec.get("chi2") is not None:
        title += f"  (chi² = {float(rec['chi2']):.4g})"
    ax1.set_title(title)
    ax1.legend()
    ax1.tick_params(axis="x", labelbottom=False)
    ax2.axhline(0.0, lw=1)
    ax2.plot(curve["q"], curve["residual"], "o", ms=3)
    ax2.set_xlabel("q")
    ax2.set_ylabel("Residual")
    plt.tight_layout()
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=160, bbox_inches="tight")
    plt.close(fig)
    buf.seek(0)
    plot_b64 = base64.b64encode(buf.read()).decode("utf-8")
    plot_html = f'<img src="data:image/png;base64,{plot_b64}" style="width:{plot_width}px; max-width:100%;">'

    html = f"""
    <div style="display:flex; flex-wrap:wrap; gap:20px; align-items:flex-start; margin-top:10px; margin-bottom:10px;">
        <div style="flex:0 0 auto;">
            <div style="font-weight:600; margin-bottom:8px;">Mixture component structures</div>
            {viewer_html}
        </div>
        <div style="flex:0 0 auto;">
            <div style="font-weight:600; margin-bottom:8px;">Approx. MultiFoXS fit</div>
            {plot_html}
        </div>
    </div>
    """
    display(HTML(html))