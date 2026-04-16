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


def visualisePrediction(directory, runNo, predNo=None, subNo=0,subRun=False):
    if not HAS_PY3DMOL:
        _warn_missing_py3dmol()
        return None
    view = py3Dmol.view(width=800, height=600)
    if subRun:
        run_dir = os.path.join(directory, f"allAtomRun{runNo}")
    else:
        run_dir =directory
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
    with open(aa_path, "r") as f:
        pdb_data_aa = f.read()
    view.addModel(pdb_data_aa, "pdb")   # model 0
    aa_chains = _chains_present_in_pdb(pdb_data_aa)

    # ---- load CA model ----
    with open(ca_path, "r") as f:
        pdb_data_ca = f.read()
    view.addModel(pdb_data_ca, "pdb")   # model 1
    ca_chains = _chains_present_in_pdb(pdb_data_ca)

    # ---- choose colors ----
    palette = ["blue", "green", "red", "yellow", "cyan", "magenta",
               "orange", "purple", "lime", "gray"]

    # Style AA chains (cartoon)
    for i, ch in enumerate(aa_chains):
        color = palette[i % len(palette)]
        view.setStyle({"model": 0, "chain": ch}, {"cartoon": {"color": color}})

    # Style CA chains (spheres)
    for i, ch in enumerate(ca_chains):
        color = palette[i % len(palette)]
        view.setStyle({"model": 1, "chain": ch}, {"sphere": {"color": color, "opacity": 0.5}})

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

    with open(aa_path, "r") as f:
        pdb_data_aa = f.read()

    with open(ca_path, "r") as f:
        pdb_data_ca = f.read()

    # --- superimpose CA model onto AA model for fair visual comparison ---
    if do_superpose:
        pdb_data_ca_aln, rmsd, nmatch = superimpose_pdb_strings_by_ca(pdb_data_aa, pdb_data_ca)
        pdb_data_ca_to_show = pdb_data_ca_aln
        print(f"Superposed model 1 onto model 0 using {nmatch} matched Cα atoms. RMSD = {rmsd:.3f} Å")
    else:
        pdb_data_ca_to_show = pdb_data_ca

    # Add models
    view.addModel(pdb_data_aa, "pdb")          # model 0 (reference)
    view.addModel(pdb_data_ca_to_show, "pdb")  # model 1 (mobile/aligned)

    aa_chains = _chains_present_in_pdb(pdb_data_aa)
    ca_chains = _chains_present_in_pdb(pdb_data_ca_to_show)

    palette = ["blue", "green", "red", "yellow", "cyan", "magenta", "orange", "purple", "lime", "gray"]

    # Style AA chains (cartoon)
    for i, ch in enumerate(aa_chains):
        color = palette[i % len(palette)]
        view.setStyle({"model": 0, "chain": ch}, {"cartoon": {"color": color}})

    # Style aligned CA model (spheres)
    for i, ch in enumerate(ca_chains):
        color = palette[i % len(palette)]
        view.setStyle({"model": 1, "chain": ch}, {"sphere": {"color": color, "opacity": 0.5}})

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
    Read Cα coordinates from PDB or mmCIF using MDTraj.

    MDTraj's traj.xyz is *normally always in nm*.
    This function can optionally auto-detect and then force Å.

    unit: "auto" | "nm" | "A"
      - "auto": infer from CA-CA spacing
      - "nm": treat coords as nm
      - "A":  treat coords as Å
    force_angstrom: if True, returns coords in Å.
    """
    traj = md.load(pdb_or_cif)
    top = traj.topology
    xyz = traj.xyz[0]  # MDTraj: nm

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
            ch = (chain.chain_id or "").strip()
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


def pairwise_structure_metrics(pdb_files, compare_func):
    """
    Compute pairwise RMSD / TM / GDT for all unique structure pairs.

    Parameters
    ----------
    pdb_files : list of str
        Paths to PDB files.
    compare_func : callable
        Function like compare_structures_vals(pdb1, pdb2).

    Returns
    -------
    results : list of dict
        Each entry contains indices, filenames, and metrics.
    """
    n = len(pdb_files)
    n_pairs = n * (n - 1) // 2

    results = []

    with tqdm(total=n_pairs, desc="Pairwise comparisons", unit="pair") as pbar:
        for i in range(n):
            for j in range(i + 1, n):
                p1 = pdb_files[i]
                p2 = pdb_files[j]

                rmsd, tm, gdt = compare_func(p1, p2)

                results.append({
                    "i": i,
                    "j": j,
                    "pdb1": p1,
                    "pdb2": p2,
                    "rmsd": float(rmsd),
                    "tm": float(tm),
                    "gdt_ts": float(gdt),
                })

                pbar.update(1)

    return results

def structure_metrics_vs_reference(pdb_files, ref_pdb, compare_func):
    """
    Compute RMSD / TM / GDT for each structure in `pdb_files`
    against a single reference structure `ref_pdb`.

    Parameters
    ----------
    pdb_files : list of str
        Paths to PDB files to compare.
    ref_pdb : str
        Path to the reference PDB file.
    compare_func : callable
        Function like compare_structures_vals(pdb1, pdb2),
        returning (rmsd, tm, gdt).

    Returns
    -------
    results : list of dict
        Each entry contains index, filenames, and metrics.
    """
    results = []

    with tqdm(total=len(pdb_files), desc="Comparisons vs reference", unit="pdb") as pbar:
        for i, pdb in enumerate(pdb_files):
            rmsd, tm, gdt = compare_func(pdb, ref_pdb)

            results.append({
                "i": i,
                "pdb": pdb,
                "ref_pdb": ref_pdb,
                "rmsd": float(rmsd),
                "tm": float(tm),
                "gdt_ts": float(gdt),
            })

            pbar.update(1)

    return results

def structure_metrics_vs_carbonara(pdb_files, carbonara_dir, compare_func):
    """
    Compare each PDB in `pdb_files` against Carbonara coordinate files
    found in `carbonara_dir`.

    Parameters
    ----------
    pdb_files : list of str
        Paths to PDB files to compare.
    carbonara_dir : str
        Directory containing coordinates*.dat files.
    compare_func : callable
        Function like compare_structures_vals_carbonara(pdb, dat).

    Returns
    -------
    results : list of dict
    """

    coord_files = sorted(glob.glob(os.path.join(carbonara_dir, "coordinates*.dat")))

    if not coord_files:
        raise ValueError(f"No coordinates*.dat files found in {carbonara_dir}")

    total = len(pdb_files) * len(coord_files)
    results = []

    with tqdm(total=total, desc="Comparisons vs Carbonara", unit="pair") as pbar:

        for dat in coord_files:
            for i, pdb in enumerate(pdb_files):

                rmsd, tm, gdt = compare_func(pdb, dat)

                results.append({
                    "i": i,
                    "pdb": pdb,
                    "carbonara_dat": dat,
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

def pairwise_structure_metrics_mp(pdb_files, nprocs=None, chunksize=20, start_method=None):
    """
    Parallel pairwise RMSD/TM/GDT over all unique pairs.

    Parameters
    ----------
    pdb_files : list[str]
    nprocs : int | None
        Defaults to os.cpu_count().
    chunksize : int
        Increase to reduce overhead; 20–200 often good.
    start_method : {"fork","spawn","forkserver", None}
        On Linux, "fork" is usually fastest. If None, use "fork" on Linux else default.

    Returns
    -------
    results : list[dict]
    """
    n = len(pdb_files)
    n_pairs = n * (n - 1) // 2
    if n_pairs == 0:
        return []

    if nprocs is None:
        nprocs = os.cpu_count() or 1

    if start_method is None:
        # Best default on Linux for this kind of numeric work
        start_method = "fork" if os.name == "posix" else "spawn"

    ctx = mp.get_context(start_method)

    # Important: pass pdb_files once (as part of each task args). For 200 files this is fine.
    tasks = ((i, j, pdb_files) for i, j in _pair_indices(n))

    results = []
    with ctx.Pool(processes=nprocs) as pool:
        it = pool.imap_unordered(_worker_compare_pair, tasks, chunksize=chunksize)
        for i, j, rmsd, tm, gdt in tqdm(it, total=n_pairs, desc="Pairwise comparisons", unit="pair"):
            results.append({
                "i": i, "j": j,
                "pdb1": pdb_files[i],
                "pdb2": pdb_files[j],
                "rmsd": rmsd,
                "tm": tm,
                "gdt_ts": gdt,
            })

    # Optional: sort for deterministic order
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
):
    """
    Calculate radius of gyration values from a flat list of PDB files.

    Parameters
    ----------
    pdb_files : list of str
        Flat list of PDB paths. These may include multiple sub-structures
        per prediction, e.g.
            mol5_sub_0_step_10__AA.pdb
            mol5_sub_1_step_10__AA.pdb
            mol5_sub_2_step_10__AA.pdb
    rg_func : callable
        Function like radius_of_gyration(pdb_path) -> float
    weighted : bool, optional
        If False, return one result per pdb file.
        If True, group pdbs by (runNo, predTag) and return one weighted
        result per prediction.
    mixtures : dict or list, optional
        Required if weighted=True.

        Supported forms:
        - dict keyed by (runNo, predTag)
        - dict keyed by runNo
        - list aligned with sorted grouped predictions

    keep_components : bool, optional
        If weighted=True, include component pdbs, component Rg values,
        and weights in the returned records.

    Returns
    -------
    results : list of dict
        Unweighted mode:
            {
                "i": i,
                "pdb": pdb,
                "rg": float(...)
            }

        Weighted mode:
            {
                "i": i,
                "runNo": runNo,
                "predTag": predTag,
                "rg": float(weighted_rg),
                ...
            }

        In both cases, the plottable quantity is always under key "rg".
    """

    pat_step = re.compile(r"mol(\d+)_sub_(\d+)_step_(\d+)__AA\.pdb$")
    pat_end  = re.compile(r"mol(\d+)_sub_(\d+)_end__AA\.pdb$")

    def parse_pdb_name(path):
        fname = os.path.basename(path)

        m = pat_step.match(fname)
        if m:
            runNo = int(m.group(1))
            subNo = int(m.group(2))
            predTag = f"step_{int(m.group(3))}"
            return runNo, subNo, predTag

        m = pat_end.match(fname)
        if m:
            runNo = int(m.group(1))
            subNo = int(m.group(2))
            predTag = "end"
            return runNo, subNo, predTag

        raise ValueError(f"Filename does not match expected pattern: {fname}")

    results = []

    # --------------------------------------------------
    # Unweighted mode: one output per pdb
    # --------------------------------------------------
    if not weighted:
        with tqdm(total=len(pdb_files), desc="Calculating Rg", unit="pdb") as pbar:
            for i, pdb in enumerate(pdb_files):
                rg = float(rg_func(pdb))

                results.append({
                    "i": i,
                    "pdb": pdb,
                    "rg": rg,
                })

                pbar.update(1)

        return results

    # --------------------------------------------------
    # Weighted mode: one output per grouped prediction
    # --------------------------------------------------
    if mixtures is None:
        raise ValueError("mixtures must be provided when weighted=True")

    grouped = defaultdict(list)
    for pdb in pdb_files:
        runNo, subNo, predTag = parse_pdb_name(pdb)
        grouped[(runNo, predTag)].append((subNo, pdb))

    grouped_sorted = {}
    for key, vals in grouped.items():
        grouped_sorted[key] = [pdb for subNo, pdb in sorted(vals, key=lambda x: x[0])]

    group_keys = sorted(grouped_sorted.keys(), key=lambda x: (x[0], x[1]))

    with tqdm(total=len(group_keys), desc="Calculating weighted Rg", unit="pred") as pbar:
        for i, key in enumerate(group_keys):
            runNo, predTag = key
            pdb_group = grouped_sorted[key]

            rg_components = [float(rg_func(pdb)) for pdb in pdb_group]

            if isinstance(mixtures, dict):
                if key in mixtures:
                    weights = mixtures[key]
                elif runNo in mixtures:
                    weights = mixtures[runNo]
                else:
                    raise KeyError(f"No mixture weights found for group {key}")
            else:
                try:
                    weights = mixtures[i]
                except IndexError:
                    raise IndexError(f"No mixture weights supplied for group {key}")

            if len(weights) != len(rg_components):
                raise ValueError(
                    f"Weight length mismatch for group {key}: "
                    f"{len(weights)} weights but {len(rg_components)} pdbs"
                )

            weighted_rg = float(sum(r * w for r, w in zip(rg_components, weights)))

            rec = {
                "i": i,
                "runNo": runNo,
                "predTag": predTag,
                "rg": weighted_rg,   # <- crucial: histogram can use "rg"
            }

            if keep_components:
                rec.update({
                    "pdbs": pdb_group,
                    "rg_components": rg_components,
                    "weights": list(weights),
                })

            results.append(rec)
            pbar.update(1)

    return results


#########################################################
#
#   For the live update fornt end this will grab all the "good" predictions (user specified chi squared)
#
#########################################################
    

def collect_good_prediction_files(
    fitdata_dir: str | Path,
    chi2_threshold: float,
    require_exists: bool = True,
    sort_by_chi2: bool = True,
) -> List[Tuple[Path, float]]:
    """
    Collect AA PDB files whose FoXS chi^2 is <= chi2_threshold.

    Parameters
    ----------
    fitdata_dir : str | Path
        Path to the fitdata directory containing allAtomRun*/foxs_results.txt
    chi2_threshold : float
        Maximum chi^2 to accept.
    require_exists : bool
        If True, only return PDB paths that currently exist on disk.
    sort_by_chi2 : bool
        If True, sort results by increasing chi^2.

    Returns
    -------
    List[Tuple[Path, float]]
        List of (pdb_path, chi2) tuples.
    """
    fitdata_dir = Path(fitdata_dir)
    good = []

    for run_dir in sorted(fitdata_dir.glob("allAtomRun*")):
        if not run_dir.is_dir():
            continue

        foxs_file = run_dir / "foxs_results.txt"
        if not foxs_file.exists():
            continue

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

            try:
                chi2 = float(chi_str)
            except ValueError:
                continue

            if chi2 <= chi2_threshold:
                pdb_path = Path(pdb_path_str)
                if (not require_exists) or pdb_path.exists():
                    good.append((pdb_path, chi2))

    if sort_by_chi2:
        good.sort(key=lambda x: x[1])

    return good




############################################

## Visulaisation routine for a specific path, used in live notebook

############################################


def visualisePredictionIndividual(aa_path):
    if not HAS_PY3DMOL:
        _warn_missing_py3dmol()
        return None
    view = py3Dmol.view(width=800, height=600)
    
    

    # ---- load AA model ----
    with open(aa_path, "r") as f:
        pdb_data_aa = f.read()
    view.addModel(pdb_data_aa, "pdb")   # model 0
    aa_chains = _chains_present_in_pdb(pdb_data_aa)


    # ---- choose colors ----
    palette = ["blue", "green", "red", "yellow", "cyan", "magenta",
               "orange", "purple", "lime", "gray"]

    # Style AA chains (cartoon)
    for i, ch in enumerate(aa_chains):
        color = palette[i % len(palette)]
        view.setStyle({"model": 0, "chain": ch}, {"cartoon": {"color": color}})

    view.zoomTo()
    view.show()


def visualisePredictionComp(pdb1, pdb2, do_superpose=True):
    if not HAS_PY3DMOL:
        _warn_missing_py3dmol()
        return None

    view = py3Dmol.view(width=800, height=600)

    with open(pdb1, "r") as f:
        pdb_data_aa = f.read()

    with open(pdb2, "r") as f:
        pdb_data_ca = f.read()

    # --- superimpose CA model onto AA model for fair visual comparison ---
    if do_superpose:
        pdb_data_ca_aln, rmsd, nmatch = superimpose_pdb_strings_by_ca(pdb_data_aa, pdb_data_ca)
        pdb_data_ca_to_show = pdb_data_ca_aln
        print(f"Superposed model 1 onto model 0 using {nmatch} matched Cα atoms. RMSD = {rmsd:.3f} Å")
    else:
        pdb_data_ca_to_show = pdb_data_ca

    # Add models
    view.addModel(pdb_data_aa, "pdb")          # model 0 (reference)
    view.addModel(pdb_data_ca_to_show, "pdb")  # model 1 (mobile/aligned)

    aa_chains = _chains_present_in_pdb(pdb_data_aa)
    ca_chains = _chains_present_in_pdb(pdb_data_ca_to_show)

    palette = ["blue", "green", "red", "yellow", "cyan", "magenta",
               "orange", "purple", "lime", "gray"]

    # Style model 0
    for i, ch in enumerate(aa_chains):
        color = palette[i % len(palette)]
        view.setStyle({"model": 0, "chain": ch}, {"cartoon": {"color": color}})

    # Style model 1
    for i, ch in enumerate(ca_chains):
        color = palette[(i + 1) % len(palette)]
        view.setStyle({"model": 1, "chain": ch}, {"cartoon": {"color": color}})

    view.zoomTo()
    view.show()
    return view

