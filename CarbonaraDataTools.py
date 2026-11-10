'''
------------------------------------------------------------------------------------------------
Weclome to Carbonara's Data Tools (CDT)!
------------------------------------------------------------------------------------------------

This package provides processing for Carbonara specific data input/outputs
and generating useful inputs for external programs 

------------------------------------------------------------------------------------------------
'''



# -----------------------------------------------------------------------------
# Assemble a mixture from independently prepared Carbonara runs
# -----------------------------------------------------------------------------

def merge_carbonara_setups(
    source_dirs,
    name,
    root_dir=None,
    *,
    make_runme=True,
    max_mixture_combos=30,
    mixture_step=0.1,
    mixture_snap=0.05,
    mixture_dirichlet_alpha=1.0,
    run_options=None,
):
    """Merge existing, independently prepared Carbonara components into a new run.

    Each entry of ``source_dirs`` is either a folder containing a prepared
    Carbonara run, or ``(folder, component_number)`` to choose another numbered
    component from an existing multi-structure run.  Nothing is reprocessed:
    fingerprints, coordinates, flexible-section selections and restraints are
    copied byte-for-byte and only the filename indices change.

    The output is ``root_dir/carbonara_runs/name`` (root_dir defaults to cwd).
    If ``make_runme=True``, reuse ``setup_carbonara_allAtom.write_runme`` to make
    a launcher alongside ``carbonara_runs``.  Because that existing generator
    uses the current directory as its data root, ``root_dir`` must then equal
    the current working directory.  Use ``make_runme=False`` to prepare data
    under another root without generating the launcher.

    Components need not contain the same protein or the same number of chains.
    If at least one supplies fixed distances, empty constraint files represent
    components without restraints (the manifest distinguishes these).
    One experimental SAXS file is used for the fit; differences are warned about.

    Example::

        merge_carbonara_setups(
            ['carbonara_runs/monomer', 'carbonara_runs/dimer'],
            'monomer_dimer',
            run_options={'fit_n_times': 10, 'max_q': 0.2},
        )

    Returns the Path to the new run directory.  Refuses to overwrite a run.
    """
    import hashlib
    import importlib.util
    import warnings

    if not isinstance(name, str) or not name or name in {'.', '..'} or Path(name).name != name:
        raise ValueError('name must be a nonempty run-folder name, without path separators')
    source_dirs = list(source_dirs)
    if not source_dirs:
        raise ValueError('Supply at least two independently prepared Carbonara components')
    if len(source_dirs) < 2:
        raise ValueError('Mixture merging requires at least two components')

    root = Path(root_dir if root_dir is not None else os.getcwd()).expanduser().resolve()
    target = root / 'carbonara_runs' / name
    launch_path = root / ('RunMe_' + name + '.sh')
    if target.exists():
        raise FileExistsError('Refusing to overwrite existing run: ' + str(target))
    if make_runme:
        if root != Path.cwd().resolve():
            raise ValueError('For make_runme=True, run from root_dir (or use make_runme=False)')
        if launch_path.exists():
            raise FileExistsError('Refusing to overwrite existing launcher: ' + str(launch_path))

    # Preflight every component before creating anything.  In particular, a
    # missing numbered file must never silently fall back to component 1.
    components = []
    for chosen in source_dirs:
        if isinstance(chosen, (tuple, list)) and len(chosen) == 2:
            folder, original_idx = chosen
        else:
            folder, original_idx = chosen, 1
        folder = Path(folder).expanduser().resolve()
        original_idx = int(original_idx)
        if original_idx < 1 or not folder.is_dir():
            raise ValueError('Invalid source directory or component number: ' + str(chosen))
        files = {key: folder / (prefix + str(original_idx) + '.dat')
                 for key, prefix in (
                     ('coordinates', 'coordinates'),
                     ('fingerprint', 'fingerPrint'),
                     ('flexibility', 'varyingSectionSecondary'),
                 )}
        files['constraints'] = folder / ('fixedDistanceConstraints' + str(original_idx) + '.dat')
        files['saxs'] = folder / 'Saxs.dat'
        files['chain_lengths'] = (
            folder / ('chainLengths' + str(original_idx) + '.dat')
            if (folder / ('chainLengths' + str(original_idx) + '.dat')).is_file()
            else folder / 'chainLengths.dat'
        )
        for key in ('coordinates', 'fingerprint', 'flexibility', 'saxs', 'chain_lengths'):
            if not files[key].is_file():
                raise FileNotFoundError('Missing ' + key + ' input: ' + str(files[key]))
        if files['coordinates'].stat().st_size == 0 or files['fingerprint'].stat().st_size == 0:
            raise ValueError('Empty structural input in ' + str(folder))
        components.append((folder, original_idx, files))

    # For mixtures, Carbonara uses one experimental profile.  Check every
    # candidate against the first by numeric content, not only file bytes.
    def _saxs_array(path):
        array = np.loadtxt(path, comments='#', ndmin=2)
        if array.ndim != 2 or array.shape[1] not in (2, 3) or not np.isfinite(array).all():
            raise ValueError('Expected finite 2/3-column SAXS data in ' + str(path))
        return array

    reference_saxs = _saxs_array(components[0][2]['saxs'])
    saxs_checks = []
    for i, (_folder, _idx, files) in enumerate(components, start=1):
        current_saxs = _saxs_array(files['saxs'])
        same = (current_saxs.shape == reference_saxs.shape and
                np.allclose(current_saxs, reference_saxs, rtol=1e-7, atol=1e-10))
        saxs_checks.append(bool(same))
        if not same:
            warnings.warn(
                'SAXS data for component ' + str(i) + ' differ from component 1; '
                'the merged fit will use the first component\'s Saxs.dat',
                RuntimeWarning, stacklevel=2,
            )

    # Reuse the established run-script writer.  No DSSP, automatic linker
    # selection or setup regeneration is performed by this import.
    setup = None
    if make_runme:
        setup_path = Path(__file__).resolve().with_name('setup_carbonara_allAtom.py')
        if not setup_path.is_file():
            raise FileNotFoundError(
                'For RunMe generation, setup_carbonara_allAtom.py must be alongside '
                'CarbonaraDataTools.py: ' + str(setup_path)
            )
        spec = importlib.util.spec_from_file_location('_carbonara_setup_for_merge', setup_path)
        setup = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(setup)
        if not callable(getattr(setup, 'write_runme', None)):
            raise RuntimeError('Installed setup script has no write_runme() function')

    has_constraints = any(files['constraints'].is_file() and
                          bool(files['constraints'].read_text().strip())
                          for _folder, _idx, files in components)
    if not np.isfinite(mixture_step) or mixture_step <= 0 or mixture_step > 1:
        raise ValueError('mixture_step must be in (0, 1]')
    if not np.isfinite(mixture_dirichlet_alpha) or mixture_dirichlet_alpha <= 0:
        raise ValueError('mixture_dirichlet_alpha must be > 0')
    if not np.isfinite(mixture_snap) or mixture_snap < 0 or mixture_snap > 1:
        raise ValueError('mixture_snap must be in [0, 1]')
    if int(max_mixture_combos) < 1:
        raise ValueError('max_mixture_combos must be >= 1')

    options = dict(run_options or {})
    controlled = {'working_path', 'fit_name', 'no_structures', 'pairedQ'}
    if controlled.intersection(options):
        raise ValueError('The merger sets working_path, fit_name, no_structures and pairedQ itself')
    defaults = {'fit_n_times': 20, 'min_q': 0.01, 'max_q': 0.2,
                'max_q_start': 0.2, 'max_fit_steps': 10000}
    defaults.update(options)

    # Make the same weight-list convention used by setup_carbonara_allAtom.py.
    # The unnormalised component SAXS profiles are combined with these weights;
    # weights are not scattering-intensity fractions for unequal-size species.
    n = len(components)
    def _snap(row):
        row = np.maximum(np.asarray(row, dtype=float), 0)
        if mixture_snap > 0:
            row = np.round(row / mixture_snap) * mixture_snap
        total = row.sum()
        if total == 0:
            row = np.full(n, 1 / n)
        else:
            row = row / total
        return row

    rows = []
    if int(max_mixture_combos) == 1:
        rows.append(np.full(n, 1 / n))
    elif n == 2:
        steps = int(round(1 / mixture_step))
        rows.extend([np.array([min(1.0, k * mixture_step),
                               max(0.0, 1 - k * mixture_step)]) for k in range(steps + 1)])
        rows.extend([np.array([0., 1.]), np.array([1., 0.])])
    else:
        rows.append(np.full(n, 1 / n))
        rows.extend(np.eye(n))
        rng = np.random.default_rng(0)
        for _ in range(int(max_mixture_combos) * 20):
            if len(rows) >= int(max_mixture_combos) * 3:
                break
            rows.append(rng.dirichlet(np.full(n, mixture_dirichlet_alpha)))
    unique, seen = [], set()
    for row in rows:
        row = _snap(row)
        key = tuple(np.round(row, 8))
        if key not in seen:
            unique.append(row)
            seen.add(key)
        if len(unique) >= int(max_mixture_combos):
            break

    # All checks are done; make an entirely new output, never mutate sources.
    target.mkdir(parents=True, exist_ok=False)
    try:
        shutil.copy2(components[0][2]['saxs'], target / 'Saxs.dat')
        (target / 'fitdata').mkdir()
        (target / 'redundant').touch()
        manifest = {'name': name, 'mixture_n': n,
                    'saxs_reference': str(components[0][2]['saxs']),
                    'saxs_files_match_reference': saxs_checks, 'components': []}
        for new_idx, (folder, original_idx, files) in enumerate(components, start=1):
            for key, prefix in (('coordinates', 'coordinates'),
                                ('fingerprint', 'fingerPrint'),
                                ('flexibility', 'varyingSectionSecondary')):
                shutil.copy2(files[key], target / (prefix + str(new_idx) + '.dat'))
            original_constraints = files['constraints'].is_file() and bool(files['constraints'].read_text().strip())
            if has_constraints:
                outfile = target / ('fixedDistanceConstraints' + str(new_idx) + '.dat')
                if files['constraints'].is_file():
                    shutil.copy2(files['constraints'], outfile)
                else:
                    outfile.touch()  # this component has zero restraints
            shutil.copy2(files['chain_lengths'], target / ('chainLengths' + str(new_idx) + '.dat'))
            if new_idx == 1:
                # Legacy name used by existing notebook/backmapping code.  Store
                # component-specific files too, without combining chain counts.
                shutil.copy2(files['chain_lengths'], target / 'chainLengths.dat')
            manifest['components'].append({
                'new_index': new_idx,
                'source_directory': str(folder),
                'source_index': original_idx,
                'has_distance_constraints': bool(original_constraints),
                'fingerprint_sha256': hashlib.sha256(files['fingerprint'].read_bytes()).hexdigest(),
            })
        np.savetxt(target / 'mixtureFile.dat', np.vstack(unique), fmt='%.3f')
        (target / 'mixture_manifest.json').write_text(
            json.dumps(manifest, indent=2) + '\n', encoding='utf-8'
        )
        if make_runme:
            setup.write_runme(working_path=str(target), fit_name=name,
                              no_structures=n, pairedQ=has_constraints, **defaults)
    except Exception:
        # Do not leave a misleading half-prepared folder.  If the launcher was
        # created, remove only the new launcher owned by this merge operation.
        if make_runme and launch_path.exists():
            launch_path.unlink()
        shutil.rmtree(target)
        raise

    print('Merged ' + str(n) + ' Carbonara components into ' + str(target))
    if make_runme:
        print('Run script: ' + str(launch_path))
    return target

import pandas as pd
import numpy as np

from scipy.spatial.distance import cdist

import os
# Prefer the CPU platform for any PDBFixer/OpenMM clean-up performed during
# Carbonara input preparation. This avoids broken OpenCL drivers crashing setup.
os.environ.setdefault("OPENMM_DEFAULT_PLATFORM", "CPU")
import sys
import subprocess
import shutil
#import math
#import rmsd
from tqdm import tqdm

from typing import List, Set

import shlex
from pathlib import Path

#from Bio.PDB import PDBParser
#from Bio.PDB.DSSP import DSSP
#from Bio.PDB.DSSP import make_dssp_dict
#from DSSPparser import parseDSSP

from pdbfixer import PDBFixer
from openmm.app import PDBFile
from tempfile import NamedTemporaryFile
import matplotlib.pyplot as plt

import biobox as bb

import shutil
import re

from scipy import interpolate

from plotly.subplots import make_subplots
import plotly.graph_objects as go
import json

# Set Plotly to use the correct renderer for Colab

from glob import glob


def _carbonara_openmm_cpu_platform_or_none():
    """Return the OpenMM CPU platform if available, otherwise None."""
    try:
        from openmm import Platform
        return Platform.getPlatformByName("CPU")
    except Exception:
        return None


def _carbonara_pdbfixer(filename):
    """Construct PDBFixer while explicitly preferring CPU over OpenCL/CUDA."""
    cpu_platform = _carbonara_openmm_cpu_platform_or_none()
    if cpu_platform is not None:
        try:
            return PDBFixer(filename=filename, platform=cpu_platform)
        except TypeError:
            pass
    return PDBFixer(filename=filename)


def _numbered_indices_for_carbonara_inputs(directory, prefix):
    """
    Return the numbered Carbonara component indices present in ``directory``.

    The primary signal is the matching file family itself, e.g.
    fingerPrint1.dat, fingerPrint2.dat, ... .  We also look at
    coordinates1.dat, coordinates2.dat, ... because in multi-structure runs the
    coordinate files define how many component copies should exist.
    """
    directory = Path(directory)
    indices = set()

    for pattern_prefix in (prefix, "coordinates"):
        for p in directory.glob(f"{pattern_prefix}*.dat"):
            m = re.fullmatch(rf"{re.escape(pattern_prefix)}(\d+)\.dat", p.name)
            if m:
                indices.add(int(m.group(1)))

    return sorted(i for i in indices if i >= 1)


def _replicate_numbered_file_from_1(filename, prefix):
    """
    After writing ``prefix1.dat``, copy it to the matching numbered files.

    This is used by the notebook editing cells.  The user edits/writes only
    fingerPrint1.dat or varyingSectionSecondary1.dat, but in a multi-structure
    run the same edited file must also be present as fingerPrint2.dat,
    fingerPrint3.dat, ... or varyingSectionSecondary2.dat,
    varyingSectionSecondary3.dat, ... .

    For ordinary single-structure runs this is a no-op.  For arbitrary filenames
    not ending in ``prefix1.dat`` this is also a no-op, preserving old behaviour.
    """
    path = Path(filename)
    if path.name != f"{prefix}1.dat":
        return [str(path)]

    indices = _numbered_indices_for_carbonara_inputs(path.parent, prefix)
    if not indices:
        indices = [1]

    written = [str(path)]
    for idx in indices:
        if idx == 1:
            continue
        dst = path.with_name(f"{prefix}{idx}.dat")
        shutil.copy2(path, dst)
        written.append(str(dst))

    return written

#import hdbscan

import mdtraj as md

# Utility stuff

def list_nohidden(path):
    lst = []
    for f in os.listdir(path):
        if not f.startswith('.'):
            lst.append(f)
    return lst


def sort_by_creation(file_lst):

    files = list(filter(os.path.isfile, file_lst))
    files.sort(key=lambda x: os.path.getmtime(x))
    return files


# Pulling in PDB to Carbonara format

def pdb_2_biobox(pdb_file):
    M = bb.Molecule()
    ext = os.path.splitext(pdb_file)[1].lower()
    if ext == ".pdb":
        M.import_pdb(pdb_file)
    elif ext == ".cif":
        # Load with PDBFixer
        fixer = _carbonara_pdbfixer(pdb_file)
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens()

        # Write to temporary PDB file and load with MDTraj
        with NamedTemporaryFile(mode="w", suffix=".pdb", delete=False) as temp:
           PDBFile.writeFile(fixer.topology, fixer.positions, temp)
           temp_path = temp.name
        M.import_pdb(temp_path)
        os.remove(temp_path)  # clean up temporary file
    else:
        raise ValueError(f"Unsupported file extension: {ext}")
    return M



def extract_CA_coordinates(M):
    ca_idx = (M.data['name']=='CA').values
    ca_coords = M.coordinates[0][ca_idx]
   # if ca_coords.shape[0] != M.data['resid'].nunique():
    #    raise Exception("You better check your PDB... The number of CA atoms does not equal the number of ResIDs in your PDB file!")
    #else:
    return ca_coords


def histConv(aminoName):

    if(aminoName == 'HIP' or aminoName == 'HID' or aminoName == 'HSD' or aminoName == 'HSE' or aminoName == 'HIE'):
         return "HIS"
    elif(aminoName == 'GLH' or aminoName == 'GLUP'):
         return "GLU"
    elif(aminoName == 'CYX' or aminoName == 'CYM'):
         return "CYS"
    elif(aminoName == 'ASH' or aminoName == 'ASPP'):
         return "ASP"
    elif(aminoName == 'LYN' or aminoName == 'LSN'):
         return "LYS"
    else:
        return aminoName


def extract_sequence(M):

    aa_names = {
                'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU',
                'F': 'PHE', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
                'K': 'LYS', 'L': 'LEU', 'M': 'MET', 'N': 'ASN',
                'P': 'PRO', 'Q': 'GLN', 'R': 'ARG', 'S': 'SER',
                'T': 'THR', 'V': 'VAL', 'W': 'TRP', 'Y': 'TYR'
                }

    names_aa = {y: x for x, y in aa_names.items()}


    ca_idx = (M.data['name']=='CA').values

    aminoDat = M.data['resname'][ca_idx]
    
    for i in range(len(aminoDat)):
        aminoDat[i:i+1].values[0] = histConv(aminoDat[i:i+1].values[0])
    
    resnames = aminoDat.map(names_aa).values

    # resnames = M.data['resname'][ca_idx].map(names_aa).values

    #if resnames.shape[0] != M.data['resid'].nunique():
    #    raise Exception("You better check your PDB... The number of CA atoms does not equal the number of ResIDs in your PDB file!")
    #else:
    return resnames


# > From Carbonara format

def extract_coords(coords_file):

    coords = np.genfromtxt(coords_file)
    coords = coords[~np.isnan(coords).any(axis=1)]

    return coords


def read_coords(coords_file):

    coords = np.genfromtxt(coords_file)
    coords = coords[~np.isnan(coords).any(axis=1)]

    return coords


def extract_sequence_file(fingerprint_file):

    seqin = open(fingerprint_file, 'r').readlines()
    seqout=""
    for i in range(2,len(seqin),4):
         seqout= seqout+seqin[i][:-1]

    return seqout


# > returning to PDB from Carbonara

# --- Minimal patch: Carbonara_2_PDB that supports multiple chains ---

import re
import numpy as np
import pandas as pd
import biobox as bb  # assumes Matteo's biobox is installed

def extract_coords_xyz_or_carbonara(coords_file):
    """Read a 3-column XYZ-style coords file; ignore lines like 'End chain ...' (any case)."""
    raw = []
    with open(coords_file, "r") as f:
        for ln in f:
            if not ln.strip():
                continue
            if re.search(r'end\s+chain', ln, re.IGNORECASE):
                continue
            raw.append(ln.strip())
    xyz = []
    for ln in raw:
        parts = ln.split()
        if len(parts) < 3:
            raise ValueError(f"Bad coordinate line (need 3 floats): {ln!r}")
        x, y, z = map(float, parts[:3])
        xyz.append((x, y, z))
    return np.asarray(xyz, dtype=float)

def Carbonara_2_PDB_multichain(coords_file, fp_file, output_file,
                               segment_lengths, chain_ids=None):
    """
    Writes CA-only PDBs from Carbonara output with distinct chains.

    coords_file   : 3-column XYZ coords; 'End chain ...' lines are allowed and ignored
    fp_file       : Carbonara fingerprint (sequence); must match total length
    output_file   : PDB to write
    segment_lengths : list of ints, residues per chain (sum must equal #coords == len(seq))
    chain_ids     : list like ['A','B','C',...]; defaults to 'A','B','C',... as needed
    """
    # --- read coords & sequence ---
    coords = extract_coords_xyz_or_carbonara(coords_file)
    size = coords.shape[0]

    # you likely already have these helpers; keeping names from your code:
    seq = extract_sequence_file(fp_file)  # returns 1-letter string
    if len(seq) != size:
        raise ValueError(f"Sequence length {len(seq)} != number of coordinates {size}")

    if sum(segment_lengths) != size:
        raise ValueError(f"sum(segment_lengths) {sum(segment_lengths)} != #coords {size}")

    # --- chain ids ---
    if chain_ids is None:
        import string
        letters = list(string.ascii_uppercase)
        chain_ids = []
        i = 0
        while len(chain_ids) < len(segment_lengths):
            if i < 26:
                chain_ids.append(letters[i])
            else:
                n = i - 26
                a, b = divmod(n, 26)
                chain_ids.append(letters[a] + letters[b])  # AA, AB, ...
            i += 1
    if len(chain_ids) < len(segment_lengths):
        raise ValueError("Not enough chain_ids for the given segment_lengths")

    # --- 1-letter -> 3-letter ---
    aa_map = {
        'A':'ALA','C':'CYS','D':'ASP','E':'GLU','F':'PHE','G':'GLY','H':'HIS','I':'ILE',
        'K':'LYS','L':'LEU','M':'MET','N':'ASN','P':'PRO','Q':'GLN','R':'ARG','S':'SER',
        'T':'THR','V':'VAL','W':'TRP','Y':'TYR'
    }
    seq_3 = [aa_map[a] for a in seq]

    # --- build per-residue chain labels and per-chain residue numbers ---
    chains = []
    resids = []
    offset = 0
    for L, cid in zip(segment_lengths, chain_ids):
        chains.extend([cid] * L)
        # per-chain residue numbers start at 1
        resids.extend(list(range(1, L + 1)))
        offset += L
    chains = np.array(chains)
    resids = np.array(resids, dtype=int)

    # --- biobox molecule ---
    df = pd.DataFrame({
        'atom'     : ['ATOM'] * size,
        'index'    : np.arange(1, size + 1),  # atom serial (1-based)
        'name'     : ['CA'] * size,
        'resname'  : seq_3,
        'chain'    : chains,
        'resid'    : resids,
        'occupancy': [1.0] * size,
        'beta'     : [50.0] * size,
        'atomtype' : ['C'] * size,
        'radius'   : [1.7] * size,
        'charge'   : [0.0] * size,
    })

    molecule = bb.Molecule()
    molecule.data = df
    molecule.coordinates = np.expand_dims(coords, axis=0)  # (frames, atoms, 3)

    # NOTE: many PDB readers infer chain breaks from chain ID changes;
    # if you *need* explicit TER records, write directly instead of via biobox,
    # or postprocess to insert TER lines. Most tools are fine with chain IDs.
    molecule.write_pdb(output_file)

# ---- example call (your case) ----
# Carbonara_2_PDB_multichain(
#     coords_file="mol1_sub_0_step_5_xyz.dat",
#     fp_file="your_fingerprint.fp",
#     output_file="mol1_sub_0_step_5_CA_chains_from_biobox.pdb",
#     segment_lengths=[900, 900, 900, 900],
#     chain_ids=list("ABCD"),
# )


def Carbonara_2_PDB(coords_file, fp_file, output_file):

    '''
    Writes alpha carbon PDBs from Carbonara output

    Input
        coords_file      : coordinates of the carbon alpha chain
        fingerprint_file : Carbonara specific format containing secondary structure and sequence
        output_file      : define name of write output
    '''

    # read in coordinates and fingerprint
    coords = extract_coords(coords_file)
    size = coords.shape[0]
    seq = extract_sequence_file(fp_file)

    # map the AA shorthand to 3 letter abr.
    aa_map = {
            'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU',
            'F': 'PHE', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
            'K': 'LYS', 'L': 'LEU', 'M': 'MET', 'N': 'ASN',
            'P': 'PRO', 'Q': 'GLN', 'R': 'ARG', 'S': 'SER',
            'T': 'THR', 'V': 'VAL', 'W': 'TRP', 'Y': 'TYR'
        }

    # map sequence to 3 character abr.
    seq_3 = []
    for a in list(seq):
        seq_3.append(aa_map[a])

    # create dataframe for biobox molecule type
    df = pd.DataFrame({'atom':['ATOM']*size, 'index':np.arange(size), 'name':['CA']*size,
                       'resname':seq_3, 'chain':['A']*size, 'resid':np.arange(size),
                       'occupancy':[1]*size, 'beta':[50]*size, 'atomtype':['C']*size,
                       'radius':[1.7]*size, 'charge':[0]*size})

    # take full advantage of Matteo's lovely biobox library - manually 'create' a molecule
    molecule = bb.Molecule()
    molecule.data = df
    molecule.coordinates = np.expand_dims(coords, axis=0)

    # write out!
    molecule.write_pdb(output_file)


# The meeet of extracting

def clean_s_sequences(arr):
    arr = np.array(arr)
    result = arr.copy()

    i = 0
    while i < len(arr):
        if arr[i] == 'S':
            # Check for isolated single 'S'
            if (i == 0 or arr[i - 1] != 'S') and (i + 1 == len(arr) or arr[i + 1] != 'S'):
                result[i] = '-'
                i += 1
            # Check for isolated pair 'S S' not part of longer sequence
            elif (i + 1 < len(arr) and arr[i + 1] == 'S' and
                  (i == 0 or arr[i - 1] != 'S') and
                  (i + 2 == len(arr) or arr[i + 2] != 'S')):
                result[i] = '-'
                result[i + 1] = '-'
                i += 2
            else:
                # Part of a longer sequence, skip to end of sequence
                while i < len(arr) and arr[i] == 'S':
                    i += 1
        else:
            i += 1
    return result


def pull_structure_from_pdb(file_path):
    """
    Pulls the structure from a (single) PDB file using MDTraj and returns the coordinates
    and sequence of the chain(s).

    Parameters:
        pdb_file (str): The path of the PDB file.

    Returns:
        coords_chain (list): A list of numpy arrays containing the CA coordinates of the chain(s).
        sequence_chain (list): A list of numpy arrays containing the sequence of the chain(s).
        secondary_structure_chains (list): A list of predicted secondary structures for each chain.
        missing_residues_chain (list): A list of numpy arrays containing the missing residues of the chain(s).

    Raises:
        ValueError: If no chains are found in the PDB file.

    Example:
        coords_chains, sequence_chains, secondary_structure_chains, missing_residues_chains =
            pull_structure_from_pdb_mdtraj('/path/to/pdb/file.pdb')
    """
    ext = os.path.splitext(file_path)[1].lower()

    """
    To deal with alphaFold3 cif format
    """
    if ext == ".pdb":
        traj = md.load(file_path)
    elif ext == ".cif":
        # Load with PDBFixer
        fixer = _carbonara_pdbfixer(file_path)
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens()

        # Write to temporary PDB file and load with MDTraj
        with NamedTemporaryFile(mode="w", suffix=".pdb", delete=False) as temp:
           PDBFile.writeFile(fixer.topology, fixer.positions, temp)
           temp_path = temp.name
        traj = md.load(temp_path)
        os.remove(temp_path)  # clean up temporary file
    else:
        raise ValueError(f"Unsupported file extension: {ext}")

    topology = traj.topology
    three_to_one = get_residue_map()

    chains = list(topology.chains)
    if len(chains) == 0:
        raise ValueError("No chains found in structure")

    coords_chains = []
    sequence_chains = []
    secondary_structure_chains = []
    missing_residues_chains = []

    ss_pred = md.compute_dssp(traj, simplified=True)[0]
    ss_map = {'H': 'H', 'E': 'S', 'C': '-', 'NA': '-'}
    ss_pred_mapped = np.array([ss_map.get(ss, '-') for ss in ss_pred])
    ss_pred_mapped = clean_s_sequences(ss_pred_mapped)

    residue_index = 0
    for chain in chains:
        residues = list(chain.residues)
        ca_atoms_indices = []
        resids = []
        seq = []

        for res in residues:
            resids.append(res.resSeq)
            seq.append(three_to_one.get(res.name, 'X'))

            for atom in res.atoms:
                if atom.name == 'CA':
                    ca_atoms_indices.append(atom.index)
                    break

        if ca_atoms_indices:
            ca_coords = traj.xyz[0, ca_atoms_indices, :] * 10  # nm to Å
            if len(ca_coords) > 10:
                coords_chains.append(ca_coords)
                sequence_chains.append(np.array(seq))

                chain_ss = ss_pred_mapped[residue_index:residue_index + len(residues)]
                secondary_structure_chains.append(chain_ss)

                resids = np.array(resids)
                missing_residues = find_missing_residues(resids)
                missing_residues_chains.append(missing_residues)

                residue_index += len(residues)

    return coords_chains, sequence_chains, secondary_structure_chains, missing_residues_chains



def pull_structure_from_pdb(pdb_file):
    """
    Pulls the structure from a (single) PDB file using MDTraj and returns the coordinates
    and sequence of the chain(s).

    Parameters:
        pdb_file (str): The path of the PDB file.

    Returns:
        coords_chain (list): A list of numpy arrays containing the CA coordinates of the chain(s).
        sequence_chain (list): A list of numpy arrays containing the sequence of the chain(s).
        secondary_structure_chains (list): A list of predicted secondary structures for each chain.
        missing_residues_chain (list): A list of numpy arrays containing the missing residues of the chain(s).

    Raises:
        ValueError: If no chains are found in the PDB file.

    Example:
        coords_chains, sequence_chains, secondary_structure_chains, missing_residues_chains =
            pull_structure_from_pdb_mdtraj('/path/to/pdb/file.pdb')
    """
    # Load the PDB file
    traj = md.load(pdb_file)
    
    # Get topology
    topology = traj.topology
    
    # Create a mapping for three-letter to one-letter amino acid codes
    three_to_one = get_residue_map()

    
    # Get unique chains
    chains = [chain for chain in topology.chains]

    if len(chains) == 0:
        raise ValueError("No chains found in pdb file")
    
    coords_chains = []
    sequence_chains = []
    secondary_structure_chains = []
    missing_residues_chains = []
    
    # Compute secondary structure for the entire trajectory
    ss_pred = md.compute_dssp(traj, simplified=True)[0]
    ss_map = {'H': 'H', 'E': 'S', 'C': '-','NA': '-'}
   
    ss_pred_mapped = np.array([ss_map[ss] for ss in ss_pred])

    ss_pred_mapped = clean_s_sequences(ss_pred_mapped)

    # For each chain in the PDB
    residue_index = 0
    for chain in chains:
        # Get residues in this chain
        residues = list(chain.residues)
        
        # Extract CA atoms for this chain
        ca_atoms_indices = []
        resids = []
        seq = []

        for res in residues:
            resids.append(res.resSeq)
            
            # Get one letter code for the residue
            if res.name in three_to_one:
                seq.append(three_to_one[res.name])
            else:
                seq.append('X')  # Unknown amino acid
                
            # Find CA atom index
            for atom in res.atoms:
                if atom.name == 'CA':
                    ca_atoms_indices.append(atom.index)
                    break
        
        # Get coordinates of CA atoms
        if ca_atoms_indices:
            ca_coords = traj.xyz[0, ca_atoms_indices, :]*10 # << nm to A!!!
            if(len(ca_coords)>10):
                coords_chains.append(ca_coords)
                sequence_chains.append(np.array(seq))
                
                # Get secondary structure for this chain
                chain_ss = ss_pred_mapped[residue_index:residue_index + len(residues)]
                secondary_structure_chains.append(chain_ss)
                
                # Find missing residues
                resids = np.array(resids)
                missing_residues = find_missing_residues(resids)
                missing_residues_chains.append(missing_residues)
                
                residue_index += len(residues)
            
    return coords_chains, sequence_chains, secondary_structure_chains, missing_residues_chains




# dealing with unclean PDBs


def find_missing_residues(resIDs):
    
    missing_residues = []

    # missing begining residues
    for i in range(1, resIDs[0]):
        missing_residues.append(i)
    
    # looking for residues labels with gaps greater than 1
    res_diff = np.diff(resIDs)
    missing_indices = np.where( res_diff > 1 )[0]

    for m in missing_indices:
        
        number_missing = res_diff[m] - 1

        # account for missing sequential residues
        for i in range(number_missing):
            missing_residues.append( resIDs[m] + i + 1 )

    return np.asarray( missing_residues )


def get_residue_map(direction='321'):
    
    aa_names = {
                'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU',
                'F': 'PHE', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
                'K': 'LYS', 'L': 'LEU', 'M': 'MET', 'N': 'ASN',
                'P': 'PRO', 'Q': 'GLN', 'R': 'ARG', 'S': 'SER',
                'T': 'THR', 'V': 'VAL', 'W': 'TRP', 'Y': 'TYR'
                }

    names_aa = {y: x for x, y in aa_names.items()}

    if direction == '123':
        return aa_names
    else:
        return names_aa



# Geometrical check for chain breaks
def missing_ca_check(coords, threshold_dist_Å = 10):

    breaking_indices = np.where( np.linalg.norm(np.diff(coords, axis=0), axis=1) > threshold_dist_Å )[0] + 1

    return breaking_indices


# break into chains based on the breaking indices found geometrically
def break_into_chains(coords, sequence, breaking_indices):

    coords_chains = np.array_split(coords, breaking_indices)
    sequence_chains = np.array_split(sequence, breaking_indices)

    return coords_chains, sequence_chains


# > Curvature & Torsion calculation for geometrically inferring secondary structure


# > CB inferrence stuff (not yet used too much)

def CA_PDB_2_CB_PDB(CA_pdb, output_pdb):

    # Load protein into BioBox object
    protein = bb.Molecule(CA_pdb)

    # Get CA coordinates
    CA_xyz = protein.coordinates[0]

    # infer the CB positions
    CB_xyz = infer_CB_positions(CA_xyz)

    interlace_CA_CB_write(CA_xyz, CB_xyz, protein, output_pdb)



# Geometric justified chain splitting

def splitMolecule(coords,sequence,secondary):
    splitList=[]
    for i in range(0,coords.shape[0]-1):
        dist = np.linalg.norm(coords[i+1]-coords[i])
        if dist>10.0:
            splitList.append(i+1)
    if len(splitList)>0:
        coordssplit = np.array_split(coords,splitList)
        sequencesplit = np.array_split(sequence,splitList)
        secondarysplit = np.array_split(secondary,splitList)
    else:
        coordssplit = [coords]
        sequencesplit = [sequence]
        secondarysplit = [secondary]
    return coordssplit,sequencesplit,secondarysplit,splitList


def splitMoleculeNoSec(coords, sequence):
    '''
    Split a molecule into segments based on a distance threshold and return the split coordinates and sequences.

    Parameters:
    coords (numpy.ndarray): The coordinates of the molecule.
    sequence (numpy.ndarray): The sequence of the molecule.

    Returns:
    coordssplit (list): A list of numpy arrays containing the split coordinates.
    sequencesplit (list): A list of numpy arrays containing the split sequences.
    splitList (list): A list of indices where the molecule is split.
    '''
    
    splitList = []
    for i in range(0, coords.shape[0] - 1):
        dist = np.linalg.norm(coords[i + 1] - coords[i])
        if dist > 10.0:
            splitList.append(i + 1)
    if len(splitList) > 0:
        coordssplit = np.array_split(coords, splitList)
        sequencesplit = np.array_split(sequence, splitList)
    else:
        coordssplit = [coords]
        sequencesplit = [sequence]

    return coordssplit, sequencesplit, splitList


# Dealing with secondary section selection / manipulation

def section_finder(ss):

    '''Find protein sub-unit sections from the full secondary structure'''

    sections = []
    structure_change = np.diff(np.unique(ss, return_inverse=True)[1])

    for i, c in enumerate( structure_change ):

        if c!=0:
            sections.append(ss[i])

        if i==structure_change.shape[0]-1:
            sections.append(ss[i])

    sections = np.array(sections)

    return sections #, linker_indices #, structure_change


def section_finder_sub(ss):

    '''Find protein sub-unit sections from the full secondary structure and split the structure into these subsections'''

    sections = []
    subsection =[]
    for sec in ss:
        structure_change = np.diff(np.unique(sec, return_inverse=True)[1])
        for i, c in enumerate( structure_change ):
            subsection.append(sec[i])
            if c!=0:
                sections.append(subsection)
                subsection=[]

            if i==structure_change.shape[0]-1:
                subsection.append(sec[i+1])
                sections.append(subsection)
                subsection=[]

    return sections


def find_sheet_indices(sections):

    '''Find sheet sub-unit section indices'''

    sheet_indices = np.where(sections=='S')[0]
    return sheet_indices


def find_linker_indices(sections):

    '''Find linker sub-unit section indices'''

    linker_indices = np.where(sections=='-')[0]
    return linker_indices

def sheet_group_mask(ss,group,sheet_groups,endindex):

    '''Groups adjacent sheets in secondary structure file and returns a grouping mask ( 0 : not a sheet;  1+: sheet )

    Parameters
    ss (numpy array):            Secondary structure labels (array of strings)

    Returns
    sheet_groups (numpy array):  Mask of grouped sheet sections
    '''

    sheet_mask = (ss == 'S')*1

    if sheet_mask[0] == 1:
        label = True
    else:
        label = False

    for i, c in enumerate(np.diff(sheet_mask)):


        if c == 1:
            label = True

        elif c==-1:
            label=False
            group += 1

        else:
            pass

        if label == True:
            if ss[i+1] == 'S':
                sheet_groups[i+1+endindex] = group
    endindex = endindex + ss.shape[0]
    return group,endindex


def linker_group_mask(ss):

    '''Groups adjacent linkers in secondary structure file and returns a grouping mask ( 0 : not a linker;  1+: linker )

    Parameters
    ss (numpy array):             Secondary structure labels (array of strings)

    Returns
    linker_groups (numpy array):  Mask of grouped linker sections
    '''

    linker_mask = (ss == '-')*1
    linker_groups = np.zeros(ss.shape[0])
    group = 1

    # checking first index for linker
    if linker_mask[0] == 1:
        label = True
        linker_groups[0] = group
    else:
        label = False

    for i, c in enumerate(np.diff(linker_mask)):

        if c == 1:
            label = True

        elif c==-1:
            label=False
            group += 1

        else:
            pass

        if label == True:

            linker_groups[i+1] = group

    return linker_groups #, linker_mask


def get_sheet_coords(coords, sheet_groups):

    '''Finds CA coordinates of

    Parameters
    coords (numpy array):        xyz coordinates of all protein CA atoms
    sheet_groups (numpy array):  Mask of grouped sheet sections

    Returns
    sheet_coords (numpy array):  xyz coordinates of CA atoms in each sheet structure [ [...sheet 1 coords...] [...sheet 2 coords...] ... ]
    '''

    sheet_coords = []

    for g in np.unique(sheet_groups):
        if g>0:
            sheet_coords.append(coords[sheet_groups==g])

    return sheet_coords



def get_section_groupings(ss, structure_change):

    group = 0
    structural_groups = np.zeros(ss.shape)
    structural_groups[0] = group

    for i, c in enumerate(structure_change):

        if c != 0:
            group += 1

        structural_groups[i+1] = group
    return structural_groups



def get_secondary(fingerprint_file):
    fplst=np.loadtxt(fingerprint_file, str)
    fplstout= [np.asarray(list(fplst[i])) for i in range(2,len(fplst),2)]
    return fplstout
    
    
def sheet_pipe(coords_file, fingerprint_file):

    coords = read_coords(coords_file)
    ss = get_secondary(fingerprint_file)
    sheet_groups = np.zeros(coords.shape[0])
    group=1
    endindex=0
    for i in range(0,len(ss)):
        group,endindex = sheet_group_mask(ss[i],group,sheet_groups,endindex);
    sheet_coords = get_sheet_coords(coords, sheet_groups)
    return sheet_coords


# - * - Inferring flexible regions of protein - * -

def generate_random_structures(coords_file, fingerprint_file, return_failures=False):

    '''Generate random structures changing one linker section at a time

    Parameters
    coords_file:       /path/ to CA coordinates.dat file
    fingerprint_file:  /path/ to fingerprint.dat file

    Return
    Generated structures are written to ~/rand_structures/.. section_*LINKERINDEX*.dat as xyz
    Linker Indices
    '''
    secondarystruct = get_secondary(fingerprint_file)

    # Build the same global section numbering used downstream, but do not pass
    # very short linker sections into the C++ generate_structure executable.
    # The later auto selector already uses len(section) > 3; doing the same
    # here prevents length-1/2/3 terminal or boundary fragments from being
    # trial-moved during the sheet/disulfide pre-check stage.
    min_linker_len = 4
    linker_indices = []
    skipped_short_linkers = []
    skipped_chain_start_linkers = []
    currMax = 0
    for chain_i, ss in enumerate(secondarystruct):
        sections = split_into_sections(ss)
        for local_i, sec in enumerate(sections):
            if len(sec) == 0 or sec[0] != '-':
                continue
            global_i = currMax + local_i
            # An N-terminal linker has no preceding structural frame.  The
            # original selector excludes only global section 0 *after* trial
            # generation; apply the same exclusion at every chain boundary
            # *before* invoking the empirical C++ internal-linker sampler.
            if local_i == 0:
                skipped_chain_start_linkers.append((global_i, chain_i + 1, len(sec)))
                continue
            if len(sec) >= min_linker_len:
                linker_indices.append(global_i)
            else:
                skipped_short_linkers.append((global_i, chain_i + 1, local_i, len(sec)))
        currMax += len(sections)

    if skipped_chain_start_linkers:
        print('CarbonaraDataTools: excluded N-terminal linkers from internal-linker '
              'trials (global section, chain, residues): '
              + ', '.join(f'{g} (chain {c}, n={n})'
                          for g, c, n in skipped_chain_start_linkers))

    if skipped_short_linkers:
        preview = ', '.join(
            f'{g}(chain {c}, local {l}, n={n})'
            for g, c, l, n in skipped_short_linkers[:10]
        )
        if len(skipped_short_linkers) > 10:
            preview += ', ...'
        print(
            'CarbonaraDataTools: skipped '
            f'{len(skipped_short_linkers)} linker section(s) shorter than '
            f'{min_linker_len} residues before generate_structure: {preview}'
        )

    linker_indices = np.asarray(linker_indices, dtype=int)
    # The importable CarbonaraDataTools.py lives in the Carbonara installation root.
    # Unlike sys.argv[0], this is independent of notebook / launcher arguments.
    current = str(Path(__file__).resolve().parent)
    generator_executable = os.path.join(current, 'build', 'bin', 'generate_structure')
    if not os.path.isfile(generator_executable):
        raise FileNotFoundError('Carbonara trial generator not found: ' + generator_executable)

    # generate_structure requires the interpolation-data root. Its C++ entry
    # point may not set CARBONARA_ROOT, so carbonaraConfig.cpp reads this env var.
    generator_env = os.environ.copy()
    generator_env['CARBONARA_ROOT'] = current

    random = 'rand_structures'
    random_working = os.path.join(os.getcwd(), random)

    if os.path.exists(random_working) and os.path.isdir(random_working):
        shutil.rmtree(random_working)

    os.mkdir(random_working)

    # try:

    # except OSError as error:
    #     print(str(error)[11:])

    # print('Beginning random structures generation \n')

    rand_file_dict = {}
    successful_indices = []
    failed_sections = []
    for l in tqdm(linker_indices):

        outputname = random_working+'/section_'+str(l)

#         !./generate_structure {fingerprint_file} {coords_file} {outputname} {l}
        # Retry only stochastic no-output runs. An unchanged but written trial
        # remains valid, as in the existing audited beta-sheet selection.
        max_generator_attempts = 4
        for attempt in range(1, max_generator_attempts + 1):
            result = subprocess.run(
                [generator_executable, str(fingerprint_file), str(coords_file), outputname, str(l)],
                cwd=current,
                env=generator_env,
                capture_output=True,
                text=True,
            )
            # A nonzero exit is a real C++ / installation failure, not simply
            # a rejected geometric trial: preserve the existing hard error.
            if result.returncode != 0:
                raise RuntimeError(
                    f'Carbonara generate_structure failed for section {int(l)} '
                    f'(return code {result.returncode}).\n'
                    f'Command: {generator_executable} ... {int(l)}\n'
                    f'stderr:\n{result.stderr[-3000:]}\nstdout:\n{result.stdout[-3000:]}'
                )

            if any(Path(random_working).glob(f'section_{int(l)}_*')):
                successful_indices.append(int(l))
                if attempt > 1:
                    print(f'Carbonara: section {int(l)} generated trial files '
                          f'on attempt {attempt}/{max_generator_attempts}.')
                break
            if attempt < max_generator_attempts:
                print(f'Carbonara: section {int(l)} generated no trial files '
                      f'on attempt {attempt}; retrying.')
        else:
            failure = {
                'section_index': int(l),
                'attempts': max_generator_attempts,
                'reason': 'generator returned success but produced no trial files',
                'selected_automatically': False,
            }
            failed_sections.append(failure)
            print(f'WARNING: Carbonara section {int(l)} produced no trial files '
                  f'after {max_generator_attempts} attempts. Excluding it from '
                  'automatic linker selection; continuing setup.')

    # print('')
    # print('Finished generating random structures')

    successful_indices = np.asarray(successful_indices, dtype=int)
    if return_failures:
        return successful_indices, failed_sections
    return successful_indices


def sheet_pairwise_bond_number(sheet_coords, thr=5.5):

    '''Finds the number of pairs of CA atoms within some threshold between all sheet sections

    Parameters
    sheet_coords (numpy array): xyz coordinates of CA atoms in each sheet structure [ [...sheet 1 coords...] [...sheet 2 coords...] ... ]
    thr (float) {optional}:     Cutoff distance for inter-sheet bonding (default = 5.5 Å)

    Returns
    pairwise_bond_num (numpy array): Lower triangular array containing the number of individual CA bonds within threshold between each sheet pair

    '''

    number_bonds = 0

    pairwise_bond_num = np.zeros([len(sheet_coords), len(sheet_coords)])

    for i in range(1,len(sheet_coords)):

        for j in range(0,i):

            arr1, arr2 = sheet_coords[j], sheet_coords[i]
            dist_matrix = cdist(arr1, arr2)
            indices = np.where(dist_matrix < thr)

            pairwise_bond_num[i,j] = indices[0].shape[0]

            number_bonds += indices[0].shape[0]
    return pairwise_bond_num


def random_bond_finder(rand_file_dir, fingerprint_file, linker_indices):

    # grouping all random structure for each linker together

    struture_lst = list_nohidden(rand_file_dir)

    linker_file_dict = {}
    for l in linker_indices:
        tmp = []

        for file in np.sort(struture_lst):
            if str(l) == file.split('_')[1]:
                tmp.append(file)

        linker_file_dict[l] = tmp

    # Pairwise sheet bonds for each random str for each linker
    linker_bond_dict = {}

    for l in linker_indices:

        tmp = []

        for file in linker_file_dict[l]:
            coords_file = rand_file_dir+file
            sheet_coords = sheet_pipe(coords_file, fingerprint_file)
            tmp.append( sheet_pairwise_bond_number(sheet_coords) )

        linker_bond_dict[l] = tmp

    return linker_bond_dict


def find_non_varying_linkers(initial_coords_file, fingerprint_file):

    # initial_coords_file = 'Fitting/coordinates1.dat'
    # fingerprint_file = 'Fitting/fingerPrint1.dat'

    # Reference initial structure
    sheet_coords = sheet_pipe(initial_coords_file,
                              fingerprint_file)
    ref_bonds = sheet_pairwise_bond_number(sheet_coords, thr=5.5)

    # Generate the random structure changing each linker section
    linker_indices = generate_random_structures(initial_coords_file, fingerprint_file)

    # Calculate the number of inter-sheet bonds for each rand struct
    linker_bond_arr_dict = random_bond_finder('rand_structures/',
                                              fingerprint_file,
                                              linker_indices)

    # Find number of bond breaks relative to initial structure
    bond_breaks_dict = {}


    for l in linker_indices:

        bond_break_lst = []
        for bond_arr in linker_bond_arr_dict[l]:


            bond_break_lst.append( (ref_bonds > bond_arr).sum() )

        bond_breaks_dict[l] = sum(bond_break_lst)/(len(linker_bond_arr_dict[l])+1)

    # Linker indices that cause no bond breaks
    conds = np.asarray(list(bond_breaks_dict.values())) < 0.0000001


    allowed_linker = linker_indices[conds]

    if 0 in linker_indices:
        linker_indices = np.delete(linker_indices, np.where(linker_indices==0)[0].item())


    if 0 in allowed_linker:
        allowed_linker = np.delete(allowed_linker, np.where(allowed_linker==0)[0].item())

    return allowed_linker, linker_indices



def auto_select_varying_linker(coords_file, fingerprint_file):

    allowed_linker, linker_indices = find_non_varying_linkers(initial_coords_file = coords_file,
                                                                fingerprint_file = fingerprint_file)

    secondary = get_secondary(fingerprint_file)
    sections = section_finder_sub(secondary)
    varying_linker_indices = []
    for section_index in allowed_linker:
        if len(sections[section_index]) > 3:
            varying_linker_indices.append(section_index)

    # dict of linker lengths - maybe we priotise longer earlier or something?
    # can we find a way to equate overall structure impact to each linker? Have some sort of scale change approach?
    linker_length_dict = {}
    for section_index in varying_linker_indices:
        linker_length_dict[section_index] = len(sections[section_index])
    
    return varying_linker_indices


# ------ Carbonara Setup Methods ---------


def setup_fit_master_dir(root_dir=os.getcwd(), fit_master_name='newFitData'):
    """
    Creates a master directory containing all refinement directories.

    Parameters:
        root_dir (str): The root directory where the fit master directory will be created. Default is the current working directory.
        fit_master_name (str): The name of the fit master directory. Default is 'newFitData'.

    Returns:
        str: The path of the fit master directory.

    Raises:
        OSError: If the fit master directory cannot be created.

    Example:
        fit_master_dir = setup_fit_master_dir('/path/to/root', 'fitMaster')
    """
    
    fit_master_dir = os.path.join(root_dir, fit_master_name)

    if not os.path.exists(fit_master_dir):
        try:
            os.makedirs(fit_master_dir)
        except OSError as e:
            raise OSError(f"Failed to create fit master directory: {e}")

    return fit_master_dir


def setup_refinement_dir(refine_name, fit_master_dir = 'newFitData'):
    """
    Creates the refinement directory where carbonara will refine a target (in /fit master directory/).
    
    Parameters:
        refine_name (str): The name of the refinement directory.
        fit_master_dir (str): The name of the fit master directory. Default is 'newFitData'.

    Returns:
        str: The path of the refinement directory.

    Raises:
        OSError: If the refinement directory cannot be created.
    
    Example:
        refine_dir = setup_refinement_dir('refine1', 'fitMaster')
    """

    refine_dir = os.path.join(fit_master_dir, refine_name)

    if not os.path.exists(refine_dir):
        try:
            os.makedirs(refine_dir)
        except OSError as e:
            raise OSError(f"Failed to create fit refinement directory: {e}")

    return refine_dir


def setup_runs(refine_dir, number_runs = 3):
    """
    Creates a number of run directories within the refinement directory.

    Parameters:
        refine_dir (str): The path of the refinement directory.
        number_runs (int): The number of run directories to create. Default is 3.

    Returns:
        list: A list of run directories.

    Raises:
        OSError: If a run directory cannot be created.

    Example:
        run_dirs = setup_runs('path/to/refine/directory', 3)
    """

    run_dirs = []

    for i in range(1, number_runs+1):
        run_name = 'run' + str(i)
        run_dir = os.path.join(refine_dir, run_name)

        if not os.path.exists(run_dir):
            try:
                os.makedirs(run_dir)
            except OSError as e:
                raise OSError(f"Failed to create run_dir {i}: {e}")

        run_dirs.append(run_dir)

    return run_dirs



# Writing to carbonara format

def write_coordinates_file(coords, working_path, carb_index=1):

    assert type(coords).__module__ == np.__name__, 'Thats never good... the CA coordinates are not a numpy array'

    file_write_name = working_path+'/coordinates' + str(carb_index) + '.dat'
    np.savetxt(file_write_name, coords, delimiter=' ', fmt='%s',newline='\n', header='', footer='')

    return file_write_name


def write_fingerprint_file(number_chains, sequence, secondary_structure, working_path):

    assert isinstance(number_chains, int), 'Yikes... The number of chains is not int type!'

    if number_chains > 1:
        print('Are sure you have more than one chain - if not this will cause segmentation errors later! You have been warned...')

    #seq_run = ''.join(list(sequence))
    #ss_run = ''.join(list(secondary_structure))

    # if len(seq_run) != len(ss_run):
    #    raise Exception("Uh Oh... The length of sequence and secondary structure is not equal!")
    file_name_path = working_path+"/fingerPrint1.dat"
    f = open(file_name_path, "w")
    f.write(str(number_chains))
    for i in range(0,number_chains):
        seq_run =''.join(list(sequence[i]))
        ss_run = ''.join(list(secondary_structure[i]))
        if len(seq_run) != len(ss_run):
            raise Exception("Uh Oh... The length of sequence and secondary structure is not equal!")
        f.write('\n \n')
        f.write(seq_run)
        f.write('\n \n')
        f.write(ss_run)
    f.close()

    _replicate_numbered_file_from_1(file_name_path, "fingerPrint")
    return file_name_path


def write_varysections_file(varying_sections, working_path, carb_index=1):
    # auto: run beta sheet breaking code; write output sections to file
    file_write_name = working_path+"/varyingSectionSecondary" + str(carb_index) + ".dat"
    f = open(file_write_name, "w")

    for i, s in enumerate(varying_sections):
        f.write(str(s))

        if i < len(varying_sections)-1:
            f.write('\n')
    f.close()

    if int(carb_index) == 1:
        _replicate_numbered_file_from_1(file_write_name, "varyingSectionSecondary")

    return file_write_name


def write_mixture_file(working_path):

    '''NOT CURRENTLY USEFUL FOR MULTIMERS YET'''

    file_path_name = working_path+"/mixtureFile.dat"
    # if default:
    f = open(file_path_name, "w")
    f.write(str(1))

    return file_path_name

#     else:
#          copy input file

def write_saxs(SAXS_file, working_path):
    with open(SAXS_file) as oldfile, open('temp.txt', 'w') as newfile:
        for line in oldfile:
            final_list = []
            for elem in line.split():
                try:
                    float(elem)
                except ValueError:
                    final_list.append(elem)
            if len(final_list)==0:
                newfile.write(line)

    saxs_arr = read_triplets_from_file('temp.txt')

    if saxs_arr.shape[1] == 3:
        saxs_arr = saxs_arr[:,:3]

    #check if it is in angstoms, if the last value is >1 we assume its in nanometers.
    if saxs_arr[-1,0] >1:
        for i in range(0,len(saxs_arr)):
            saxs_arr[i,0]=saxs_arr[i,0]/10.0

    file_path_name = working_path+'/Saxs.dat'
    np.savetxt(file_path_name, saxs_arr, delimiter=' ', fmt='%s',newline='\n', header='', footer='')
    os.remove("temp.txt")

    return file_path_name


def get_sses(ss_file):
    '''
    From a SS fingerprint, returns the SS with the subsection length, eg. ---SSS--- becomes [[-,3],[S,3],[-,3]]]
    '''
    lines = []
    with open(ss_file,'r') as fin:
        for line in fin:
            lines+= [line.split()]
    ss_tensor=[]
    for i in range(4,4*int(lines[0][0])+1,4):
        ss = lines[i][0]
        sses = []
        count = 1
        i = 0
        while i<len(ss)-1:
            if ss[i+1] == ss[i]:
                count += 1
                i += 1
            else:
                sses.append([ss[i], count])
                count = 1
                i += 1
        sses.append([ss[-1], count])
        ss_tensor.append(sses)
    return ss_tensor


def read_triplets_from_file(filename):
    data = []
    with open(filename, 'r') as f:
        for line in f:
            # Try converting the line to 3 floats
            try:
                numbers = list(map(float, line.strip().split()))
                if len(numbers)>1:
                    data.append(numbers)
            except ValueError:
                continue  # Skip non-numeric lines
    return np.array(data)
    
def SAXS_fit_plotter(SAXS_file, fit_file, full_q=True):

    fig = make_subplots(rows=2, cols=1,row_heights=[0.7,0.3],vertical_spacing=0,shared_xaxes=True)

    SAXS = read_triplets_from_file(SAXS_file)

    fitting = np.genfromtxt(fit_file, skip_footer=1)
    fit_q = fitting[:,0]
    fit_I = fitting[:,2]

    q = SAXS[:,0]
    I = np.log(np.where(SAXS[:,1] <= 0, np.nan, SAXS[:,1]))

    min_q = fit_q.min()
    max_q = fit_q.max()

    if q[-1] >1:
        for i in range(0,len(q)):
            q[i]=q[i]/10.0


    cond = (q >= min_q) & (q <= max_q)
    q_range = q[cond]
    I_range = I[cond]

    
    tck = interpolate.splrep(fit_q, fit_I)
    spli_I = interpolate.splev(q_range,tck)

    #     q_selected = q[(q>=min_q)&(q<=max_q)]
    #     q_grey = q[(q<min_q) | (q>=max_q)]

    #     I_selected = I[(q>=min_q)&(q<=max_q)]
    #     I_grey = I[(q<min_q) | (q>=max_q)]

    # fig = go.Figure(

    residuals = spli_I - I_range

    if full_q:
        fig.add_trace( go.Scatter(x=q, y=I, mode='markers', line=dict(color="grey"), opacity=0.7, name='Data'),row=1,col=1 )

    else:
        fig.add_trace( go.Scatter(x=q_range, y=I_range, mode='markers', line=dict(color="grey"), opacity=0.7, name='Data'),row=1,col=1 )

    fig.add_trace( go.Scatter(x=fit_q, y=fit_I, mode='markers',
                        marker=dict( color='crimson', size=8),
                         name='Fit'),row=1,col=1 )

    fig.add_trace( go.Scatter(x=q_range, y=spli_I, mode='lines', line=dict(color="crimson", width=3), name='Fit'),row=1,col=1 )

    fig.add_trace(go.Scatter(x=q_range,y=residuals,mode='lines',name='Residual',showlegend=False,line=dict(color='red')),row=2,col=1)
    fig.add_trace(go.Scatter(x=q_range,y=np.zeros_like(q_range),mode='lines',showlegend=False,line=dict(color='black',dash='dash',width=1)),row=2,col=1)

    fig.update_layout(
        # title='Experiment vs Fit',
                      # yaxis_type = "log",
                    template='simple_white',
                    # width=1200, height=800,
                    font_size=18)

    fig.update_yaxes(title_text="Intensity I(q)", row=1, col=1)
    fig.update_yaxes(title_text="Residual", row=2, col=1)
    fig.update_xaxes(title_text="q", row=2, col=1)

    # max_res = max( np.abs(residuals).max(), .5)
    max_res = np.abs(residuals).max()*1.3
    fig.update_yaxes(range=[-max_res,max_res],row=2,col=1)
    fig.update_layout(margin=dict(l=0, r=0, t=0, b=0))
    fig.update_traces(showlegend=False)
    return fig

def highlightVaryingSections(MolPath,PDB_fl,varyingSections,chain=1):
    resids = getResIDs_from_structure(PDB_fl, MolPath + '/fingerPrint1.dat')
    ss = get_sses(MolPath+'/fingerPrint1.dat')[chain-1]
    cols=[]
    varcols = []
    sscoldict = {'H': 'rgb(240,0,128)', 'S':'rgb(255,255,0)', '-': 'grey'}
    fp=[]
    for i in range(len(ss)):
        for j in range(ss[i][1]):
            fp.append(ss[i][0])
            if i in varyingSections:
                varcols.append('red')
            else:
                varcols.append('black')
    sscols = [sscoldict[i] for i in fp]
    coords_chains = pull_structure_from_pdb(PDB_fl)[0]
    mol = coords_chains[chain - 1]
    for coords in coords_chains:
        breaking_indices = missing_ca_check(coords)
        if len(breaking_indices) > 0:
            coords_chains = np.array_split(coords_chains,breaking_indices)
    mol = coords_chains[chain-1]
    hover_texts = ['ResID: '+str(resids[chain-1][i]) + ', SS: ' + list(fp)[i] for i in range(len(fp))]
    fig = go.Figure()
    fig.add_trace(go.Scatter3d(
            x=mol[:,0], 
            y=mol[:,1], 
            z=mol[:,2],
            text=hover_texts,
            hoverinfo='text',
            name='Varying Sections',
            marker=dict(
                size=1,
                color=varcols,
            ),
            line=dict(
                color=varcols,
                width=10
            )
        ))
    fig.add_trace(go.Scatter3d(
            x=mol[:,0], 
            y=mol[:,1], 
            z=mol[:,2],
            text=hover_texts,
            hoverinfo='text',
            visible='legendonly',
            name='Secondary Structure',
            marker=dict(
                size=1,
                color=sscols,
            ),
            line=dict(
                color=sscols,
                width=12.5
            )
        ))
    fig.update_layout(width=1000,height=1000)
    fig.update_layout(margin=dict(l=0, r=0, t=0, b=0))
    fig.update_layout(
        showlegend=True,
        legend=dict(x=0),
        scene=dict(
            xaxis_title='',
            yaxis_title='',
            zaxis_title='',
            aspectratio = dict( x=1, y=1, z=1 ),
            aspectmode = 'manual',
            xaxis = dict(
                gridcolor="white",
                showbackground=False,
                zerolinecolor="white",
                nticks=0,
                showticklabels=False),
            yaxis = dict(
                gridcolor="white",
                showbackground=False,
                zerolinecolor="white",
                nticks=0,
                showticklabels=False),
            zaxis = dict(
                gridcolor="white",
                showbackground=False,
                zerolinecolor="white",
                nticks=0,
                showticklabels=False),),
    )
    return fig


def getResIDs(pdb_fl):
    M = pdb_2_biobox(pdb_fl)
    ca_idx = (M.data['name']=='CA').values
    resids = M.get_data(indices=ca_idx)
    coords_chains_in,sequence_chains_in,_,_= pull_structure_from_pdb(pdb_fl)
    resids = resids[:,5]
    
    # >> check for missing residues geometrically in the chain(s)
    for j in range(len(coords_chains_in)):
        breaking_indices = missing_ca_check(coords_chains_in[j])
        if len(breaking_indices) > 0:
            coords_chains_in[j] = break_into_chains(coords_chains_in[j],sequence_chains_in[j],breaking_indices)[0]
    coords_chains = []
    for i in range(len(coords_chains_in)):
        if isinstance(coords_chains_in[i],list):
            for j in range(len(coords_chains_in[i])):
                coords_chains.append(coords_chains_in[i][j])
        else:
            coords_chains.append(coords_chains_in[i])
    idx = 0
    resid_tensor=[]
    for i in range(len(coords_chains)):
        resid_tensor.append(resids[idx:idx+len(coords_chains[i])])
        idx+=len(coords_chains[i])
    return resid_tensor

def groupResIDs(pdb_fl,fp_fl,chain=1):
    resids = getResIDs(pdb_fl)
    ss = get_sses(fp_fl)
    grouped = []
    idx=0
    for i in ss[chain-1]:
        grouped.append([i[0],resids[chain-1][idx:idx+i[1]][0],resids[chain-1][idx:idx+i[1]][-1]])
        idx+=i[1]
    return(grouped)


def possibleLinkerList(pdb_fl,fp_fl,chain=1):
    grouped = groupResIDs(pdb_fl,fp_fl,chain)
    poss = []
    for i in range(len(grouped)):
        if grouped[i][0]=='-':
            poss.append([i,'ResID: ' + str(grouped[i][1]) + '-' + str(grouped[i][2])])
    return np.array(poss)


def linkerLengthCheck(resid_str):
    res_range = resid_str.split(' ')[-1]
    start = int(res_range.split('-')[0])
    end = int(res_range.split('-')[1])
    if end-start<3:
        return False
    else:
        return True

def toggle_paired_predictions(script_path):
    with open(script_path, 'r') as f:
        lines = f.readlines()

    with open(script_path, 'w') as f:
        for line in lines:
            if line.strip().startswith("pairedPredictions="):
                line = "pairedPredictions=True\n"
                f.write(line)
            else:
                f.write(line)

def toggle_startk(script_path, kmaxstart):
    with open(script_path, 'r') as f:
        lines = f.readlines()

    with open(script_path, 'w') as f:
        for line in lines:
            if line.strip().startswith("kmaxStart="):
                f.write(f"kmaxStart={kmaxstart}\n")
            else:
                f.write(line)

def toggle_kmin(script_path,kmin):
    with open(script_path, 'r') as f:
        lines = f.readlines()

    with open(script_path, 'w') as f:
        for line in lines:
            if line.strip().startswith("kmin="):
                f.write(f"kmin={kmin}\n")
            else:
                f.write(line)

def read_triplets_from_file(filename):
    data = []
    with open(filename, 'r') as f:
        for line in f:
            # Try converting the line to 3 floats
            try:
                numbers = list(map(float, line.strip().split()))
                if len(numbers) >1:
                    data.append(numbers)
            except ValueError:
                continue  # Skip non-numeric lines
    return np.array(data)

def translate_distance_constraints(contactPredsIn,coords,working_path,fixedDistList=[]):
    # shift the coordinates back one to fit [0,1, array labelling
    contactPreds =contactPredsIn
    dists= []
    for i in range(len(contactPredsIn)):
        contactPreds[i][0] = contactPredsIn[i][0]
        contactPreds[i][1] = contactPredsIn[i][1]
    contactPredNara = []
    for i in range(len(contactPreds)):
        contactPreds[i].sort()
        currIndex=0;
        ss = get_secondary(working_path+"/fingerPrint1.dat")
        sections = section_finder_sub(ss)
        currMax=len(sections[0])
        prevMax=0
        while (contactPreds[i][0]>currMax and currIndex<len(sections)):
            currIndex= currIndex+1
            currMax=currMax+len(sections[currIndex])
            prevMax = prevMax+len(sections[currIndex-1])
           # second coord of pair
        pair1 =[currIndex,contactPreds[i][0]-prevMax-1]
        currIndex=0;
        currMax=len(sections[0])
        prevMax=0
        while (contactPreds[i][1]>currMax and currIndex<len(sections)):
            currIndex= currIndex+1
            currMax=currMax+len(sections[currIndex])
            prevMax = prevMax+len(sections[currIndex-1])
        pair2 =[currIndex,contactPreds[i][1]-prevMax-1]
        if len(fixedDistList)>0:
            dist = fixedDistList[i]
        else:
            dist = np.linalg.norm(coords[contactPreds[i][1]-1]-coords[contactPreds[i][0]-1])
        # contactPredNara.append(pair1+pair2+[dist])
        contactPredNara.append(pair1+pair2+[dist]+[0.5])
        dists.append(dist)

        # now write to file
    np.savetxt(working_path+"/fixedDistanceConstraints1.dat",contactPredNara,fmt="%i %i %i %i %1.10f %1.10f")


def get_secondary(fingerprint_file):
    # Read as plain text; stable across all NumPy versions
    with open(fingerprint_file, "r") as f:
        # Drop empty / whitespace-only lines (these are what break np.loadtxt)
        fplst = [ln.strip() for ln in f if ln.strip() != ""]
    # Now fplst is a simple Python list of lines, always iterable
    # Preserve your original indexing logic: 2,4,6,... (0-based)
    fplstout = [np.asarray(list(fplst[i])) for i in range(2, len(fplst), 2)]
    return fplstout

def flatten_extend(matrix):
     flat_list = []
     for row in matrix:
         flat_list.extend(row)
     return flat_list


def read_coords_from_file(filename):
    data = []
    with open(filename, 'r') as f:
        for line in f:
            # Try converting the line to 3 floats
            try:
                numbers = list(map(float, line.strip().split()))
                if len(numbers) ==3:
                    data.append(numbers)
            except ValueError:
                continue  # Skip non-numeric lines
    return np.array(data)

def plotMolAndSAXS(RunPath,saxs_fl,mol_fl):
    '''
    Combined plot of a given mol and saxs fit.
    '''
    fig = make_subplots(rows=2, cols=1,row_heights=[0.7,0.3],vertical_spacing=0, specs=[[{'type': 'scene'}], [{'type': 'xy'}]])
    ### First plot mol
   # Load coordinates
    coords =read_coords_from_file(mol_fl)

    # Load secondary structure list
    sec = flatten_extend(get_secondary(RunPath + "fingerPrint1.dat"))

    # Ensure lengths match
    assert len(coords) == len(sec), "Length of coordinates and secondary structure must match"

    # Color and opacity logic
    def get_style(a, b):
        if 'H' in (a, b):
            return 'red', 1.0
        elif 'S' in (a, b):
            return 'green', 1.0
        else:
            return 'blue', 0.3  # transparent coil

    # Add line segments conditionally
    for i in range(len(coords) - 1):
        point_a = coords[i]
        point_b = coords[i + 1]
        distance = np.linalg.norm(point_a - point_b)
    
        if distance > 5:
            continue  # skip long segments

        color, opacity = get_style(sec[i], sec[i+1])
    
        fig.add_trace(go.Scatter3d(
            x=[point_a[0], point_b[0]],
            y=[point_a[1], point_b[1]],
            z=[point_a[2], point_b[2]],
            mode='lines',
            line=dict(color=color, width=10),
            opacity=opacity,
            showlegend=False
        ), row=1, col=1)

    # Final layout
    fig.update_layout(width=1000, height=1000)
    ### Now plot SAXS
    SAXS = read_triplets_from_file(RunPath+"Saxs.dat")
    fitting = np.genfromtxt(saxs_fl, skip_footer=1)
    fit_q = fitting[:,0]
    fit_I = fitting[:,2]

    q = SAXS[:,0]
    I = np.log(SAXS[:,1])

    min_q = fit_q.min()
    max_q = fit_q.max()



    cond = (q >= min_q) & (q <= max_q)
    q_range = q[cond]
    I_range = I[cond]
    kval  = np.min(np.array([len(fit_I)-1,5]))
    tck = interpolate.splrep(fit_q, fit_I,k=int(kval))
    spli_I = interpolate.splev(q_range,tck)

    residuals = spli_I - I_range
    residuals = residuals/max(residuals)

    fig.add_trace(go.Scatter(x=q,y=I, mode='markers', line=dict(color="grey"), opacity=0.7, name='Data'),row=2,col=1 )
    fig.add_trace(go.Scatter(x=fit_q, y=fit_I, mode='markers', marker=dict( color='crimson', size=8), name='Fit'),row=2,col=1 )
    fig.add_trace( go.Scatter(x=q_range, y=spli_I, mode='lines', line=dict(color="crimson", width=3), name='Fit'),row=2,col=1 )

    fig.update_layout(
                    template='simple_white',
                    font_size=18)

    fig.update_yaxes(title_text="Intensity I(q)", row=2, col=1)
    fig.update_xaxes(title_text="q", row=2, col=1)
    fig.update_layout(margin=dict(l=0, r=0, t=0, b=0))
    fig.update_traces(showlegend=False)
    fig.update_layout(
        showlegend=False,
        scene=dict(
            xaxis_title='',
            yaxis_title='',
            zaxis_title='',
            aspectratio = dict( x=1, y=1, z=1 ),
            aspectmode = 'manual',
            xaxis = dict(
                visible=False,
                showbackground=False,
                showticklabels=False,
                ),
            yaxis = dict(
                visible=False,
                showbackground=False,
                showticklabels=False),
            zaxis = dict(
                visible=False,
                showbackground=False,
                showticklabels=False)))
    return fig.show()


def load_pae_matrix(json_path):
    """
    Load a Predicted Aligned Error (PAE) matrix from a JSON file.
    
    Supports AlphaFold DB format (list of dicts with 'predicted_aligned_error')
    and AlphaFold v3/custom format (dict with 'pae' or 'predicted_aligned_error' key).
    Returns:
        numpy.ndarray: 2D array of PAE values.
    Raises:
        ValueError: If no supported PAE format is found in the JSON.
    """
    # Read the JSON file
    with open(json_path, 'r') as f:
        data = json.load(f)
    
    # Determine the format and extract the PAE matrix
    if isinstance(data, list):
        # Format 1: AlphaFold DB style (list containing dict with PAE data)
        if not data:
            raise ValueError("PAE JSON list is empty.")
        entry = data[0]
        if isinstance(entry, dict):
            # Check known keys for PAE matrix
            if 'predicted_aligned_error' in entry:
                matrix = entry['predicted_aligned_error']
            elif 'predicted_alignment_error' in entry:
                matrix = entry['predicted_alignment_error']
            elif 'pae' in entry:
                matrix = entry['pae']
            else:
                raise ValueError("No PAE matrix found under expected keys in the JSON list entry.")
        else:
            raise ValueError("PAE JSON list does not contain a dictionary object.")
    
    elif isinstance(data, dict):
        # Format 2: AlphaFold v3 or custom style (PAE under top-level keys in a dict)
        if 'predicted_aligned_error' in data:
            matrix = data['predicted_aligned_error']
        elif 'predicted_alignment_error' in data:
            matrix = data['predicted_alignment_error']
        elif 'pae' in data:
            matrix = data['pae']
        else:
            raise ValueError("No supported PAE key ('predicted_aligned_error' or 'pae') found in the JSON file.")
    
    else:
        # Unsupported JSON structure
        raise ValueError("Unsupported JSON format for PAE data.")
    
    # Convert the matrix to a NumPy array and return
    return np.array(matrix)

# def getFlexibleSections(file_path,pae_threshold = 1.5): 
#     pae_data = load_pae_matrix(file_path)
#     if len(pae_data)==0:
#         print("No 'pae' data found in the JSON file.")
#     else:
#         # Convert PAE data to a NumPy array for easier manipulation
#         pae_array = np.array(pae_data)
#         window_size = 2  # Define a window size
#         high_pae_residues = []
#         for i in range(len(pae_array) - window_size + 1):
#             window = pae_array[i:i + window_size, i:i + window_size]
#             avg_pae = np.mean(window)
#             print(avg_pae)
#             if avg_pae > pae_threshold:
#                 #print(f"Residue range {i}-{i + window_size} exceeds threshold with avg PAE {avg_pae}")
#                 high_pae_residues.extend(range(i, i + window_size))
#     # Remove duplicate residue indices
#         high_pae_residues = sorted(set(high_pae_residues))
#         # Plotting the PAE heatmap with inverted colors
#         plt.figure(figsize=(8, 6))
#         heatmap = plt.imshow(pae_array, cmap='Greens_r', interpolation='nearest', aspect='auto')
#         plt.colorbar(heatmap, label='Expected position error (Ångströms)')
#         # Highlight regions identified as linkers
#         for res in high_pae_residues:
#             plt.axvline(x=res, color='red', linestyle='--', alpha=0.3)
#             plt.axhline(y=res, color='red', linestyle='--', alpha=0.3)
#         # Set labels
#         plt.xlabel('Scored residue')
#         plt.ylabel('Aligned residue')
#         plt.title('Predicted Aligned Error (PAE) Heatmap')

#         # Show the plot
#         plt.show()
#         # Second plot: Highlighting selected residues
#         plt.figure(figsize=(10, 2))
#         plt.bar(range(len(pae_array)), [1 if i in high_pae_residues else 0 for i in range(len(pae_array))], color='green', alpha=0.6)
#         plt.xlabel('Residue Index')
#         plt.ylabel('')
#         plt.title('Residues Identified as Linkers or Between Domains')
#         plt.yticks([0, 1], [])
#         plt.grid(axis='x')

#         # Show the second plot
#         plt.show()
#         return np.array([1 if i in high_pae_residues else 0 for i in range(len(pae_array))])



def getFlexibleSections(file_path, pae_threshold=16.0):
    """
    Returns: np.array of 0/1, length L.

    Decision logic (not percentage-based):
      - Compute off-diagonal mean PAE score per residue
      - Flag if score exceeds:
          * absolute Å threshold if pae_threshold >= 2 (recommended: 8–12)
          * median + (pae_threshold * MAD) if pae_threshold < 2 (use ~2.0–3.0)
    """
    pae = np.asarray(load_pae_matrix(file_path), dtype=float)
    if pae.ndim != 2 or pae.shape[0] != pae.shape[1]:
        raise ValueError(f"PAE matrix is not square: {pae.shape}")
    L = pae.shape[0]

    exclude_diag = 20
    smooth_w = 9

    score = np.zeros(L, dtype=float)
    for i in range(L):
        far = np.ones(L, dtype=bool)
        far[max(0, i-exclude_diag):min(L, i+exclude_diag+1)] = False
        score[i] = float(pae[i, far].mean()) if far.any() else float(pae[i].mean())

    if smooth_w > 1:
        score = np.convolve(score, np.ones(smooth_w)/smooth_w, mode="same")

    # Threshold: absolute Å, or robust "median + k*MAD"
    if pae_threshold >= 2.0:
        thr = float(pae_threshold)
    else:
        # interpret as k for MAD; practical values are ~2.0–3.0
        k = float(pae_threshold)
        med = float(np.median(score))
        mad = float(np.median(np.abs(score - med))) + 1e-12
        thr = med + k * mad
    flags = (score >= thr).astype(int)
    return flags

import numpy as np
from typing import List, Set, Union

def find_flexible_linker_sections_multi(
    ss_chains: List[Union[str, np.ndarray]],
    pae_flags: Union[List[int], np.ndarray],
    min_linker_len: int = 3,   # <--- new rule
) -> Set[int]:
    """
    - Respects chain boundaries: linker segments cannot merge across chains.
    - Only selects '-' segments of length >= min_linker_len.
    - Returns global section indices (as before).
    """
    pae_flags = np.asarray(pae_flags).astype(int)
    if pae_flags.ndim == 0:
        raise ValueError("pae_flags is scalar; expected 1-D.")
    pae_flags = np.atleast_1d(pae_flags)

    # normalize ss_chains to list of Python strings
    ss_strings: List[str] = []
    for ss in ss_chains:
        if isinstance(ss, np.ndarray):
            ss_strings.append("".join(ss.tolist()))
        else:
            ss_strings.append(str(ss))

    total_L = sum(len(s) for s in ss_strings)
    if total_L != pae_flags.shape[0]:
        raise ValueError(f"Length mismatch: sum(ss_chains)={total_L} vs pae_flags={pae_flags.shape[0]}")

    flexible_linker_indices: Set[int] = set()

    r0 = 0
    section_index = 0

    for ss in ss_strings:
        Lc = len(ss)
        flags_c = pae_flags[r0:r0 + Lc]
        r0 += Lc

        i = 0
        while i < Lc:
            current_char = ss[i]
            start = i
            while i < Lc and ss[i] == current_char:
                i += 1
            end = i

            if current_char == '-':
                seg_len = end - start
                if seg_len >= min_linker_len and flags_c[start:end].any():
                    flexible_linker_indices.add(section_index)

            section_index += 1

    return flexible_linker_indices

def drop_zero(int_set):
    int_set.discard(0)   # removes 0 if present, does nothing otherwise
    return int_set

def getFlexibility(paeFile,fingerprint_file, abs_thr=16.0):
    # if file is in noy format convert to json   
    if paeFile.endswith('.npy'):
        # Load the .npy file
        pae_matrix = np.load('paerank_2.npy')
        
        # Convert to list of lists for JSON
        pae_list = pae_matrix.tolist()
        
        # Define the JSON structure
        pae_json = {
            "predicted_aligned_error": pae_list,
            "max_predicted_aligned_error": float(np.max(pae_matrix))
        }

        # Save as JSON
        with open('paerank_2_converted.json', 'w') as f:
            json.dump(pae_json, f)
        flexsec =getFlexibleSections('paerank_2_converted.json',pae_threshold=abs_thr)
        fingerprint = get_secondary(fingerprint_file)
        return [list(sorted(drop_zero(find_flexible_linker_sections_multi(fingerprint, flexsec))))]
    else:
        flexsec =getFlexibleSections(paeFile,pae_threshold=abs_thr)
        fingerprint = get_secondary(fingerprint_file)
        return [list(sorted(drop_zero(find_flexible_linker_sections_multi(fingerprint, flexsec))))]


    

def parse_structures_with_segments(filename):
    with open(filename, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]

    n = int(lines[0])
    sequences = lines[1::2]
    structures = lines[2::2]

    chains = []
    segment_index = 0  # global segment counter

    for seq, struct in zip(sequences, structures):
        segments = []
        for match in re.finditer(r'(-+|S+|H+)', struct):
            start, end = match.start(), match.end()
            label = match.group()
            segments.append({
                'start': start,
                'end': end,
                'type': label[0],
                'segment_num': segment_index
            })
            segment_index += 1

        chains.append({
            'sequence': seq,
            'structure': struct,
            'segments': segments
        })

    return chains

def colorize(text, color='red'):
    # ANSI red for terminal; in Colab can also use HTML span + display(HTML(...))
    return f"\033[91m{text}\033[0m"  # red

def print_structure_with_highlights(chains, highlight_segments):
    for i, chain in enumerate(chains):
        print(f"\n=== Chain {i+1} ===")
        struct_colored = list(chain['structure'])  # mutable char list

        for segment in chain['segments']:
            if segment['segment_num'] in highlight_segments:
                for pos in range(segment['start'], segment['end']):
                    struct_colored[pos] = colorize(struct_colored[pos])

        print("Structure:")
        print("".join(struct_colored))
        print("Sequence:")
        print(chain['sequence'])



def merge_chains_robust_all(chains, merge_indices):
    """
    Merge specified chains, remap segment numbers globally across all chains.
    Returns:
        new_chains: list of updated chain dicts
        segment_map: dict of old segment number -> new segment number
    """
    from collections import defaultdict, Counter

    merge_set = set(idx - 1 for idx in merge_indices)
    new_chains = []
    old_to_new_segment = {}
    merged_seq = ''
    merged_struct = ''
    merged_segment_ids = []
    merged_original_segments = []

    segment_counter = 0  # For assigning new segment numbers
    segment_offset = 0   # How many segments we've added so far (to shift future chains)

    # Phase 1: Process and build merged chain
    for i, chain in enumerate(chains):
        if i in merge_set:
            struct = chain['structure']
            segment_ids = [None] * len(struct)

            for seg in chain['segments']:
                for j in range(seg['start'], seg['end']):
                    segment_ids[j] = seg['segment_num']

            merged_seq += chain['sequence']
            merged_struct += struct
            merged_segment_ids.extend(segment_ids)
            merged_original_segments.extend(chain['segments'])

    # Recompute segments for merged chain
    merged_segments = []
    char_to_new_segment = [None] * len(merged_struct)
    for match in re.finditer(r'(-+|S+|H+)', merged_struct):
        start, end = match.start(), match.end()
        type_char = match.group()[0]
        merged_segments.append({
            'start': start,
            'end': end,
            'type': type_char,
            'segment_num': segment_counter
        })
        for i in range(start, end):
            char_to_new_segment[i] = segment_counter
        segment_counter += 1

    # Map old segment numbers from merged chains
    seg_votes = defaultdict(list)
    for old_id, new_id in zip(merged_segment_ids, char_to_new_segment):
        if old_id is not None:
            seg_votes[old_id].append(new_id)

    for old_id, new_ids in seg_votes.items():
        most_common = Counter(new_ids).most_common(1)[0][0]
        old_to_new_segment[old_id] = most_common

    # Add merged chain
    new_chains.append({
        'sequence': merged_seq,
        'structure': merged_struct,
        'segments': merged_segments
    })

    segment_offset = segment_counter  # update offset for next chains

    # Phase 2: Process unmerged chains and shift their segment numbers
    for i, chain in enumerate(chains):
        if i in merge_set:
            continue  # already handled

        new_segments = []
        for seg in chain['segments']:
            new_seg = seg.copy()
            new_seg['segment_num'] = seg['segment_num'] + segment_offset
            new_segments.append(new_seg)
            old_to_new_segment[seg['segment_num']] = new_seg['segment_num']

        new_chains.append({
            'sequence': chain['sequence'],
            'structure': chain['structure'],
            'segments': new_segments
        })

        segment_offset += len(new_segments)  # accumulate for next chain

    boundary_segments = set()

    # Identify boundary segments before any remapping
    for chain in chains:
        if chain['segments']:
            boundary_segments.add(chain['segments'][0]['segment_num'])      # first
            boundary_segments.add(chain['segments'][-1]['segment_num'])     # last

    # ⬅ rest of the function stays the same
    return new_chains, old_to_new_segment, boundary_segments
    
def filter_boundary_segments(updated_segments, boundary_segments, reverse_map=None):
    """
    Removes any segments that were originally the first or last in their chains.
    
    Args:
        updated_segments (set of int): segments after merging.
        boundary_segments (set of int): original boundary segment numbers.
        reverse_map (dict): optional, maps new segment → original segment(s)
        
    Returns:
        kept_segments (set): filtered updated_segments with boundaries removed
        removal_reasons (list of str): explanation for each removed segment
    """
    kept = set()
    reasons = []

    for seg in updated_segments:
        # If reverse_map is provided, check if *any* of the originals were boundaries
        originals = reverse_map.get(seg, [seg]) if reverse_map else [seg]
        if any(orig in boundary_segments for orig in originals):
            reasons.append(
                f"Removed segment {seg} (from original segment(s) {originals}) — at chain boundary"
            )
        else:
            kept.add(seg)

    return kept, reasons


def update_modified_segments(original_segment_list, mapping):
    return {mapping.get(s, s) for s in original_segment_list}

# If you are happy with this merging then set it in place

def export_chains_to_file(chains, filename):
    """
    Export a list of chains (with 'sequence' and 'structure') to a file
    in the original alternating format.

    If this writes fingerPrint1.dat inside a multi-structure Carbonara run,
    the edited file is copied to fingerPrint2.dat, fingerPrint3.dat, ... so
    all structure copies remain consistent with the notebook edits.
    """
    with open(filename, 'w') as f:
        f.write(f"{len(chains)}\n\n")  # Write number of chains

        for chain in chains:
            f.write(f"{chain['sequence']}\n\n")
            f.write(f"{chain['structure']}\n\n")

    _replicate_numbered_file_from_1(filename, "fingerPrint")

def export_segment_list(segment_set, filename):
    """
    Write a set of segment numbers to a file, one per line, without trailing newline at the end.

    If this writes varyingSectionSecondary1.dat inside a multi-structure
    Carbonara run, the edited file is copied to varyingSectionSecondary2.dat,
    varyingSectionSecondary3.dat, ... so all structure copies use the same
    updated varying-section labels.
    """
    segments = sorted(segment_set)
    with open(filename, 'w') as f:
        for i, seg in enumerate(segments):
            if i < len(segments) - 1:
                f.write(f"{seg}\n")
            else:
                f.write(f"{seg}")  # last line, no newline

    _replicate_numbered_file_from_1(filename, "varyingSectionSecondary")

def _canonical_structure_for_fingerprint(pdb_fl, structure_file):
    """Return the structure that was actually used to create a run fingerprint.

    ``setup_carbonara_allAtom.py`` writes ``carbonara_structure_path.txt`` beside
    ``fingerPrint1.dat``.  This is necessary because setup is often run in a
    subprocess/imported module, so changing a ``pdb_name`` global inside setup
    cannot change an already-existing Jupyter variable in the caller.

    If no marker exists, preserve the historical behaviour and use ``pdb_fl``.
    """
    from pathlib import Path

    marker = Path(structure_file).expanduser().resolve().parent / "carbonara_structure_path.txt"
    if not marker.exists():
        return str(pdb_fl)

    recorded = marker.read_text().strip()
    if not recorded:
        raise ValueError(f"Canonical structure marker is empty: {marker}")

    canonical = Path(recorded).expanduser()
    if not canonical.is_absolute():
        canonical = (marker.parent / canonical).resolve()

    if not canonical.exists():
        raise FileNotFoundError(
            f"Canonical structure recorded for this Carbonara run does not exist: {canonical} "
            f"(marker: {marker})"
        )

    return str(canonical)


def canonical_structure_for_run(run_name, fallback_pdb=None):
    """Public helper for notebooks: return the run's canonical/repaired PDB path."""
    from pathlib import Path

    marker = Path("carbonara_runs") / str(run_name) / "carbonara_structure_path.txt"
    if marker.exists():
        recorded = marker.read_text().strip()
        if not recorded:
            raise ValueError(f"Canonical structure marker is empty: {marker}")
        canonical = Path(recorded).expanduser()
        if not canonical.is_absolute():
            canonical = (marker.parent / canonical).resolve()
        if not canonical.exists():
            raise FileNotFoundError(
                f"Canonical structure recorded for run {run_name!r} does not exist: {canonical}"
            )
        return str(canonical)

    if fallback_pdb is None:
        raise FileNotFoundError(
            f"No carbonara_structure_path.txt found for run {run_name!r}, and no fallback PDB was supplied."
        )
    return str(fallback_pdb)


def getResIDs_from_structure(pdb_fl, structure_file):
    """
    Returns resid_tensor split to match the chains in the structure file,
    ignoring original PDB chain IDs.

    For run fingerprints, automatically use the canonical/repaired PDB recorded
    by setup rather than a stale original ``pdb_name`` from a notebook.
    """
    import re

    pdb_fl = _canonical_structure_for_fingerprint(pdb_fl, structure_file)

    # Get all CA atoms
    M = pdb_2_biobox(pdb_fl)
    ca_idx = (M.data['name'] == 'CA').values
    resids = M.get_data(indices=ca_idx)[:, 5]  # residue numbers

    # Read structure file and get lengths of each structure line
    with open(structure_file, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]
    n = int(lines[0])
    structures = lines[2::2]
    lengths = [len(struct) for struct in structures]

    if len(structures) != n:
        raise ValueError(
            f"Fingerprint says {n} chains, but contains {len(structures)} structure lines: {structure_file}"
        )

    expected = int(sum(lengths))
    observed = int(len(resids))
    if observed != expected:
        raise ValueError(
            "Fingerprint/PDB residue-count mismatch: "
            f"{structure_file} describes {expected} residues, but {pdb_fl} contains "
            f"{observed} CA atoms. This would otherwise produce an out-of-bounds linker index. "
            "Check that setup created carbonara_structure_path.txt and that it points to the "
            "repaired run-local structure."
        )

    # Split resids according to structure lengths
    resid_tensor = []
    idx = 0
    for length in lengths:
        resid_tensor.append(resids[idx:idx+length])
        idx += length

    return resid_tensor


def possibleLinkerList(fp_fl, pdb_fl, chain=1):
    """
    Reads a secondary structure file (fp_fl) and PDB (pdb_fl),
    and returns a list of dash-only ('-') segments for the given chain.

    Works whether the structure file is merged or not, assuming correct chain index is used.

    Returns:
        numpy array of [segment_number, 'ResID: start-end']
    """

    # Parse residue IDs from the PDB
    resid_tensor = getResIDs_from_structure(pdb_fl,fp_fl)
    resids = resid_tensor[chain - 1]

    # Read and parse the structure file
    with open(fp_fl, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]
    n_chains = int(lines[0])
    structures = lines[2::2]

    if chain > n_chains:
        raise IndexError(f"Chain {chain} requested but only {n_chains} chains in file.")

    structure = structures[chain - 1]
    dash_segments = []
    idx = 0
    segment_number = 0

    for match in re.finditer(r'(-+|S+|H+)', structure):
        length = match.end() - match.start()
        if match.group()[0] == '-':
            res_start = resids[idx]
            res_end = resids[idx + length - 1]
            dash_segments.append([segment_number, f"ResID: {res_start}-{res_end}"])
        idx += length
        segment_number += 1

    return np.array(dash_segments, dtype=object)

def getCoordsMatchingStructure(pdb_fl, structure_file):
    """
    Returns list of coordinate chains split to match structure file chains.
    Ignores PDB chain labels and assumes CA atoms in order.
    """
    M = pdb_2_biobox(pdb_fl)
    ca_idx = (M.data['name'] == 'CA').values
    coords = M.get_coord(indices=ca_idx)

    with open(structure_file, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]
    n = int(lines[0])
    structures = lines[2::2]
    lengths = [len(s) for s in structures]

    coords_chains = []
    idx = 0
    for length in lengths:
        coords_chains.append(coords[idx:idx+length])
        idx += length

    return coords_chains

def mapFixedConstraints(filename, chain_lengths):
    """
    Parse a constraints file and convert local residue indices into global indices.

    Args:
        filename (str): Path to the input file. Each line contains:
                        Chain1 Residue1 Chain2 Residue2 Value
        chain_lengths (dict): Mapping from chain letters (e.g., 'A', 'B') to chain lengths.

    Returns:
        constraint_pairs: List of [global_index1, global_index2]
        constraint_values: List of associated values (last number on each line)
    """
    chain_order = sorted(chain_lengths.keys())  # Assume alphabetical order
    chain_offsets = {}
    current_offset = 0

    print(chain_order)
    for ch in chain_order:
        chain_offsets[ch] = current_offset
        current_offset += chain_lengths[ch]

    constraint_pairs = []
    constraint_values = []

    with open(filename, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) == 5:
                res1, ch1, res2, ch2, value = parts
                global1 = chain_offsets[ch1] + int(res1)
                global2 = chain_offsets[ch2] + int(res2)
                constraint_pairs.append([global1, global2])
                constraint_values.append(int(value))
            elif len(parts) == 4:
                res1, ch1, res2, ch2 = parts
                global1 = chain_offsets[ch1] + int(res1)
                global2 = chain_offsets[ch2] + int(res2)
                constraint_pairs.append([global1, global2])

    return constraint_pairs, constraint_values
def merge_chains_only_clean_consistent_segments(chains, merge_pair):
    """
    Merge two chains (1-based indices) by concatenating their sequence and structure.
    If boundary segments are both '-', merge them into one.
    Recalculate segment numbers, start, and end positions globally after all adjustments.
    """
    import copy

    i, j = merge_pair[0] - 1, merge_pair[1] - 1
    new_chains = []

    for idx, chain in enumerate(chains):
        if idx == i:
            merged_chain = copy.deepcopy(chains[i])
            other_chain = chains[j]

            merged_sequence = merged_chain["sequence"] + other_chain["sequence"]
            merged_structure = merged_chain["structure"] + other_chain["structure"]

            merged_segments = copy.deepcopy(merged_chain["segments"])
            segs_j = copy.deepcopy(other_chain["segments"])

            if merged_segments[-1]["type"] == segs_j[0]["type"] == "-":
                merged_segments[-1]["end"] = merged_segments[-1]["end"] + (segs_j[0]["end"] - segs_j[0]["start"])
                segs_j = segs_j[1:]

            merged_segments.extend(segs_j)

            new_chains.append({
                "sequence": merged_sequence,
                "structure": merged_structure,
                "segments": merged_segments
            })
        elif idx == j:
            continue
        else:
            new_chains.append(copy.deepcopy(chain))

    # Recalculate start, end, and segment_num globally across all chains
    segment_num = 0
    for chain in new_chains:
        seq_idx = 0
        new_segments = []
        while seq_idx < len(chain["structure"]):
            seg_type = chain["structure"][seq_idx]
            start = seq_idx
            while seq_idx < len(chain["structure"]) and chain["structure"][seq_idx] == seg_type:
                seq_idx += 1
            new_segments.append({
                "start": start,
                "end": seq_idx,
                "type": seg_type,
                "segment_num": segment_num
            })
            segment_num += 1
        chain["segments"] = new_segments

    return new_chains


import copy

def create_segment_label_arrays_with_merge_v4(chains, highlighted_segments, merge_pair):
    """
    Create original and editable segment label arrays per chain, and remap highlighted segments.
    The chains themselves are not merged; only the labels and mapping are updated.

    Args:
        chains: list of chain dicts
        highlighted_segments: numpy array of segment numbers to highlight
        merge_pair: tuple (i, j) for 1-based chain indices to merge

    Returns:
        original_label_arrays: list of np.arrays of original segment labels per chain
        editable_label_arrays: list of np.arrays of updated segment labels per chain
        remapped_highlighted_segments: numpy array of updated highlighted segment numbers
    """
    import numpy as np
    import copy

    i, j = merge_pair[0] - 1, merge_pair[1] - 1  # Convert to 0-based indexing

    # Step 1: Create original label arrays
    original_label_arrays = []
    for ch in chains:
        segment_nums = [seg['segment_num'] for seg in ch['segments']]
        original_label_arrays.append(np.array(segment_nums, dtype=int))

    editable_label_arrays = copy.deepcopy(original_label_arrays)
    updated_highlighted_segments = set(highlighted_segments.tolist())

    # Step 2: Check if boundary segment should be removed
    last_seg_i = original_label_arrays[i][-1]
    first_seg_j = original_label_arrays[j][0]

    last_type_i = chains[i]['segments'][-1]['type']
    first_type_j = chains[j]['segments'][0]['type']

    removed_boundary = False
    if last_type_i == first_type_j == '-':
        if last_seg_i in updated_highlighted_segments:
            updated_highlighted_segments.remove(last_seg_i)
            removed_boundary = True

    # Step 3: Update labels in chain j (shifted to follow chain i)
    li = editable_label_arrays[i][-1]
    len_j = len(editable_label_arrays[j])
    editable_label_arrays[j] = np.arange(li + 1, li + 1 + len_j)

    # Step 4: Shift chains between i and j by len_j
    for k in range(i + 1, j):
        editable_label_arrays[k] += len_j

    # Step 5: If boundary removed, shift all chains after i down by 1
    if removed_boundary:
        for k in range(i + 1, len(editable_label_arrays)):
            editable_label_arrays[k] -= 1

    # Step 6: Build old-to-new segment map and remap highlighted segments
    old2new_segment_map = {}
    for orig_arr, edit_arr in zip(original_label_arrays, editable_label_arrays):
        for old, new in zip(orig_arr, edit_arr):
            old2new_segment_map[old] = new

    remapped_highlighted_segments = np.array(
        sorted(old2new_segment_map[s] for s in updated_highlighted_segments if s in old2new_segment_map)
    )

    return original_label_arrays, editable_label_arrays, remapped_highlighted_segments

from typing import List, Tuple

def read_secondary_structures_from_file(filepath: str) -> List[str]:
    """
    Reads a file with the custom format:
    <num_chains>
    <empty line>
    <sequence1>
    <empty line>
    <secondary1>
    ...
    Returns a list of secondary structure strings, one per chain.
    """
    with open(filepath, 'r') as f:
        lines = [line.strip() for line in f if line.strip() != '']

    num_chains = int(lines[0])
    ss_chains = []

    # Each chain consists of 2 lines (sequence, secondary), starting from index 1
    for i in range(num_chains):
        seq_index = 1 + i * 2
        ss_index = seq_index + 1
        ss_chains.append(lines[ss_index])

    return ss_chains

def parse_secondary_structure_chainwise(ss_chains: List[str]) -> List[Tuple[int, str, int, int]]:
    """
    Parses multiple secondary structure strings into segments across chains.
    Each segment is (segment_id, symbol, start, end) where [start:end] is the local position in the chain.
    Segment IDs are globally unique and increment continuously across chains.
    """
    segments = []
    current_id = 0
    for chain_ss in ss_chains:
        i = 0
        while i < len(chain_ss):
            symbol = chain_ss[i]
            start = i
            while i < len(chain_ss) and chain_ss[i] == symbol:
                i += 1
            end = i
            segments.append((current_id, symbol, start, end))
            current_id += 1
    return segments

def get_segment_lengths_from_file(filepath: str, target_segments: List[int], min_length: int = 3) -> List[int]:
    """
    Reads secondary structures from a file, computes segment lengths, and filters by min_length.
    Returns a list of segment IDs (from target_segments) that are >= min_length.
    """
    ss_chains = read_secondary_structures_from_file(filepath)
    segments = parse_secondary_structure_chainwise(ss_chains)
    id_to_length = {seg_id: end - start for seg_id, symbol, start, end in segments}
    return [seg_id for seg_id in target_segments if id_to_length.get(seg_id, 0) >= min_length]

def split_into_sections(ss):
    """
    Split a secondary structure sequence into contiguous sections.
    Example: ['H','H','-','-','S'] -> [['H','H'], ['-','-'], ['S']]
    """
    sections = []
    if len(ss) == 0:
        return sections

    current = [ss[0]]
    for c in ss[1:]:
        if c == current[-1]:
            current.append(c)
        else:
            sections.append(current)
            current = [c]
    sections.append(current)
    return sections


def build_global_sections(secondary):
    """
    Take list of secondary-structure arrays (one per chain) and flatten
    into a single global sections list across all chains.
    """
    global_sections = []
    for sec in secondary:
        global_sections.extend(split_into_sections(sec))
    return global_sections


def find_non_varying_linkers(initial_coords_file, fingerprint_file):
    """Original beta-sheet selection with an auditable, fail-closed trial check.

    The empirical decision criterion is unchanged: a linker is allowed only if
    no evaluated trial reduces the per-strand-pair CA contact count. HDX does not
    participate in this function.

    Previously zero matched trials produced zero broken contacts, so ALL linkers
    passed despite no useful structural evidence. Missing or invalid trial sets
    still trigger an error. Valid unchanged trials count as zero contact loss.
    """
    sheet_coords = sheet_pipe(initial_coords_file, fingerprint_file)
    ref_bonds = sheet_pairwise_bond_number(sheet_coords, thr=5.5)
    ref_coords = read_coords(initial_coords_file)
    total_sheet_contacts = int(ref_bonds.sum())
    n_residues = sum(len(ss) for ss in get_secondary(fingerprint_file))
    if ref_coords.shape != (n_residues, 3):
        raise RuntimeError(
            f'Beta-sheet selector: {ref_coords.shape} reference-coordinate shape '
            f'does not match fingerprint with {n_residues} CA positions.')
    if len(sheet_coords) >= 2 and total_sheet_contacts == 0:
        raise RuntimeError(
            f'Beta-sheet selector: {len(sheet_coords)} strand sections, but zero '
            'inter-strand contacts. Check fingerprint and coordinate alignment.')

    linker_indices, skipped_generator_sections = generate_random_structures(
        initial_coords_file, fingerprint_file, return_failures=True
    )
    linker_bond_arr_dict = random_bond_finder(
        'rand_structures/', fingerprint_file, linker_indices
    )

    # These are precisely the filenames matched by random_bond_finder().
    # Require valid evaluated trial files, but allow unchanged structures:
    # a failed regeneration can leave the original contacts intact.
    # Do not change the scoring rule.
    trial_dir = Path('rand_structures')
    names = [p for p in trial_dir.iterdir() if p.is_file() and not p.name.startswith('.')]
    bond_breaks_dict = {}
    audit = []
    errors = []
    for l0 in linker_indices:
        l = int(l0)
        matches = sorted(p for p in names if len(p.name.split('_')) > 1
                         and p.name.split('_')[1] == str(l))
        matrices = linker_bond_arr_dict.get(l0, [])
        if len(matrices) != len(matches):
            errors.append(f'section {l}: {len(matches)} matching files but '
                          f'{len(matrices)} contact matrices')
        changed = 0
        invalid = []
        for p in matches:
            try:
                xyz = read_coords(str(p))
                if xyz.shape != ref_coords.shape or not np.isfinite(xyz).all():
                    invalid.append(p.name + ' (wrong size or nonfinite)')
                elif np.max(np.abs(xyz - ref_coords)) > 1e-5:
                    changed += 1
            except Exception as exc:
                invalid.append(p.name + f' ({exc})')
        unchanged = len(matches) - changed - len(invalid)
        if not matrices or invalid:
            errors.append(f'section {l}: {len(matches)} matched files, '
                          f'{len(matrices)} evaluated matrices, '
                          f'{changed} genuinely changed structures, '
                          f'{len(invalid)} invalid files')
        if any(np.shape(a) != ref_bonds.shape for a in matrices):
            errors.append(f'section {l}: trial contact matrix shape mismatch')
        # Retain the original decision calculation exactly for valid trials.
        breaks = [int((ref_bonds > a).sum()) for a in matrices
                  if np.shape(a) == ref_bonds.shape]
        score = sum(breaks) / (len(matrices) + 1)
        bond_breaks_dict[l] = score
        audit.append({
            'section_index': l,
            'matching_trial_files': len(matches),
            'evaluated_trial_matrices': len(matrices),
            'changed_trial_files': changed,
            'unchanged_trial_files': unchanged,
            'invalid_trial_files': invalid,
            'trials_breaking_at_least_one_sheet_pair': sum(x > 0 for x in breaks),
            'total_broken_sheet_pair_counts': sum(breaks),
            'original_score': score,
            'selected_by_original_sheet_rule': bool(score < 1e-7),
        })

    report = {
        'reference_structure': str(initial_coords_file),
        'fingerprint': str(fingerprint_file),
        'strand_sections': len(sheet_coords),
        'reference_interstrand_ca_contacts_at_5_5A': total_sheet_contacts,
        'eligible_linker_count': len(linker_indices),
        'attempted_linker_count': len(linker_indices) + len(skipped_generator_sections),
        'generator_excluded_linker_count': len(skipped_generator_sections),
        'generator_excluded_linkers': skipped_generator_sections,
        'beta_sheet_approved_count': sum(x['selected_by_original_sheet_rule'] for x in audit),
        'errors': errors,
        'linkers': audit,
    }
    report_path = Path(fingerprint_file).parent / 'beta_sheet_selection_audit.json'
    report_path.write_text(json.dumps(report, indent=2) + '\n')
    print(f'Beta-sheet selector: {len(sheet_coords)} strand sections, '
          f'{total_sheet_contacts} reference CA contacts, '
          f'{len(linker_indices)} eligible linkers, '
          f'{report["beta_sheet_approved_count"]} passed original contact-loss rule. '
          f'Audit: {report_path}')
    if errors:
        raise RuntimeError(
            'Beta-sheet selection cannot be trusted; see ' + str(report_path)
            + '. First problems: ' + '; '.join(errors[:5])
            + '. No automatic linker selection was written.')

    conds = np.asarray([bond_breaks_dict[int(l)] < 1e-7 for l in linker_indices], dtype=bool)
    allowed_linker = linker_indices[conds]
    secondary = get_secondary(fingerprint_file)
    global_sections = build_global_sections(secondary)
    global_linker_indices = [i for i, sec in enumerate(global_sections) if sec[0] == '-']
    allowed_linker = np.array([i for i in allowed_linker if i in global_linker_indices], dtype=int)
    global_linker_indices = np.array([i for i in global_linker_indices if i != 0], dtype=int)
    return allowed_linker, global_linker_indices


def _structure_text_handle(path):
    """Open PDB/mmCIF text, allowing plain files and .gz files."""
    path = str(path)
    if path.lower().endswith(".gz"):
        import gzip
        return gzip.open(path, "rt", errors="replace")
    return open(path, "rt", errors="replace")


def _is_missing_cif_value(value):
    return value is None or value in ("", "?", ".")


def _as_int_resseq(value):
    """Convert common PDB/mmCIF residue sequence fields to int when possible."""
    if _is_missing_cif_value(value):
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        # Some auth_seq_id fields can arrive as "23.0" or occasionally "23A".
        try:
            return int(float(value))
        except (TypeError, ValueError):
            m = re.match(r"^\s*(-?\d+)", str(value))
            if m:
                return int(m.group(1))
    return None


def _cys_name(resname):
    """Treat common oxidised/protonation CYS names as cysteine-like."""
    return str(resname).strip().upper() in {"CYS", "CYX", "CYM"}


def _atom_choice_score(alt_id, occupancy):
    """
    Score alternate atom records.  Prefer blank/default altlocs, then A, then
    other altlocs; within that, prefer higher occupancy.
    """
    alt = "" if alt_id is None else str(alt_id).strip()
    if alt in ("", ".", "?"):
        alt_rank = 3
    elif alt.upper() == "A":
        alt_rank = 2
    else:
        alt_rank = 1
    try:
        occ = float(occupancy)
    except (TypeError, ValueError):
        occ = 0.0
    return (alt_rank, occ)


def _store_best_atom(atom_dict, key, xyz, alt_id=None, occupancy=None):
    score = _atom_choice_score(alt_id, occupancy)
    old = atom_dict.get(key)
    if old is None or score > old[1]:
        atom_dict[key] = (xyz, score)


def _read_cys_atom_coords_from_pdb_records(structure_path, atom_name):
    """
    Read CYS-like atom coordinates from fixed-width PDB/PDB-like ATOM records.

    Returns
    -------
    dict
        {(chain_id, resseq): np.array([x, y, z])}
    """
    atom_name = str(atom_name).strip().upper()
    atoms = {}
    current_model = None

    with _structure_text_handle(structure_path) as fh:
        for line in fh:
            rec = line[0:6].strip()
            if rec == "MODEL":
                try:
                    current_model = int(line[10:14])
                except ValueError:
                    current_model = 1
                continue
            if rec == "ENDMDL" and current_model == 1:
                break
            if current_model not in (None, 1):
                continue
            if rec not in ("ATOM", "HETATM"):
                continue
            if line[12:16].strip().upper() != atom_name:
                continue
            if not _cys_name(line[17:20]):
                continue

            chain = line[21:22].strip() or " "
            resseq = _as_int_resseq(line[22:26])
            if resseq is None:
                continue
            try:
                xyz = np.array([
                    float(line[30:38]),
                    float(line[38:46]),
                    float(line[46:54]),
                ])
            except ValueError:
                continue

            alt_id = line[16:17].strip()
            occupancy = line[54:60].strip()
            _store_best_atom(atoms, (chain, resseq), xyz, alt_id, occupancy)

    return {key: value[0] for key, value in atoms.items()}


def _split_cif_tokens(line):
    """Small mmCIF token splitter suitable for atom_site rows."""
    try:
        return shlex.split(line, comments=False, posix=True)
    except ValueError:
        # Fall back to whitespace splitting for slightly malformed atom_site rows.
        return line.split()


def _read_cys_atom_coords_from_cif(structure_path, atom_name):
    """
    Read CYS-like atom coordinates from an mmCIF _atom_site loop.

    The returned key uses author identifiers where available:
        (auth_asym_id, int(auth_seq_id))
    falling back to label_asym_id / label_seq_id.  This keeps the key space as
    close as possible to the PDB-style chain/residue labels used elsewhere.
    """
    atom_name = str(atom_name).strip().upper()
    atoms = {}

    with _structure_text_handle(structure_path) as fh:
        lines = list(fh)

    i = 0
    n = len(lines)
    while i < n:
        line = lines[i].strip()
        if line.lower() != "loop_":
            i += 1
            continue

        i += 1
        columns = []
        while i < n and lines[i].strip().startswith("_"):
            col = lines[i].strip().split()[0]
            columns.append(col)
            i += 1

        if not columns or not columns[0].lower().startswith("_atom_site."):
            # Skip non-atom_site loop data.
            while i < n:
                s = lines[i].strip()
                if s.lower() == "loop_" or s.startswith("_"):
                    break
                i += 1
            continue

        col_index = {c.lower(): k for k, c in enumerate(columns)}
        row_tokens = []

        def get(row, *names):
            for name in names:
                idx = col_index.get(name.lower())
                if idx is not None and idx < len(row):
                    value = row[idx]
                    if not _is_missing_cif_value(value):
                        return value
            return None

        while i < n:
            s = lines[i].strip()
            if not s or s == "#":
                i += 1
                if row_tokens:
                    row_tokens = []
                continue
            if s.lower() == "loop_" or s.startswith("_"):
                break
            if s.startswith(";"):
                # Multiline values should not occur in atom_site coordinate rows.
                # Skip the block safely rather than trying to interpret it as atoms.
                i += 1
                while i < n and not lines[i].startswith(";"):
                    i += 1
                if i < n:
                    i += 1
                continue

            row_tokens.extend(_split_cif_tokens(s))
            while len(row_tokens) >= len(columns):
                row = row_tokens[:len(columns)]
                row_tokens = row_tokens[len(columns):]

                group = get(row, "_atom_site.group_pdb")
                if group is not None and str(group).upper() not in {"ATOM", "HETATM"}:
                    continue

                model = get(row, "_atom_site.pdbx_pdb_model_num")
                if model is not None and str(model) not in {"1", "1.0"}:
                    continue

                this_atom = get(row, "_atom_site.auth_atom_id", "_atom_site.label_atom_id")
                if this_atom is None or str(this_atom).strip().upper() != atom_name:
                    continue

                resname = get(row, "_atom_site.auth_comp_id", "_atom_site.label_comp_id")
                if not _cys_name(resname):
                    continue

                chain = get(row, "_atom_site.auth_asym_id", "_atom_site.label_asym_id")
                if _is_missing_cif_value(chain):
                    chain = " "

                resseq = _as_int_resseq(
                    get(row, "_atom_site.auth_seq_id", "_atom_site.label_seq_id")
                )
                if resseq is None:
                    continue

                try:
                    xyz = np.array([
                        float(get(row, "_atom_site.cartn_x")),
                        float(get(row, "_atom_site.cartn_y")),
                        float(get(row, "_atom_site.cartn_z")),
                    ])
                except (TypeError, ValueError):
                    continue

                alt_id = get(row, "_atom_site.label_alt_id", "_atom_site.pdbx_pdb_ins_code")
                occupancy = get(row, "_atom_site.occupancy")
                _store_best_atom(atoms, (str(chain), resseq), xyz, alt_id, occupancy)

            i += 1

    return {key: value[0] for key, value in atoms.items()}


def _read_cys_atom_coords(structure_path, atom_name):
    """
    Format-safe CYS atom coordinate reader for PDB and mmCIF input.

    The old disulfide helper parsed raw fixed-width PDB records directly, which
    fails on AlphaFold / PDBx mmCIF files.  This dispatcher keeps the PDB path
    lightweight and adds a minimal mmCIF _atom_site parser for .cif/.mmcif.
    """
    suffixes = [s.lower() for s in Path(str(structure_path)).suffixes]
    if ".cif" in suffixes or ".mmcif" in suffixes:
        return _read_cys_atom_coords_from_cif(structure_path, atom_name)
    return _read_cys_atom_coords_from_pdb_records(structure_path, atom_name)


def find_disulfide_bonds_from_pdb(pdb_path, thr_min=1.8, thr_max=2.5):
    """
    Find Cys SG-SG pairs at real disulfide-bond distance.

    Despite the historical function name, ``pdb_path`` may now be either PDB
    or mmCIF/PDBx CIF.  PDB files are read from fixed-width ATOM/HETATM records;
    CIF files are read from the _atom_site loop.  The function returns
    ``[((chain, resSeq), (chain, resSeq)), ...]`` using author chain/residue
    identifiers where available.

    The search is intentionally performed on the raw input structure before
    Carbonara's own sanitizing/renumbering/chain-splitting.  Callers map these
    labels onto Carbonara's post-split coordinate arrays via the matching logic
    in ``disulfide_safe_linkers``.
    """
    sg = _read_cys_atom_coords(pdb_path, "SG")

    keys = list(sg.keys())
    bonds = []
    for i in range(len(keys)):
        for j in range(i + 1, len(keys)):
            d = np.linalg.norm(sg[keys[i]] - sg[keys[j]])
            if thr_min <= d <= thr_max:
                bonds.append((keys[i], keys[j]))
    return bonds


def _read_cys_ca_from_pdb(pdb_path):
    """
    Cys CA coordinates keyed by (chain, resSeq), in the same key space as
    find_disulfide_bonds_from_pdb.

    Despite the historical function name, ``pdb_path`` may be either PDB or
    mmCIF/PDBx CIF.  These CA coordinates are used to relocate a bonded Cys
    inside Carbonara's post-split coordinate arrays, since SG atoms do not
    survive into the CA-only coarse-grained representation.
    """
    return _read_cys_atom_coords(pdb_path, "CA")


def disulfide_safe_linkers(varying_indices, pdb_path, coords_file, fingerprint_file,
                            rand_dir='rand_structures', thr=1.5):
    """
    Friend to find_non_varying_linkers! 

    returns the input varying_indices if there are no disulfide bonds, or if pdb_path/rand_dir can't be read (maybe we should give a warning?)

    
    
    """
    bonds = find_disulfide_bonds_from_pdb(pdb_path)
    if not bonds:
        return list(varying_indices)

    secondary = get_secondary(fingerprint_file)
    chain_lengths = [len(ss) for ss in secondary]

    try:
        ref_coords_full = read_coords_from_file(coords_file)
    except OSError:
        return list(varying_indices)
    ref_frag_coords, offset = [], 0
    for n in chain_lengths:
        ref_frag_coords.append(ref_coords_full[offset:offset + n])
        offset += n

    def find_pos(xyz):
        for chain_i, arr in enumerate(ref_frag_coords):
            if len(arr) == 0:
                continue
            d2 = np.sum((arr - xyz) ** 2, axis=1)
            idx = int(np.argmin(d2))
            if d2[idx] < 0.01:  # same atom, not just nearby -- sanitizing never moves coordinates
                return chain_i, idx
        return None

    ca = _read_cys_ca_from_pdb(pdb_path)
    flat_offsets = np.cumsum([0] + chain_lengths)
    pairs = []  # (flat_index_a, flat_index_b, reference_distance)
    for key1, key2 in bonds:
        if key1 not in ca or key2 not in ca:
            continue
        loc1, loc2 = find_pos(ca[key1]), find_pos(ca[key2])
        if loc1 is None or loc2 is None:
            continue
        flat1 = flat_offsets[loc1[0]] + loc1[1]
        flat2 = flat_offsets[loc2[0]] + loc2[1]
        ref_dist = np.linalg.norm(ref_coords_full[flat1] - ref_coords_full[flat2])
        pairs.append((flat1, flat2, ref_dist))
    if not pairs:
        return list(varying_indices)

    try:
        rand_files = os.listdir(rand_dir)
    except OSError:
        return list(varying_indices)

    safe = []
    for l in varying_indices:
        li = int(l)
        sample_files = [f for f in rand_files if len(f.split('_')) > 1 and f.split('_')[1] == str(li)]
        stretched = False
        for fname in sample_files:
            try:
                sample_coords = read_coords_from_file(os.path.join(rand_dir, fname))
            except OSError:
                continue
            if len(sample_coords) != len(ref_coords_full):
                continue  # a repeat that failed the CA-CA check writes nothing; a length
                          # mismatch here means something else is off -- skip, don't guess
            for flat1, flat2, ref_dist in pairs:
                d = np.linalg.norm(sample_coords[flat1] - sample_coords[flat2])
                if d > ref_dist + thr:
                    stretched = True
                    break
            if stretched:
                break
        if not stretched:
            safe.append(l)
    return safe


def auto_select_varying_linker(coords_file, fingerprint_file, pdb_path=None):
    """
    Select varying linkers (longer coil regions that can vary without breaking sheets).

    pdb_path: optional. If given, linkers are also filtered through
    disulfide_safe_linkers -- any linker whose reshaping stretches a real
    disulfide bond (detected from the PDB via find_disulfide_bonds_from_pdb)
    beyond a small tolerance is dropped from the varying set. Omit (default
    None) to get the old, disulfide-unaware behaviour unchanged.
    """
    allowed_linker, linker_indices = find_non_varying_linkers(coords_file, fingerprint_file)

    secondary = get_secondary(fingerprint_file)
    sections = build_global_sections(secondary)

    varying_linker_indices = []
    for section_index in allowed_linker:
        sec = sections[section_index]
        if sec[0] == '-' and len(sec) > 3:
            varying_linker_indices.append(section_index)

    if pdb_path is not None:
        varying_linker_indices = disulfide_safe_linkers(
            varying_linker_indices, pdb_path, coords_file, fingerprint_file
        )

    return varying_linker_indices

def choose_sections_by_number(fullPoss, selected_ids):
    """
    fullPoss: list like ['ResID: 1-6', 'ResID: 14-17', ...]
    selected_ids: 1-based list like [1, 5, 7, 26, 33]

    Returns
    -------
    new_selected_secs : list[str]
    """
    n = len(fullPoss)
    new_selected_secs = []

    for i in selected_ids:
        if not isinstance(i, int):
            print(f"Warning: {i!r} is not an integer, skipping")
            continue

        if 1 <= i <= n:
            new_selected_secs.append(fullPoss[i - 1])
        else:
            print(f"Warning: index {i} out of range 1-{n}, skipping")

    # remove duplicates, preserve order
    seen = set()
    new_selected_secs = [x for x in new_selected_secs if not (x in seen or seen.add(x))]

    return new_selected_secs

def choose_sections_by_number_index(fullPoss, selected_ids):
    return [fullPoss[ss-1] for ss in selected_ids]


def write_varysections_file_frontend(selected_secs,working_path):
    ss_len_tensor = [len(i) for i in get_sses(working_path+"/fingerPrint1.dat")]
    # Calculate cumulative lengths
    cumulative_lengths = np.cumsum(ss_len_tensor)

    # Initialize an empty list to hold the chunked groups
    chunked_groups = [[] for _ in range(len(cumulative_lengths))]
    
    for index in selected_secs:
        for i, cum_length in enumerate(cumulative_lengths):
            if index < cum_length:
                chunked_groups[i].append(index)
                break
    flattened = [item for sublist in chunked_groups for item in sublist]
    write_varysections_file(flattened, working_path)


def view_selected_sections_sequence(run_name):
    filename = "carbonara_runs/"+run_name+"/fingerPrint1.dat"
    highlight_segments = np.loadtxt("carbonara_runs/"+run_name+"/varyingSectionSecondary1.dat",dtype=int)
    chains = parse_structures_with_segments(filename)
    print_structure_with_highlights(chains,highlight_segments)

def view_newly_selected_sections_sequence(run_name,highlight_segments):
    filename = "carbonara_runs/"+run_name+"/fingerPrint1.dat"
    chains = parse_structures_with_segments(filename)
    print_structure_with_highlights(chains,highlight_segments)

#The segments highlighted in red below are those selected for changing

def possibleLinkerList_all_chains(fp_fl, pdb_fl):
    """
    Reads a secondary structure file (fp_fl) and structure file (pdb_fl),
    and returns dash-only ('-') segments for all chains.

    Returns
    -------
    numpy.ndarray
        Array of shape (n_linkers, 2) with rows:
            [global_segment_number, 'Chain X ResID: start-end']

    Here global_segment_number counts all segments ('-', 'S', 'H')
    across all chains in order, without restarting at each chain.
    """

    resid_tensor = getResIDs_from_structure(pdb_fl, fp_fl)

    with open(fp_fl, "r") as f:
        lines = [line.strip() for line in f if line.strip()]

    n_chains = int(lines[0])
    structures = lines[2::2]

    if len(structures) != n_chains:
        raise ValueError(
            f"Fingerprint file says {n_chains} chains, but found {len(structures)} structure entries."
        )

    dash_segments = []
    global_segment_number = 0

    for chain in range(1, n_chains + 1):
        resids = resid_tensor[chain - 1]
        structure = structures[chain - 1]

        idx = 0

        for match in re.finditer(r"(-+|S+|H+)", structure):
            length = match.end() - match.start()

            if match.group()[0] == "-":
                res_start = resids[idx]
                res_end = resids[idx + length - 1]
                label = f"Chain {chain} ResID: {res_start}-{res_end}"
                dash_segments.append([global_segment_number, label])

            idx += length
            global_segment_number += 1

    return np.array(dash_segments, dtype=object)


def parse_linker_label(label):
    m = re.search(r"Chain\s+(\d+)\s+ResID:\s*(\d+)-(\d+)", label)
    if m is None:
        raise ValueError(f"Could not parse linker label: {label}")
    return int(m.group(1)), int(m.group(2)), int(m.group(3))


def ranges_overlap(a1, a2, b1, b2):
    return not (a2 < b1 or b2 < a1)

def linker_ids_from_ranges_new(fullPoss, user_ranges):
    """
    Returns a numpy array with same structure as fullPoss:
        [segment_id, label]
    """

    selected_rows = []

    for row in fullPoss:
        seg_id = int(row[0])
        label = row[1]

        chain, lstart, lend = parse_linker_label(label)

        if chain not in user_ranges:
            continue

        for rstart, rend in user_ranges[chain]:
            if ranges_overlap(lstart, lend, rstart, rend):
                selected_rows.append([seg_id, label])
                break

    # remove duplicates but preserve order
    seen = set()
    final_rows = []
    for seg_id, label in selected_rows:
        key = (seg_id, label)
        if key not in seen:
            seen.add(key)
            final_rows.append([seg_id, label])

    return np.array(final_rows, dtype=object)


def convert_user_ranges_to_structure(user_ranges, resid_tensor, warn=True):
    """
    Convert chain-local residue ranges to structure residue IDs.

    Adds warnings if user ranges fall outside the chain bounds.
    """
    converted = {}

    for chain, ranges in user_ranges.items():
        resids = resid_tensor[chain - 1]

        min_res = 1
        max_res = len(resids)

        offset = int(resids[0]) - 1

        new_ranges = []

        for start, end in ranges:
            # --- sanity checks ---
            if warn:
                if start < min_res or end > max_res:
                    print(
                        f"Warning: Chain {chain} range {start}-{end} "
                        f"is outside valid range {min_res}-{max_res}"
                    )
                elif start > end:
                    print(
                        f"Warning: Chain {chain} range {start}-{end} "
                        f"has start > end"
                    )

            # still convert (do not silently drop)
            new_ranges.append((start + offset, end + offset))

        converted[chain] = new_ranges

    return converted
    
def linker_ids_from_ranges(run_name,pdb_name,user_ranges):
    fullPoss = possibleLinkerList_all_chains("carbonara_runs/" + run_name + "/fingerPrint1.dat",pdb_name)
    resid_tensor = getResIDs_from_structure(pdb_name, "carbonara_runs/" + run_name + "/fingerPrint1.dat")
    user_ranges_mod =convert_user_ranges_to_structure(user_ranges, resid_tensor)
    selected_secs = linker_ids_from_ranges_new(fullPoss, user_ranges_mod)
    return selected_secs[:,0].tolist(),selected_secs[:,1].tolist()

def initial_linker_ids(run_name,pdb_name):
    allowed_linker = np.loadtxt("carbonara_runs/"+run_name+"/varyingSectionSecondary1.dat",dtype=int).tolist()
    fullPoss = possibleLinkerList_all_chains("carbonara_runs/"+run_name+"/fingerPrint1.dat",pdb_name)
    selected_secs = [i[1] for i in fullPoss if int(i[0]) in allowed_linker and linkerLengthCheck(i[1])]
    fullPossIndicies = fullPoss[:,0].tolist()
    fullPoss = fullPoss[:,1].tolist()
    return selected_secs,fullPossIndicies,fullPoss


def _normalise_cmd(cmd):
    if isinstance(cmd, str):
        return shlex.split(cmd)
    if isinstance(cmd, (list, tuple)):
        return list(cmd)
    raise TypeError("foxs_cmd must be a string or a list/tuple")


def load_numeric_table_loose(path, min_cols=2):
    """
    Load only numeric rows from a whitespace-delimited file.
    Skips headers and other non-numeric lines.
    """
    rows = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            try:
                vals = [float(x) for x in parts]
            except ValueError:
                continue
            if len(vals) >= min_cols:
                rows.append(vals)

    if not rows:
        raise ValueError(f"No numeric rows found in {path}")

    ncols = min(len(r) for r in rows)
    return np.array([r[:ncols] for r in rows], dtype=float)


def _find_foxs_fit_file(pdb_path, saxs_path):
    """
    Try to locate a likely FoXS fit/output file.
    Searches near the input files and current working directory.
    """
    search_dirs = []
    for d in [Path.cwd(), pdb_path.parent, saxs_path.parent]:
        if d not in search_dirs:
            search_dirs.append(d)

    pdb_stem = pdb_path.stem.lower()
    saxs_stem = saxs_path.stem.lower()

    candidates = []

    for d in search_dirs:
        for f in d.iterdir():
            if not f.is_file():
                continue

            name = f.name.lower()
            score = 0

            if "fit" in name:
                score += 5
            if "foxs" in name:
                score += 3
            if "profile" in name:
                score += 2
            if pdb_stem in name:
                score += 2
            if saxs_stem in name:
                score += 2
            if f.suffix.lower() in {".fit", ".dat", ".txt"}:
                score += 1

            if score > 0:
                candidates.append((score, f.stat().st_mtime, f))

    if not candidates:
        return None

    candidates.sort(key=lambda x: (x[0], x[1]), reverse=True)
    return candidates[0][2]


def _resolve_saxs_for_initial_foxs_check(saxs_name, prepared_saxs_name=None, run_dir=None):
    """
    Choose the SAXS file used by the initial FoXS sanity check.

    Prefer the Carbonara-prepared/trimmed Saxs.dat when the caller supplies
    either ``run_dir`` or ``prepared_saxs_name``.  Otherwise keep the legacy
    behaviour and use ``saxs_name`` exactly as supplied.
    """
    original_saxs = Path(saxs_name).resolve()

    if run_dir is not None:
        candidate = Path(run_dir).resolve() / "Saxs.dat"
        if not candidate.exists():
            raise FileNotFoundError(
                f"run_dir was supplied for the initial FoXS check, but no Saxs.dat was found: {candidate}"
            )
        return candidate, original_saxs, "run_dir/Saxs.dat"

    if prepared_saxs_name is not None:
        candidate = Path(prepared_saxs_name).resolve()
        if candidate.is_dir():
            candidate = candidate / "Saxs.dat"
        if not candidate.exists():
            raise FileNotFoundError(
                f"Prepared SAXS file for the initial FoXS check was not found: {candidate}"
            )
        return candidate, original_saxs, "prepared_saxs_name"

    if original_saxs.is_dir():
        candidate = original_saxs / "Saxs.dat"
        if not candidate.exists():
            raise FileNotFoundError(
                f"SAXS argument is a directory, but it does not contain Saxs.dat: {candidate}"
            )
        return candidate, original_saxs, "saxs_name_directory/Saxs.dat"

    if not original_saxs.exists():
        raise FileNotFoundError(f"SAXS file not found: {original_saxs}")

    return original_saxs, original_saxs, "saxs_name"


def _write_protein_only_pdb_for_foxs(pdb_path):
    """
    Write a temporary PDB containing only the protein ATOM records used by
    Carbonara-style backmapped structures.

    This deliberately drops HETATM records such as waters, ligands and metal
    ions.  pyFoXS can fail on some metal element labels, for example cobalt
    written as ``CO`` rather than ``Co``.  Carbonara's generated all-atom
    structures are protein-only, so the initial FoXS check should use the same
    convention.
    """
    pdb_path = Path(pdb_path)
    tmp = NamedTemporaryFile(
        suffix="_protein_only_for_foxs.pdb",
        prefix=pdb_path.stem + "_",
        mode="w",
        delete=False,
    )

    kept_atoms = 0
    dropped_hetatm = 0
    dropped_other_atomlike = 0
    chains = set()
    wrote_end = False
    in_first_model = True
    saw_model = False
    saw_any_atom = False

    try:
        with open(pdb_path, "r") as fin, tmp:
            for line in fin:
                rec = line[:6]

                if rec.startswith("MODEL"):
                    if saw_model and saw_any_atom:
                        break
                    saw_model = True
                    in_first_model = True
                    continue

                if rec.startswith("ENDMDL"):
                    if saw_any_atom:
                        break
                    in_first_model = False
                    continue

                if not in_first_model:
                    continue

                if rec == "ATOM  ":
                    tmp.write(line)
                    kept_atoms += 1
                    saw_any_atom = True
                    if len(line) > 21 and line[21].strip():
                        chains.add(line[21])
                elif rec == "HETATM":
                    dropped_hetatm += 1
                elif rec in {"ANISOU", "SIGATM", "SIGUIJ"}:
                    dropped_other_atomlike += 1
                elif rec.startswith("TER") and saw_any_atom:
                    tmp.write(line)
                elif rec.startswith("END"):
                    wrote_end = True
                    # Write a clean END below, after filtering.
                    break

            if kept_atoms == 0:
                raise ValueError(
                    f"No ATOM records were found after making a protein-only FoXS input from {pdb_path}."
                )

            if not wrote_end:
                tmp.write("END\n")
            else:
                tmp.write("END\n")

    except Exception:
        try:
            os.remove(tmp.name)
        except OSError:
            pass
        raise

    return Path(tmp.name), {
        "kept_atom_records": kept_atoms,
        "dropped_hetatm_records": dropped_hetatm,
        "dropped_other_atomlike_records": dropped_other_atomlike,
        "chains": sorted(chains),
    }

def run_initial_foxs_check(pdb_name, saxs_name, foxs_cmd="pyfoxs", max_q=None, protein_only=True, prepared_saxs_name=None, run_dir=None, keep_temp_inputs=False):
    """
    Run pyFoXS, extract chi^2, locate the fit file, and plot
    the fit with a FoXS-style residual panel.

    Parameters
    ----------
    pdb_name : str
        Path to structure file.  By default the FoXS check uses a temporary
        protein-only PDB generated from this file, matching Carbonara's
        backmapped AA outputs.
    saxs_name : str
        Path to SAXS data file.  Prefer passing Carbonara's prepared/trimmed
        Saxs.dat, or pass ``run_dir``/``prepared_saxs_name`` below.
    foxs_cmd : str or list
        Examples:
            "pyfoxs"
            "python3 /path/to/foxs.py"
            ["python3", "/path/to/foxs.py"]
    max_q : float or None
        Optional maximum q-value to pass to pyFoXS.
    protein_only : bool
        If True, run FoXS on a temporary PDB containing only ATOM records.
        This drops HETATM waters, ligands and metal ions, matching the
        protein-only structures produced by Carbonara backmapping and avoiding
        pyFoXS parser failures on unusual element labels.
    prepared_saxs_name : str or None
        Optional explicit path to the Carbonara-prepared SAXS file.  If this is
        a directory, ``Saxs.dat`` inside it is used.
    run_dir : str or None
        Optional Carbonara run directory.  If supplied, ``run_dir/Saxs.dat`` is
        used for FoXS, so the already normalised and Guinier-trimmed data are
        used.
    keep_temp_inputs : bool
        Keep temporary FoXS PDB inputs for debugging.

    Returns
    -------
    dict
        Keys:
            chi2, stdout, stderr, fit_file
    """
    pdb_path = Path(pdb_name).resolve()
    if not pdb_path.exists():
        raise FileNotFoundError(f"Structure file not found: {pdb_path}")

    saxs_path, original_saxs_path, saxs_source = _resolve_saxs_for_initial_foxs_check(
        saxs_name, prepared_saxs_name=prepared_saxs_name, run_dir=run_dir
    )

    temp_pdbs_to_clean = []
    protein_only_info = None

    try:
        pdb_for_foxs, temp_pdb_to_clean = _convert_cif_to_pdb_for_foxs(pdb_path)
        pdb_for_foxs = Path(pdb_for_foxs)
        if temp_pdb_to_clean is not None:
            temp_pdbs_to_clean.append(Path(temp_pdb_to_clean))

        if protein_only:
            protein_pdb, protein_only_info = _write_protein_only_pdb_for_foxs(pdb_for_foxs)
            temp_pdbs_to_clean.append(protein_pdb)
            pdb_for_foxs = protein_pdb
            print(
                "Initial FoXS check: using protein-only temporary PDB "
                f"({protein_only_info['kept_atom_records']} ATOM records; "
                f"dropped {protein_only_info['dropped_hetatm_records']} HETATM records; "
                f"chains={','.join(protein_only_info['chains']) or 'unknown'})."
            )

        if saxs_source != "saxs_name":
            print(f"Initial FoXS check: using Carbonara-prepared SAXS data from {saxs_path}")
        else:
            print(f"Initial FoXS check: using SAXS data from {saxs_path}")

        base_cmd = _normalise_cmd(foxs_cmd)
        cmd = base_cmd + [str(pdb_for_foxs), str(saxs_path)]

        if max_q is not None:
            cmd += ["--max_q", str(max_q)]

        print(cmd)
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

        fit_file = _find_foxs_fit_file(pdb_for_foxs, saxs_path)

        if chi2 is not None:
            print(f"Initial FoXS chi^2: {chi2:.4g}")
        else:
            print("Initial FoXS chi^2: not parsed from output")

        if fit_file is None:
            print("Could not identify a fit file automatically.")
            return {
                "chi2": chi2,
                "stdout": stdout,
                "stderr": stderr,
                "fit_file": None,
                "pdb_for_foxs": str(pdb_for_foxs),
                "saxs_for_foxs": str(saxs_path),
                "original_saxs": str(original_saxs_path),
                "saxs_source": saxs_source,
                "protein_only": protein_only,
                "protein_only_info": protein_only_info,
            }

        #print(f"Using fit file: {fit_file}")

        try:
            fit = load_numeric_table_loose(fit_file, min_cols=2)

            q = fit[:, 0]

            if fit.shape[1] >= 4:
                i_exp = fit[:, 1]
                sigma = fit[:, 2]
                i_fit = fit[:, 3]
                residual = (i_exp - i_fit) / sigma
                residual_label = r"$(I_{\rm exp}-I_{\rm fit})/\sigma$"
            elif fit.shape[1] == 3:
                i_exp = fit[:, 1]
                sigma = None
                i_fit = fit[:, 2]
                residual = i_exp - i_fit
                residual_label = r"$I_{\rm exp}-I_{\rm fit}$"
            else:
                raise ValueError(
                    "Fit file must have at least 3 columns for residual plotting."
                )

            fig = plt.figure(figsize=(7, 7))
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
            # ax2.set_ylabel(residual_label)

            plt.tight_layout()
            plt.show()

        except Exception as e:
            print(f"FoXS ran, but plotting failed: {e}")

        return {
            "chi2": chi2,
            "stdout": stdout,
            "stderr": stderr,
            "fit_file": fit_file,
            "pdb_for_foxs": str(pdb_for_foxs),
            "saxs_for_foxs": str(saxs_path),
            "original_saxs": str(original_saxs_path),
            "saxs_source": saxs_source,
            "protein_only": protein_only,
            "protein_only_info": protein_only_info,
        }

    finally:
        if not keep_temp_inputs:
            for tmp_pdb in temp_pdbs_to_clean:
                try:
                    os.remove(tmp_pdb)
                except OSError:
                    pass
        elif temp_pdbs_to_clean:
            print("Initial FoXS check: kept temporary input(s): " + ", ".join(str(p) for p in temp_pdbs_to_clean))


from Bio.PDB import MMCIFParser, PDBIO
import tempfile

def _convert_cif_to_pdb_for_foxs(structure_path):
    structure_path = Path(structure_path)
    suffix = structure_path.suffix.lower()

    if suffix not in [".cif", ".mmcif"]:
        return str(structure_path), None

    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("model", str(structure_path))

    tmp = tempfile.NamedTemporaryFile(suffix=".pdb", delete=False)
    tmp.close()

    io = PDBIO()
    io.set_structure(structure)
    io.save(tmp.name)

    return tmp.name, tmp.name

# -----------------------------------------------------------------------------
# Carbonara PDB/mmCIF robustness patch
# Paste this at the END of CarbonaraDataTools.py so it overrides the earlier
# definitions. It can then be removed cleanly if needed.
#
# What this overrides/adds:
#   - sanitize_pdb_for_carbonara(...)
#   - pdb_2_biobox(...)
#   - find_missing_residues(...)
#   - pull_structure_from_pdb(...)
#
# Defaults are chosen to preserve Carbonara's legacy global input numbering
# convention while fixing duplicate/alternate atom records and odd author
# residue numbering.
# -----------------------------------------------------------------------------

from collections import OrderedDict
from tempfile import NamedTemporaryFile
import os
import numpy as np
import mdtraj as md
import biobox as bb


_CARBONARA_AA3 = {
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
    'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
    'SEC', 'PYL',
    'HIP', 'HID', 'HSD', 'HSE', 'HIE',
    'GLH', 'GLUP',
    'CYX', 'CYM',
    'ASH', 'ASPP',
    'LYN', 'LSN',
    'MSE',
}


def _carbonara_normalise_resname(resname):
    """Map common residue variants onto names understood by get_residue_map()."""
    resname = str(resname).strip()
    if resname in ('HIP', 'HID', 'HSD', 'HSE', 'HIE'):
        return 'HIS'
    if resname in ('GLH', 'GLUP'):
        return 'GLU'
    if resname in ('CYX', 'CYM'):
        return 'CYS'
    if resname in ('ASH', 'ASPP'):
        return 'ASP'
    if resname in ('LYN', 'LSN'):
        return 'LYS'
    if resname == 'MSE':
        return 'MET'
    return resname


def _carbonara_safe_occupancy(line):
    """Read PDB occupancy; malformed/missing occupancy is treated as 0."""
    try:
        return float(line[54:60])
    except Exception:
        return 0.0


def _carbonara_replace_pdb_fields(line, serial=None, resseq=None, altloc=' '):
    """Replace selected fixed-width PDB fields while preserving coordinates etc."""
    line = line.rstrip('\n')
    if len(line) < 80:
        line = line.ljust(80)
    chars = list(line)

    if serial is not None:
        chars[6:11] = f"{int(serial):5d}"

    if altloc is not None:
        chars[16] = altloc

    if resseq is not None:
        resseq = int(resseq)
        if not (1 <= resseq <= 9999):
            raise ValueError(
                f"PDB residue number {resseq} is outside the fixed-width PDB range 1..9999. "
                "Use mmCIF output or set renumber_residues=False for this structure."
            )
        chars[22:26] = f"{resseq:4d}"
        chars[26] = ' '

    return ''.join(chars).rstrip() + '\n'


def _carbonara_make_ter_line(serial, last_atom_line):
    """Create a simple TER line matching the previous written atom."""
    if last_atom_line is None:
        return 'TER\n'

    line = last_atom_line.rstrip('\n')
    if len(line) < 80:
        line = line.ljust(80)

    ter = f"TER   {int(serial):5d}      {line[17:20]} {line[21]}{line[22:26]}{line[26]}"
    return ter.ljust(80).rstrip() + '\n'


def sanitize_pdb_for_carbonara(
    pdb_in,
    pdb_out=None,
    keep_hetatm=True,
    renumber_residues=True,
    renumber_mode='global',
    first_model_only=True,
    preserve_ter=True,
):
    """
    Write a Carbonara-safe PDB and return (path, report).

    Fixes:
      1. Alternate-location atoms and duplicate atom records written as if normal.
      2. Odd residue numbering, e.g. chain A starting at 7 or chain B at 0.
      3. Non-protein HETATM records being mistaken for CA atoms.
      4. Blank chain IDs where chains are defined only by TER records.

    renumber_mode:
      'global'    -> one global residue axis across all chains/TER segments
      'per_chain' -> each chain/TER segment starts at residue 1
    """
    if pdb_out is None:
        tmp = NamedTemporaryFile(mode='w', suffix='.carbonara_clean.pdb', delete=False)
        pdb_out = tmp.name
        tmp.close()

    chosen = OrderedDict()
    events = []
    input_atom_records = 0
    seen_model = False
    segment_index = 0

    with open(pdb_in, errors='replace') as handle:
        for original_i, line in enumerate(handle):
            if line.startswith('MODEL'):
                if seen_model and first_model_only:
                    break
                seen_model = True
                continue

            if line.startswith('ENDMDL') and first_model_only:
                break

            if preserve_ter and line.startswith('TER'):
                events.append((original_i, 'TER', segment_index, None))
                segment_index += 1
                continue

            is_atom = line.startswith('ATOM')
            is_hetatm = line.startswith('HETATM')
            if not is_atom and not (keep_hetatm and is_hetatm):
                continue

            resname = line[17:20].strip()
            if resname not in _CARBONARA_AA3:
                continue

            input_atom_records += 1

            atom_key = (
                segment_index,
                line[21],
                line[22:26],
                line[26],
                resname,
                line[12:16],
                line[76:78] if len(line) >= 78 else '',
            )
            events.append((original_i, 'ATOM', segment_index, atom_key))

            occ = _carbonara_safe_occupancy(line)
            if atom_key not in chosen:
                chosen[atom_key] = (original_i, occ, line)
            else:
                old_i, old_occ, old_line = chosen[atom_key]
                if occ > old_occ:
                    chosen[atom_key] = (original_i, occ, line)

    selected_by_original_i = {item[0]: item[2] for item in chosen.values()}

    residue_map = OrderedDict()
    chain_counters = {}
    next_global_resseq = 1
    out_lines = []
    serial = 1
    atoms_since_ter = False
    last_atom_line = None
    output_atom_records = 0

    for original_i, event_type, segment_index, atom_key in events:
        if event_type == 'TER':
            if preserve_ter and atoms_since_ter:
                out_lines.append(_carbonara_make_ter_line(serial, last_atom_line))
                serial += 1
                atoms_since_ter = False
                last_atom_line = None
            continue

        line = selected_by_original_i.get(original_i)
        if line is None:
            continue

        new_resseq = None
        if renumber_residues:
            residue_key = (segment_index, line[21], line[22:26], line[26], line[17:20].strip())
            if residue_key not in residue_map:
                if renumber_mode == 'global':
                    residue_map[residue_key] = next_global_resseq
                    next_global_resseq += 1
                elif renumber_mode == 'per_chain':
                    chain_key = (segment_index, line[21])
                    chain_counters[chain_key] = chain_counters.get(chain_key, 0) + 1
                    residue_map[residue_key] = chain_counters[chain_key]
                else:
                    raise ValueError("renumber_mode must be either 'global' or 'per_chain'")
            new_resseq = residue_map[residue_key]

        out_line = _carbonara_replace_pdb_fields(
            line,
            serial=serial,
            resseq=new_resseq,
            altloc=' ',
        )
        out_lines.append(out_line)
        serial += 1
        output_atom_records += 1
        atoms_since_ter = True
        last_atom_line = out_line

    if preserve_ter and atoms_since_ter:
        out_lines.append(_carbonara_make_ter_line(serial, last_atom_line))

    out_lines.append('END\n')
    with open(pdb_out, 'w') as handle:
        handle.writelines(out_lines)

    report = {
        'input_atom_records': input_atom_records,
        'output_atom_records': output_atom_records,
        'removed_duplicate_atom_records': input_atom_records - output_atom_records,
        'output_residues': len(residue_map) if renumber_residues else None,
        'renumber_residues': renumber_residues,
        'renumber_mode': renumber_mode if renumber_residues else None,
        'preserve_ter': preserve_ter,
    }
    return pdb_out, report


def _carbonara_cif_to_pdb(cif_file):
    """
    Convert mmCIF to a temporary PDB for BioBox and the PDB sanitizer.

    This imports PDBFixer/OpenMM lazily so ordinary PDB loading is unaffected if
    those optional dependencies are unavailable. It does not add missing residues,
    missing atoms, or hydrogens; it only rewrites the existing topology/positions.
    """
    try:
        from pdbfixer import PDBFixer
        from openmm.app import PDBFile
    except Exception as exc:
        raise ImportError(
            "Reading .cif files through this Carbonara patch requires pdbfixer and openmm. "
            "Install them or convert the mmCIF to PDB before calling Carbonara."
        ) from exc

    fixer = _carbonara_pdbfixer(cif_file)
    tmp = NamedTemporaryFile(mode='w', suffix='.pdb', delete=False)
    tmp_path = tmp.name
    tmp.close()

    with open(tmp_path, 'w') as handle:
        try:
            PDBFile.writeFile(fixer.topology, fixer.positions, handle, keepIds=True)
        except TypeError:
            PDBFile.writeFile(fixer.topology, fixer.positions, handle)

    return tmp_path


def pdb_2_biobox(
    pdb_file,
    sanitize=True,
    renumber_residues=True,
    renumber_mode='global',
    preserve_ter=True,
):
    """
    Load a PDB/mmCIF into BioBox after applying the same cleaning used by MDTraj.

    This intentionally overrides the earlier pdb_2_biobox definition when pasted
    at the end of CarbonaraDataTools.py. possibleLinkerList therefore sees the
    cleaned/global residue numbering too.
    """
    M = bb.Molecule()
    ext = os.path.splitext(pdb_file)[1].lower()
    temp_paths = []

    try:
        if ext == '.pdb':
            load_path = pdb_file
            if sanitize:
                load_path, _ = sanitize_pdb_for_carbonara(
                    pdb_file,
                    renumber_residues=renumber_residues,
                    renumber_mode=renumber_mode,
                    preserve_ter=preserve_ter,
                )
                temp_paths.append(load_path)
            M.import_pdb(load_path)

        elif ext == '.cif':
            pdb_tmp = _carbonara_cif_to_pdb(pdb_file)
            temp_paths.append(pdb_tmp)
            load_path = pdb_tmp
            if sanitize:
                clean_tmp, _ = sanitize_pdb_for_carbonara(
                    pdb_tmp,
                    renumber_residues=renumber_residues,
                    renumber_mode=renumber_mode,
                    preserve_ter=preserve_ter,
                )
                temp_paths.append(clean_tmp)
                load_path = clean_tmp
            M.import_pdb(load_path)

        else:
            raise ValueError(f"Unsupported structure file extension: {ext}")

        return M

    finally:
        for path in temp_paths:
            try:
                os.remove(path)
            except OSError:
                pass


def find_missing_residues(resIDs, include_leading=True):
    """
    Return missing residue numbers caused by gaps in residue numbering.

    include_leading=True preserves legacy Carbonara behaviour: if a chain starts
    at residue 242, residues 1..241 are reported as missing. This can look odd,
    but it keeps the old behaviour for globally numbered multichain inputs.

    The sanitizer normally renumbers single-chain inputs to 1..N, so excised
    domains like 325..870 no longer falsely report 1..324 as missing.
    """
    resIDs = np.asarray(resIDs, dtype=int)
    if resIDs.size == 0:
        return np.asarray([], dtype=int)

    missing_residues = []

    if include_leading and resIDs[0] > 1:
        missing_residues.extend(range(1, int(resIDs[0])))

    for left, right in zip(resIDs[:-1], resIDs[1:]):
        left = int(left)
        right = int(right)
        if right > left + 1:
            missing_residues.extend(range(left + 1, right))

    return np.asarray(missing_residues, dtype=int)


def pull_structure_from_pdb(
    file_path,
    sanitize=True,
    renumber_residues=True,
    renumber_mode='global',
    preserve_ter=True,
    missing_include_leading=True,
    verbose=False,
):
    """
    Robust replacement for pull_structure_from_pdb().

    Returns the same four objects as before:
      coords_chains, sequence_chains, secondary_structure_chains, missing_residues_chains

    For clean legacy inputs this should leave coordinates, sequence and secondary
    structure unchanged. It only intervenes before MDTraj/BioBox loading to remove
    duplicate/alternate atom records and normalise problematic residue numbering.
    """
    ext = os.path.splitext(file_path)[1].lower()
    temp_paths = []

    try:
        if ext == '.pdb':
            load_path = file_path
            report = None
            if sanitize:
                load_path, report = sanitize_pdb_for_carbonara(
                    file_path,
                    renumber_residues=renumber_residues,
                    renumber_mode=renumber_mode,
                    preserve_ter=preserve_ter,
                )
                temp_paths.append(load_path)

            if verbose and report is not None:
                print(
                    "Carbonara PDB clean-up: "
                    f"input atoms={report['input_atom_records']}, "
                    f"output atoms={report['output_atom_records']}, "
                    f"removed duplicates/alternates={report['removed_duplicate_atom_records']}, "
                    f"renumber={report['renumber_residues']}:{report['renumber_mode']}"
                )

            traj = md.load(load_path)

        elif ext == '.cif':
            pdb_tmp = _carbonara_cif_to_pdb(file_path)
            temp_paths.append(pdb_tmp)
            load_path = pdb_tmp
            report = None

            if sanitize:
                clean_tmp, report = sanitize_pdb_for_carbonara(
                    pdb_tmp,
                    renumber_residues=renumber_residues,
                    renumber_mode=renumber_mode,
                    preserve_ter=preserve_ter,
                )
                temp_paths.append(clean_tmp)
                load_path = clean_tmp

            if verbose and report is not None:
                print(
                    "Carbonara mmCIF clean-up: "
                    f"input atoms={report['input_atom_records']}, "
                    f"output atoms={report['output_atom_records']}, "
                    f"removed duplicates/alternates={report['removed_duplicate_atom_records']}, "
                    f"renumber={report['renumber_residues']}:{report['renumber_mode']}"
                )

            traj = md.load(load_path)

        else:
            raise ValueError(f"Unsupported structure file extension: {ext}")

        topology = traj.topology
        chains = list(topology.chains)
        if len(chains) == 0:
            raise ValueError("No chains found in structure")

        three_to_one = get_residue_map()
        coords_chains = []
        sequence_chains = []
        secondary_structure_chains = []
        missing_residues_chains = []

        ss_pred = md.compute_dssp(traj, simplified=True)[0]
        ss_map = {'H': 'H', 'E': 'S', 'C': '-', 'NA': '-'}
        ss_pred_mapped = np.array([ss_map.get(str(ss), '-') for ss in ss_pred])
        ss_pred_mapped = clean_s_sequences(ss_pred_mapped)

        residue_index = 0
        for chain in chains:
            residues = list(chain.residues)
            chain_ss_all = ss_pred_mapped[residue_index:residue_index + len(residues)]

            ca_atom_indices = []
            resids = []
            seq = []
            chain_ss = []

            for local_i, res in enumerate(residues):
                ca_atom = next((atom for atom in res.atoms if atom.name == 'CA'), None)
                if ca_atom is None:
                    continue

                ca_atom_indices.append(ca_atom.index)
                resids.append(res.resSeq)

                resname = _carbonara_normalise_resname(res.name)
                seq.append(three_to_one.get(resname, 'X'))
                chain_ss.append(chain_ss_all[local_i])

            if ca_atom_indices:
                ca_coords = traj.xyz[0, ca_atom_indices, :] * 10.0

                if len(ca_coords) > 10:
                    if not (len(ca_coords) == len(seq) == len(chain_ss)):
                        raise ValueError(
                            "Internal structure parsing mismatch after clean-up: "
                            f"coords={len(ca_coords)}, sequence={len(seq)}, ss={len(chain_ss)}"
                        )

                    coords_chains.append(ca_coords)
                    sequence_chains.append(np.asarray(seq))
                    secondary_structure_chains.append(np.asarray(chain_ss))
                    missing_residues_chains.append(
                        find_missing_residues(
                            np.asarray(resids),
                            include_leading=missing_include_leading,
                        )
                    )

            residue_index += len(residues)

        return coords_chains, sequence_chains, secondary_structure_chains, missing_residues_chains

    finally:
        for path in temp_paths:
            try:
                os.remove(path)
            except OSError:
                pass


# Update the generated RunMe script to use a specific disulfide constraint file.
def set_disulfide_constraints_for_run(run_name, disulfide_file):
    runme = Path(f"RunMe_{run_name}.sh")
    disulfide_file = Path(disulfide_file).resolve()

    if not runme.exists():
        raise FileNotFoundError(f"Could not find generated run script: {runme}")

    if not disulfide_file.exists():
        raise FileNotFoundError(f"Could not find disulfide constraint file: {disulfide_file}")

    text = runme.read_text()

    replacement = f'DISULFIDE_CONSTRAINTS_FILE="{disulfide_file}"'
    text_new, n = re.subn(
        r'(?m)^DISULFIDE_CONSTRAINTS_FILE=.*$',
        replacement,
        text,
    )

    if n == 0:
        raise RuntimeError(
            "Could not find DISULFIDE_CONSTRAINTS_FILE line in "
            f"{runme}. Regenerate the RunMe script or add the variable manually."
        )

    runme.write_text(text_new)
    print(f"Updated {runme}:")
    print(replacement)


# -----------------------------------------------------------------------------
# HDX-MS-guided flexibility for Carbonara setups.
# HDX is an exchange/protection proxy, not a direct measure of backbone mobility.
# Positions may refer to the concatenated fingerprint chains, one chain, or repeated subunits.
# Apply any helical override BEFORE writing/numbering the fingerprint.
# -----------------------------------------------------------------------------
def prepare_hdx_flexibility_from_csv(
    csv_file,
    sequence_chains,
    secondary_structure_chains,
    *,
    chain=None,
    mapping='auto',
    include_medium=False,
    override_helices=False,
    position_offset=0,
    min_helix_run=4,
    fast_min=1.0,
    medium_min=0.1,
    return_measured_masks=False,
):
    """Map sparse per-residue HDX measurements to Carbonara's prepared chains.

    ``position`` is ALWAYS retained: a missing ``k_est`` does not shift any
    subsequent residue. The four mapping schemes are:

    * global: CSV positions refer to the concatenation of fingerprint chains,
      starting at chain 1, residue 1. Suitable for a 900-residue tetramer
      whose *total* length is 900.
    * repeat: the same measured 1..L profile applies separately to every
      chain. Suitable for a homotetramer with L residues PER subunit.
    * chain: positions refer only to ``chain`` (1-based), optionally using
      ``position_offset`` to convert original numbering to chain-local.
    * auto: choose global if the full CSV length equals total chain length;
      choose repeat if every chain matches that length; single-chain cases
      choose chain. Otherwise require an explicit mapping.

    Missing values are unknown, not evidence of stability or flexibility.
    ``aa`` sequence validation is required for helix override and recommended
    otherwise. ``H`` -> ``-`` is performed only on consecutive selected helical
    runs of ``min_helix_run`` or more residues, before fingerprint numbering.
    Strands and coordinates are never touched.

    Returns (secondary_structure_chains, selected_masks, diagnostic_report),
    or those three values plus measured_masks when return_measured_masks=True.
    The optional fourth return preserves the existing API for other callers.
    """
    from pathlib import Path as _Path
    csv_file = _Path(csv_file).expanduser().resolve()
    if not csv_file.is_file():
        raise FileNotFoundError(f'HDX CSV not found: {csv_file}')
    if not (0 < medium_min < fast_min) or int(min_helix_run) < 4:
        raise ValueError('Require 0 < medium_min < fast_min and min_helix_run >= 4')
    if mapping not in ('auto', 'global', 'repeat', 'chain'):
        raise ValueError('HDX mapping must be auto, global, repeat or chain')
    if len(sequence_chains) != len(secondary_structure_chains) or len(sequence_chains) == 0:
        raise ValueError('Sequence/secondary-structure chain counts differ or are empty')
    if chain is not None and (int(chain) < 1 or int(chain) > len(sequence_chains)):
        raise ValueError(f'HDX chain {chain} invalid; fingerprint has {len(sequence_chains)} chains')
    if chain is not None and mapping in ('global', 'repeat'):
        raise ValueError('--hdx-chain cannot be combined with global/repeat mapping')

    seqs = [''.join(str(x) for x in s).upper() for s in sequence_chains]
    ss_chains = [''.join(str(x) for x in s) for s in secondary_structure_chains]
    lengths = [len(s) for s in seqs]
    if any(len(seq) != len(ss) for seq, ss in zip(seqs, ss_chains)):
        raise ValueError('HDX sequence and secondary structure lengths differ')

    df = pd.read_csv(csv_file)
    missing_cols = {'position', 'k_est'} - set(df.columns)
    if missing_cols:
        raise ValueError(f'HDX CSV missing columns: {sorted(missing_cols)}')
    pos = pd.to_numeric(df['position'], errors='coerce')
    if pos.isna().any() or ((pos % 1) != 0).any() or (pos <= 0).any():
        raise ValueError('HDX positions must be positive integer residue indices')
    pos = pos.astype(int)
    if pos.duplicated().any():
        repeated = pos[pos.duplicated()].head(6).tolist()
        raise ValueError(f'Duplicate HDX residue positions: {repeated}')
    rates = pd.to_numeric(df['k_est'], errors='coerce').to_numpy(dtype=float)
    if np.any(np.isinf(rates)) or np.any(rates[np.isfinite(rates)] < 0):
        raise ValueError('HDX k_est contains infinite or negative exchange rates')
    # Proline has no backbone amide proton. In conventional backbone HDX,
    # k_est=0 for proline is NOT evidence of protection or rigidity; retaining
    # it as a measured low-exchange restraint can spuriously preserve an H1
    # island between otherwise fast-exchanging linkers. Exclude it from the
    # HDX evidence masks, without removing rows or changing residue numbering.
    # (All non-proline finite exchange values retain their original meaning.)
    proline_uninformative = np.zeros(len(rates), dtype=bool)
    if 'aa' in df.columns:
        aa_letters = df['aa'].fillna('').astype(str).str.strip().str.upper().to_numpy()
        proline_uninformative = (aa_letters == 'P') & np.isfinite(rates)
    measured = np.isfinite(rates) & ~proline_uninformative
    offset = int(position_offset)
    positions = pos.to_numpy(dtype=int) - offset
    nrows = len(df)
    # Determine whether this is a whole-complex or repeated-subunit profile
    # from its ORIGINAL numbering, before any offset is applied. In particular,
    # an offset of one must not stop a 900-row profile being recognised as a
    # repeated profile for 900-residue subunits.
    csv_length = int(pos.max())
    full_length = sum(lengths)
    if mapping == 'auto':
        if chain is not None or len(lengths) == 1:
            mapping = 'chain'
        elif csv_length == full_length:
            mapping = 'global'
        elif all(L == csv_length for L in lengths):
            mapping = 'repeat'
        else:
            raise ValueError(
                f'Cannot safely map {len(df)} HDX CSV positions (maximum {csv_length}) '
                f'to {len(lengths)} Carbonara chains of lengths {lengths}. '
                'Specify --hdx-mapping global, repeat, or chain with --hdx-chain. '
                'Missing measurements must not be compressed or shifted.'
            )
    if mapping == 'chain' and chain is None:
        if len(lengths) != 1:
            raise ValueError('HDX chain mapping requires --hdx-chain for multichain structures')
        chain = 1

    # Translate a CSV row to zero or more (chain-index-0, residue-index-0)
    # targets.  A 900-row CSV may map to 900 positions across four chains,
    # or (repeat mode) to 900 in each of four identical subunits.
    boundaries = np.cumsum([0] + lengths)
    targets_per_row = []
    in_bounds_measured_rows = 0
    for p, known in zip(positions, measured):
        targets = []
        if mapping == 'global':
            if 1 <= p <= full_length:
                ci = int(np.searchsorted(boundaries, p, side='left') - 1)
                # At a chain boundary p == boundaries[j], this residue
                # belongs to the PREVIOUS chain.
                ci = min(max(ci, 0), len(lengths) - 1)
                local = int(p - boundaries[ci] - 1)
                targets = [(ci, local)]
        elif mapping == 'repeat':
            targets = [(ci, int(p - 1)) for ci, L in enumerate(lengths) if 1 <= p <= L]
        else:
            ci = int(chain) - 1
            if 1 <= p <= lengths[ci]:
                targets = [(ci, int(p - 1))]
        targets_per_row.append(targets)
        if known and targets:
            in_bounds_measured_rows += 1

    # Unmeasured rows are placeholders, not restraints. An offset may place
    # empty terminal placeholders outside the structure, which is harmless.
    # Measured rows must NEVER be silently discarded.
    out_of_bounds_measured = measured & np.asarray([not t for t in targets_per_row])
    if mapping in ('global', 'repeat') and np.any(out_of_bounds_measured):
        bad = pos[out_of_bounds_measured].head(8).tolist()
        raise ValueError(
            f'Measured HDX positions outside {mapping} mapping for prepared '
            f'chains {lengths}, with position_offset={offset}: {bad}. '
            'Check numbering or structure length; no measured HDX rows were dropped.'
        )
    if in_bounds_measured_rows == 0:
        raise ValueError('No measured HDX positions map to selected Carbonara chain(s)')

    comparisons = []
    mismatches = []
    if 'aa' in df.columns:
        for csv_p, value, known, targets in zip(pos, df['aa'], measured, targets_per_row):
            if not known or pd.isna(value):
                continue
            letter = str(value).strip().upper()
            if len(letter) != 1 or not letter.isalpha():
                continue
            for ci, local in targets:
                expected = seqs[ci][local]
                if letter == 'X' or expected == 'X':
                    continue
                comparisons.append((int(csv_p), int(ci + 1), int(local + 1), letter, expected))
                if letter != expected:
                    mismatches.append((int(csv_p), int(ci + 1), int(local + 1), letter, expected))
        if comparisons and len(mismatches) > max(2, int(np.ceil(0.05 * len(comparisons)))):
            example = ', '.join(
                f'CSV {p}->chain {c}:{l} {a}!={b}' for p,c,l,a,b in mismatches[:6]
            )
            # A single missing N-terminal residue or a numbering convention
            # change shifts all comparisons. Suggest an offset only when the
            # FULL observed sequence gives strong evidence for it. We do not
            # silently remap the user's experimental data.
            observed_rows = []
            for raw_p, letter, known in zip(pos.to_numpy(dtype=int), df['aa'], measured):
                if not known or pd.isna(letter):
                    continue
                letter = str(letter).strip().upper()
                if len(letter) == 1 and letter.isalpha() and letter != 'X':
                    observed_rows.append((int(raw_p), letter))

            def offset_fit(candidate_offset):
                matches, compared = 0, 0
                for raw_p, letter in observed_rows:
                    p = raw_p - candidate_offset
                    if mapping == 'global':
                        if not (1 <= p <= full_length):
                            continue
                        ci = int(np.searchsorted(boundaries, p, side='left') - 1)
                        ci = min(max(ci, 0), len(lengths) - 1)
                        pairs = [(ci, int(p - boundaries[ci] - 1))]
                    elif mapping == 'repeat':
                        pairs = [(ci, p - 1) for ci, L in enumerate(lengths) if 1 <= p <= L]
                    else:
                        ci = int(chain) - 1
                        pairs = [(ci, p - 1)] if 1 <= p <= lengths[ci] else []
                    for ci, local in pairs:
                        expected = seqs[ci][local]
                        if expected != 'X':
                            matches += (letter == expected)
                            compared += 1
                return matches, compared

            fit_candidates = []
            for candidate in range(offset - 100, offset + 101):
                hits, total = offset_fit(candidate)
                if total >= max(30, int(0.75 * len(comparisons))):
                    fit_candidates.append((hits / total, total, candidate))
            best = max(fit_candidates, default=None)
            hint = ''
            if best is not None:
                frac, n_checked, best_offset = best
                if best_offset != offset and frac >= 0.95:
                    hint = (f' Best constant-offset alignment: '
                            f'--hdx-position-offset {best_offset} gives '
                            f'{100 * frac:.1f}% amino-acid identity across '
                            f'{n_checked} comparisons. Try that offset; '
                            'the mapping is NOT changed automatically.')
                else:
                    hint = (f' Best tested constant offset was {best_offset} '
                            f'({100 * frac:.1f}% identity); no convincing '
                            'simple offset correction was found. Confirm the '
                            'PDB and HDX refer to the same sequence/construct.')
            raise ValueError(
                f'HDX CSV does not align with prepared sequence: '
                f'{len(mismatches)}/{len(comparisons)} known amino acids differ ({example}). '
                'Check whether the CSV covers the complete tetramer, '
                'one repeated subunit, or a single chain.' + hint
            )
    elif override_helices:
        raise ValueError('Helix override requires HDX CSV aa column for sequence validation')
    if override_helices and len(comparisons) < 10:
        raise ValueError('Helix override requires at least 10 measured amino-acid comparisons')

    # Residue-by-residue mapping, without compressing positions with NaN rates.
    # Do not mark mismatched amino acids as flexible even if their rate is high.
    bad_residues = {(c - 1, l - 1) for _, c, l, _, _ in mismatches}
    selected = [np.zeros(L, dtype=bool) for L in lengths]
    measured_masks = [np.zeros(L, dtype=bool) for L in lengths]
    for targets, k, informative in zip(targets_per_row, rates, measured):
        if not informative:
            continue
        for ci, local in targets:
            measured_masks[ci][local] = True
            if k >= (medium_min if include_medium else fast_min) and (ci, local) not in bad_residues:
                selected[ci][local] = True

    changed = []
    all_converted = []
    all_skipped = []
    per_chain = []
    for ci, (ss, mask, observed) in enumerate(zip(ss_chains, selected, measured_masks)):
        before = np.asarray(list(ss), dtype='<U1')
        after = before.copy()
        converted, skipped = [], []
        if override_helices:
            moving = (before == 'H') & mask
            i = 0
            while i < len(moving):
                if not moving[i]:
                    i += 1
                    continue
                start = i
                while i < len(moving) and moving[i]:
                    i += 1
                end = i
                run = [int(start + 1), int(end)]
                if end - start >= min_helix_run:
                    after[start:end] = '-'
                    converted.append(run)
                else:
                    skipped.append(run)
        changed.append(after)
        all_converted.extend([{'chain': ci + 1, 'start': a, 'end': b} for a,b in converted])
        all_skipped.extend([{'chain': ci + 1, 'start': a, 'end': b} for a,b in skipped])
        per_chain.append({
            'chain': ci + 1,
            'length': len(ss),
            'measured_residues': int(observed.sum()),
            'unmeasured_residues': int((~observed).sum()),
            'selected_residues': int(mask.sum()),
            'helix_residues_relabelled': int(np.sum(before != after)),
            'helix_to_coil_runs_chain_local_inclusive': converted,
        })

    report = {
        'source_csv': str(csv_file),
        'source_rows': nrows,
        'csv_position_min': int(pos.min()),
        'csv_position_max': int(pos.max()),
        'mapping': mapping,
        'prepared_chain_lengths': lengths,
        'target_chain_1based': chain if mapping == 'chain' else None,
        'position_offset': offset,
        'numbering': 'CSV positions are 1-based; subtract position_offset before mapping',
        'mode': 'fast_and_medium' if include_medium else 'fast_only',
        'thresholds': {'fast_k_min': float(fast_min), 'medium_k_min': float(medium_min)},
        'measured_rows': int(measured.sum()),
        'unmeasured_rows': int((~measured).sum()),
        'source_finite_k_est_rows': int(np.isfinite(rates).sum()),
        'nonexchangeable_proline_rate_rows_ignored': int(proline_uninformative.sum()),
        'nonexchangeable_proline_positions_1based_csv': pos[proline_uninformative].astype(int).tolist(),
        'fast_rows_total': int(np.sum(measured & (rates >= fast_min))),
        'medium_rows_total': int(np.sum(measured & (rates >= medium_min) & (rates < fast_min))),
        'slow_rows_total': int(np.sum(measured & (rates < medium_min))),
        'measured_positions_mapped_in_csv': in_bounds_measured_rows,
        'measured_positions_outside_chain': int(np.sum(measured & np.asarray([not t for t in targets_per_row]))),
        'amino_acids_compared': len(comparisons),
        'amino_acid_mismatches': [list(x) for x in mismatches],
        'selected_positions_count': int(sum(int(mask.sum()) for mask in selected)),
        'selected_chain_local_positions_1based': [
            {'chain': ci + 1, 'residues': (np.flatnonzero(mask) + 1).tolist()}
            for ci, mask in enumerate(selected)
        ],
        'override_helices': bool(override_helices),
        'helix_override_min_run': int(min_helix_run),
        'helix_to_coil_runs_chain_local_inclusive': all_converted,
        'short_helix_runs_not_overridden': all_skipped,
        'helix_residues_relabelled': int(sum(c['helix_residues_relabelled'] for c in per_chain)),
        'strands_relabelled': 0,
        'per_chain': per_chain,
        'note': 'Missing k_est and proline k_est are uninformative for backbone HDX (not evidence of protection); no rows are deleted or shifted. HDX exchange/protection is not a direct mobility measurement.',
    }
    result = (np.array(changed, dtype=object), selected, report)
    if return_measured_masks:
        return (*result, measured_masks)
    return result

def select_hdx_linkers_from_fingerprint(fingerprint_file, selected_masks, *, min_length=4):
    """Select entire internal '-' sections intersecting chosen HDX residues.

    Segment IDs are global zero-based Carbonara indices. The fingerprint must
    already contain any HDX helix overrides and long-linker splitting.
    """
    chains = parse_structures_with_segments(str(fingerprint_file))
    if len(chains) != len(selected_masks):
        raise ValueError('HDX mask count does not match fingerprint chain count')
    selected_ids = []
    selected_info = []
    for chain_i, (chain, mask) in enumerate(zip(chains, selected_masks), start=1):
        mask = np.asarray(mask, dtype=bool)
        if len(mask) != len(chain['sequence']) or len(mask) != len(chain['structure']):
            raise ValueError(f'HDX mask/sequence mismatch in chain {chain_i}')
        for segment in chain['segments']:
            start, end = segment['start'], segment['end']
            if (segment['type'] == '-' and end - start >= min_length
                    and start > 0 and end < len(mask) and mask[start:end].any()):
                selected_ids.append(int(segment['segment_num']))
                selected_info.append({
                    'segment_id_0based_global': int(segment['segment_num']),
                    'chain': chain_i, 'start_1based': start + 1,
                    'end_1based_inclusive': end,
                    'length': end - start,
                    'hdx_selected_residues_in_segment': int(mask[start:end].sum()),
                })
    return selected_ids, selected_info


def select_hdx_linkers_with_automatic_fallback(
    fingerprint_file, selected_masks, measured_masks, automatic_linkers,
    *, min_length=4, hdx_only=False,
):
    """Combine HDX evidence with Carbonara's ordinary linker selection.

    The fingerprint already includes optional H->'-' HDX overrides and
    long-linker splitting. We NEVER change it here.

    For each internal whole coil section of >= min_length residues:
      1. At least one HDX-selected position -> select via HDX, even if
         automatic selection rejects this linker (experimental guidance).
      2. No HDX measurements anywhere in section -> inherit the ordinary
         automatic selection, unless hdx_only=True.
      3. Some measured positions, none selected -> do not select from auto.
         This makes measured data the authority over the entire section.

    A section containing both measured and unmeasured positions is explicitly
    recorded as mixed coverage; it is not mislabelled 'unmeasured'.
    All eligible sections are reported, including ones not selected.

    Returns (selected_global_section_ids, per_section_report).
    """
    chains = parse_structures_with_segments(str(fingerprint_file))
    if not (len(chains) == len(selected_masks) == len(measured_masks)):
        raise ValueError('HDX selection, coverage and fingerprint chain counts differ')

    auto_ids = set(int(x) for x in automatic_linkers)
    selected_ids = []
    section_report = []
    for chain_i, (chain, positives, observed) in enumerate(
        zip(chains, selected_masks, measured_masks), start=1
    ):
        positives = np.asarray(positives, dtype=bool)
        observed = np.asarray(observed, dtype=bool)
        chain_length = len(chain['sequence'])
        if (len(positives) != chain_length or len(observed) != chain_length
                or len(chain['structure']) != chain_length):
            raise ValueError(f'HDX masks/sequence mismatch in chain {chain_i}')
        if np.any(positives & ~observed):
            raise ValueError(f'HDX-selected residues lack measurements in chain {chain_i}')
        for section in chain['segments']:
            start, end = int(section['start']), int(section['end'])
            if not (section['type'] == '-' and end - start >= min_length
                    and start > 0 and end < chain_length):
                continue
            sid = int(section['segment_num'])
            n_positive = int(positives[start:end].sum())
            n_measured = int(observed[start:end].sum())
            n_unknown = end - start - n_measured
            auto_selected = sid in auto_ids
            if n_positive:
                reason = 'hdx_selected'
                keep = True
            elif n_measured:
                reason = 'hdx_measured_below_threshold'
                keep = False
            elif hdx_only:
                reason = 'unmeasured_hdx_only'
                keep = False
            elif auto_selected:
                reason = 'automatic_unmeasured'
                keep = True
            else:
                reason = 'automatic_unmeasured_rejected'
                keep = False
            if keep:
                selected_ids.append(sid)
            section_report.append({
                'segment_id_0based_global': sid,
                'chain': chain_i,
                'start_1based': start + 1,
                'end_1based_inclusive': end,
                'length': end - start,
                'hdx_selected_residues_in_segment': n_positive,
                'measured_residues_in_segment': n_measured,
                'unmeasured_residues_in_segment': n_unknown,
                'automatically_selected': auto_selected,
                'selected': keep,
                'selection_reason': reason,
                'mixed_coverage': n_measured > 0 and n_unknown > 0,
            })
    return selected_ids, section_report


# Carbonara requires at least three C-alpha points for most transition frames.
# This is deliberately an input-geometry policy, not a change to the empirical
# curvature/torsion distributions in randomMolGen.cpp.
def _carbonara_structure_runs(labels):
    """Return (type, start, exclusive_end) for each contiguous section."""
    labels = list(labels)
    if not labels:
        return []
    out = []
    start = 0
    for i in range(1, len(labels) + 1):
        if i == len(labels) or labels[i] != labels[start]:
            out.append((str(labels[start]), start, i))
            start = i
    return out


def regularize_secondary_structure_for_carbonara(
    secondary_structure_chains, *, measured_masks=None, flexible_masks=None,
    min_frame_length=3, min_sampled_linker_length=4,
):
    """Conservatively regularise short sections before fingerprint numbering.

    The input is a chainwise series of H, S, - labels, after any HDX helices
    have been overridden, and before long-linker splitting. Coordinates and
    sequence remain untouched. Chain boundaries are strict.

    Rules:
      * Isolated 1-2 residue helices (often left behind by HDX overrides) may
        become coil only if adjacent to coil AND every such residue is either
        HDX-fast or unmeasured. Slow/medium-measured helices are protected.
      * A three-residue H island entirely HDX-fast and flanked by two coils
        may be absorbed to bridge those coils (a special HDX-guided exception).
      * Very short coils (<3), or fast-measured internal 3-residue coils, can
        acquire neighbouring helix residues if the remaining helix stays >=3,
        and those residues aren't HDX-measured below the fast threshold.
      * No strand is automatically relabelled. Unrepairable short H/S/coil
        sections are reported as unsafe and rejected. We do not invent a
        helix by converting an unsupported coil into H.
      * HDX-fast coils must be long enough to be selectable (>=4) if internal.

    Important: this only enforces structural frame-length conditions. It does
    not prove all stochastic coordinate-generation attempts will succeed.

    Returns (modified_secondary_chains, JSON-serialisable report).
    """
    if min_frame_length < 3 or min_sampled_linker_length < 4:
        raise ValueError('Carbonara frame minimum must be >=3 and sampled linker minimum >=4')
    chains = list(secondary_structure_chains)
    if measured_masks is not None and len(measured_masks) != len(chains):
        raise ValueError('HDX measured masks do not match fingerprint chain count')
    if flexible_masks is not None and len(flexible_masks) != len(chains):
        raise ValueError('HDX selected masks do not match fingerprint chain count')
    output = []
    report = {
        'policy': 'conservative_frame_compatibility',
        'minimum_frame_section_length': int(min_frame_length),
        'minimum_hdx_selected_internal_coil_length': int(min_sampled_linker_length),
        'strand_relabelling_allowed': False,
        'measured_below_threshold_relabelling_allowed': False,
        'edits': [], 'unresolved': [], 'chains': [],
    }
    def sections(arr):
        return _carbonara_structure_runs(arr)

    for cidx, ss in enumerate(chains):
        arr = np.asarray(list(ss) if isinstance(ss, str) else ss, dtype='<U1').copy().reshape(-1)
        n = len(arr)
        if n == 0 or not set(arr.tolist()).issubset({'H','S','-'}):
            raise ValueError(f'Chain {cidx+1} has empty or invalid secondary structure labels')
        measured = (np.zeros(n, dtype=bool) if measured_masks is None
                    else np.asarray(measured_masks[cidx], dtype=bool).reshape(-1))
        positive = (np.zeros(n, dtype=bool) if flexible_masks is None
                    else np.asarray(flexible_masks[cidx], dtype=bool).reshape(-1))
        if len(measured) != n or len(positive) != n or np.any(positive & ~measured):
            raise ValueError(f'Chain {cidx+1}: invalid HDX mask length or HDX-selected unmeasured residue')
        protected = measured & ~positive
        chain_edits = []
        sections_before = len(sections(arr))

        def log_edit(kind, a, b, previous, new, justification):
            edit = {'chain': cidx+1, 'start_1based': int(a+1),
                    'end_1based_inclusive': int(b), 'length': int(b-a),
                    'old_label': previous, 'new_label': new,
                    'action': kind, 'reason': justification}
            chain_edits.append(edit)
            report['edits'].append(edit)

        # Pass 1: absorb truncated helices. Re-run because merging can change
        # neighbouring runs, but edits only move H -> '-' (cannot oscillate).
        while True:
            change = False
            for j, (typ, start, end) in enumerate(sections(arr)):
                if typ != 'H' or end-start >= min_frame_length:
                    continue
                adjacent_coil = (start > 0 and arr[start-1] == '-') or (end < n and arr[end] == '-')
                if adjacent_coil and not protected[start:end].any():
                    arr[start:end] = '-'
                    log_edit('absorb_short_helix', start, end, 'H', '-',
                             'Too short to define a transition frame; adjoining coil, no protected HDX measurements')
                    change = True
                    break
            if not change:
                break

        # Pass 1b: a three-residue H island wholly supported by FAST HDX
        # evidence, flanked on both sides by coil, is not a robust helix to
        # preserve. Converting it to coil merges two otherwise short linkers.
        # Do NOT generalise this to other H3 sections: these have the minimum
        # three points needed for the empirical transition frame.
        while True:
            change = False
            for typ, start, end in sections(arr):
                if (typ == 'H' and end-start == min_frame_length
                        and start > 0 and end < n
                        and arr[start-1] == '-' and arr[end] == '-'
                        and positive[start:end].all()):
                    arr[start:end] = '-'
                    log_edit('bridge_fast_hdx_helix_island', start, end, 'H', '-',
                             'Three-residue helix between coils, all residues fast-exchanging in HDX; '
                             'merges short linkers without changing protected structural assignments')
                    change = True
                    break
            if not change:
                break

        # Pass 2: extend only coil sections that cannot supply a 3-point frame,
        # or contain fast HDX and would otherwise be too short for selection.
        # Borrow edge residues only from an H run that remains >=3 afterwards.
        # A high-confidence helix is not demolished to repair a short linker.
        while True:
            change = False
            for typ, start, end in sections(arr):
                if typ != '-':
                    continue
                target = min_frame_length
                if start > 0 and end < n and positive[start:end].any():
                    target = max(target, min_sampled_linker_length)
                if end-start >= target:
                    continue
                need = target - (end-start)
                for _ in range(need):
                    candidates = []
                    for neighbor in (start-1, end):
                        if neighbor < 0 or neighbor >= n or arr[neighbor] != 'H' or protected[neighbor]:
                            continue
                        for t,a,b in sections(arr):
                            if a <= neighbor < b:
                                if b-a-1 >= min_frame_length:
                                    # Prefer experimentally fast residues; then unmeasured.
                                    candidates.append((0 if positive[neighbor] else 1, neighbor))
                                break
                    if not candidates:
                        break
                    _, idx = min(candidates, key=lambda p: (p[0], p[1]))
                    arr[idx] = '-'
                    log_edit('extend_short_coil', idx, idx+1, 'H', '-',
                             'Needed for transition frame/HDX linker length; retained >=3 neighbouring helix')
                    change = True
                    # Expansion changes the boundaries; derive them again below.
                    if idx == start-1: start -= 1
                    else: end += 1
                if change:
                    break
            if not change:
                break

        # Validation is deliberately strict for ALL short H/S sections. We do
        # not silently change a strand or protected helix to make C++ run.
        chain_issues = []
        for typ,start,end in sections(arr):
            ln = end-start
            requirement = min_frame_length
            if typ == '-' and start > 0 and end < n and positive[start:end].any():
                requirement = max(requirement, min_sampled_linker_length)
            if ln < requirement:
                reason = ('short_beta_strand_not_relabelled' if typ == 'S' else
                          'short_protected_or_isolated_helix' if typ == 'H' else
                          'short_coil_cannot_extend_without_compromising_helix_or_HDX')
                issue = {'chain': cidx+1, 'type': typ, 'start_1based': int(start+1),
                         'end_1based_inclusive': int(end), 'length': int(ln),
                         'required_length': int(requirement), 'reason': reason,
                         'measured_below_threshold': int(protected[start:end].sum()),
                         'hdx_selected': int(positive[start:end].sum())}
                chain_issues.append(issue)
                report['unresolved'].append(issue)

        output.append(arr)
        report['chains'].append({
            'chain': cidx+1, 'length': int(n), 'edited_residues': len(chain_edits),
            'n_edits': len(chain_edits), 'unresolved': len(chain_issues),
            'sections_before': sections_before, 'sections_after': len(sections(arr)),
        })

    report['n_edits'] = len(report['edits'])
    report['n_unresolved'] = len(report['unresolved'])
    return np.array(output, dtype=object), report


def validate_carbonara_secondary_sections(
    secondary_structure_chains, *, min_frame_length=3, min_flexible_length=4,
    flexible_masks=None,
):
    """Verify section frame lengths after long-linker splitting (no edits)."""
    chains = list(secondary_structure_chains)
    if flexible_masks is not None and len(flexible_masks) != len(chains):
        raise ValueError('Mismatch between flexible masks and chain count')
    issues = []
    for ci, ss in enumerate(chains):
        arr = np.asarray(list(ss) if isinstance(ss, str) else ss, dtype='<U1').reshape(-1)
        positive = (np.zeros(len(arr), dtype=bool) if flexible_masks is None
                    else np.asarray(flexible_masks[ci], dtype=bool).reshape(-1))
        if len(arr) != len(positive):
            raise ValueError(f'Chain {ci+1}: flexibility mask length mismatch')
        for typ,a,b in _carbonara_structure_runs(arr):
            required = min_frame_length
            if typ == '-' and a > 0 and b < len(arr) and positive[a:b].any():
                required = max(required, min_flexible_length)
            if b-a < required:
                issues.append({'chain': ci+1, 'type': typ, 'start_1based': a+1,
                               'end_1based_inclusive': b, 'length': b-a,
                               'required_length': required})
    return issues
