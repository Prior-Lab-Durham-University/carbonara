#!/usr/bin/env python3

import argparse
import os
# Force PDBFixer/OpenMM setup-time repair onto the CPU by default.
# Some HPC desktop sessions expose a broken OpenCL platform, which can segfault
# before Python can raise a useful exception. Users can override this externally.
os.environ.setdefault("OPENMM_DEFAULT_PLATFORM", "CPU")
import shutil
import sys
import pickle
import shlex
import json
import subprocess
import re
from pathlib import Path
from string import ascii_uppercase
from typing import Optional, List

import CarbonaraDataTools as cdt
import numpy as np


def _openmm_cpu_platform_or_none():
    """Return the OpenMM CPU platform if available, otherwise None."""
    try:
        from openmm import Platform
        return Platform.getPlatformByName("CPU")
    except Exception:
        return None


def _pdbfixer_with_safe_platform(PDBFixerClass, filename):
    """Construct PDBFixer while explicitly preferring CPU over OpenCL/CUDA.

    PDBFixer versions that do not accept a platform keyword fall back to the
    normal constructor, but OPENMM_DEFAULT_PLATFORM is still set above.
    """
    cpu_platform = _openmm_cpu_platform_or_none()
    if cpu_platform is not None:
        try:
            return PDBFixerClass(filename=filename, platform=cpu_platform)
        except TypeError:
            pass
    return PDBFixerClass(filename=filename)


def _is_cif_path(path: str) -> bool:
    """Return True for mmCIF/PDBx CIF filenames, including .gz variants."""
    suffixes = [s.lower() for s in Path(str(path)).suffixes]
    return ".cif" in suffixes or ".mmcif" in suffixes


def _open_text_maybe_gz(path: str):
    """Open plain or gzipped text files without forcing callers to care."""
    path = str(path)
    if path.lower().endswith(".gz"):
        import gzip
        return gzip.open(path, "rt", errors="replace")
    return open(path, "rt", errors="replace")


def preflight_mmcif_atom_site_table(cif_path: str) -> None:
    """
    Fail early with a clear message if an mmCIF _atom_site loop has truncated
    ATOM/HETATM rows.

    OpenMM/PDBFixer can throw a cryptic IndexError when a CIF row has fewer
    fields than the _atom_site loop declares.  This preflight catches the common
    failure mode before Carbonara calls cdt.pull_structure_from_pdb().

    It is deliberately conservative: valid PDB files and non-CIF files are left
    untouched, and a well-formed CIF passes silently.
    """
    if not _is_cif_path(cif_path):
        return

    in_loop = False
    in_atom_site_loop = False
    atom_site_columns = []
    checked_rows = 0

    with _open_text_maybe_gz(cif_path) as fh:
        for line_no, raw in enumerate(fh, start=1):
            line = raw.strip()
            if not line:
                continue

            if line == "#":
                in_loop = False
                in_atom_site_loop = False
                atom_site_columns = []
                continue

            low = line.lower()
            if low.startswith("loop_"):
                in_loop = True
                in_atom_site_loop = False
                atom_site_columns = []
                continue

            if not in_loop:
                continue

            if low.startswith("_atom_site."):
                in_atom_site_loop = True
                atom_site_columns.append(line.split()[0])
                continue

            if low.startswith("_") or low.startswith("data_") or low.startswith("save_"):
                in_loop = False
                in_atom_site_loop = False
                atom_site_columns = []
                continue

            if not in_atom_site_loop or not atom_site_columns:
                continue

            # We only need to guard atom records. Other CIF loop rows can have
            # non-atom content and are irrelevant here.
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue

            try:
                fields = shlex.split(line, posix=True)
            except ValueError as exc:
                raise ValueError(
                    f"Malformed CIF: could not parse _atom_site row at line {line_no} "
                    f"of {cif_path}: {exc}. Row preview: {line[:120]!r}"
                ) from exc

            checked_rows += 1
            expected = len(atom_site_columns)
            found = len(fields)
            if found < expected:
                missing = atom_site_columns[found:]
                preview = line[:160]
                raise ValueError(
                    "Malformed CIF: incomplete _atom_site ATOM/HETATM row at "
                    f"line {line_no} of {cif_path}. The loop declares {expected} "
                    f"columns but this row has only {found}. Missing columns start "
                    f"with: {missing[:6]}. Row preview: {preview!r}. "
                    "Re-download the CIF or provide a PDB file."
                )

    if checked_rows == 0:
        raise ValueError(
            f"Malformed CIF: no ATOM/HETATM rows were found in the _atom_site loop of {cif_path}."
        )



def _is_pdb_path(path: str) -> bool:
    suffixes = [s.lower() for s in Path(str(path)).suffixes]
    return ".pdb" in suffixes or ".ent" in suffixes


def _looks_like_mmcif_text(path: str, max_lines: int = 80) -> bool:
    """Conservative content sniff for a CIF/ModelCIF file saved with the wrong extension."""
    try:
        with _open_text_maybe_gz(path) as fh:
            seen_nonempty = 0
            for raw in fh:
                line = raw.strip()
                if not line:
                    continue
                seen_nonempty += 1
                low = line.lower()
                if low.startswith("data_"):
                    return True
                if low.startswith("_entry.") or low.startswith("_atom_site.") or low.startswith("_audit_conform."):
                    return True
                if low.startswith("loop_"):
                    # A PDB file should not have mmCIF loop_ syntax near the top.
                    return True
                if low.startswith(("header", "atom", "hetatm", "model", "remark", "title", "compnd", "source", "seqres")):
                    return False
                if seen_nonempty >= max_lines:
                    break
    except OSError:
        return False
    return False


def preflight_structure_extension_matches_content(structure_path: str) -> None:
    """Fail early if an mmCIF/ModelCIF text file has been named .pdb."""
    if _is_pdb_path(structure_path) and _looks_like_mmcif_text(structure_path):
        raise ValueError(
            "Structure file appears to be mmCIF/ModelCIF text but has a PDB extension: "
            f"{structure_path}. Carbonara chooses the structure parser from the filename "
            "extension; parsing mmCIF text as fixed-width PDB can produce cryptic MDTraj "
            "errors. Rename the file to .cif or .mmcif and rerun."
        )


def _safe_int_residue_id(value):
    """Best-effort conversion of a PDBFixer/OpenMM residue id to an integer resSeq."""
    if value is None:
        return None
    text_value = str(value).strip()
    match = re.match(r"^-?\d+", text_value)
    if not match:
        return None
    try:
        return int(match.group(0))
    except ValueError:
        return None


def _filter_missing_residue_entries(fixer, entries, max_gap: int, internal_only: bool = True):
    """Filter PDBFixer missingResidues entries to internal, reasonably sized gaps."""
    chains = list(fixer.topology.chains())
    kept = {}
    skipped = []
    for key, names in list(entries.items()):
        try:
            chain_index, residue_index = int(key[0]), int(key[1])
        except Exception:
            skipped.append({"key": repr(key), "n_residues": len(names), "reason": "unrecognised_key"})
            continue

        residues = list(chains[chain_index].residues()) if 0 <= chain_index < len(chains) else []
        n_residues_in_chain = len(residues)
        names = [str(name).upper() for name in names]
        reason = None
        if internal_only and (residue_index <= 0 or residue_index >= n_residues_in_chain):
            reason = "terminal_gap"
        elif len(names) == 0:
            reason = "empty_gap"
        elif max_gap is not None and max_gap > 0 and len(names) > max_gap:
            reason = f"gap_longer_than_{max_gap}"

        if reason is not None:
            skipped.append({
                "chain_index": chain_index,
                "insert_before_residue_index": residue_index,
                "n_residues": len(names),
                "reason": reason,
            })
            continue
        kept[(chain_index, residue_index)] = names
    return kept, skipped


def _infer_missing_residues_from_number_gaps(fixer, residue_name: str = "GLY", max_gap: int = 80):
    """Infer internal residue gaps from adjacent residue ids in the loaded topology."""
    residue_name = residue_name.upper()
    inferred = {}
    gaps = []
    skipped = []
    for chain_index, chain in enumerate(fixer.topology.chains()):
        residues = list(chain.residues())
        chain_id = getattr(chain, "id", str(chain_index))
        for i in range(len(residues) - 1):
            left = residues[i]
            right = residues[i + 1]
            left_id = _safe_int_residue_id(getattr(left, "id", None))
            right_id = _safe_int_residue_id(getattr(right, "id", None))
            if left_id is None or right_id is None:
                continue
            missing_count = right_id - left_id - 1
            if missing_count <= 0:
                continue
            gap_info = {
                "source": "residue_number_gap",
                "chain_index": chain_index,
                "chain_id": chain_id,
                "left_residue_id": left_id,
                "right_residue_id": right_id,
                "insert_before_residue_index": i + 1,
                "n_residues": missing_count,
                "residue_name": residue_name,
            }
            if max_gap is not None and max_gap > 0 and missing_count > max_gap:
                gap_info["reason"] = f"gap_longer_than_{max_gap}"
                skipped.append(gap_info)
                continue
            key = (chain_index, i + 1)
            # Do not overwrite a SEQRES-derived entry if PDBFixer already has one.
            inferred.setdefault(key, [residue_name] * missing_count)
            gaps.append(gap_info)
    return inferred, gaps, skipped


def pdbfixer_prepare_structure_before_carbonara(
    structure_path: str,
    refine_dir: str,
    enabled: bool = True,
    residue_name: str = "GLY",
    max_gap: int = 80,
    internal_only: bool = True,
) -> str:
    """
    Use PDBFixer before Carbonara reads the structure to build internal missing-residue gaps.

    This deliberately writes a repaired PDB and then lets the ordinary Carbonara reader run on
    that file.  We do not alter CarbonaraDataTools section numbering/selection logic after the
    fact, because that is exactly where label drift can occur.
    """
    report_path = os.path.join(refine_dir, "missing_residue_pdbfixer_report.json")
    report = {
        "enabled": bool(enabled),
        "input_structure": str(structure_path),
        "output_structure": None,
        "residue_name_for_number_gap_inference": str(residue_name).upper(),
        "max_gap": int(max_gap),
        "internal_only": bool(internal_only),
        "status": "disabled" if not enabled else "not_run",
        "openmm_default_platform": os.environ.get("OPENMM_DEFAULT_PLATFORM"),
        "pdbfixer_missing_residues": [],
        "number_gap_inferred_missing_residues": [],
        "skipped_missing_residues": [],
        "total_inserted_residues": 0,
    }

    def write_report():
        with open(report_path, "w") as fh:
            json.dump(report, fh, indent=2, sort_keys=True)

    if not enabled:
        write_report()
        return structure_path

    residue_name = str(residue_name).upper()
    standard_residues = {
        "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
        "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    }
    if residue_name not in standard_residues:
        raise ValueError(
            f"--fix_missing_residue_name must be a standard 3-letter amino-acid code; got {residue_name!r}"
        )

    try:
        from pdbfixer import PDBFixer
        from openmm.app import PDBFile
    except Exception as exc:
        report["status"] = "failed_import_pdbfixer"
        report["error"] = str(exc)
        write_report()
        raise ImportError(
            "--fix_missing_residues requires pdbfixer and openmm, which Carbonara already uses for CIF/PDB repair."
        ) from exc

    try:
        fixer = _pdbfixer_with_safe_platform(PDBFixer, structure_path)

        # Let PDBFixer use its normal SEQRES/mmCIF machinery where available.
        try:
            fixer.findMissingResidues()
            pdbfixer_entries_raw = dict(fixer.missingResidues)
        except Exception:
            pdbfixer_entries_raw = {}

        pdbfixer_entries, skipped_seqres = _filter_missing_residue_entries(
            fixer, pdbfixer_entries_raw, max_gap=max_gap, internal_only=internal_only
        )

        for key, names in sorted(pdbfixer_entries.items()):
            report["pdbfixer_missing_residues"].append({
                "chain_index": int(key[0]),
                "insert_before_residue_index": int(key[1]),
                "n_residues": len(names),
                "residue_names": list(names),
            })

        inferred_entries, inferred_gaps, skipped_inferred = _infer_missing_residues_from_number_gaps(
            fixer, residue_name=residue_name, max_gap=max_gap
        )
        report["number_gap_inferred_missing_residues"] = inferred_gaps
        report["skipped_missing_residues"].extend(skipped_seqres)
        report["skipped_missing_residues"].extend(skipped_inferred)

        # Merge entries.  Prefer PDBFixer's sequence-derived residue identities where present;
        # otherwise use residue-number gaps filled with residue_name (default GLY).
        final_entries = dict(pdbfixer_entries)
        for key, names in inferred_entries.items():
            final_entries.setdefault(key, names)

        if not final_entries:
            report["status"] = "no_internal_missing_residues_detected"
            write_report()
            return structure_path

        fixer.missingResidues.clear()
        fixer.missingResidues.update(final_entries)
        report["total_inserted_residues"] = int(sum(len(v) for v in final_entries.values()))

        fixer.findMissingAtoms()
        try:
            fixer.addMissingAtoms(seed=0)
        except TypeError:
            fixer.addMissingAtoms()

        out_path = os.path.join(refine_dir, Path(str(structure_path)).stem + "_pdbfixer_missing_residues.pdb")
        with open(out_path, "w") as handle:
            try:
                PDBFile.writeFile(fixer.topology, fixer.positions, handle, keepIds=True)
            except TypeError:
                PDBFile.writeFile(fixer.topology, fixer.positions, handle)

        report["status"] = "fixed_missing_residues"
        report["output_structure"] = out_path
        write_report()
        print(
            "PDBFixer missing-residue repair: inserted "
            f"{report['total_inserted_residues']} residue(s) across {len(final_entries)} internal gap(s); "
            f"wrote {out_path}"
        )
        return out_path

    except Exception as exc:
        report["status"] = "failed"
        report["error"] = repr(exc)
        write_report()
        raise



def _pdb_atom_line(atom):
    """Write a minimal standard fixed-width PDB ATOM record."""
    return (
        f"ATOM  {atom['serial']:5d} {atom['name']:<4s}{atom.get('altloc',' '):1s}"
        f"{atom['resname']:>3s} {atom['chain'][:1]:1s}"
        f"{int(atom['resseq']):4d}{atom.get('icode',' '):1s}   "
        f"{float(atom['x']):8.3f}{float(atom['y']):8.3f}{float(atom['z']):8.3f}"
        f"{float(atom.get('occupancy',1.0)):6.2f}{float(atom.get('bfactor',50.0)):6.2f}          "
        f"{atom.get('element', atom['name'].strip()[0]).strip()[:2]:>2s}\n"
    )


def _parse_pdb_atom_record_for_bridge(line: str):
    if not line.startswith(('ATOM  ', 'HETATM')):
        return None
    try:
        resseq = int(line[22:26])
        x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
    except Exception:
        return None
    name = line[12:16].strip()
    element = line[76:78].strip() if len(line) >= 78 else ''
    if not element:
        element = ''.join(ch for ch in name if ch.isalpha())[:1] or 'C'
    return {
        'record': line[0:6].strip(),
        'name': name,
        'altloc': line[16:17] if len(line) > 16 else ' ',
        'resname': line[17:20].strip(),
        'chain': (line[21:22] if len(line) > 21 else ' ') or ' ',
        'resseq': resseq,
        'icode': line[26:27] if len(line) > 26 else ' ',
        'x': x,
        'y': y,
        'z': z,
        'occupancy': float(line[54:60]) if len(line) >= 60 and line[54:60].strip() else 1.0,
        'bfactor': float(line[60:66]) if len(line) >= 66 and line[60:66].strip() else 50.0,
        'element': element,
    }


def _bridge_unit_perpendicular(direction):
    import numpy as _np
    direction = _np.asarray(direction, dtype=float)
    norm = float(_np.linalg.norm(direction))
    if norm <= 1e-12:
        direction = _np.array([1.0, 0.0, 0.0])
    else:
        direction = direction / norm
    axis = _np.array([0.0, 0.0, 1.0])
    if abs(float(_np.dot(direction, axis))) > 0.9:
        axis = _np.array([0.0, 1.0, 0.0])
    perp = _np.cross(direction, axis)
    pnorm = float(_np.linalg.norm(perp))
    return perp / pnorm if pnorm > 1e-12 else _np.array([0.0, 1.0, 0.0])


def _make_ca_arc_points(p0, p1, n_missing: int, ca_spacing: float = 3.8):
    """Return missing CA positions on a smooth circular bridge.

    Adjacent CA positions, including the two observed endpoint CAs, are placed at
    approximately ``ca_spacing`` (exactly so to floating-point precision whenever
    the endpoint geometry permits).  This is much safer for Carbonara than asking
    PDBFixer to choose a compact loop geometry and, unlike a straight interpolation,
    it can accommodate a contour length much larger than the endpoint separation.
    """
    import numpy as _np

    p0 = _np.asarray(p0, dtype=float)
    p1 = _np.asarray(p1, dtype=float)
    n_missing = int(n_missing)
    if n_missing <= 0:
        return []

    n_segments = n_missing + 1
    spacing = float(ca_spacing)
    if spacing <= 0:
        raise ValueError("ca_spacing must be > 0")

    chord_vec = p1 - p0
    chord = float(_np.linalg.norm(chord_vec))
    max_reach = n_segments * spacing
    if chord > max_reach + 1e-5:
        raise ValueError(
            f"Cannot bridge {n_missing} missing residues: endpoint CA distance "
            f"{chord:.3f} A exceeds the available {n_segments} x {spacing:.3f} A "
            f"CA contour ({max_reach:.3f} A)."
        )

    # Nearly fully extended: avoid a numerically huge circle radius.
    if chord >= max_reach - 1e-6:
        pts = _np.linspace(p0, p1, n_segments + 1)
        return [pts[i].copy() for i in range(1, n_missing + 1)]

    if chord <= 1e-8:
        # Degenerate case: build a closed regular polygon in an arbitrary plane.
        u = _np.array([1.0, 0.0, 0.0])
        v = _np.array([0.0, 1.0, 0.0])
        delta = 2.0 * _np.pi / n_segments
        radius = spacing / (2.0 * _np.sin(delta / 2.0))
        centre = p0 + radius * v
        r0 = p0 - centre
        normal = _np.array([0.0, 0.0, 1.0])
    else:
        u = chord_vec / chord
        v = _bridge_unit_perpendicular(u)

        # For N equal circular-arc chords of length s, solve
        #   D/s = sin(N*delta/2) / sin(delta/2)
        # on 0 < delta < 2*pi/N.  This branch is monotone from N to 0.
        target = chord / spacing
        lo = 1e-10
        hi = 2.0 * _np.pi / n_segments - 1e-10

        def _ratio(delta):
            return _np.sin(n_segments * delta / 2.0) / _np.sin(delta / 2.0)

        for _ in range(80):
            mid = 0.5 * (lo + hi)
            if _ratio(mid) > target:
                lo = mid
            else:
                hi = mid
        delta = 0.5 * (lo + hi)
        radius = spacing / (2.0 * _np.sin(delta / 2.0))
        total_angle = n_segments * delta
        midpoint = 0.5 * (p0 + p1)
        centre = midpoint + radius * _np.cos(total_angle / 2.0) * v
        r0 = p0 - centre
        normal = _np.cross(u, v)
        nrm = float(_np.linalg.norm(normal))
        if nrm <= 1e-12:
            normal = _np.array([0.0, 0.0, 1.0])
        else:
            normal = normal / nrm

    def _rotate(vec, angle):
        return (
            vec * _np.cos(angle)
            + _np.cross(normal, vec) * _np.sin(angle)
            + normal * float(_np.dot(normal, vec)) * (1.0 - _np.cos(angle))
        )

    # Select the rotation direction that lands on p1.
    total = n_segments * delta
    plus_end = centre + _rotate(r0, total)
    minus_end = centre + _rotate(r0, -total)
    sign = 1.0 if _np.linalg.norm(plus_end - p1) <= _np.linalg.norm(minus_end - p1) else -1.0

    pts = [centre + _rotate(r0, sign * k * delta) for k in range(n_segments + 1)]
    pts[0] = p0.copy()
    pts[-1] = p1.copy()
    return [_np.asarray(pts[i], dtype=float).copy() for i in range(1, n_missing + 1)]

def _make_fake_gly_backbone(ca_points, chain, resseqs, serial_start, bfactor=80.0):
    """Create minimal N/CA/C/O atoms for fake GLY residues around CA arc points."""
    import numpy as _np
    atoms = []
    if not ca_points:
        return atoms, serial_start
    pts = [_np.asarray(p, dtype=float) for p in ca_points]
    perp = _bridge_unit_perpendicular((pts[-1] - pts[0]) if len(pts) > 1 else _np.array([1.0, 0.0, 0.0]))
    serial = int(serial_start)
    for i, ca in enumerate(pts):
        if len(pts) == 1:
            tangent = _np.array([1.0, 0.0, 0.0])
        elif i == 0:
            tangent = pts[1] - pts[0]
        elif i == len(pts) - 1:
            tangent = pts[-1] - pts[-2]
        else:
            tangent = pts[i + 1] - pts[i - 1]
        tnorm = float(_np.linalg.norm(tangent))
        tangent = tangent / tnorm if tnorm > 1e-12 else _np.array([1.0, 0.0, 0.0])
        coords = {
            'N':  ca - 1.30 * tangent,
            'CA': ca,
            'C':  ca + 1.30 * tangent,
            'O':  ca + 1.30 * tangent + 1.00 * perp,
        }
        for name in ('N', 'CA', 'C', 'O'):
            xyz = coords[name]
            atoms.append({
                'serial': serial,
                'name': name,
                'altloc': ' ',
                'resname': 'GLY',
                'chain': chain,
                'resseq': int(resseqs[i]),
                'icode': ' ',
                'x': float(xyz[0]), 'y': float(xyz[1]), 'z': float(xyz[2]),
                'occupancy': 1.0,
                'bfactor': float(bfactor),
                'element': 'N' if name == 'N' else ('O' if name == 'O' else 'C'),
            })
            serial += 1
    return atoms, serial


def bridge_missing_residue_gaps_to_pdb_before_carbonara(
    structure_path: str,
    refine_dir: str,
    enabled: bool = True,
    residue_name: str = 'GLY',
    max_gap: int = 80,
    ca_spacing: float = 3.8,
) -> str:
    """Create a run-local repaired PDB before Carbonara reads the structure.

    The user's input file is *never* changed.  Every original PDB record is copied
    verbatim to the run-local repaired file; only new minimal N/CA/C/O records for
    genuinely absent internal residues are inserted.  This preserves SSBOND, HELIX,
    SHEET, HETATM, CONECT, REMARK and other records from the original file.

    Gap detection uses observed ATOM residues with CA atoms and therefore does not
    mistake modified protein residues such as MSE for missing residues.
    """
    report_path = os.path.join(refine_dir, 'missing_residue_bridge_report.json')
    report = {
        'enabled': bool(enabled),
        'method': 'pre_carbonara_equal_ca_spacing_bridge',
        'input_structure': str(structure_path),
        'output_structure': None,
        'status': 'disabled' if not enabled else 'not_run',
        'residue_name': str(residue_name).upper(),
        'max_gap': int(max_gap),
        'ca_spacing': float(ca_spacing),
        'gaps': [],
        'total_inserted_residues': 0,
        'original_input_modified': False,
    }

    def write_report():
        with open(report_path, 'w') as fh:
            json.dump(report, fh, indent=2, sort_keys=True)

    if not enabled:
        write_report()
        return structure_path

    if str(residue_name).upper() != 'GLY':
        print('WARNING: pre-Carbonara gap bridge currently writes GLY residues; ignoring --fix_missing_residue_name')

    if not _is_pdb_path(structure_path):
        report['status'] = 'not_pdb_no_bridge_applied'
        write_report()
        return structure_path

    with open(structure_path, 'r', errors='replace') as fh:
        original_lines = fh.readlines()

    # Existing atom serials are retained exactly.  New fake atoms get fresh serials
    # above the existing maximum, so CONECT/SSBOND/header records remain untouched.
    max_serial = 0
    for line in original_lines:
        if line.startswith(('ATOM  ', 'HETATM')):
            try:
                max_serial = max(max_serial, int(line[6:11]))
            except Exception:
                pass
    next_serial = max_serial + 1

    # Build residue blocks from ATOM records in file order.  Do not filter on the
    # 3-letter residue name: modified amino acids (e.g. MSE) are real observed
    # residues and must not be mistaken for gaps.
    residues = []
    current = None
    for line_index, line in enumerate(original_lines):
        atom = _parse_pdb_atom_record_for_bridge(line)
        if atom is None or atom['record'] != 'ATOM':
            continue
        key = (atom['chain'], int(atom['resseq']), atom['icode'])
        if current is None or key != current['key']:
            if current is not None:
                residues.append(current)
            current = {
                'key': key,
                'resname': atom['resname'],
                'atoms': [],
                'first_line_index': int(line_index),
                'last_line_index': int(line_index),
            }
        current['atoms'].append(atom)
        current['last_line_index'] = int(line_index)
    if current is not None:
        residues.append(current)

    if not residues:
        report['status'] = 'no_protein_atom_records_no_bridge_applied'
        write_report()
        return structure_path

    def ca_of(residue):
        ca_atoms = [a for a in residue['atoms'] if a['name'].strip().upper() == 'CA']
        if not ca_atoms:
            return None
        # Prefer the highest-occupancy CA if alternate locations are present.
        ca = max(ca_atoms, key=lambda a: float(a.get('occupancy', 0.0)))
        return np.array([ca['x'], ca['y'], ca['z']], dtype=float)

    insert_before = {}
    inserted_total = 0
    inserted_gap_count = 0

    for left, right in zip(residues[:-1], residues[1:]):
        left_chain, left_resseq, left_icode = left['key']
        right_chain, right_resseq, right_icode = right['key']
        if right_chain != left_chain:
            continue

        # Never bridge across an explicit TER record, even if chain IDs happen to
        # be reused afterwards.
        between = original_lines[left['last_line_index'] + 1:right['first_line_index']]
        if any(line.startswith('TER') for line in between):
            continue

        missing_count = int(right_resseq) - int(left_resseq) - 1
        if missing_count <= 0:
            continue

        gap_info = {
            'chain': left_chain,
            'left_resseq': int(left_resseq),
            'right_resseq': int(right_resseq),
            'n_missing': int(missing_count),
        }
        if max_gap is not None and int(max_gap) > 0 and missing_count > int(max_gap):
            gap_info['status'] = f'skipped_gap_longer_than_{int(max_gap)}'
            report['gaps'].append(gap_info)
            continue

        p0 = ca_of(left)
        p1 = ca_of(right)
        if p0 is None or p1 is None:
            gap_info['status'] = 'skipped_missing_endpoint_CA'
            report['gaps'].append(gap_info)
            continue

        try:
            ca_points = _make_ca_arc_points(p0, p1, missing_count, ca_spacing=ca_spacing)
        except ValueError as exc:
            gap_info['status'] = 'skipped_geometry_impossible'
            gap_info['error'] = str(exc)
            report['gaps'].append(gap_info)
            continue

        fake_resseqs = list(range(int(left_resseq) + 1, int(right_resseq)))
        fake_atoms, next_serial = _make_fake_gly_backbone(
            ca_points, left_chain, fake_resseqs, next_serial, bfactor=80.0
        )
        fake_lines = [_pdb_atom_line(atom) for atom in fake_atoms]
        insert_before.setdefault(right['first_line_index'], []).extend(fake_lines)

        full_ca_path = [p0] + [np.asarray(p, dtype=float) for p in ca_points] + [p1]
        ca_steps = [float(np.linalg.norm(b - a)) for a, b in zip(full_ca_path[:-1], full_ca_path[1:])]
        gap_info.update({
            'status': 'inserted',
            'residue_numbers': fake_resseqs,
            'endpoint_ca_distance': float(np.linalg.norm(p1 - p0)),
            'ca_step_min': float(min(ca_steps)),
            'ca_step_max': float(max(ca_steps)),
            'insert_before_pdb_line_1based': int(right['first_line_index'] + 1),
        })
        report['gaps'].append(gap_info)
        inserted_total += missing_count
        inserted_gap_count += 1

    if inserted_total <= 0:
        report['status'] = 'no_internal_missing_residue_gaps_detected'
        write_report()
        return structure_path

    out_path = os.path.join(
        refine_dir,
        Path(str(structure_path)).stem + '_precarbonara_missing_linkers.pdb'
    )

    # Preserve every original line verbatim and insert only the fake residue records.
    with open(out_path, 'w', newline='') as fh:
        for idx, line in enumerate(original_lines):
            for added in insert_before.get(idx, []):
                fh.write(added)
            fh.write(line)

    report['total_inserted_residues'] = int(inserted_total)
    report['output_structure'] = out_path
    report['status'] = 'fixed_missing_residues'
    write_report()
    print(
        'Pre-Carbonara missing-linker bridge: inserted '
        f'{inserted_total} residue(s) across {inserted_gap_count} internal gap(s); '
        f'wrote {out_path}'
    )
    return out_path

def prepare_structure_before_carbonara(
    structure_path: str,
    refine_dir: str,
    enabled: bool = True,
    residue_name: str = 'GLY',
    max_gap: int = 80,
    method: str = 'bridge',
) -> str:
    """Select the pre-Carbonara missing-residue repair strategy."""
    method = str(method).lower()
    if not enabled:
        return structure_path
    if method == 'pdbfixer':
        return pdbfixer_prepare_structure_before_carbonara(
            structure_path,
            refine_dir,
            enabled=True,
            residue_name=residue_name,
            max_gap=max_gap,
            internal_only=True,
        )
    if method == 'bridge':
        # The explicit bridge writer is intentionally PDB-only.  For CIF/mmCIF,
        # retain automatic repair by falling back to PDBFixer rather than silently
        # returning an unrepaired structure.
        if _is_pdb_path(structure_path):
            return bridge_missing_residue_gaps_to_pdb_before_carbonara(
                structure_path,
                refine_dir,
                enabled=True,
                residue_name=residue_name,
                max_gap=max_gap,
            )
        # CIF/mmCIF is checked with PDBFixer because the explicit bridge writer
        # is PDB-only.  If no genuine internal missing residues are detected,
        # PDBFixer returns the original CIF/mmCIF path unchanged and prints no
        # misleading "bridge" message.
        return pdbfixer_prepare_structure_before_carbonara(
            structure_path,
            refine_dir,
            enabled=True,
            residue_name=residue_name,
            max_gap=max_gap,
            internal_only=True,
        )
    raise ValueError(f"Unknown --fix_missing_residue_method {method!r}; expected 'bridge' or 'pdbfixer'")

def split_long_linkers_in_secondary_structure(
    secondary_structure_chains,
    enabled: bool = True,
    max_linker_len: int = 25,
    fake_helix_len: int = 3,
    report_path: str = None,
):
    """
    Break very long flexible ('-') secondary-structure runs into smaller linker
    sections before the Carbonara fingerprint is written.

    The sequence and coordinates are untouched.  We only replace a few residues
    inside an over-long '-' run by a short artificial helical marker (default HHH).
    This gives Carbonara several manageable flexible sections instead of one very
    large linker, while keeping all section labels internally consistent because
    the modification happens before write_fingerprint_file().

    A linker of exactly max_linker_len residues is left unchanged; only runs with
    length > max_linker_len are split.
    """
    max_linker_len = int(max_linker_len)
    fake_helix_len = int(fake_helix_len)
    if max_linker_len < 4:
        raise ValueError("max_linker_len must be >= 4")
    if fake_helix_len < 1:
        raise ValueError("fake_helix_len must be >= 1")

    report = {
        "enabled": bool(enabled),
        "max_linker_len": max_linker_len,
        "fake_helix_len": fake_helix_len,
        "modified_linker_count": 0,
        "chains": [],
    }

    output = []
    for chain_index, ss in enumerate(secondary_structure_chains):
        arr = np.asarray(ss, dtype="<U1").copy()
        chain_report = {"chain_index": int(chain_index), "modified_linkers": []}

        if enabled and len(arr):
            i = 0
            while i < len(arr):
                if arr[i] != '-':
                    i += 1
                    continue
                start = i
                while i < len(arr) and arr[i] == '-':
                    i += 1
                end = i
                run_len = end - start
                if run_len <= max_linker_len:
                    continue

                # Use the smallest number of H...H separators that leaves every
                # remaining flexible chunk <= max_linker_len.  Because separators
                # replace residues rather than insert residues, sequence/coordinate
                # indexing is unchanged.
                n_breaks = 1
                while True:
                    flexible_residues = run_len - n_breaks * fake_helix_len
                    n_chunks = n_breaks + 1
                    if flexible_residues >= n_chunks:
                        largest_chunk = (flexible_residues + n_chunks - 1) // n_chunks
                        if largest_chunk <= max_linker_len:
                            break
                    n_breaks += 1
                    if n_breaks * fake_helix_len >= run_len:
                        raise ValueError(
                            f"Cannot split linker of length {run_len} with fake_helix_len={fake_helix_len}"
                        )

                base, remainder = divmod(flexible_residues, n_chunks)
                chunk_sizes = [base + (1 if k < remainder else 0) for k in range(n_chunks)]

                pos = start
                separator_ranges = []
                for k, chunk_len in enumerate(chunk_sizes):
                    pos += chunk_len
                    if k < n_breaks:
                        sep_start = pos
                        sep_end = pos + fake_helix_len
                        arr[sep_start:sep_end] = 'H'
                        separator_ranges.append([int(sep_start), int(sep_end)])
                        pos = sep_end

                chain_report["modified_linkers"].append({
                    "start_index_0based": int(start),
                    "end_index_exclusive_0based": int(end),
                    "start_residue_1based": int(start + 1),
                    "end_residue_1based": int(end),
                    "original_length": int(run_len),
                    "flexible_chunk_lengths": [int(x) for x in chunk_sizes],
                    "separator_ranges_0based_halfopen": separator_ranges,
                })
                report["modified_linker_count"] += 1

        output.append(arr)
        report["chains"].append(chain_report)

    if report_path is not None:
        with open(report_path, "w") as fh:
            json.dump(report, fh, indent=2, sort_keys=True)

    if enabled:
        if report["modified_linker_count"]:
            print(
                "Long-linker split: broke "
                f"{report['modified_linker_count']} linker section(s) longer than {max_linker_len} residues "
                f"using {fake_helix_len}-residue H separators."
            )
        else:
            print(f"Long-linker split: no linker sections longer than {max_linker_len} residues found.")

    return np.array(output, dtype=object), report

def _snap_simplex_row(w: np.ndarray, snap: float) -> np.ndarray:
    """Snap weights to multiples of `snap` and renormalize to sum to 1."""
    if snap is None or snap <= 0:
        w = np.clip(w, 0.0, 1.0)
        s = float(w.sum())
        return w / s if s > 0 else np.ones_like(w) / len(w)

    w = np.maximum(w, 0.0)
    w = np.round(w / snap) * snap
    s = float(w.sum())
    if s <= 0:
        w[:] = 1.0 / len(w)
    else:
        w = w / s

    # final clean-up
    w = np.clip(w, 0.0, 1.0)
    w = w / float(w.sum())
    return w


def write_sparse_mixture_file(
    outpath: str,
    n: int,
    max_combos: int = 30,
    step_2: float = 0.1,
    dirichlet_alpha: float = 1.0,
    snap: float = 0.05,
    seed: int = 0,
    include_uniform: bool = True,
    include_corners: bool = True,
    decimals: int = 3,
) -> str:
    """
    Write mixture weights to `outpath`.

    Special case:
      - If max_combos == 1 and include_uniform is True, write ONLY the uniform mixture:
        [1/n, 1/n, ..., 1/n]
    """
    if n < 1:
        raise ValueError("mixture_n must be >= 1")
    if max_combos < 1:
        raise ValueError("max_mixture_combos must be >= 1")

    # --- FIX: enforce "one combo => uniform" (when requested) ---
    if max_combos == 1 and include_uniform:
        row = np.ones(n, dtype=float) / float(n)
        os.makedirs(os.path.dirname(outpath) or ".", exist_ok=True)
        fmt = f"{{:.{decimals}f}}"
        with open(outpath, "w") as f:
            f.write(" ".join(fmt.format(x) for x in row) + "\n")
        return outpath
    # -----------------------------------------------------------

    rows: List[np.ndarray] = []

    if n == 1:
        rows = [np.array([1.0], dtype=float)]

    elif n == 2:
        if step_2 <= 0 or step_2 > 1:
            raise ValueError("mixture_step must be in (0, 1]")
        k = int(round(1.0 / step_2))
        for i in range(k + 1):
            a = i * step_2
            b = 1.0 - a
            rows.append(np.array([a, b], dtype=float))

        # force exact endpoints if step doesn't divide 1 cleanly
        if rows[0][0] != 0.0 or rows[0][1] != 1.0:
            rows.insert(0, np.array([0.0, 1.0], dtype=float))
        if rows[-1][0] != 1.0 or rows[-1][1] != 0.0:
            rows.append(np.array([1.0, 0.0], dtype=float))

    else:
        rng = np.random.default_rng(seed)

        # If you prefer uniform to appear first even when max_combos>1,
        # you can append it before corners. (Optional.)
        if include_uniform:
            rows.append(np.ones(n, dtype=float) / float(n))

        if include_corners:
            for i in range(n):
                e = np.zeros(n, dtype=float)
                e[i] = 1.0
                rows.append(e)

        alpha_vec = np.full(n, float(dirichlet_alpha), dtype=float)
        if np.any(alpha_vec <= 0):
            raise ValueError("mixture_dirichlet_alpha must be > 0")

        attempts = 0
        max_attempts = max(200, 20 * max_combos)
        while len(rows) < max_combos and attempts < max_attempts:
            w = rng.dirichlet(alpha_vec)
            w = _snap_simplex_row(w, snap)
            rows.append(w)
            attempts += 1

    # snap + de-duplicate preserving order
    snapped: List[np.ndarray] = []
    seen = set()
    for r in rows:
        r = _snap_simplex_row(np.array(r, dtype=float), snap)
        key = tuple(np.round(r, 8))
        if key not in seen:
            seen.add(key)
            snapped.append(r)

    snapped = snapped[:max_combos]

    os.makedirs(os.path.dirname(outpath) or ".", exist_ok=True)
    fmt = f"{{:.{decimals}f}}"
    with open(outpath, "w") as f:
        for r in snapped:
            f.write(" ".join(fmt.format(x) for x in r) + "\n")

    return outpath

def replicate_numbered_files(refine_dir: str, n: int) -> None:
    """
    Replicate per-structure inputs for ensemble refinement:
      - coordinates1.dat -> coordinates{i}.dat
      - fingerPrint1.dat -> fingerPrint{i}.dat
      - varyingSectionSecondary1.dat -> varyingSectionSecondary{i}.dat
    Additionally replicate fixedDistanceConstraints{i}.dat ONLY if fixedDistanceConstraints1.dat exists.
    """
    if n <= 1:
        return

    def cp(src: str, dst: str) -> None:
        shutil.copy2(os.path.join(refine_dir, src), os.path.join(refine_dir, dst))

    # Always replicate these
    for i in range(2, n + 1):
        cp("coordinates1.dat", f"coordinates{i}.dat")
        cp("fingerPrint1.dat", f"fingerPrint{i}.dat")
        cp("varyingSectionSecondary1.dat", f"varyingSectionSecondary{i}.dat")

    # Only replicate fixed distance constraints if they exist
    fd1 = os.path.join(refine_dir, "fixedDistanceConstraints1.dat")
    if os.path.exists(fd1):
        for i in range(2, n + 1):
            cp("fixedDistanceConstraints1.dat", f"fixedDistanceConstraints{i}.dat")


def write_runme(
    working_path,
    fit_name,
    fit_n_times,
    min_q,
    max_q,
    max_q_start,
    max_fit_steps,
    no_structures: int = 1,
    pairedQ: bool = False,
    rotation: bool = False,
    max_backmap: int = 3,
    defer_backmap_seconds: int = 600,
    do_foxs: bool = True,
    backend: str = "modeller",
    cg2all_exec: str = "./bin/micromamba run -p /root/micromamba/envs/cg2all convert_cg2all",
    disulfide_constraints_file: str = "",
    foxs_cmd_default: str = "pyfoxs",
    python_exe: Optional[str] = None,
    terminate_on_foxs: bool = False,
    terminate_threshold: float = 2.5,
    terminate_confirmation_count: int = 1,
):
    curr = os.getcwd()
    script_name = "RunMe_" + str(fit_name) + ".sh"
    run_file = os.path.join(curr, script_name)
    
    carbonaradir = os.path.dirname(os.path.realpath(sys.argv[0]))

    # Use the Python that ran setup by default. The generated shell script still
    # allows PYTHON_EXE to override this at runtime, which is useful for notebooks,
    # HPC modules, or advanced users.
    python_exe = python_exe or sys.executable

    # Path to the data directory (relative to ROOT)
    data_path = f"carbonara_runs/{fit_name}"

    # Create new directories
    new_data_dir = os.path.join(curr, "carbonara_runs", fit_name)
    fitdata_dir = os.path.join(new_data_dir, "fitdata")

    os.makedirs(new_data_dir, exist_ok=True)
    os.makedirs(fitdata_dir, exist_ok=True)

    # Copy necessary files from refine_dir to new_data_dir
    try:
        files_to_copy = ["Saxs.dat", "mixtureFile.dat", "chainLengths.dat"]

        # Copy numbered per-structure files 1..no_structures
        for i in range(1, int(no_structures) + 1):
            files_to_copy.extend(
                [
                    f"fingerPrint{i}.dat",
                    f"varyingSectionSecondary{i}.dat",
                    f"coordinates{i}.dat",
                ]
            )

        # Only copy fixedDistanceConstraints{i}.dat if constraints1 exists
        if os.path.exists(os.path.join(working_path, "fixedDistanceConstraints1.dat")):
            for i in range(1, int(no_structures) + 1):
                files_to_copy.append(f"fixedDistanceConstraints{i}.dat")

        for file in files_to_copy:
            source = os.path.join(working_path, file)
            destination = os.path.join(new_data_dir, file)

            if os.path.exists(source):
                shutil.copy2(source, destination)
                print(f"Copied {file} to {destination}")
            else:
                print(f"Warning: Source file {source} not found")

        # Create an empty redundant file
        with open(os.path.join(new_data_dir, "redundant"), "w") as f:
            f.write("")

    except Exception as e:
        print(f"Warning: Could not copy files: {e}")

    lines = [
        "#!/bin/bash",
        "echo $0",
        "echo $SHELL",
        "set -euo pipefail",
        "set +m   # ensure background jobs stay in same job-control context",
        "",
        "# Determine the root directory based on the script location",
        'ROOT=$(dirname "$(readlink -f "$0")")',
        'CARBONARADIR='+str(carbonaradir),
        'export PATH=$PATH:'+str(carbonaradir),
        "",
        "# Optional first argument: FoXS command",
        f'FOXS_CMD="${{1:-{foxs_cmd_default}}}"',
        "",
        "# Python executable for watcher/backmapping.",
        "# Default is the Python used to generate this script; PYTHON_EXE can override it.",
        'if [[ -z "${PYTHON_EXE:-}" ]]; then',
        f'    PYTHON_EXE={shlex.quote(str(python_exe))}',
        "fi",
        'echo "Using Python for watcher/backmapping: $PYTHON_EXE"',
        "",
        "# Directory to clear before running",
        f'CLEAR_DIR="$ROOT/{data_path}/fitdata"',
        "",
        "# ========= NEW: bookkeeping for clean shutdown =========",
        "PIDS=()",
        'WATCHER_PID=""',
        "cleanup() {",
        '    echo',
        '    echo ">>> Stopping Carbonara frontend..."',
        "    # Stop watcher first",
        '    if [[ -n "${WATCHER_PID}" ]]; then',
        '        kill -INT "$WATCHER_PID" 2>/dev/null || true',
        "        sleep 0.5",
        '        kill -TERM "$WATCHER_PID" 2>/dev/null || true',
        "        sleep 0.5",
        '        kill -KILL "$WATCHER_PID" 2>/dev/null || true',
        "    fi",
        "    # Stop all predictor processes explicitly",
        '    if ((${#PIDS[@]})); then',
        '        echo ">>> Stopping predictor processes (${#PIDS[@]})"',
        '        kill -INT  "${PIDS[@]}" 2>/dev/null || true',
        "        sleep 1",
        '        kill -TERM "${PIDS[@]}" 2>/dev/null || true',
        "        sleep 1",
        '        kill -KILL "${PIDS[@]}" 2>/dev/null || true',
        "    fi",
        '    echo ">>> Stopped."',
        "    exit 0",
        "}",
        "trap cleanup SIGINT SIGTERM",
        "# ======================================================",
        "",
        "# Clear the directory",
        'echo "Clearing directory: $CLEAR_DIR"',
        'rm -rf "$CLEAR_DIR"/*',
        'mkdir -p "$CLEAR_DIR"',
        "",
        "### argv[ 1] scattering data file",
        f'ScatterFile="$ROOT/{data_path}/Saxs.dat"',
        "",
        "### argv[ 2] sequence file location",
        f'fileLocs="$ROOT/{data_path}/"',
        "",
        "### argv[ 3] restart tag (use to start from existing prediction)",
        'initialCoordsFile="frompdb"',
        "",
        "### argv[ 4] paired distances file (can be empty)",
        f'pairedPredictions="$ROOT/{data_path}/fixedDistanceConstraints1.dat"' if pairedQ else 'pairedPredictions="False"',
        "",
        "### argv[ 5] fixed sections file (again can be empty)",
        f'fixedsections="$ROOT/{data_path}/varyingSectionSecondary1.dat"',
        "",
        "### argv[ 6] number of structures",
        f"noStructures={int(no_structures)}",
        "",
        "### argv[ 7] request to apply hydrophobic covering WITHIN monomers will be a list of sections on which to apply it -- Currently not used",
        'withinMonomerHydroCover="none"',
        "",
        "### argv[ 8] kmin",
        f"kmin={min_q}",
        "",
        "### argv[ 9] kmax",
        f"kmax={max_q}",
        "",
        "### argv[ 10] kmax Start",
        f"kmaxStart={max_q_start}",
        "",
        "### argv[11] Max number of fitting steps",
        f"maxNoFitSteps={max_fit_steps}",
        "",
        "### argv[12] prediction file - mol[i] in the fitting folder",
        f'predictionFile="$ROOT/{data_path}/fitdata"',
        "",
        "### argv[13] scattering output file",
        f'scatterOut="$ROOT/{data_path}/fitdata"',
        "",
        "### argv[14] mixture list file",
        f'mixtureFile="$ROOT/{data_path}/mixtureFile.dat"',
        "",
        "### argv[15] previous fit string",
        f'prevFitStr="$ROOT/{data_path}/redundant"',
        "",
        "### argv[16] log file location",
        f'logLoc="$ROOT/{data_path}/fitdata"',
        "",
        "### argv[17] last line of the previous fit log",
        'endLinePrevLog="null"',
        "",
        "### argv[18] apply affine rotations",
        'affineTrans="True"' if rotation else 'affineTrans="False"',
        "",
        "### argv[19] use errors in scattering calculation",
        'useErrors="True"',
        "",
        "# ========= NEW: backmapping backend =========",
        f'BACKMAP_BACKEND="{backend}"',
        f'CG2ALL_EXEC="{cg2all_exec}"',
        f'DISULFIDE_CONSTRAINTS_FILE="{disulfide_constraints_file}"',
        "# ===========================================",
        "",
        "# ========= FoXS-based early termination mode =========",
        f'TERMINATE_ON_FOXS="{"True" if terminate_on_foxs else "False"}"',
        f'TERMINATE_FOXS_THRESHOLD="{float(terminate_threshold)}"',
        f'TERMINATE_CONFIRMATION_COUNT="{int(terminate_confirmation_count)}"',
        "# Default mode is maximal exploration: leave all runs to finish.",
        "# If TERMINATE_ON_FOXS=True, watcher may stop individual runs once FoXS is good enough.",
        "# =====================================================",
        "",
        "# ========= NEW: start watcher (background) =========",
        'WATCHER_SCRIPT="$CARBONARADIR/watch_and_backmap.py"',
        'BACKMAP_SCRIPT="$CARBONARADIR/backmap_cli.py"',
        'WATCHER_LOG="$predictionFile/watcher.out"',
        "",
        'if [[ ! -f "$WATCHER_SCRIPT" ]]; then',
        '    echo "ERROR: watcher script not found: $WATCHER_SCRIPT"',
        "    exit 1",
        "fi",
        'if [[ ! -f "$BACKMAP_SCRIPT" ]]; then',
        '    echo "ERROR: backmap script not found: $BACKMAP_SCRIPT"',
        "    exit 1",
        "fi",
        "",
        "# Check that the Python used by the watcher/backmapper has core dependencies.",
        'if ! "$PYTHON_EXE" -c "import sys; print(sys.executable); import Bio; from Bio.PDB import PDBParser, PDBIO" >/dev/null 2>&1; then',
        '    echo "ERROR: chosen PYTHON_EXE cannot import Biopython/Bio.PDB: $PYTHON_EXE"',
        '    echo "       Activate the correct environment, install biopython, or set PYTHON_EXE=/path/to/python."',
        "    exit 1",
        "fi",
        'if [[ "$BACKMAP_BACKEND" == "modeller" ]]; then',
        '    if ! "$PYTHON_EXE" -c "from modeller import environ; env=environ()" >/dev/null 2>&1; then',
        '        echo "ERROR: BACKMAP_BACKEND=modeller but chosen PYTHON_EXE cannot import/configure MODELLER: $PYTHON_EXE"',
        '        echo "       Install/configure MODELLER in this Python, choose --backend cg2all, or set PYTHON_EXE=/path/to/python."',
        "        exit 1",
        "    fi",
        "fi",
        "",
        "WATCHER_ARGS=(",
        '    --watch-dir "$predictionFile"',
        f'    --scenario-root "$ROOT/{data_path}"',
        '    --backmap-script "$BACKMAP_SCRIPT"',
        f'    --max-backmap {int(max_backmap)}',
        '    --no-structures "$noStructures"',
        f'    --defer-backmap-seconds {int(defer_backmap_seconds)}',
        '    --backend "$BACKMAP_BACKEND"',
    ]

    if do_foxs:
        lines.extend([
            '    --do-foxs',
            '    --foxs-py "$FOXS_CMD"',
            '    --saxs "$ScatterFile"',
            '    --max-q "$kmax"',
        ])

    if terminate_on_foxs:
        lines.extend([
            '    --terminate-on-foxs',
            '    --terminate-threshold "$TERMINATE_FOXS_THRESHOLD"',
            '    --terminate-confirmation-count "$TERMINATE_CONFIRMATION_COUNT"',
        ])

    lines.extend([
        ")",
        "",
        'if [[ "$BACKMAP_BACKEND" == "cg2all" ]]; then',
        '    WATCHER_ARGS+=(--cg2all-exec "$CG2ALL_EXEC")',
        "fi",
        "",
        'if [[ -n "$DISULFIDE_CONSTRAINTS_FILE" ]]; then',
        '    WATCHER_ARGS+=(--disulfide-file "$DISULFIDE_CONSTRAINTS_FILE")',
        "fi",
        "",
        '"$PYTHON_EXE" "$WATCHER_SCRIPT" "${WATCHER_ARGS[@]}" > "$WATCHER_LOG" 2>&1 &',
        'WATCHER_PID=$!',
        'echo "Watcher started (PID=$WATCHER_PID)"',
        "# ==================================================",
        "",
        f"for i in {{1..{fit_n_times}}}",
        "do",
        '    echo ""',
        '    echo " >> Run number : $i "',
        '    echo ""',
        '    echo "Max number of fitting steps: " $maxNoFitSteps',
        '    echo ""',
        "",
        '    stdbuf -oL -eL \\',
        '    "$CARBONARADIR/build/bin/predictStructureQvary" \\',
        '        "$ScatterFile" \\',
        '        "$fileLocs" \\',
        '        "$initialCoordsFile" \\',
        '        "$pairedPredictions" \\',
        '        "$fixedsections" \\',
        '        "$noStructures" \\',
        '        "$withinMonomerHydroCover" \\',
        '        "$kmin" \\',
        '        "$kmax" \\',
        '        "$kmaxStart" \\',
        '        "$maxNoFitSteps" \\',
        '        "$predictionFile/mol$i" \\',
        '        "$scatterOut/scatter$i.dat" \\',
        '        "$mixtureFile" \\',
        '        "$prevFitStr" \\',
        '        "$logLoc/fitLog$i.dat" \\',
        '        "$endLinePrevLog" \\',
        '        "$affineTrans" \\',
        '        "$useErrors" \\',
        '        > "$predictionFile/run$i.out" 2> "$predictionFile/run$i.err" &',
        "",
        '    pid=$!',
        '    PIDS+=("$pid")',
        '    echo "$pid" > "$predictionFile/run${i}.pid"',
        '    echo "Predictor run $i started (PID=$pid)"',
        "done",
        "",
        'echo',
        'echo ">>> All runs launched"',
        'echo ">>> Press Ctrl+C to stop everything"',
        'echo',
        "",
        "set +e",
        'for pid in "${PIDS[@]}"; do',
        '    wait "$pid" || true',
        "done",
        "set -e",
        'echo ">>> Predictor processes have finished or been stopped."',
        'echo ">>> Watcher remains active for any final backmapping/FoXS scoring."',
        'echo ">>> Press Ctrl+C to stop the watcher when finished."',
        'wait "$WATCHER_PID" || true',
        "",
    ])

    with open(run_file, "w", newline="\n") as fout:
        fout.write("\n".join(lines))
    os.chmod(run_file, 0o755)
    return run_file


def parse_structure_lengths(filename: str) -> dict:
    with open(filename, "r") as f:
        lines = f.read().splitlines()

    # Filter out lines that look like secondary structure (contain only '-', 'S', 'H')
    structure_lines = [
        line
        for line in lines
        if set(line.strip()).issubset({"-", "S", "H"}) and len(line) > 2
    ]

    # Label them A, B, C, ... and count their lengths
    result = {
        label: len(structure)
        for label, structure in zip(ascii_uppercase, structure_lines)
    }

    return result


def _load_numeric_saxs_table(path: str) -> np.ndarray:
    """Load a Carbonara-style SAXS table and keep only finite numeric rows."""
    arr = np.loadtxt(path)
    arr = np.atleast_2d(np.asarray(arr, dtype=float))
    if arr.shape[1] < 2:
        raise ValueError(f"SAXS file {path!r} must contain at least q and I columns")
    finite = np.isfinite(arr[:, 0]) & np.isfinite(arr[:, 1])
    if arr.shape[1] >= 3:
        finite &= np.isfinite(arr[:, 2])
    arr = arr[finite]
    if arr.shape[0] < 3:
        raise ValueError(f"SAXS file {path!r} has too few finite numeric rows")
    return arr


def _weighted_linear_fit(x: np.ndarray, y: np.ndarray, sigma_y: Optional[np.ndarray] = None):
    """Weighted fit y = intercept + slope*x. Returns dict or None."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if sigma_y is None:
        w = np.ones_like(y)
    else:
        sigma_y = np.asarray(sigma_y, dtype=float)
        good = np.isfinite(sigma_y) & (sigma_y > 0)
        if not np.any(good):
            w = np.ones_like(y)
        else:
            med = float(np.nanmedian(sigma_y[good]))
            sigma_y = np.where(good, sigma_y, med)
            w = 1.0 / np.maximum(sigma_y, 1e-12) ** 2

    X = np.column_stack([np.ones_like(x), x])
    sw = np.sqrt(w)
    try:
        beta, *_ = np.linalg.lstsq(X * sw[:, None], y * sw, rcond=None)
    except np.linalg.LinAlgError:
        return None

    intercept, slope = float(beta[0]), float(beta[1])
    yhat = intercept + slope * x
    resid = y - yhat
    ybar = float(np.average(y, weights=w))
    ss_res = float(np.sum(w * resid ** 2))
    ss_tot = float(np.sum(w * (y - ybar) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 1.0

    if sigma_y is None:
        scale = 1.4826 * np.median(np.abs(resid - np.median(resid)))
        if not np.isfinite(scale) or scale <= 1e-12:
            scale = float(np.sqrt(np.mean(resid ** 2))) if len(resid) else 1.0
        z = resid / max(scale, 1e-12)
    else:
        z = resid / np.maximum(sigma_y, 1e-12)

    return {
        "intercept": intercept,
        "slope": slope,
        "r2": float(r2),
        "residuals": resid,
        "max_abs_z": float(np.max(np.abs(z))) if len(z) else 0.0,
        "yhat": yhat,
    }


def _guinier_fit_for_window(q: np.ndarray, intensity: np.ndarray,
                            sigma_i: Optional[np.ndarray], start: int, end: int):
    """Fit ln(I) versus q^2 on q[start:end]."""
    qwin = q[start:end]
    iwin = intensity[start:end]
    if np.any(qwin <= 0) or np.any(iwin <= 0):
        return None
    x = qwin ** 2
    y = np.log(iwin)
    sigma_y = None
    if sigma_i is not None:
        swin = sigma_i[start:end]
        good = np.isfinite(swin) & (swin > 0) & np.isfinite(iwin) & (iwin > 0)
        if np.any(good):
            sigma_y = np.where(good, swin / iwin, np.nan)

    fit = _weighted_linear_fit(x, y, sigma_y=sigma_y)
    if fit is None:
        return None
    if fit["slope"] >= 0:
        return None

    rg = float(np.sqrt(-3.0 * fit["slope"]))  # q in A^-1 gives Rg in A
    fit["rg"] = rg
    fit["qrg_max"] = float(np.max(qwin) * rg)
    fit["start"] = int(start)
    fit["end"] = int(end)
    fit["n_points"] = int(end - start)
    return fit



def _robust_mad_scale(values: np.ndarray, floor: float = 1e-12) -> float:
    """Return a robust residual scale using MAD, with an RMS fallback."""
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return 1.0
    med = float(np.median(values))
    mad = float(np.median(np.abs(values - med)))
    scale = 1.4826 * mad
    if not np.isfinite(scale) or scale <= floor:
        scale = float(np.sqrt(np.mean((values - med) ** 2)))
    if not np.isfinite(scale) or scale <= floor:
        scale = floor
    return scale


def _score_internal_guinier_fit(fit: dict, q: np.ndarray, intensity: np.ndarray,
                                sigma_i: Optional[np.ndarray], start: int, end: int,
                                r2_min: float, max_abs_z: float,
                                qrg_max: float, qrg_min: float) -> Optional[dict]:
    """Attach AutoRg-style diagnostic and quality fields to a Guinier fit."""
    n = int(end - start)
    qwin = q[start:end]
    iwin = intensity[start:end]
    resid = np.asarray(fit["residuals"], dtype=float)

    robust_scale = _robust_mad_scale(resid)
    robust_z = resid / robust_scale
    max_robust_z = float(np.max(np.abs(robust_z))) if robust_z.size else 0.0
    rms_resid = float(np.sqrt(np.mean(resid ** 2))) if resid.size else float("inf")

    sigma_chi2 = None
    if sigma_i is not None:
        swin = sigma_i[start:end]
        sigma_ln = np.where((swin > 0) & np.isfinite(swin) & (iwin > 0), swin / iwin, np.nan)
        good = np.isfinite(sigma_ln) & (sigma_ln > 0)
        if np.count_nonzero(good) >= max(3, n // 2):
            z = resid[good] / sigma_ln[good]
            sigma_chi2 = float(np.mean(z ** 2))

    qrg_lo = float(np.min(qwin) * fit["rg"])
    qrg_hi = float(np.max(qwin) * fit["rg"])
    if not np.isfinite(qrg_hi) or qrg_hi <= 0:
        return None
    if qrg_hi > float(qrg_max):
        return None
    # Very small qRg intervals can appear linear simply because there is too
    # little leverage to measure an Rg.  We keep this as a soft criterion in
    # the score but reject extremely low-leverage intervals.
    if qrg_hi < max(0.2, 0.5 * float(qrg_min)):
        return None

    r2 = float(fit["r2"])
    # Treat R^2 and residual-z as the main linearity tests.  Experimental
    # sigma values in merged SAXS files can be optimistic, so sigma_chi2 is
    # diagnostic rather than an absolute rejection criterion.
    if r2 < float(r2_min):
        return None
    if max_robust_z > float(max_abs_z):
        return None

    r2_score = min(1.0, max(0.0, (r2 - r2_min) / max(1e-12, 1.0 - r2_min)))
    z_score = min(1.0, max(0.0, 1.0 - max_robust_z / max(1e-12, max_abs_z)))
    qrg_score = min(1.0, qrg_hi / max(1e-12, qrg_min))
    length_score = min(1.0, n / max(1.0, 2.0 * 8.0))
    # Slightly prefer intervals that extend farther in Guinier space but remain
    # inside the accepted qRg limit.
    range_score = 0.5 + 0.5 * min(1.0, qrg_hi / max(1e-12, qrg_max))
    quality = float((0.45 * r2_score + 0.35 * z_score + 0.20 * qrg_score) * length_score * range_score)

    out = dict(fit)
    out.update({
        "start": int(start),
        "end": int(end),
        "first_point": int(start + 1),   # 1-indexed, AUTORG-style
        "last_point": int(end),          # 1-indexed inclusive
        "n_points": int(n),
        "qrg_min": qrg_lo,
        "qrg_max": qrg_hi,
        "rms_log_residual": rms_resid,
        "robust_residual_scale": float(robust_scale),
        "max_abs_robust_z": max_robust_z,
        "sigma_chi2": sigma_chi2,
        "quality": quality,
    })
    return out


def _select_internal_guinier_interval(
    arr: np.ndarray,
    min_points: int = 8,
    window_points: int = 80,
    qrg_max: float = 1.3,
    r2_min: float = 0.98,
    max_abs_z: float = 4.0,
    max_trim_fraction: float = 0.25,
) -> dict:
    """
    Self-contained AutoRg-style Guinier interval selection.

    The protocol mirrors the standard manual/AUTORG idea without requiring an
    ATSAS install:
      1. work in Guinier coordinates, ln(I) versus q^2;
      2. scan many contiguous low-q intervals longer than a minimum length;
      3. keep intervals with negative slope, qRg inside the Guinier range,
         and good linear residuals;
      4. require the fitted Rg to be stable across neighbouring intervals;
      5. choose the closest-to-origin stable interval and trim to its first
         point.
    """
    q = np.asarray(arr[:, 0], dtype=float)
    intensity = np.asarray(arr[:, 1], dtype=float)
    sigma_i = np.asarray(arr[:, 2], dtype=float) if arr.shape[1] >= 3 else None
    n = int(arr.shape[0])
    min_points = int(min_points)
    max_width = int(max(window_points, min_points))
    # Need enough candidate width to reach qRg~1 for large molecules.  Cap for
    # speed, but still search a generous low-q window.
    max_width = min(n, max(max_width, min_points + 20))
    max_start = min(n - min_points, int(np.floor(float(max_trim_fraction) * n)))
    if max_start < 0:
        raise RuntimeError("not enough points for Guinier interval selection")

    qrg_min = min(1.0, max(0.45, 0.55 * float(qrg_max)))
    by_start = {}
    all_fits = []

    for start in range(max_start + 1):
        end_max = min(n, start + max_width)
        for end in range(start + min_points, end_max + 1):
            fit = _guinier_fit_for_window(q, intensity, sigma_i, start, end)
            if fit is None:
                continue
            scored = _score_internal_guinier_fit(
                fit, q, intensity, sigma_i, start, end,
                r2_min=r2_min,
                max_abs_z=max_abs_z,
                qrg_max=qrg_max,
                qrg_min=qrg_min,
            )
            if scored is None:
                continue
            by_start.setdefault(start, []).append(scored)
            all_fits.append(scored)

    if not all_fits:
        raise RuntimeError(
            "No acceptable Guinier interval found by internal AutoRg-style scan. "
            f"Tried starts 1-{max_start + 1}, min_points={min_points}, "
            f"max_width={max_width}, qRg_max={qrg_max}, r2_min={r2_min}, "
            f"max_abs_z={max_abs_z}."
        )

    accepted_starts = []
    for start, fits in sorted(by_start.items()):
        if len(fits) < 2:
            continue
        rgs = np.asarray([f["rg"] for f in fits], dtype=float)
        rg_med = float(np.median(rgs))
        rg_cv = float(np.std(rgs) / max(abs(rg_med), 1e-12))
        best = max(fits, key=lambda f: (f["quality"], f["n_points"], -f["rms_log_residual"]))
        # Stability is deliberately a little permissive: real SAXS files can be
        # noisy/merged, and this trim is only choosing the first reliable point.
        if rg_cv <= 0.20 and best["quality"] >= 0.12:
            best = dict(best)
            best["rg_stability_cv_for_same_start"] = rg_cv
            best["candidate_intervals_same_start"] = int(len(fits))
            accepted_starts.append(best)

    if accepted_starts:
        selected = sorted(
            accepted_starts,
            key=lambda f: (f["start"], -f["quality"], -f["n_points"]),
        )[0]
        selection_status = "accepted_closest_to_origin_stable_interval"
    else:
        # Fallback: choose the best global interval, but report that the stability
        # criterion was not fully satisfied.
        selected = max(all_fits, key=lambda f: (f["quality"], -f["start"], f["n_points"]))
        selected = dict(selected)
        selected["rg_stability_cv_for_same_start"] = None
        selected["candidate_intervals_same_start"] = int(len(by_start.get(selected["start"], [])))
        selection_status = "best_linear_interval_no_stable_start_cluster"

    selected["selection_status"] = selection_status
    selected["candidate_interval_count"] = int(len(all_fits))
    selected["searched_start_points"] = int(max_start + 1)
    selected["max_interval_width_points"] = int(max_width)
    selected["qrg_target_min_used"] = float(qrg_min)
    return selected


def guinier_trim_saxs_file_inplace(
    saxs_path: str,
    min_points: int = 8,
    window_points: int = 80,
    qrg_max: float = 1.3,
    r2_min: float = 0.98,
    max_abs_z: float = 4.0,
    max_trim_fraction: float = 0.25,
    autorg_cmd: str = None,
    require_success: bool = False,
) -> dict:
    """
    Remove the low-q prefix before the first internally selected Guinier point.

    This is self-contained and does not require ATSAS/AUTORG.  It follows the
    standard Guinier protocol: identify a linear interval in ln(I) versus q^2,
    check qRg, then discard points before that interval so Carbonara reads only
    the usable low-q region onward.
    """
    saxs_path = str(saxs_path)
    backup_path = saxs_path + ".pre_guinier_trim"

    arr0 = _load_numeric_saxs_table(saxs_path)
    try:
        shutil.copy2(saxs_path, backup_path)
    except OSError as exc:
        print(f"WARNING: could not write Guinier pre-trim backup {backup_path}: {exc}")
        backup_path = None

    positive = (arr0[:, 0] > 0) & (arr0[:, 1] > 0)
    dropped_nonpositive = int(np.count_nonzero(~positive))
    arr = arr0[positive]

    if arr.shape[0] < max(3, int(min_points)):
        report = {
            "enabled": True,
            "method": "internal_autorg_like",
            "status": "not_enough_positive_points_no_trim_applied",
            "saxs_path": saxs_path,
            "backup_path": backup_path,
            "points_before": int(arr0.shape[0]),
            "points_after": int(arr.shape[0]),
            "trimmed_lowq_points": 0,
            "dropped_nonpositive_or_invalid_rows": dropped_nonpositive,
        }
        _write_guinier_trim_report(saxs_path, report)
        msg = "WARNING: Guinier trim skipped: not enough positive SAXS points"
        if require_success:
            raise RuntimeError(msg)
        print(msg)
        return report

    order = np.argsort(arr[:, 0])
    arr = arr[order]
    n = int(arr.shape[0])

    try:
        selected = _select_internal_guinier_interval(
            arr,
            min_points=min_points,
            window_points=window_points,
            qrg_max=qrg_max,
            r2_min=r2_min,
            max_abs_z=max_abs_z,
            max_trim_fraction=max_trim_fraction,
        )
    except Exception as exc:
        report = {
            "enabled": True,
            "method": "internal_autorg_like",
            "status": "no_acceptable_guinier_window_found_no_trim_applied",
            "saxs_path": saxs_path,
            "backup_path": backup_path,
            "points_before": int(arr0.shape[0]),
            "points_after": int(n),
            "trimmed_lowq_points": 0,
            "dropped_nonpositive_or_invalid_rows": dropped_nonpositive,
            "error": str(exc),
            "settings": {
                "min_points": int(min_points),
                "window_points": int(window_points),
                "qrg_max": float(qrg_max),
                "r2_min": float(r2_min),
                "max_abs_z": float(max_abs_z),
                "max_trim_fraction": float(max_trim_fraction),
                "autorg_cmd_ignored": autorg_cmd,
            },
        }
        _write_guinier_trim_report(saxs_path, report)
        msg = "WARNING: Guinier trim found no acceptable interval; leaving SAXS data unchanged. See guinier_trim_report.json"
        if require_success:
            raise RuntimeError(msg + "\n" + str(exc))
        print(msg)
        return report

    trim_start = int(selected["start"])
    max_trim_points = int(np.floor(float(max_trim_fraction) * n))
    if trim_start > max_trim_points:
        report = {
            "enabled": True,
            "method": "internal_autorg_like",
            "status": "selected_interval_exceeds_max_trim_fraction_no_trim_applied",
            "saxs_path": saxs_path,
            "backup_path": backup_path,
            "points_before": int(arr0.shape[0]),
            "points_after": int(n),
            "trimmed_lowq_points": 0,
            "requested_trim_lowq_points": int(trim_start),
            "dropped_nonpositive_or_invalid_rows": dropped_nonpositive,
            "selected_interval": selected,
            "settings": {
                "min_points": int(min_points),
                "window_points": int(window_points),
                "qrg_max": float(qrg_max),
                "r2_min": float(r2_min),
                "max_abs_z": float(max_abs_z),
                "max_trim_fraction": float(max_trim_fraction),
                "autorg_cmd_ignored": autorg_cmd,
            },
        }
        _write_guinier_trim_report(saxs_path, report)
        msg = (
            "WARNING: selected Guinier interval starts after the allowed trim cap "
            f"--guinier_trim_max_fraction={max_trim_fraction}; leaving SAXS data unchanged"
        )
        if require_success:
            raise RuntimeError(msg)
        print(msg)
        return report

    trimmed = arr[trim_start:]
    if trim_start > 0 or dropped_nonpositive > 0:
        np.savetxt(saxs_path, trimmed, delimiter=" ", fmt="%.10g")

    report = {
        "enabled": True,
        "method": "internal_autorg_like",
        "status": "trim_applied" if (trim_start > 0 or dropped_nonpositive > 0) else "already_good_no_trim_needed",
        "saxs_path": saxs_path,
        "backup_path": backup_path,
        "points_before": int(arr0.shape[0]),
        "points_after": int(trimmed.shape[0]),
        "trimmed_lowq_points": int(trim_start),
        "dropped_nonpositive_or_invalid_rows": dropped_nonpositive,
        "first_q_before": float(arr[0, 0]),
        "first_q_after": float(trimmed[0, 0]),
        "selected_first_point_1_indexed": int(selected["first_point"]),
        "selected_last_point_1_indexed_before_trim": int(selected["last_point"]),
        "selected_last_point_1_indexed_after_trim": int(selected["last_point"] - trim_start),
        "selected_rg": float(selected["rg"]),
        "selected_i0": float(np.exp(selected["intercept"])),
        "selected_r2": float(selected["r2"]),
        "selected_quality": float(selected["quality"]),
        "selected_qrg_min": float(selected["qrg_min"]),
        "selected_qrg_max": float(selected["qrg_max"]),
        "selected_max_abs_robust_z": float(selected["max_abs_robust_z"]),
        "selected_sigma_chi2": selected.get("sigma_chi2"),
        "selection_status": selected.get("selection_status"),
        "candidate_interval_count": int(selected.get("candidate_interval_count", 0)),
        "searched_start_points": int(selected.get("searched_start_points", 0)),
        "settings": {
            "min_points": int(min_points),
            "window_points": int(window_points),
            "qrg_max": float(qrg_max),
            "r2_min": float(r2_min),
            "max_abs_z": float(max_abs_z),
            "max_trim_fraction": float(max_trim_fraction),
            "autorg_cmd_ignored": autorg_cmd,
        },
    }
    _write_guinier_trim_report(saxs_path, report)
    print(
        "Guinier trim: "
        f"status={report['status']}; removed {report['trimmed_lowq_points']} leading points; "
        f"selected points {report['selected_first_point_1_indexed']}-{report['selected_last_point_1_indexed_before_trim']}; "
        f"first q {report['first_q_before']:.5g} -> {report['first_q_after']:.5g}; "
        f"Rg={report['selected_rg']:.5g}; R2={report['selected_r2']:.5g}; "
        f"qRg={report['selected_qrg_min']:.3g}-{report['selected_qrg_max']:.3g}; "
        f"quality={report['selected_quality']:.3g}"
    )
    return report


def _write_guinier_trim_report(saxs_path: str, report: dict) -> None:
    """Write a small JSON report next to Saxs.dat."""
    report_path = os.path.join(os.path.dirname(str(saxs_path)), "guinier_trim_report.json")
    try:
        import json
        with open(report_path, "w", encoding="utf-8") as fh:
            json.dump(report, fh, indent=2, sort_keys=True)
    except OSError as exc:
        print(f"WARNING: could not write Guinier trim report {report_path}: {exc}")


def main():
    parser = argparse.ArgumentParser(description="Setup Carbonara processing pipeline")
    parser.add_argument("-p", "--pdb", required=True, help="Path to input PDB file")
    parser.add_argument("-s", "--saxs", required=True, help="Path to input SAXS data file")
    parser.add_argument("-n", "--name", required=True, help="Name for this protein/refinement")
    parser.add_argument(
        "-f",
        "--pae",
        required=False,
        help="PAE file for this protein, used to specify flexibility",
    )
    parser.add_argument(
        "-d", "--dir", default=os.getcwd(), help="Base directory (default: current directory)"
    )

    # Additional parameters for write_runme
    parser.add_argument(
        "--fit_n_times", type=int, default=20, help="Number of times to run the fit (default: 20)"
    )
    parser.add_argument("--min_q", type=float, default=0.01, help="Minimum q-value (default: 0.01)")
    parser.add_argument("--max_q", type=float, default=0.2, help="Maximum q-value (default: 0.2)")
    parser.add_argument(
        "--max_q_start",
        type=float,
        default=0.2,
        help="Maximum q-value to start fitting to (default: 0.2)",
    )
    parser.add_argument(
        "--max_fit_steps",
        type=int,
        default=10000,
        help="Maximum number of fitting steps (default: 10000)",
    )
    parser.add_argument("--pairedQ", action="store_true", help="Use paired predictions")
    parser.add_argument("--rotation", action="store_true", help="Apply affine rotations")
    parser.add_argument(
        "--alphaFoldFlex",
        action="store_true",
        help="Use an alphaFold pae file to specify the flexibility of the molecule",
    )
    parser.add_argument('--pae_flex_threshold', type=float, default=16.0,
                    help="Absolute Å threshold if --pae_flex_mode=absolute (default: 16).")

    # Ensemble / mixture controls (replicating setup)
    parser.add_argument(
        "--mixture_n",
        type=int,
        default=1,
        help="Number of structures/species in mixture/ensemble (default: 1)",
    )
    parser.add_argument(
        "--max_mixture_combos",
        type=int,
        default=30,
        help="Max number of mixture combinations to write when mixture_n>1 (default: 30)",
    )
    parser.add_argument(
        "--mixture_step",
        type=float,
        default=0.1,
        help="Step for n=2 mixture grid (default: 0.1)",
    )
    parser.add_argument(
        "--mixture_dirichlet_alpha",
        type=float,
        default=1.0,
        help="Dirichlet alpha for n>=3 mixture sampling (default: 1.0)",
    )
    parser.add_argument(
        "--mixture_snap",
        type=float,
        default=0.05,
        help="Snap sampled mixtures to multiples of this (0 disables). Default: 0.05",
    )
    parser.add_argument("--max_backmap", type=int, default=3,
                    help="Maximum number of concurrent backmapping jobs (default: 3)")
    parser.add_argument("--defer_backmap_seconds", type=int, default=600,
                    help="Ignore early structures for this many seconds before backmapping starts (default: 600)")
    parser.add_argument("--no_foxs", action="store_true",
                    help="Do not enable FoXS in the generated run script")
    parser.add_argument("--guinier_trim", "--guinier-trim", dest="guinier_trim", action="store_true", default=True,
                    help="Apply setup-time Guinier low-q trimming before Carbonara reads Saxs.dat (default: on)")
    parser.add_argument("--no_guinier_trim", "--no-guinier-trim", dest="guinier_trim", action="store_false",
                    help="Disable setup-time Guinier low-q trimming and keep Saxs.dat unchanged after cdt.write_saxs")
    parser.add_argument("--guinier_trim_min_points", type=int, default=8,
                    help="Minimum number of points in the accepted Guinier window (default: 8)")
    parser.add_argument("--guinier_trim_window_points", type=int, default=80,
                    help="Maximum interval width, in points, for the internal AutoRg-style Guinier scan (default: 80)")
    parser.add_argument("--guinier_trim_qrg_max", type=float, default=1.3,
                    help="Maximum q*Rg allowed in the accepted Guinier window (default: 1.3)")
    parser.add_argument("--guinier_trim_r2_min", type=float, default=0.98,
                    help="Minimum R^2 for internally selected Guinier-linearity intervals (default: 0.98)")
    parser.add_argument("--guinier_trim_max_abs_z", type=float, default=4.0,
                    help="Maximum robust standardized residual in internally selected Guinier intervals (default: 4.0)")
    parser.add_argument("--guinier_trim_max_fraction", type=float, default=0.25,
                    help="Safety cap on fraction of leading points that Guinier trim may remove (default: 0.25)")
    parser.add_argument("--guinier_trim_autorg_cmd", default=None,
                    help="Deprecated and ignored: Guinier trim is now self-contained and does not require ATSAS/autorg")
    parser.add_argument("--guinier_trim_require_success", action="store_true",
                    help="Fail setup if the internal Guinier trim cannot select an acceptable interval")
    parser.add_argument("--backend", choices=["modeller", "cg2all"], default="modeller",
                    help="Backmapping backend for the generated RunMe script")
    parser.add_argument("--disulfide_constraints_file", default="",
                    help="Optional constraint file to treat as disulfides during backmapping")
    parser.add_argument(
        "--no-disulfide-linker-check", "--no_disulfide_linker_check",
        dest="no_disulfide_linker_check",
        action="store_true",
        help=(
            "Disable only the setup-time disulfide-aware automatic linker filter. "
            "This does not disable --disulfide_constraints_file or backmapping disulfide enforcement."
        ),
    )
    parser.add_argument("--fix_missing_residues", "--fix-missing-residues", dest="fix_missing_residues",
                    action="store_true", default=True,
                    help=(
                        "Repair internal missing-residue gaps before Carbonara reads the structure "
                        "(default: on; Carbonara-friendly bridge method for PDB inputs)."
                    ))
    parser.add_argument("--no_fix_missing_residues", "--no-fix-missing-residues", dest="fix_missing_residues",
                    action="store_false",
                    help="Disable setup-time repair of internal missing-residue gaps")
    parser.add_argument("--fix_missing_residue_name", default="GLY",
                    help="3-letter residue name used for residue-number gaps when no SEQRES identity is available (default: GLY)")
    parser.add_argument("--fix_missing_residue_max_gap", type=int, default=80,
                    help="Do not auto-build any individual missing-residue gap longer than this many residues (default: 80)")
    parser.add_argument("--fix_missing_residue_method", choices=["bridge", "pdbfixer"], default="bridge",
                    help=("Missing-residue repair method. 'bridge' (default) writes a run-local, "
                          "Carbonara-friendly GLY backbone bridge before any CDT reading; "
                          "'pdbfixer' uses PDBFixer addMissingAtoms."))
    parser.add_argument("--split_long_linkers", "--split-long-linkers", dest="split_long_linkers",
                    action="store_true", default=True,
                    help="Break flexible linker sections longer than --max_linker_len before writing the fingerprint (default: on)")
    parser.add_argument("--no_split_long_linkers", "--no-split-long-linkers", dest="split_long_linkers",
                    action="store_false",
                    help="Disable automatic splitting of very long flexible linker sections")
    parser.add_argument("--max_linker_len", "--max-linker-len", type=int, default=25,
                    help="Split flexible '-' sections only when their length is greater than this value (default: 25)")
    parser.add_argument("--long_linker_fake_helix_len", "--long-linker-fake-helix-len", type=int, default=3,
                    help="Length of each artificial H separator used to break a long linker (default: 3)")
    parser.add_argument("--foxs_cmd_default", default="pyfoxs",
                    help="Default FoXS command for the generated RunMe script; user can still override as first shell arg")
    parser.add_argument("--python_exe", default=None,
                    help="Python executable for watcher/backmapping in the generated RunMe script (default: the Python running setup). Runtime override: PYTHON_EXE=/path/to/python ./RunMe_name.sh")
    parser.add_argument("--cg2all_exec",default=None,help="Override cg2all executable command string (advanced users only)")
    parser.add_argument("--terminate-on-foxs", action="store_true",
                    help="Opt-in batch mode: terminate individual predictStructureQvary runs once FoXS chi^2 is good enough")
    parser.add_argument("--terminate-threshold", type=float, default=2.5,
                    help="FoXS chi^2 threshold for --terminate-on-foxs (default: 2.5)")
    parser.add_argument("--terminate-confirmation-count", type=int, default=1,
                    help="Number of qualifying FoXS scores required before stopping a run (default: 1)")
    args = parser.parse_args()

    if args.terminate_confirmation_count < 1:
        parser.error("--terminate-confirmation-count must be >= 1")
    if args.terminate_threshold <= 0:
        parser.error("--terminate-threshold must be > 0")
    if args.terminate_on_foxs and args.no_foxs:
        parser.error("--terminate-on-foxs requires FoXS; remove --no_foxs")
    
    if args.guinier_trim:
        if args.guinier_trim_min_points < 4:
            parser.error("--guinier_trim_min_points must be >= 4")
        if args.guinier_trim_window_points < args.guinier_trim_min_points:
            parser.error("--guinier_trim_window_points must be >= --guinier_trim_min_points")
        if args.guinier_trim_qrg_max <= 0:
            parser.error("--guinier_trim_qrg_max must be > 0")
        if not (0.0 <= args.guinier_trim_r2_min <= 1.0):
            parser.error("--guinier_trim_r2_min must be between 0 and 1")
        if args.guinier_trim_max_abs_z <= 0:
            parser.error("--guinier_trim_max_abs_z must be > 0")
        if not (0.0 <= args.guinier_trim_max_fraction <= 1.0):
            parser.error("--guinier_trim_max_fraction must be between 0 and 1")

    if args.fix_missing_residues:
        if args.fix_missing_residue_max_gap < 1:
            parser.error("--fix_missing_residue_max_gap must be >= 1")
        args.fix_missing_residue_name = str(args.fix_missing_residue_name).upper()

    if args.split_long_linkers:
        if args.max_linker_len < 4:
            parser.error("--max_linker_len must be >= 4")
        if args.long_linker_fake_helix_len < 1:
            parser.error("--long_linker_fake_helix_len must be >= 1")

    print("DIR: "+str(args.dir))
    print(os.getcwd())
    
    print("new")
    def detect_cg2all_exec(user_value):
        if user_value:
            return user_value

        if os.path.exists("/content/bin/micromamba"):
            return "/content/bin/micromamba run -p /root/micromamba/envs/cg2all convert_cg2all"

        if os.path.exists("./bin/micromamba"):
            return "./bin/micromamba run -p /root/micromamba/envs/cg2all convert_cg2all"

        return "convert_cg2all"

    if args.backend == "cg2all":
        args.cg2all_exec = detect_cg2all_exec(args.cg2all_exec)
        print(f"Using CG2ALL executable: {args.cg2all_exec}")
    else:
        args.cg2all_exec = ""
    try:
        if args.mixture_n < 1:
            raise ValueError("--mixture_n must be >= 1")

        # Setup master directory
        fit_master_dir = cdt.setup_fit_master_dir(root_dir=args.dir, fit_master_name="carbonara_runs")

        # Setup refinement directory
        refine_dir = cdt.setup_refinement_dir(args.name, fit_master_dir)
        print(f"Created directory structure in: {refine_dir}")

        # Process PDB/CIF and extract structure information.
        # First perform any repairs that must happen before Carbonara sees the structure.
        # In particular, build internal missing-residue gaps as a run-local physical
        # structure before Carbonara reads it, rather than patching section labels afterwards.
        preflight_structure_extension_matches_content(args.pdb)
        preflight_mmcif_atom_site_table(args.pdb)

        pdb_for_carbonara = prepare_structure_before_carbonara(
            args.pdb,
            refine_dir,
            enabled=args.fix_missing_residues,
            residue_name=args.fix_missing_residue_name,
            max_gap=args.fix_missing_residue_max_gap,
            method=args.fix_missing_residue_method,
        )

        # From this point onward the repaired copy (when one was needed) is the
        # canonical structure for the run.  The user's original input file is
        # never renamed, moved, overwritten, or otherwise modified.
        #
        # Do NOT rely on assigning a module global to update a Jupyter variable:
        # setup is commonly run as an imported module or subprocess, in which case
        # its globals live in a different namespace.  Persist the canonical path
        # beside fingerPrint1.dat so CarbonaraDataTools can recover it later even
        # when a notebook still holds the original pdb_name.
        canonical_structure_file = os.path.join(refine_dir, "carbonara_structure_path.txt")
        canonical_structure_path = os.path.abspath(str(pdb_for_carbonara))
        with open(canonical_structure_file, "w") as fh:
            fh.write(canonical_structure_path + "\n")
        print(f"Canonical Carbonara structure: {canonical_structure_path}")

        # Retain the module-level assignment as a convenience for direct/%run use,
        # but correctness no longer depends on it.
        global pdb_name
        pdb_name = canonical_structure_path
        args.pdb = pdb_name

        coords_chains, sequence_chains, secondary_structure_chains, missing_residues_chains = (
            cdt.pull_structure_from_pdb(pdb_name)
        )

        print("The number of chains is ", len(coords_chains))
        new_coords_chains = []
        new_sequence_chains = []
        new_secondary_structure_chains = []

        for i in range(len(coords_chains)):
            breaking_indices = cdt.missing_ca_check(coords_chains[i], threshold_dist_Å=7)
            if len(breaking_indices) > 0:
                print("Warning: Missing segments of chain found: ", len(breaking_indices), breaking_indices)

            # Always split — even if indices is empty, returns [full array]
            split_coords = np.array_split(coords_chains[i], breaking_indices)
            split_seq = np.array_split(sequence_chains[i], breaking_indices)
            split_ss = np.array_split(secondary_structure_chains[i], breaking_indices)

            new_coords_chains.extend(split_coords)
            new_sequence_chains.extend(split_seq)
            new_secondary_structure_chains.extend(split_ss)

        coords_chains = np.array(new_coords_chains, dtype=object)
        sequence_chains = np.array(new_sequence_chains, dtype=object)
        secondary_structure_chains = np.array(new_secondary_structure_chains, dtype=object)

        # Break very long flexible regions before any Carbonara section numbering is
        # created.  This deliberately modifies only the secondary-structure labels;
        # sequence and coordinates retain exactly the same length/indexing.
        secondary_structure_chains, _long_linker_report = split_long_linkers_in_secondary_structure(
            secondary_structure_chains,
            enabled=args.split_long_linkers,
            max_linker_len=args.max_linker_len,
            fake_helix_len=args.long_linker_fake_helix_len,
            report_path=os.path.join(refine_dir, "long_linker_break_report.json"),
        )

        # collapse coordinates file into one chain
        coords_full = None
        for i, coords in enumerate(coords_chains):
            coords_full = coords if i == 0 else np.concatenate((coords_full, coords), axis=0)

        # write this to file
        coords_files = []
        coords_files.append(cdt.write_coordinates_file(coords_full, working_path=refine_dir, carb_index=1))

        # Write fingerprint file
        number_of_chains = len(coords_chains)
        fingerprint_file = cdt.write_fingerprint_file(
            number_chains=number_of_chains,
            sequence=sequence_chains,
            secondary_structure=secondary_structure_chains,
            working_path=refine_dir,
        )

        # Copy SAXS file to Saxs.dat (this is the file that Carbonara will use).
        # Guinier trimming is applied here by default, before Carbonara reads Saxs.dat.
        # It is self-contained: no ATSAS/autorg installation is required.
        # Use --no_guinier_trim / --no-guinier-trim to keep the normalized SAXS file unchanged.
        saxs_path = cdt.write_saxs(args.saxs, refine_dir)
        if args.guinier_trim:
            guinier_trim_saxs_file_inplace(
                saxs_path,
                min_points=args.guinier_trim_min_points,
                window_points=args.guinier_trim_window_points,
                qrg_max=args.guinier_trim_qrg_max,
                r2_min=args.guinier_trim_r2_min,
                max_abs_z=args.guinier_trim_max_abs_z,
                max_trim_fraction=args.guinier_trim_max_fraction,
                autorg_cmd=args.guinier_trim_autorg_cmd,
                require_success=args.guinier_trim_require_success,
            )

        # check the requested min q is not less than the minimum value in the saxs file
        qmin = np.max([np.loadtxt(refine_dir + "/Saxs.dat")[0][0], args.min_q])

        # set alphaFold flexibility
        varying_linker_chains = []
        if args.alphaFoldFlex:
            # use pae scores to specify flexibility
            varying_linker_chains = cdt.getFlexibility(
                args.pae, fingerprint_file,
                abs_thr=args.pae_flex_threshold
            )
        else:
            # auto select flexible linker chains that dont break inter-beta sheets
            # Default: pass args.pdb so CarbonaraDataTools can protect disulfide-linked
            # regions when selecting flexible linkers.  The flag below disables only
            # this setup-time linker filter; it does not affect disulfide constraints
            # passed to the backmapper/watcher.
            linker_check_pdb = None if args.no_disulfide_linker_check else pdb_name
            if args.no_disulfide_linker_check:
                print(
                    "WARNING: setup-time disulfide linker check disabled; "
                    "backmapping disulfide constraints are unchanged."
                )
            for coord_file in coords_files:
                varying_linker_chains.append(
                    cdt.auto_select_varying_linker(coord_file, fingerprint_file, linker_check_pdb)
                )

        # write flexible linkers to files (varysections1.dat, varysections2.dat, etc [each file is for a different chain])
        varying_section_files = []
        for varying_linkers in varying_linker_chains:
            varying_section_files.append(cdt.write_varysections_file(varying_linkers, refine_dir))

        # Final guard for varying sections.  This must run whether the file is
        # empty or not: section IDs passed to the C++ code must be real internal
        # linker sections of length >= 4.
        filepath = refine_dir + "/fingerPrint1.dat"
        vs_path = refine_dir + "/varyingSectionSecondary1.dat"

        try:
            target_segments = np.loadtxt(vs_path, dtype=int, ndmin=1)
        except ValueError:
            # happens if the file is empty / whitespace
            target_segments = np.array([], dtype=int)

        target_segments = np.atleast_1d(np.asarray(target_segments, dtype=int))
        if target_segments.size == 0:
            open(vs_path, "w").close()
        else:
            if hasattr(cdt, "validate_varying_linker_indices"):
                filtered_segments = cdt.validate_varying_linker_indices(
                    target_segments, filepath, min_length=4, exclude_terminal=True,
                    label="setup final varyingSectionSecondary filter"
                )
            else:
                filtered_segments = cdt.get_segment_lengths_from_file(
                    filepath, target_segments, min_length=4
                )
            filtered_segments = np.asarray(filtered_segments, dtype=int)
            np.savetxt(vs_path, filtered_segments, fmt="%i")

        
         # store chain lengths
        chain_lengths = parse_structure_lengths(refine_dir + "/fingerPrint1.dat")
        with open(refine_dir + "/chainLengths.dat", "wb") as f:
            pickle.dump(chain_lengths, f)
            
        # Mixture file
        # - keep original behavior for mixture_n=1 (whatever cdt.write_mixture_file does)
        # - for mixture_n>1 write a sparse mixture list (<= max_mixture_combos)
        if args.mixture_n <= 1:
            mixture_file = cdt.write_mixture_file(working_path=refine_dir)
        else:
            mixture_file = write_sparse_mixture_file(
                os.path.join(refine_dir, "mixtureFile.dat"),
                n=args.mixture_n,
                max_combos=args.max_mixture_combos,
                step_2=args.mixture_step,
                dirichlet_alpha=args.mixture_dirichlet_alpha,
                snap=args.mixture_snap,
                seed=0,
            )

        # Replicate required numbered files for ensemble refinement (mixture_n > 1)
        replicate_numbered_files(refine_dir, args.mixture_n)

        # Write the RunMe_<name>.sh script
        run_script = write_runme(
            working_path=refine_dir,
            fit_name=args.name,
            fit_n_times=args.fit_n_times,
            min_q=args.min_q,
            max_q=args.max_q,
            max_q_start=args.max_q_start,
            max_fit_steps=args.max_fit_steps,
            no_structures=args.mixture_n,
            pairedQ=args.pairedQ,
            rotation=args.rotation,
            max_backmap=args.max_backmap,
            defer_backmap_seconds=args.defer_backmap_seconds,
            do_foxs=(not args.no_foxs),
            backend=args.backend,
            cg2all_exec=args.cg2all_exec,
            disulfide_constraints_file=args.disulfide_constraints_file,
            foxs_cmd_default=args.foxs_cmd_default,
            python_exe=args.python_exe or sys.executable,
            terminate_on_foxs=args.terminate_on_foxs,
            terminate_threshold=args.terminate_threshold,
            terminate_confirmation_count=args.terminate_confirmation_count,
        )


        new_data_dir = os.path.join(os.getcwd(), "carbonara_runs", args.name)
        print("\nSetup completed successfully!")
        print(f"Initial files were created in: {refine_dir}")
        print(f"Files for Carbonara were copied to: {new_data_dir}")
        print(f"Run script created at: {run_script}")
        print(f"Watcher/backmapping Python: {args.python_exe or sys.executable}")
        if args.terminate_on_foxs:
            print(f"Run mode: terminate-on-FoXS, threshold={args.terminate_threshold}, confirmation_count={args.terminate_confirmation_count}")
        else:
            print("Run mode: maximal exploration (no individual predictor processes are terminated early)")
        print("\nTo run the refinement, execute:")
        print(f"cd {os.path.dirname(run_script)} && ./RunMe_" + str(args.name) + ".sh")

    except Exception as e:
        raise e
        print(f"Error during setup: {str(e)}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
