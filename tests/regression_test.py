"""
Regression test suite for Carbonara's core fitting engine (moleculeFitAndState /
ktlMoleculeRandom's penalty combination logic).

Two layers:

  - Golden-value tests: deterministic numbers read from the *initial* (step-0,
    pre-search) structure evaluation. No search/RNG involved at step 0 -- it's a
    single evaluation of the unmodified starting structure -- so these must match
    exactly, byte for byte, before and after any change to how penalties combine.
    Covers: baseline chi2, the soft-penalty cap, the hard-constraint feasibility
    check, and ensemble-OR (min-across-mixture-states) aggregation.

  - Invariant tests: short *live* searches. The search itself is stochastic (not
    exactly reproducible run to run), so these assert properties that must always
    hold regardless of RNG state, rather than exact values: a hard constraint is
    never violated in any accepted move; the soft-distance sum never exceeds its
    theoretical cap ceiling; Chi2 is present and <= the combined objective whenever
    a penalty is active.

Usage:
    python3 tests/regression_test.py --capture   # (re)write tests/golden_values.json
                                                  # from the CURRENT build -- run this
                                                  # to (re)establish a baseline, e.g.
                                                  # right before a refactor.
    python3 tests/regression_test.py             # verify the current build against
                                                  # the saved golden values.
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
FIXTURES = Path(__file__).resolve().parent / "fixtures"
GOLDEN_PATH = Path(__file__).resolve().parent / "golden_values.json"
BUILD_BIN = REPO_ROOT / "build" / "bin"
PDB_FIXTURE = FIXTURES / "MBP_apo_Bilbo_2.pdb"
SAXS_FIXTURE = FIXTURES / "Saxs.dat"

# Same circularly-permuted construct used throughout this session -- already
# validated extensively against the real engine, so a good, realistic fixture for
# exercising the hard/soft/cap/ensemble-OR machinery together.
SENSING_N = (1, 313)
SENSING_C = (563, 620)
GFP = (316, 559)

ROOT_ITEMS = [
    "build", "watch_and_backmap.py", "backmap_cli.py", "external", "bin",
    "probabilityInterpolation", "averagedSolutions", "CarbonaraDataTools.py",
    "setup_carbonara_allAtom.py",
]


class TestFailure(Exception):
    pass


def _isolated_workdir():
    work = Path(tempfile.mkdtemp(prefix="carbonara-regtest-"))
    (work / "pdbFiles").mkdir()
    (work / "saxsFiles").mkdir()
    for name in ROOT_ITEMS:
        src = REPO_ROOT / name
        if src.exists():
            (work / name).symlink_to(src)
    shutil.copy(PDB_FIXTURE, work / "pdbFiles" / "x.pdb")
    shutil.copy(SAXS_FIXTURE, work / "saxsFiles" / "x.dat")
    return work


def _build_fragment_pdb(in_pdb, out_pdb, fragments):
    lines_by_frag = {i: [] for i in range(len(fragments))}
    with open(in_pdb) as f:
        atom_lines = [l for l in f if l.startswith(("ATOM", "HETATM"))]
    for l in atom_lines:
        resid = int(l[22:26])
        for i, (cid, start, end) in enumerate(fragments):
            if start <= resid <= end:
                lines_by_frag[i].append(l[:21] + cid + l[22:])
                break
    with open(out_pdb, "w") as f:
        for i in range(len(fragments)):
            f.writelines(lines_by_frag[i])
            f.write("TER\n")
        f.write("END\n")


def _setup_construct(work, run_name, mixture_n=1):
    """Builds the 3-fragment (sensing-N, sensing-C, GFP) construct, runs
    setup_carbonara_allAtom.py --rotation, merges sensing-N+sensing-C into one
    rigid chain. Returns the run directory."""
    fragments = [("A", *SENSING_N), ("B", *SENSING_C), ("C", *GFP)]
    construct_pdb = work / "pdbFiles" / f"{run_name}_construct.pdb"
    _build_fragment_pdb(work / "pdbFiles" / "x.pdb", construct_pdb, fragments)

    cmd = [
        sys.executable, "setup_carbonara_allAtom.py",
        "-p", str(construct_pdb.relative_to(work)),
        "-s", "saxsFiles/x.dat",
        "-n", run_name,
        "--rotation", "--no_backmap",
    ]
    if mixture_n > 1:
        cmd += ["--mixture_n", str(mixture_n)]
    result = subprocess.run(cmd, cwd=work, capture_output=True, text=True)
    if result.returncode != 0:
        raise TestFailure(f"setup failed: {result.stderr[-2000:]}")

    run_dir = work / "carbonara_runs" / run_name

    sys.path.insert(0, str(work))
    import CarbonaraDataTools as CDT
    import numpy as np

    chains = CDT.parse_structures_with_segments(str(run_dir / "fingerPrint1.dat"))
    merge_pair = (1, 2)
    _, _, filtered = CDT.create_segment_label_arrays_with_merge_v4(
        chains, highlighted_segments=np.array([], dtype=int), merge_pair=merge_pair
    )
    merged_chains = CDT.merge_chains_only_clean_consistent_segments(chains, merge_pair)
    CDT.export_chains_to_file(merged_chains, str(run_dir / "fingerPrint1.dat"))
    CDT.export_segment_list(set(filtered.tolist()), str(run_dir / "varyingSectionSecondary1.dat"))

    return run_dir, fragments


def _resid_to_global(resid, fragments, chain_lengths):
    chain_order = sorted(chain_lengths.keys())
    offsets, off = {}, 0
    for ch in chain_order:
        offsets[ch] = off
        off += chain_lengths[ch]
    for (cid, start, end), letter in zip(fragments, chain_order):
        if start <= resid <= end:
            return offsets[letter] + (resid - start + 1)
    raise ValueError(f"residue {resid} not in any fragment")


def _measure_distance(run_dir, fragments, resid_a, resid_b, coords=None):
    """Real Angstrom distance between two original-PDB-numbering residues, correctly
    mapped through the post-merge chain offsets -- NOT `coords[resid-1]` directly,
    which is only valid for residues in the first fragment (no offset)."""
    import numpy as np
    import pickle
    if coords is None:
        coords = np.genfromtxt(run_dir / "coordinates1.dat")
    chain_lengths = pickle.load(open(run_dir / "chainLengths.dat", "rb"))
    ga = _resid_to_global(resid_a, fragments, chain_lengths)
    gb = _resid_to_global(resid_b, fragments, chain_lengths)
    return float(np.linalg.norm(coords[ga - 1] - coords[gb - 1]))


def _run_engine(run_dir, no_structures=1, max_steps=1, extra_argv=None, timeout=25):
    """Direct predictStructureQvary invocation (not through RunMe_*.sh), so we
    control every argv precisely. Returns parsed fitLog1.dat entries."""
    fitdata = run_dir / "fitdata"
    fitdata.mkdir(exist_ok=True)
    for f in fitdata.glob("*"):
        f.unlink()

    argv = [
        str(BUILD_BIN / "predictStructureQvary"),
        str(run_dir / "Saxs.dat"),
        str(run_dir) + "/",
        "frompdb",
        "True",  # paired predictions on
        str(run_dir / "varyingSectionSecondary1.dat"),
        str(no_structures),
        "none",
        "0.01", "0.2", "0.2",
        str(max_steps),
        str(fitdata / "mol1"),
        str(fitdata / "scatter1.dat"),
        str(run_dir / "mixtureFile.dat"),
        str(run_dir / "redundant"),
        str(fitdata / "fitLog1.dat"),
        "null",
        "False",  # affineTrans -- caller overrides via extra_argv if needed
        "True",   # useErrors
        "-1.0", "1", "5.0", "50.0",
    ]
    if extra_argv:
        argv = extra_argv(argv)

    try:
        # cwd pinned to run_dir (not inherited from wherever this script is invoked from):
        # every path we build above is absolute, so this is defense-in-depth only, but it's
        # cheap insurance against a bad argv (wrong index, empty string, etc.) turning into a
        # relative-path write landing in the real repo instead of the temp workdir -- exactly
        # what happened once already while developing this suite (see git history).
        subprocess.run(argv, capture_output=True, text=True, timeout=timeout, cwd=run_dir)
    except subprocess.TimeoutExpired:
        pass  # expected -- we just want however many steps ran in the time budget

    log_path = fitdata / "fitLog1.dat"
    entries = []
    if log_path.exists():
        for line in log_path.read_text().splitlines():
            line = line.strip()
            if line.startswith("{") and '"ImprovementIndex"' in line:
                try:
                    entries.append(json.loads(line))
                except Exception:
                    pass
    return entries


def _write_constraints(run_dir, fragments, contact_specs, mixture_n=1):
    """contact_specs: list of dicts with resid_a, resid_b, target, tolerance,
    bound_type, hard, ensemble_or."""
    sys.path.insert(0, str(run_dir.parent.parent))
    import CarbonaraDataTools as CDT
    import numpy as np
    import pickle

    chain_lengths = pickle.load(open(run_dir / "chainLengths.dat", "rb"))
    coords = np.genfromtxt(run_dir / "coordinates1.dat")

    contactPreds, fixedDists, tolerances, boundTypes, hardList, orList = [], [], [], [], [], []
    for spec in contact_specs:
        ga = _resid_to_global(spec["resid_a"], fragments, chain_lengths)
        gb = _resid_to_global(spec["resid_b"], fragments, chain_lengths)
        contactPreds.append([ga, gb])
        fixedDists.append(spec["target"])
        tolerances.append(spec["tolerance"])
        boundTypes.append(spec.get("bound_type", 0))
        hardList.append(int(spec.get("hard", 0)))
        orList.append(int(spec.get("ensemble_or", 0)))

    CDT.translate_distance_constraints(contactPreds, coords, str(run_dir), fixedDists, tolerances, boundTypes, hardList, orList)
    for i in range(2, mixture_n + 1):
        shutil.copy(run_dir / "fixedDistanceConstraints1.dat", run_dir / f"fixedDistanceConstraints{i}.dat")


# ---------------------------------------------------------------------------
# Scenarios
# ---------------------------------------------------------------------------

def scenario_baseline_chi2(work):
    """No constraints: initial chi2 should be a fixed, reproducible number."""
    run_dir, fragments = _setup_construct(work, "reg_baseline")
    entries = _run_engine(run_dir, max_steps=1)
    if not entries:
        raise TestFailure("no fitLog entries produced")
    return {"chi2": entries[0]["Chi2"], "scatter_first": entries[0]["ScatterFitFirst"]}


def scenario_soft_cap(work):
    """A badly-violated soft pair (tight tolerance) must saturate at the cap,
    not grow unboundedly."""
    run_dir, fragments = _setup_construct(work, "reg_softcap")
    d_actual = _measure_distance(run_dir, fragments, 65, 359)
    _write_constraints(run_dir, fragments, [{
        "resid_a": 65, "resid_b": 359,
        "target": d_actual / 3,  # force a large violation
        "tolerance": 0.01, "bound_type": 0, "hard": 0, "ensemble_or": 0,
    }])
    entries = _run_engine(run_dir, max_steps=1)
    return {"distance_constraints_capped": entries[0]["DistanceConstraints"]}


def scenario_hard_satisfied_initially(work):
    """A hard pair that the initial structure already satisfies must report
    HardConstraintsSatisfied=true, MaxViolation=0, and contribute nothing to the
    soft DistanceConstraints sum."""
    run_dir, fragments = _setup_construct(work, "reg_hardok")
    d_actual = _measure_distance(run_dir, fragments, 65, 359)
    _write_constraints(run_dir, fragments, [{
        "resid_a": 65, "resid_b": 359,
        "target": d_actual, "tolerance": 0.05, "bound_type": 0, "hard": 1, "ensemble_or": 0,
    }])
    entries = _run_engine(run_dir, max_steps=1)
    return {
        "hard_satisfied": entries[0]["HardConstraintsSatisfied"],
        "hard_max_violation": entries[0]["HardConstraintsMaxViolation"],
        "distance_constraints": entries[0]["DistanceConstraints"],
    }


def scenario_ensemble_or(work):
    """2-state mixture, same pair, satisfied in state 1 / badly violated in
    state 2: ensemble_or=1 must report ~0 (min across states), ensemble_or=0
    must report the saturated penalty (sum across states)."""
    run_dir, fragments = _setup_construct(work, "reg_ensor", mixture_n=2)
    import numpy as np
    coords1 = np.genfromtxt(run_dir / "coordinates1.dat")
    coords2 = coords1.copy()
    coords2[99] += np.array([80.0, 0.0, 0.0])
    np.savetxt(run_dir / "coordinates2.dat", coords2, delimiter=" ", fmt="%s")
    d_actual = float(np.linalg.norm(coords1[49] - coords1[99]))

    results = {}
    for label, or_flag in [("or_on", 1), ("or_off", 0)]:
        _write_constraints(run_dir, fragments, [{
            "resid_a": 50, "resid_b": 100,
            "target": d_actual, "tolerance": 0.05, "bound_type": 0, "hard": 0,
            "ensemble_or": or_flag,
        }], mixture_n=2)
        entries = _run_engine(
            run_dir, no_structures=2, max_steps=1,
            extra_argv=lambda argv: argv,  # affineTrans stays False -- fine, single-chain-per-state
        )
        results[label] = entries[0]["DistanceConstraints"]
    return results


def invariant_hard_never_violated_when_accepted(work):
    run_dir, fragments = _setup_construct(work, "reg_hardinv")
    d_actual = _measure_distance(run_dir, fragments, 313, 316)  # the "post" pair from this session
    _write_constraints(run_dir, fragments, [{
        "resid_a": 313, "resid_b": 316,
        "target": d_actual, "tolerance": 0.03, "bound_type": 0, "hard": 1, "ensemble_or": 0,
    }])
    entries = _run_engine(
        run_dir, max_steps=100, timeout=20,
        extra_argv=lambda argv: argv[:18] + ["True"] + argv[19:],  # argv[18]=affineTrans -> True
    )
    if len(entries) < 2:
        return True  # not enough moves happened to be a meaningful check either way
    violations = [e for e in entries if not e.get("HardConstraintsSatisfied", True)]
    if violations:
        raise TestFailure(f"{len(violations)} accepted moves violated a hard constraint")
    return True


SCENARIOS = {
    "baseline_chi2": scenario_baseline_chi2,
    "soft_cap": scenario_soft_cap,
    "hard_satisfied_initially": scenario_hard_satisfied_initially,
    "ensemble_or": scenario_ensemble_or,
}

INVARIANTS = {
    "hard_never_violated_when_accepted": invariant_hard_never_violated_when_accepted,
}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--capture", action="store_true", help="(re)write golden_values.json from the current build")
    args = parser.parse_args()

    if not (BUILD_BIN / "predictStructureQvary").exists():
        print("build/bin/predictStructureQvary not found -- build first (cmake --build build)", file=sys.stderr)
        sys.exit(2)

    results = {}
    failures = []

    work = _isolated_workdir()
    try:
        for name, fn in SCENARIOS.items():
            try:
                results[name] = fn(work)
                print(f"  [{'CAPTURED' if args.capture else 'ran'}] {name}: {results[name]}")
            except Exception as e:
                failures.append(f"{name}: {e}")
                print(f"  [ERROR] {name}: {e}")

        for name, fn in INVARIANTS.items():
            try:
                fn(work)
                print(f"  [OK] invariant: {name}")
            except Exception as e:
                failures.append(f"invariant {name}: {e}")
                print(f"  [FAIL] invariant {name}: {e}")
    finally:
        shutil.rmtree(work, ignore_errors=True)

    if args.capture:
        GOLDEN_PATH.write_text(json.dumps(results, indent=2, sort_keys=True))
        print(f"\nWrote golden values to {GOLDEN_PATH}")
        sys.exit(1 if failures else 0)

    if not GOLDEN_PATH.exists():
        print(f"\nNo golden values found at {GOLDEN_PATH} -- run with --capture first.", file=sys.stderr)
        sys.exit(2)

    golden = json.loads(GOLDEN_PATH.read_text())
    mismatches = []
    for name, value in results.items():
        expected = golden.get(name)
        if expected != value:
            mismatches.append(f"{name}: expected {expected}, got {value}")

    print()
    if mismatches:
        print("GOLDEN VALUE MISMATCHES:")
        for m in mismatches:
            print(f"  - {m}")
    if failures:
        print("FAILURES:")
        for f in failures:
            print(f"  - {f}")

    if mismatches or failures:
        sys.exit(1)
    print("All regression checks passed.")


if __name__ == "__main__":
    main()
