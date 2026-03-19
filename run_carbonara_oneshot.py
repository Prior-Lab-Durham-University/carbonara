#!/usr/bin/env python3

import argparse
import os
import shutil
import sys
import pickle
from string import ascii_uppercase
from typing import Optional, List
import subprocess

import CarbonaraDataTools as cdt
import numpy as np


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
):
    curr = os.getcwd()
    script_name = "RunMe_" + str(fit_name) + ".sh"
    run_file = os.path.join(curr, script_name)

    # Path to the data directory (relative to ROOT)
    data_path = f"carbonara_runs/{fit_name}"

    # Create new directories
    new_data_dir = os.path.join(curr, "carbonara_runs", fit_name)
    fitdata_dir = os.path.join(new_data_dir, "fitdata")

    os.makedirs(new_data_dir, exist_ok=True)
    os.makedirs(fitdata_dir, exist_ok=True)

    # Copy necessary files from refine_dir to new_data_dir
    try:
        files_to_copy = ["Saxs.dat", "mixtureFile.dat"]

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

    # Write the script
    with open(run_file, "w+") as fout:
        fout.write("#!/bin/bash\n")
        fout.write("# Determine the root directory based on the script location\n")
        fout.write('ROOT=$(dirname "$(readlink -f "$0")")\n\n')

        fout.write("# Directory to clear before running\n")
        fout.write(f'CLEAR_DIR="$ROOT/{data_path}/fitdata"\n\n')

        fout.write("# Clear the directory\n")
        fout.write('echo "Clearing directory: $CLEAR_DIR"\n')
        fout.write('rm -rf "$CLEAR_DIR"/*\n')
        fout.write('mkdir -p "$CLEAR_DIR"\n\n')

        fout.write("### argv[ 1] scattering data file\n")
        fout.write(f"ScatterFile=$ROOT/{data_path}/Saxs.dat\n")

        fout.write("### argv[ 2] sequence file location\n")
        fout.write(f"fileLocs=$ROOT/{data_path}/\n")

        fout.write("### argv[ 3] restart tag (use to start from existing prediction)\n")
        fout.write("initialCoordsFile=frompdb\n")

        fout.write("### argv[ 4] paired distances file (can be empty)\n")
        if pairedQ:
            fout.write(f"pairedPredictions=$ROOT/{data_path}/fixedDistanceConstraints1.dat\n")
        else:
            fout.write("pairedPredictions=False\n")

        fout.write("### argv[ 5] fixed sections file (again can be empty)\n")
        fout.write(f"fixedsections=$ROOT/{data_path}/varyingSectionSecondary1.dat\n")

        fout.write("### argv[ 6] number of structures\n")
        fout.write(f"noStructures={int(no_structures)}\n")

        fout.write(
            "### argv[ 7] request to apply hydrophobic covering WITHIN monomers will be a list of sections on which to apply it -- Currently not used\n"
        )
        fout.write("withinMonomerHydroCover=none\n")

        fout.write("### argv[ 8] kmin\n")
        fout.write(f"kmin={min_q}\n")

        fout.write("### argv[ 9] kmax\n")
        fout.write(f"kmax={max_q}\n")

        fout.write("### argv[ 10] kmax Start\n")
        fout.write(f"kmaxStart={max_q_start}\n")

        fout.write("### argv[11] Max number of fitting steps\n")
        fout.write(f"maxNoFitSteps={max_fit_steps}\n")

        fout.write("### argv[12] prediction file - mol[i] in the fitting folder\n")
        fout.write(f"predictionFile=$ROOT/{data_path}/fitdata\n")

        fout.write("### argv[13] scattering output file\n")
        fout.write(f"scatterOut=$ROOT/{data_path}/fitdata\n")

        fout.write(
            "### argv[14] mixture list file, alist of sets of numbers indicatig the allowed set of mixture percentages of each species (e.g. dimer 20 monomer 80)\n"
        )
        fout.write(f"mixtureFile=$ROOT/{data_path}/mixtureFile.dat\n")

        fout.write(
            "### argv[15] previous fit string in form fitname/mol6Substep_10_1.dat+fitname/mol6Substep_10_2.dat\n"
        )
        fout.write(f"prevFitStr=$ROOT/{data_path}/redundant\n")

        fout.write("### argv[16] log file location\n")
        fout.write(f"logLoc=$ROOT/{data_path}/fitdata\n")

        fout.write(
            "### argv[17] last line of the previous fit log, this is only used for a restart if argv[3] = True\n"
        )
        fout.write("endLinePrevLog=null\n")

        fout.write(
            "### argv[18] is true if we want to apply affine rotations,false if not.\n"
        )
        fout.write("affineTrans=True\n" if rotation else "affineTrans=False\n")

        fout.write(
            "### argv[19] is true if we want to use errors in the scattering calculation false if not.\n"
        )
        fout.write("useErrors=True\n")

        fout.write(f"for i in {{1..{fit_n_times}}}\n")
        fout.write("do\n")
        fout.write('    echo "\\n"\n')
        fout.write('    echo " >> Run number : $i "\n')
        fout.write('    echo "\\n"\n')
        fout.write('    echo "Max number of fitting steps: " $maxNoFitSteps\n')
        fout.write('    echo "\\n"\n')
        fout.write(
            "    $ROOT/build/bin/predictStructureQvary "
            "$ScatterFile $fileLocs $initialCoordsFile $pairedPredictions $fixedsections "
            "$noStructures $withinMonomerHydroCover $kmin $kmax $kmaxStart $maxNoFitSteps "
            "$predictionFile/mol$i $scatterOut/scatter$i.dat $mixtureFile $prevFitStr "
            "$logLoc/fitLog$i.dat $endLinePrevLog $affineTrans $useErrors &\n"
        )
        fout.write("done\n")

    # Make the script executable
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
    result = {label: len(structure) for label, structure in zip(ascii_uppercase, structure_lines)}
    return result


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
    

    args = parser.parse_args()

    try:
        if args.mixture_n < 1:
            raise ValueError("--mixture_n must be >= 1")

        # Setup master directory
        fit_master_dir = cdt.setup_fit_master_dir(root_dir=args.dir, fit_master_name="carbonara_runs")

        # Setup refinement directory
        refine_dir = cdt.setup_refinement_dir(args.name, fit_master_dir)
        print(f"Created directory structure in: {refine_dir}")

        # Process PDB and extract structure information
        coords_chains, sequence_chains, secondary_structure_chains, missing_residues_chains = (
            cdt.pull_structure_from_pdb(args.pdb)
        )

        print("number of chains is ", len(coords_chains))
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

        # Copy SAXS file to Saxs.dat (this is the file that Carbonara will use)
        cdt.write_saxs(args.saxs, refine_dir)

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
            for coord_file in coords_files:
                varying_linker_chains.append(cdt.auto_select_varying_linker(coord_file, fingerprint_file))

        # write flexible linkers to files (varysections1.dat, varysections2.dat, etc [each file is for a different chain])
        varying_section_files = []
        for varying_linkers in varying_linker_chains:
            varying_section_files.append(cdt.write_varysections_file(varying_linkers, refine_dir))

        # check for length 2 varying sections and filter them out
        filepath = refine_dir + "/fingerPrint1.dat"
        vs_path = refine_dir + "/varyingSectionSecondary1.dat"

        # load robustly: always get a 1-D array (even if file has 1 int)
        try:
            target_segments = np.loadtxt(vs_path, dtype=int, ndmin=1)
        except ValueError:
            # happens if the file is empty / whitespace
            target_segments = np.array([], dtype=int)
            
            target_segments = np.atleast_1d(target_segments)
            
            # If nothing to filter, keep file empty and move on
            if target_segments.size == 0:
                # optional: ensure empty file exists
                open(vs_path, "w").close()
            else:
                filtered_segments = cdt.get_segment_lengths_from_file(filepath, target_segments)
                filtered_segments = np.atleast_1d(np.asarray(filtered_segments, dtype=int))
                np.savetxt(vs_path, filtered_segments, fmt="%i")


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
        )

        # store chain lengths
        chain_lengths = parse_structure_lengths(refine_dir + "/fingerPrint1.dat")
        with open(refine_dir + "/chainLengths.dat", "wb") as f:
            pickle.dump(chain_lengths, f)

        new_data_dir = os.path.join(os.getcwd(), "carbonara_runs", args.name)
        print("\nSetup completed successfully!")
        print(f"Initial files were created in: {refine_dir}")
        print(f"Files for Carbonara were copied to: {new_data_dir}")
        print(f"Run script created at: {run_script}")
        print("\nLaunching refinement now...")
        print(f"Working directory: {os.path.dirname(run_script)}")

        # finally run the script without waiting for background jobs
        subprocess.Popen(
            ["bash", run_script],
            cwd=os.path.dirname(run_script)
        )
        
    except Exception as e:
        print(f"Error during setup: {str(e)}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
