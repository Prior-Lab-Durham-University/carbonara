#!/usr/bin/env python3

import argparse
import os
import shutil
import sys
import CarbonaraDataTools as cdt
import numpy as np

def write_runme(working_path, fit_name, fit_n_times, min_q, max_q,max_q_start, max_fit_steps, pairedQ=False, rotation=False):

    curr = os.getcwd()
    script_name = 'RunMe_' + str(fit_name) + '.sh'
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
        # List of files to copy
        files_to_copy = [
            'Saxs.dat',
            'mixtureFile.dat',
            'fingerPrint1.dat',
            'varyingSectionSecondary1.dat',
            'coordinates1.dat'
        ]
        
        for file in files_to_copy:
            source = os.path.join(working_path, file)
            destination = os.path.join(new_data_dir, file)
            
            if os.path.exists(source):
                shutil.copy2(source, destination)
                print(f"Copied {file} to {destination}")
            else:
                print(f"Warning: Source file {source} not found")
        
        # Create an empty redundant file
        with open(os.path.join(new_data_dir, 'redundant'), 'w') as f:
            f.write('')
        
    except Exception as e:
        print(f"Warning: Could not copy files: {e}")
    
    # Write the script
    with open(run_file, 'w+') as fout:
        fout.write('#!/bin/bash\n')
        fout.write('# Determine the root directory based on the script location\n')
        fout.write('ROOT=$(dirname "$(readlink -f "$0")")\n\n')
        
        fout.write('# Directory to clear before running\n')
        fout.write(f'CLEAR_DIR="$ROOT/{data_path}/fitdata"\n\n')
        
        fout.write('# Clear the directory\n')
        fout.write('echo "Clearing directory: $CLEAR_DIR"\n')
        fout.write('rm -r "$CLEAR_DIR"/*\n')
        fout.write('mkdir -p "$CLEAR_DIR"\n\n')
        
        fout.write('### argv[ 1] scattering data file\n')
        fout.write(f'ScatterFile=$ROOT/{data_path}/Saxs.dat\n')
        
        fout.write('### argv[ 2] sequence file location\n')
        fout.write(f'fileLocs=$ROOT/{data_path}/\n')
        
        fout.write('### argv[ 3] restart tag (use to start from existing prediction)\n')
        fout.write('initialCoordsFile=frompdb\n')
        
        fout.write('### argv[ 4] paired distances file (can be empty)\n')
        if pairedQ:
            fout.write(f'pairedPredictions=$ROOT/{data_path}/fixedDistanceConstraints1.dat\n')
        else:
            fout.write('pairedPredictions=False\n')
        
        fout.write('### argv[ 5] fixed sections file (again can be empty)\n')
        fout.write(f'fixedsections=$ROOT/{data_path}/varyingSectionSecondary1.dat\n')
        
        fout.write('### argv[ 6] number of structures\n')
        fout.write('noStructures=1\n')
        
        fout.write('### argv[ 7] request to apply hydrophobic covering WITHIN monomers will be a list of sections on which to apply it -- Currently not used\n')
        fout.write('withinMonomerHydroCover=none\n')
        
        fout.write('### argv[ 8] kmin\n')
        fout.write(f'kmin={min_q}\n')
        
        fout.write('### argv[ 9] kmax\n')
        fout.write(f'kmax={max_q}\n')

        fout.write('### argv[ 10] kmax Start\n')
        fout.write(f'kmaxStart={max_q_start}\n')
        
        fout.write('### argv[11] Max number of fitting steps\n')
        fout.write(f'maxNoFitSteps={max_fit_steps}\n')
        
        fout.write('### argv[12] prediction file - mol[i] in the fitting folder\n')
        fout.write(f'predictionFile=$ROOT/{data_path}/fitdata\n')
        
        fout.write('### argv[13] scattering output file\n')
        fout.write(f'scatterOut=$ROOT/{data_path}/fitdata\n')
        
        fout.write('### argv[14] mixture list file, alist of sets of numbers indicatig the allowed set of mixture percentages of each species (e.g. dimer 20 monomer 80)\n')
        fout.write(f'mixtureFile=$ROOT/{data_path}/mixtureFile.dat\n')
        
        fout.write('### argv[15] previous fit string in form fitname/mol6Substep_10_1.dat+fitname/mol6Substep_10_2.dat\n')
        fout.write(f'prevFitStr=$ROOT/{data_path}/redundant\n')
        
        fout.write('### argv[16] log file location\n')
        fout.write(f'logLoc=$ROOT/{data_path}/fitdata\n')
        
        fout.write('### argv[17] last line of the previous fit log, this is only used for a restart if argv[3] = True\n')
        fout.write('endLinePrevLog=null\n')
        
        fout.write('### argv[18] is true if we want to apply affine rotations,false if not.\n')
        if rotation:
            fout.write('affineTrans=True\n')
        else:
            fout.write('affineTrans=False\n')
        
        fout.write(f'for i in {{1..{fit_n_times}}}\n')
        fout.write('do\n')
        fout.write('    echo "\\n"\n')
        fout.write('    echo " >> Run number : $i "\n')
        fout.write('    echo "\\n"\n')
        fout.write('    echo "Max number of fitting steps: " $maxNoFitSteps\n')
        fout.write('    echo "\\n"\n')
        fout.write('    $ROOT/build/bin/predictStructureQvary $ScatterFile $fileLocs $initialCoordsFile $pairedPredictions $fixedsections $noStructures $withinMonomerHydroCover $kmin $kmax $kmaxStart $maxNoFitSteps $predictionFile/mol$i $scatterOut/scatter$i.dat $mixtureFile $prevFitStr $logLoc/fitLog$i.dat $endLinePrevLog $affineTrans\n')
        fout.write('done\n')
    
    # Make the script executable
    os.chmod(run_file, 0o755)
    
    return run_file

def main():
    parser = argparse.ArgumentParser(description='Setup Carbonara processing pipeline')
    parser.add_argument('-p', '--pdb', required=True, help='Path to input PDB file')
    parser.add_argument('-s', '--saxs', required=True, help='Path to input SAXS data file')
    parser.add_argument('-n', '--name', required=True, help='Name for this protein/refinement')
    parser.add_argument('-d', '--dir', default=os.getcwd(), help='Base directory (default: current directory)')
    
    # Additional parameters for write_runme
    parser.add_argument('--fit_n_times', type=int, default=20, help='Number of times to run the fit (default: 5)')
    parser.add_argument('--min_q', type=float, default=0.01, help='Minimum q-value (default: 0.01)')
    parser.add_argument('--max_q', type=float, default=0.2, help='Maximum q-value (default: 0.2)')
    parser.add_argument('--max_q_start', type=float, default=0.1, help='Maximum q-value to start fitting to (default: 0.1)')
    parser.add_argument('--max_fit_steps', type=int, default=10000, help='Maximum number of fitting steps (default: 1000)')
    parser.add_argument('--pairedQ', action='store_true', help='Use paired predictions')
    parser.add_argument('--rotation', action='store_true', help='Apply affine rotations')
    
    args = parser.parse_args()
    
    try:
        # Setup master directory
        fit_master_dir = cdt.setup_fit_master_dir(root_dir=args.dir, fit_master_name='carbonara_runs')
        
        # Setup refinement directory
        refine_dir = cdt.setup_refinement_dir(args.name, fit_master_dir)
        print(f"Created directory structure in: {refine_dir}")
        
        # Process PDB and extract structure information
        coords_chains, sequence_chains, secondary_structure_chains, missing_residues_chains = cdt.pull_structure_from_pdb(args.pdb)
        # Give warning if missing residues are found
        print(len(coords_chains))
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

            # Convert to array of objects if needed
        coords_chains = np.array(new_coords_chains, dtype=object)
        sequence_chains = np.array(new_sequence_chains, dtype=object)
        secondary_structure_chains = np.array(new_secondary_structure_chains, dtype=object)

        # collapse coordinates file into one chain
        
        for i, coords in enumerate(coords_chains):
            if i==0:
                coords_full=coords
            else:
                coords_full = np.concatenate((coords_full,coords),axis=0)

        # write this to file
        
        coords_files = []
        coords_files.append(cdt.write_coordinates_file(coords_full, working_path=refine_dir, carb_index=1))
        
        # Write fingerprint file
        number_of_chains = len(coords_chains)
        fingerprint_file = cdt.write_fingerprint_file(
            number_chains=number_of_chains,
            sequence=sequence_chains,
            secondary_structure=secondary_structure_chains,
            working_path=refine_dir
        )

        # write mixture file - used for ensemble refinement, currently not used - writes 1 to mixture file
        mixture_file = cdt.write_mixture_file(working_path=refine_dir)
        
        # Copy SAXS file to Saxs.dat (this is the file that Carbonara will use)
        shutil.copy2(args.saxs, os.path.join(refine_dir, 'Saxs.dat'))
        
        # auto select flexible linker chains that dont break inter-beta sheets
        varying_linker_chains = []
        index =0
        for coord_file in coords_files:
            varying_linker_chains.append(cdt.auto_select_varying_linker(coord_file, fingerprint_file))
            index = index +1
            
        # write flexible linkers to files (varysections1.dat, varysections2.dat, etc [each file is for a different chain])
        print(varying_linker_chains)
        varying_section_files = []
        for varying_linkers in varying_linker_chains:
            varying_section_files.append(cdt.write_varysections_file(varying_linkers, refine_dir))
        
        # Write the RunMe_bsa.sh script
        run_script = write_runme(
            working_path=refine_dir,
            fit_name=args.name,
            fit_n_times=args.fit_n_times,
            min_q=args.min_q,
            max_q=args.max_q,
            max_q_start=args.max_q_start,
            max_fit_steps=args.max_fit_steps,
            pairedQ=args.pairedQ,
            rotation=args.rotation
        )
        
        # Updated output message

        new_data_dir = os.path.join(os.getcwd(), "carbonara_runs", args.name)
        print("\nSetup completed successfully!")
        print(f"Initial files were created in: {refine_dir}")
        print(f"Files for Carbonara were copied to: {new_data_dir}")
        print(f"Run script created at: {run_script}")
        print("\nTo run the refinement, execute:")
        print(f"cd {os.path.dirname(run_script)} && ./RunMe_" +str(args.name)+ ".sh")
        
    except Exception as e:
        print(f"Error during setup: {str(e)}", file=sys.stderr)
        sys.exit(1)

if __name__ == '__main__':
    main()
