'''
------------------------------------------------------------------------------------------------
Weclome to Carbonara's Data Tools (CDT)!
------------------------------------------------------------------------------------------------

This package provides processing for Carbonara specific data input/outputs
and generating useful inputs for external programs 

------------------------------------------------------------------------------------------------
'''

import pandas as pd
import numpy as np

from scipy.spatial.distance import cdist

import os
import subprocess
import shutil
#import math
#import rmsd
from tqdm import tqdm

from typing import List, Set


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
        fixer = PDBFixer(filename=pdb_file)
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
        fixer = PDBFixer(filename=file_path)
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

def generate_random_structures(coords_file, fingerprint_file):

    '''Generate random structures changing one linker section at a time

    Parameters
    coords_file:       /path/ to CA coordinates.dat file
    fingerprint_file:  /path/ to fingerprint.dat file

    Return
    Generated structures are written to ~/rand_structures/.. section_*LINKERINDEX*.dat as xyz
    Linker Indices
    '''
    secondarystruct = get_secondary(fingerprint_file)

    linker_indices_sep = [find_linker_indices( section_finder(i)) for i in secondarystruct]

    linker_indices =[]


    currMax=0
    for i in range(0,len(linker_indices_sep)):
        for j in range(0,len(list(linker_indices_sep[i]))):
            linker_indices.append(list(linker_indices_sep[i])[j]+currMax)
        currMax = currMax +list(linker_indices_sep[i])[-1]+1

    #print(linker_indices)
    linker_indices =np.asarray(linker_indices)
    current = os.getcwd()
    random = 'rand_structures'
    random_working = os.path.join(current, random)

    if os.path.exists(random_working) and os.path.isdir(random_working):
        shutil.rmtree(random_working)

    os.mkdir(random_working)

    # try:

    # except OSError as error:
    #     print(str(error)[11:])

    # print('Beginning random structures generation \n')

    rand_file_dict = {}
    for l in tqdm(linker_indices):

        outputname = random_working+'/section_'+str(l)

#         !./generate_structure {fingerprint_file} {coords_file} {outputname} {l}
        result = subprocess.run([current+'/build/bin/generate_structure', fingerprint_file, coords_file, outputname, str(l)], capture_output=True, text=True)

    # print('')
    # print('Finished generating random structures')

    return linker_indices


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



# Writing to carbonara format

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
    fplst=np.loadtxt(fingerprint_file, str)
    fplstout= [np.asarray(list(fplst[i])) for i in range(2,len(fplst),2)]
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


import os
import json
import numpy as np

def load_pae_matrix(path):
    """
    Load a Predicted Aligned Error (PAE) matrix from a .json or .npy file.

    - .npy: Assumes the file contains a 2D NumPy array.
    - .json: Supports AlphaFold DB format and AlphaFold v3/custom formats.

    Args:
        path (str): Path to the .npy or .json file.

    Returns:
        numpy.ndarray: 2D array of PAE values.

    Raises:
        ValueError: If the file is unreadable or not a supported format.
    """
    ext = os.path.splitext(path)[-1].lower()

    if ext == '.npy':
        matrix = np.load(path)
        if matrix.ndim != 2:
            raise ValueError("PAE .npy file must contain a 2D array.")
        return matrix

    elif ext == '.json':
        with open(path, 'r') as f:
            data = json.load(f)

        # JSON Format 1: list of dicts (AlphaFold DB style)
        if isinstance(data, list):
            if not data:
                raise ValueError("PAE JSON list is empty.")
            entry = data[0]
            if not isinstance(entry, dict):
                raise ValueError("PAE JSON list does not contain a dictionary.")
            matrix = (entry.get('predicted_aligned_error') or
                      entry.get('predicted_alignment_error') or
                      entry.get('pae'))
            if matrix is None:
                raise ValueError("No PAE matrix found in JSON list entry.")

        # JSON Format 2: top-level dict
        elif isinstance(data, dict):
            matrix = (data.get('predicted_aligned_error') or
                      data.get('predicted_alignment_error') or
                      data.get('pae'))
            if matrix is None:
                raise ValueError("No supported PAE key found in JSON dictionary.")

        else:
            raise ValueError("Unsupported JSON structure for PAE data.")

        return np.array(matrix)

    else:
        raise ValueError(f"Unsupported file type: {ext}. Must be '.json' or '.npy'.")

def getFlexibleSections(file_path,pae_threshold = 1.0):
    pae_data = load_pae_matrix(file_path)
    if len(pae_data)==0:
        print("No 'pae' data found in the JSON file.")
    else:
        # Convert PAE data to a NumPy array for easier manipulation
        pae_array = np.array(pae_data)
        window_size = 2  # Define a window size
        high_pae_residues = []
        for i in range(len(pae_array) - window_size + 1):
            window = pae_array[i:i + window_size, i:i + window_size]
            avg_pae = np.mean(window)
            if avg_pae > pae_threshold:
                #print(f"Residue range {i}-{i + window_size} exceeds threshold with avg PAE {avg_pae}")
                high_pae_residues.extend(range(i, i + window_size))
    # Remove duplicate residue indices
        high_pae_residues = sorted(set(high_pae_residues))
        # Plotting the PAE heatmap with inverted colors
        plt.figure(figsize=(8, 6))
        heatmap = plt.imshow(pae_array, cmap='Greens_r', interpolation='nearest', aspect='auto')
        plt.colorbar(heatmap, label='Expected position error (Ångströms)')
        # Highlight regions identified as linkers
        for res in high_pae_residues:
            plt.axvline(x=res, color='red', linestyle='--', alpha=0.3)
            plt.axhline(y=res, color='red', linestyle='--', alpha=0.3)
        # Set labels
        plt.xlabel('Scored residue')
        plt.ylabel('Aligned residue')
        plt.title('Predicted Aligned Error (PAE) Heatmap')

        # Show the plot
        plt.show()
        # Second plot: Highlighting selected residues
        plt.figure(figsize=(10, 2))
        plt.bar(range(len(pae_array)), [1 if i in high_pae_residues else 0 for i in range(len(pae_array))], color='green', alpha=0.6)
        plt.xlabel('Residue Index')
        plt.ylabel('')
        plt.title('Residues Identified as Linkers or Between Domains')
        plt.yticks([0, 1], [])
        plt.grid(axis='x')

        # Show the second plot
        plt.show()
        return np.array([1 if i in high_pae_residues else 0 for i in range(len(pae_array))])

from typing import List, Set

def find_flexible_linker_sections(ss_string: str, pae_flags: List[int]) -> Set[int]:
    assert len(ss_string) == len(pae_flags), "Length of sequence and PAE list must match"
    
    flexible_linker_indices = set()
    i = 0
    section_index = 0

    while i < len(ss_string):
        current_char = ss_string[i]
        start = i

        # Move i to the end of the current segment
        while i < len(ss_string) and ss_string[i] == current_char:
            i += 1
        end = i

        # Check only if this section is a linker
        if current_char == '-' and any(pae_flags[start:end]):
            flexible_linker_indices.add(section_index)

        section_index += 1

    return flexible_linker_indices

def getFlexibility(paeFile,fingerprint_file):
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
        flexsec =getFlexibleSections('paerank_2_converted.json',pae_threshold = 0.7)
        fingerprint = get_secondary(fingerprint_file)
        return [list(find_flexible_linker_sections(np.concatenate(fingerprint), flexsec))]
    else:
        flexsec =getFlexibleSections(paeFile,pae_threshold = 0.7)
        fingerprint = get_secondary(fingerprint_file)
        return [list(find_flexible_linker_sections(np.concatenate(fingerprint), flexsec))]



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
    """
    with open(filename, 'w') as f:
        f.write(f"{len(chains)}\n\n")  # Write number of chains

        for chain in chains:
            f.write(f"{chain['sequence']}\n\n")
            f.write(f"{chain['structure']}\n\n")

def export_segment_list(segment_set, filename):
    """
    Write a set of segment numbers to a file, one per line, without trailing newline at the end.
    """
    segments = sorted(segment_set)
    with open(filename, 'w') as f:
        for i, seg in enumerate(segments):
            if i < len(segments) - 1:
                f.write(f"{seg}\n")
            else:
                f.write(f"{seg}")  # last line, no newline
def getResIDs_from_structure(pdb_fl, structure_file):
    """
    Returns resid_tensor split to match the chains in the structure file,
    ignoring original PDB chain IDs.
    """
    import re

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
