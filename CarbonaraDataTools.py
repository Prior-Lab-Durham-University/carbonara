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
import math
import rmsd
from tqdm import tqdm

from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
from Bio.PDB.DSSP import make_dssp_dict
from DSSPparser import parseDSSP

import biobox as bb

import shutil
import re

from scipy import interpolate

from plotly.subplots import make_subplots
import plotly.graph_objects as go
import json

from glob import glob

import hdbscan


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
    M.import_pdb(pdb_file)
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


# The meet of extracting 

def pull_structure_from_pdb(pdb_file):
    """
    Pulls the structure from a (single) PDB file and returns the coordinates and sequence of the chain(s).

    Parameters:
        pdb_file (str): The path of the PDB file.

    Returns:
        coords_chain (list): A list of numpy arrays containing the coordinates of the chain(s).
        sequence_chain (list): A list of numpy arrays containing the sequence of the chain(s).
        missing_residues_chain (list): A list of numpy arrays containing the missing residues of the chain(s).

    Raises:
        ValueError: If no chains are found in the PDB file.

    Example:
        coords_chains, sequence_chains, missing_residues_chains = pull_structure_from_pdb('/path/to/pdb/file.pdb') 
    """
    
    M = pdb_2_biobox(pdb_file)

    coords_chains = []
    sequence_chains = []
    secondary_structure_chains = []

    missing_residues_chains = []

    chains = M.data['chain'].unique()
    ca_cond = M.data['name']=='CA'
    

    if len(chains) > 1:

        coords_chains = []
        sequence_chains = []

        for c in chains:

            # indices within chain that are CA
            cond = (M.data['chain']==c) & (ca_cond)

            # append coordinates of CA atoms in chain
            coords = M.coordinates[0][cond]
            coords_chains.append( coords )

            # append sequence of chain
            sequence_chains.append( M.data['resname'][cond].map( get_residue_map() ).values )

            # find secondary structure of chain using geometric method
            secondary_structure_chains.append( geometric_secondary_structure(coords, threshold=0.2) )

            # find missing residues in chain
            resIDs = M.data['resid'][cond].values
            missing_residues = find_missing_residues(resIDs)
            missing_residues_chains.append(missing_residues)

    elif len(chains) == 1:

        coords = M.coordinates[0][ca_cond]
        coords_chains = [coords]
        
        sequence_chains = [M.data['resname'][ca_cond].map( get_residue_map() ).values]
        
        secondary_structure_chains = [geometric_secondary_structure(coords, threshold=0.2)]

        resIDs = M.data['resid'][ca_cond].values
        missing_residues_chains = [find_missing_residues(resIDs)]

    else:
        raise ValueError("No chains found in pdb file")

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


# Secondary structure

# > DSSP for secondary structure

def DSSP_structure_extractor(pdb_file):
    
    '''
    Use DSSP for predict secondary structure from a PDB file
    
    returns a numpy array of residue secondary structure labels
    '''
    
    p = PDBParser()
    structure = p.get_structure("PDB_file", pdb_file)
    model = structure[0]

    # **************** LOOK HERE! THE MKDSSP LOC NEEDS TO BE CHANGED FOR PLATFORM ********************
    dssp = DSSP(model, pdb_file)

    simplify_dict = {'H' : 'H', 'P' : 'H', 'B': 'S', 'E': 'S', 'G': 'H', 'I': 'H', 'T': '-', 'S': '-', '-': '-', ' ': '-'}
    secondary_struct = []

    for key in list(dssp.keys()):
        secondary_struct.append( simplify_dict[ dssp[key][2] ] )
        
    return np.asarray(secondary_struct)  


def DSSP_structure_extractor_from_file(dssp_file):

    parser = parseDSSP(dssp_file)
    parser.parse()
    pddict = parser.dictTodataframe()


    simplify_dict = {'  H' : 'H', '  P' : 'H', '  B': 'S', '  E': 'S', '  G': 'H', '  I': 'H', '  T': '-', '  S': '-', '-': '-', ' ': '-','   ': '-','*  ':'-'}
    secondary_struct = []

    for key in pddict["struct"]:
         secondary_struct.append( simplify_dict[key] )

    return np.asarray(secondary_struct)
   # return np.asarray(secondary_struct)


def read_dssp_file(dssp_filename):

    simplify_dict = {'H': 'H', 'B': 'S', 'E': 'S', 'G': 'H', 'I': 'H', 'T': '-', 'S': '-', '-': '-', ' ': '-'}

    lines=[]
    with open(dssp_filename) as input_data:
        # Skips text before the beginning of the interesting block:
        for line in input_data:
            if line.strip() == '#  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    N-H-->O    O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA':
                break
        # Reads text until the end of the block:
        for line in input_data:  # This keeps reading the file
            lines.append(simplify_dict[line[16]])
    return ''.join(lines)


def simplify_secondary(dssp_struct):

    simplify_dict = {'H': 'H', 'B': 'S', 'E': 'S', 'G': 'H', 'I': 'H', 'T': '-', 'S': '-', '-': '-', ' ': '-'}

    secondary_structure = []

    for s in dssp_struct:

        if s not in list(simplify_dict.keys()):
            print('>>> ', s, ' <<<')
            raise Exception('Secondary structure not recognised!')

        secondary_structure.append(simplify_dict[s])

    return secondary_structure


def get_secondary(fingerprint_file):
    fplst=np.loadtxt(fingerprint_file, str)
    fplstout= [np.asarray(list(fplst[i])) for i in range(2,len(fplst),2)]
    return fplstout


# > Curvature & Torsion calculation for geometrically inferring secondary structure 

def dotproduct(v1, v2):
  return sum((a*b) for a, b in zip(v1, v2))


def length(v):
  return math.sqrt(dotproduct(v, v))


def nvec(v):
    len = length(v)
    if(len>0.000000001):
        return v/len
    else:
        return v


def angle(v1, v2):
  return math.acos(dotproduct(v1, v2) / (length(v1) * length(v2)))


def flatten(l):
    return [item for sublist in l for item in sublist]


def kapTau(subset):
    midPoint1 = 0.5*(subset[0]+subset[1])
    midPoint2 = 0.5*(subset[1]+subset[2])
    midPoint3 = 0.5*(subset[2]+subset[3])
    dif13 = midPoint3-midPoint1
    dif23 = midPoint3-midPoint2
    dif12 = midPoint2-midPoint1
    difcross = np.cross(dif13,dif23)
    kap_denom = length(dif12)
    kap_num = length(difcross)/length(dif13)/length(dif23)
    kap = 2.0*kap_num/kap_denom
    taudif1 = subset[1]-subset[0]
    taudif2 = subset[2]-subset[1]
    taudif3 = subset[3]-subset[2]
    N1 = np.cross(taudif1,taudif2)
    N2 = np.cross(taudif2,taudif3)
    n1 = nvec(N1);n2= nvec(N2)
    l = 0.333*(length(taudif1)+length(taudif2)+length(taudif3))
    ang = angle(n1,n2)
    tsign = np.sign(np.dot(taudif2,np.cross(N1,N2)))
    tau = tsign*2.0*np.sin(ang/2.0)/l
    return [kap,tau]


def getKapTauList(chain):
    return [kapTau(chain[i:i+4]) for i in range(len(chain)-3)]


def findSimilar(ktlist,tol):
    difList = [np.sqrt((ktlist[i][0]-ktlist[i+1][0])**2+(ktlist[i][1]-ktlist[i+1][1])**2) for i in range(len(ktlist)-1)]
    splitList = []
    splitListkt = []
    currsplit = [0]
    currsplitkt = [ktlist[0][0]]
    # find "helical geometries ", those whose curvature and torsion differs only a little
    for i in range(len(difList)):
        if difList[i]>tol:
            splitList.append(currsplit)
            splitListkt.append(currsplitkt)
            currsplit = [i+1]
            currsplitkt = [ktlist[i+1][0]]
        else:
            if(i==len(difList)-1):
                splitList.append(currsplit)
                splitListkt.append(currsplitkt)
            else: 
                currsplit.append(i+1)
                currsplitkt.append(ktlist[i+1][0])
            # collect Linker Sections
            
    sections = []
    sectionskt = []
    sublen =2
    linkerList = [splitListkt[0],splitListkt[0]]
    for i in range(len(splitList)):
        if len(splitList[i])==1:
            sublen = sublen +1
            linkerList.append(splitListkt[i])
        else:
            if sublen>0:
                sections.append(["-" for j in range(sublen)]) 
                sectionskt.append(linkerList)
            sections.append(["HL" for j in range(len(splitList[i]))])
            sectionskt.append(splitListkt[i])
            sublen = 0
    if sublen>0:
         sections.append(["-" for  j in range(sublen)])
         sectionskt.append(linkerList)
    for i in range(len(sections)):
        if sections[i][0]=="HL":            
            #print("here", sectionskt[i],np.mean(sectionskt[i]))
            if np.mean(sectionskt[i])>0.34:
                sections[i] = ["H" for j in range(len(sections[i]))]
            else:
                sections[i] = ["S" for j in range(len(sections[i]))]
    #upgrade HH to HHH 
    sectionsAlt =[]
    nextBlank =False
    for i in range(len(sections)):
        if (sections[i][0]=="H" and len(sections[i])==2 and i<(len(sections)-1)):
                if len(sections[i+1])>1:
                    sectionsAlt.append(["H" for j in range(3)])
                    sectionsAlt.append(sections[i+1][1:])
                    nextBlank=True;
                else:
                    sectionsAlt.append(["H" for j in range(3)])
                    nextBlank=True
        else:
            if nextBlank==False:
                sectionsAlt.append(sections[i])            
            else:
                nextBlank=False;
    #remove SS only, I dobn't rate them as real secondary structures  
    sections = sectionsAlt
    for i in range(len(sections)):
         if (sections[i][0]=="S" and len(sections[i])==2):
                sections[i] = ["-","-"]
    #finally recollect linkers 
    linkerList=[]
    secLen=0
    sectionsFin=[]
    for i in range(len(sections)-1):
        if sections[i][0]=="-" and sections[i+1][0]=="-":
            secLen = secLen + len(sections[i])
        elif sections[i][0]=="-" and sections[i+1][0]!="-":
            secLen = secLen + len(sections[i])
            sectionsFin.append(["-" for j in range(secLen)])
            secLen=0
        else:     
            sectionsFin.append(sections[i])
    #deal with the last element
    if  sections[len(sections)-1][0]=="-":
        secLen = secLen + len(sections[len(sections)-1])
        sectionsFin.append(["-" for j in range(secLen)])
    else:
        if secLen>0:
            sectionsFin.append(["-" for j in range(secLen)])
        sectionsFin.append(sections[len(sections)-1])
    #add last two as linkers
    sectionsFin = flatten(sectionsFin)
    sectionsFin.append("-")
    sectionsFin.append("-")
    return np.array(sectionsFin)


def geometric_secondary_structure(coords, threshold=0.2):

    kaptau = getKapTauList(coords)
    secondary = findSimilar(kaptau, threshold)
    return secondary


# Carbonara specific XYZ constraint setup

def translate_distance_constraints(contactPredsIn,coords,splitList,working_path,fixedDistList=[]):
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
        ss = get_secondary(working_path+"fingerPrint1.dat")
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
        contactPredNara.append(pair1+pair2+[dist]+[0.1])
        dists.append(dist)

        # now write to file
    print("distances in  current structure ",dists)
    np.savetxt(working_path+"fixedDistanceConstraints1.dat",contactPredNara,fmt="%i %i %i %i %1.10f %1.10f")


def findFlexibleSection(aminoIndexIn,working_path):
    # shift to [0,1,2 indexing
    aminoIndex = aminoIndexIn;
    currIndex=0;
    ss = get_secondary(working_path+"fingerPrint1.dat")
    sections= section_finder_sub(ss);
    currMax=len(sections[0])
    prevMax=0
    while aminoIndex>currMax:
        currIndex= currIndex+1
        currMax=currMax+len(sections[currIndex])
        prevMax = prevMax+len(sections[currIndex-1])
    # find the length of previous chains
    print(aminoIndex,sections[currIndex])
    return currIndex,(aminoIndex-prevMax)-1


# > CB inferrence stuff (not yet used too much)

def infer_CB_positions(CA_xyz):

    '''
    Returns CB positions (excluding end of chains)
    '''

    CA_vecs = np.diff(CA_xyz, axis=0)
    normals = np.diff(CA_vecs, axis=0)
    normals = normals/np.linalg.norm(normals, axis=1)[:,None]

    av_bond_len = 3.8

    normals = normals*av_bond_len

    CB_xyz = CA_xyz[1:-1] - normals # minus as to face outwards

    return CB_xyz


def interlace_CA_CB_write(CA_xyz, CB_xyz, protein, output_pdb):

    # positions of CB - wont add CB to the residues
    gly_idx = np.where(protein.data['resname'].values== 'GLY')[0]

    # atom names
    atom_name_lst = []

    # interlaced CA + CB coordinates
    coordinates = CA_xyz[0]

    # number of CA
    size = CA_xyz.shape[0]

    # creating index for resid
    idx_counter = 0
    idx_lst = []

    # new interlaced CA CB sequence
    new_seq = []

    # extract CA sequence
    ca_seq = protein.data['resname'].values

    for i in range(size):

        # ...CA bit...
        idx_lst.append(idx_counter)
        atom_name_lst.append('CA')
        new_seq.append(ca_seq[i])

        if i > 0:
                coordinates = np.vstack([coordinates, CA_xyz[i]])

        # ...CB bit...
        if i not in gly_idx:

            if (i > 0) & (i < size-1):
                idx_lst.append(idx_counter)
                atom_name_lst.append('CB')
                coordinates = np.vstack([coordinates, CB_xyz[i-1]])
                new_seq.append(ca_seq[i])

        idx_counter += 1


    tot_size = int(CA_xyz.shape[0]+CB_xyz.shape[0]-gly_idx.shape[0])

    if coordinates.shape[0] == tot_size:

        # create dataframe for biobox molecule type
        df = pd.DataFrame({'atom':['ATOM']*tot_size, 'index':np.arange(tot_size),
                           'name':atom_name_lst, 'resname':new_seq,
                           'chain':['A']*tot_size, 'resid':idx_lst,
                           'occupancy':[1]*tot_size, 'beta':[50]*tot_size,
                           'atomtype':['C']*tot_size, 'radius':[1.7]*tot_size,
                           'charge':[0]*tot_size})

    else:
        raise ValueError('Total number of CA + CB - no. GLY res does not equal the coordinate size!')


    molecule = bb.Molecule()
    molecule.data = df
    molecule.coordinates = np.expand_dims(coordinates, axis=0)

    molecule.write_pdb(output_pdb)


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
        result = subprocess.run(['./generate_structure', fingerprint_file, coords_file, outputname, str(l)], capture_output=True, text=True)

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
    conds = np.asarray(list(bond_breaks_dict.values())) < 1.0


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


# def setup_working_directory(run_name='Fitting'):

#     current = os.getcwd()
#     working = 'newFitData/'+run_name+'/'
#     working_path = os.path.join(current, working)

#     if os.path.exists(working_path):
#         shutil.rmtree(working_path)
#         print('Removing existing working directory')

#     os.makedirs(working_path)
#     os.mkdir(working_path+'/fitdata')

#     print('Complete')
#     return working_path


# def setup_molecule_directory(molecule_name='Test_Molecule', ignore_overwrite=True):

#     current = os.getcwd()
#     working = 'newFitData/'+molecule_name+'/'
#     working_path = os.path.join(current, working)
#     skip = False

#     if ignore_overwrite:
#         if os.path.exists(working_path):
#             shutil.rmtree(working_path)
#             print('Removing existing working directory')

#         os.makedirs(working_path)

#     else:
#         if os.path.exists(working_path):
#             skip = True
#         else:
#             os.makedirs(working_path)

#     # os.mkdir(working_path+'/fitdata')


    # return working_path

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

    saxs_arr = np.genfromtxt('temp.txt')

    if saxs_arr.shape[1] == 3:
        saxs_arr = saxs_arr[:,:2]

    #check if it is in angstoms, if the last value is >1 we assume its in nanometers.

    if saxs_arr[-1,0] >1:
        for i in range(0,len(saxs_arr)):
            saxs_arr[i,0]=saxs_arr[i,0]/10.0

    file_path_name = working_path+'/Saxs.dat'
    np.savetxt(file_path_name, saxs_arr, delimiter=' ', fmt='%s',newline='\n', header='', footer='')
    os.remove("temp.txt")

    return file_path_name


# BOTH WRITE SH METHODS ARE DEPRECIATED
def write_sh_file(working_path, fit_n_times, min_q, max_q, max_fit_steps):

    curr = os.getcwd()
    run_file = curr + '/RunMe.sh'

    with open(run_file, 'w+') as fout:
        fout.write('#!/bin/bash')


        saxs_file = working_path+'Saxs.dat'
        FP_file = working_path+"fingerPrint1.dat"
        coords_file = working_path+'coordinates1.dat'
        varying_file = working_path+"varyingSectionSecondary1.dat"
        mixture_file = working_path+"mixtureFile.dat"

        # Auto assign min / max q from SAXS profile
        # saxs_arr = np.genfromtxt(saxs_file)
        # min_q = np.round(saxs_arr[:,0].min(),2)
        # max_q = np.round(saxs_arr[:,0].max(),2)

        fout.write('\nfor i in {1..'+str(fit_n_times)+'}')

        fout.write('\n\ndo')
        fout.write('\n\n   echo " Run number : $i "')
        fout.write('\n\n   ./predictStructureQvary ' + saxs_file + ' ' + working_path+'/' + ' ' + coords_file + ' ' + 'none' + ' ' + varying_file + ' ' + '1' + ' ' + 'none' + \
                   ' ' + 'none' + ' ' + str(min_q) + ' ' + str(max_q) + ' ' + str(max_fit_steps) + ' ' + working_path+'/fitdata/fitmolecule$i' + ' ' + working_path+'/fitdata/scatter$i.dat' + ' ' + mixture_file + ' ' +'1')

        fout.write('\n\ndone')

    print('Successfully written bash script to: ', run_file)


def write_sh_qvary_file(working_path, mol_name, fit_name, fit_n_times, min_q, max_q,
                        max_fit_steps, pairedQ=False,rotation=False):

    curr = os.getcwd()
    script_name = '/RunMe_'+ mol_name + '_' + fit_name + '.sh'
    run_file = curr + script_name

    with open(run_file, 'w+') as fout:
        fout.write('#!/bin/bash')

        # argv[ 1] scattering data file
        saxs_file = working_path+'Saxs.dat'
        # argv[ 2] sequence file location
        FP_file = working_path+"fingerPrint1.dat"
        # argv[ 3] restart tag (use to start from existing prediction)
        coords_file = working_path+'coordinates1.dat'
        # argv[ 4] paired distances file (can be empty)

        # argv[ 5] fixed sections file (again can be empty)
        # argv[ 6] number of structures
        # argv[ 7] request to apply hydrophobic covering WITHIN monomers will be a list of sections on which to apply it. Will say none if not. -- Currently not used
        # argv[ 8] request to apply hydrophobic covering BETWEEN monomers will be a list of pairs to try to hydropobically pair. Will say none if not. -- currently not used
        # argv[ 9] kmin
        # argv[10] kmax
        # argv[11] Max number of fitting steps
        # argv[12] prediction file
        # argv[13] scattering output file
        # argv[14] mixture list file, alist of sets of numbers indicatig the allowed set of mixture percentages of each species (e.g. dimer 20 monomer 80)
        # argv[15] previous fit string in form fitname/mol6Substep_10_1.dat+fitname/mol6Substep_10_2.dat
        # argv[16] log file location
        # argv[17] last line of the previous fit log, this is only used for a restart if argv[3] = True
        # argv[18] is true if we want to apply affine rotations,false if not.




        varying_file = working_path+"varyingSectionSecondary1.dat"
        mixture_file = working_path+"mixtureFile.dat"

        # $1 is the moelcule name $2 the output file $3 the restart file if sied
        fout.write('\n ScatterFile='+saxs_file)
        fout.write('\n fileLocs='+working_path)
        fout.write('\n initialCoordsFile=frompdb')
        fout.write('\n noStructures=1')
        if pairedQ==False:
            fout.write('\n pairedPredictions=False')
        else:
            fout.write('\n pairedPredictions=True')

        fout.write('\n fixedsections='+working_path+'varyingSectionSecondary1.dat')
        fout.write('\n withinMonomerHydroCover=none')
        fout.write('\n betweenMonomerHydroCover=none')
        fout.write('\n kmin='+str(min_q))
        fout.write('\n kmax='+str(max_q))
        fout.write('\n maxNoFitSteps='+str(max_fit_steps))
        if rotation==False:
            fout.write('\n affineTrans=False')
        else:
            fout.write('\n affineTrans=True')

        # Auto assign min / max q from SAXS profile
        # saxs_arr = np.genfromtxt(saxs_file)
        # min_q = np.round(saxs_arr[:,0].min(),2)
        # max_q = np.round(saxs_arr[:,0].max(),2)

        fout.write('\nfor i in {1..'+str(fit_n_times)+'}')

        fout.write('\n\ndo')
        fout.write('\n\n   echo " Run number : $i "')
        fout.write('\n\n   ./predictStructureQvary $ScatterFile $fileLocs $initialCoordsFile $pairedPredictions $fixedsections $noStructures $withinMonomerHydroCover $betweenMonomerHydroCover $kmin $kmax $maxNoFitSteps '+working_path+fit_name+'/mol$i '+working_path+fit_name+'/scatter$i.dat '+working_path+'mixtureFile.dat '+working_path+'redundant '+working_path+fit_name+'/fitLog$i.dat '+'null '+'$affineTrans')

        fout.write('\n\ndone')

    # return script_name[1:]
    # print('Successfully written bash script to: ', run_file)
        fout.close()
        # now write the run folder if its empty

        working = working_path+fit_name
        if os.path.exists(working):
            skip = True
        else:
            os.makedirs(working)

    # os.mkdir(working_path+'/fitdata')


# ------ SAXS Profile Visulation ------

def SAXS_selection_plotter(SAXS_file, min_q, max_q):

    SAXS = np.genfromtxt(SAXS_file)

    q = SAXS[:,0]
    I = SAXS[:,1]

    q_selected = q[(q>=min_q)&(q<=max_q)]
    q_grey = q[(q<min_q) | (q>=max_q)]

    I_selected = I[(q>=min_q)&(q<=max_q)]
    I_grey = I[(q<min_q) | (q>=max_q)]

    fig = go.Figure( data=[go.Scatter(x=q_grey, y=I_grey, mode='markers', line=dict(color="grey"), opacity=0.7, name='Excluded')])
    fig.add_trace( go.Scatter(x=q_selected, y=I_selected, mode='markers', line=dict(color="crimson"), name='Selected') )
    fig.update_layout(
                    # title='Selected q range',
                     yaxis_type = "log", template='plotly_white',
                    width=800, height=700, font_size=28)
    fig.update_xaxes(title='q')
    fig.update_yaxes(title='I')
    fig.update_layout(margin=dict(l=0, r=0, t=0, b=0))


    return fig


def get_minmax_q(SAXS_file):

    SAXS = np.genfromtxt(SAXS_file)

    q = SAXS[:,0]

    q_exp_min = float(np.round(q.min(),2))
    q_exp_max = float(np.round(q.max(),2))

    q_spread = q_exp_max - q_exp_min

    q_Q1 = float(np.round(0.00*q_spread,2))
    q_Q3 = float(np.round(0.45*q_spread,2))

    return q_exp_min, q_exp_max, q_Q1, q_Q3


def SAXS_fit_plotter(SAXS_file, fit_file, full_q=True):

    fig = make_subplots(rows=2, cols=1,row_heights=[0.7,0.3],vertical_spacing=0,shared_xaxes=True)

    SAXS = np.genfromtxt(SAXS_file)

    fitting = np.genfromtxt(fit_file, skip_footer=1)
    fit_q = fitting[:,0]
    fit_I = fitting[:,2]

    q = SAXS[:,0]
    I = np.log(SAXS[:,1])

    min_q = fit_q.min()
    max_q = fit_q.max()



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


# ------ Protein Chain Visulation ------

def smooth_me(coords, ss, oversample=5):

    color_dic = {'-':'yellow', 'H':'red', 'S':'blue'}

    structure_change = np.diff(np.unique(ss, return_inverse=True)[1])
    sections = section_finder(ss)
    structural_groups = get_section_groupings(ss, structure_change)

    color_lst = []

    x_lst = []
    y_lst = []
    z_lst = []

    for i, sec in enumerate(sections):


        tmp = coords[np.where(structural_groups==i)]

        if tmp.shape[0]>3:

            if sec == 'H':
                tck, u = interpolate.splprep([tmp[:,0], tmp[:,1], tmp[:,2]], s=.25)

            else:
                tck, u = interpolate.splprep([tmp[:,0], tmp[:,1], tmp[:,2]], s=10)
            res_size = oversample*tmp.shape[0]
            u_fine = np.linspace(0,1,res_size)
            new_points = interpolate.splev(u_fine, tck)

            x_lst = x_lst + list(new_points[0])
            y_lst = y_lst + list(new_points[1])
            z_lst = z_lst + list(new_points[2])

            color_lst = color_lst + [color_dic[sec]] * res_size

        else:

            x_lst = x_lst + list(tmp[:,0])
            y_lst = y_lst + list(tmp[:,1])
            z_lst = z_lst + list(tmp[:,2])

            color_lst = color_lst + [color_dic[sec]] * tmp.shape[0]

    return x_lst, y_lst, z_lst, color_lst


def smooth_me_varying(coords, ss, vary_sections, oversample=5):

    structure_change = np.diff(np.unique(ss, return_inverse=True)[1])
    sections = section_finder(ss)
    structural_groups = get_section_groupings(ss, structure_change)

    color_lst = []

    x_lst = []
    y_lst = []
    z_lst = []

    for i, sec in enumerate(sections):


        tmp = coords[np.where(structural_groups==i)]

        if tmp.shape[0]>3:

            if sec == 'H':
                tck, u = interpolate.splprep([tmp[:,0], tmp[:,1], tmp[:,2]], s=1)

            else:
                tck, u = interpolate.splprep([tmp[:,0], tmp[:,1], tmp[:,2]], s=100)
            res_size = oversample*tmp.shape[0]
            u_fine = np.linspace(0,1,res_size)
            new_points = interpolate.splev(u_fine, tck)

            x_lst = x_lst + list(new_points[0])
            y_lst = y_lst + list(new_points[1])
            z_lst = z_lst + list(new_points[2])

            if i in vary_sections:
                color_lst = color_lst + ['red'] * res_size
            else:
                color_lst = color_lst + ['grey'] * res_size
        else:

            x_lst = x_lst + list(tmp[:,0])
            y_lst = y_lst + list(tmp[:,1])
            z_lst = z_lst + list(tmp[:,2])

            if i in vary_sections:
                color_lst = color_lst + ['red'] * tmp.shape[0]
            else:
                color_lst = color_lst + ['grey'] * tmp.shape[0]


    return x_lst, y_lst, z_lst, color_lst


def line_plotly(x_lst, y_lst, z_lst, color_lst, outline=False):

    fig = go.Figure()

    if outline:
        fig.add_trace(go.Scatter3d(x=x_lst, y=y_lst, z=z_lst,
                                mode='lines',
                                line=dict( width=21, color='black')))

    fig.add_trace(go.Scatter3d(x=x_lst, y=y_lst, z=z_lst,
                               mode='lines',
                               line=dict( width=12, color=color_lst)))

    fig['layout']['showlegend'] = False
    fig.update_layout( scene=dict(
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

    # fig.update_layout(width=800, height=800)
    fig.update_layout(margin=dict(l=0, r=0, t=0, b=0))
    fig.update_traces(showlegend=False)
    return fig


# Fancy tube plotting things

def data_for_cylinder_along_arb(center,tandirec,height_z,radius=0.8):
    z = np.linspace(0, height_z, int(height_z))
    theta = np.linspace(0, 2*np.pi,25)
    theta_grid, z_grid=np.meshgrid(theta, z)
    tperp = np.sqrt(tandirec[0]*tandirec[0] + tandirec[1]*tandirec[1])
    tperp2 =np.sqrt(4.0*tandirec[0]*tandirec[0]*tandirec[1]*tandirec[1] + tperp*tperp*tandirec[2]*tandirec[2])
    x_grid = -radius*tandirec[0]*np.cos(theta_grid)/tperp - radius*tandirec[1]*tandirec[2]*np.sin(theta_grid)/tperp2+ tandirec[0]*z_grid + center[0]
    y_grid = radius*tandirec[1]*np.cos(theta_grid)/tperp - radius*tandirec[0]*tandirec[2]*np.sin(theta_grid)/tperp2+ tandirec[1]*z_grid + center[1]
    z_grid = 2.0*radius*tandirec[0]*tandirec[1]*np.sin(theta_grid)/tperp2+ tandirec[2]*z_grid + center[2]
    return x_grid,y_grid,z_grid


def plot_molecule_tube(fl_path):
    times = 25
    mol = np.genfromtxt(fl_path)
    Xc,Yc,Zc = data_for_cylinder_along_arb(mol[0],(mol[1]-mol[0])/np.linalg.norm(mol[1]-mol[0]),np.linalg.norm(mol[1]-mol[0]))
    for i in range(2,len(mol)):
        normalised_tangent =  (mol[i]-mol[i-1])/np.linalg.norm(mol[i]-mol[i-1])
        prev_normalised_tangent = (mol[i-1]-mol[i-2])/np.linalg.norm(mol[i-1]-mol[i-2])
        for j in range(1,times):
            cap_data = data_for_cylinder_along_arb(mol[i-1],
                                                   ((times-j)/(times-1))*prev_normalised_tangent + (j/(times-1))*normalised_tangent,
                                                   1,
                                                  )
            Xc = np.concatenate((Xc,cap_data[0]),axis=0)
            Yc = np.concatenate((Yc,cap_data[1]),axis=0)
            Zc = np.concatenate((Zc,cap_data[2]),axis=0)
        edge_data = data_for_cylinder_along_arb(mol[i-1],
                                                normalised_tangent,
                                                int(np.linalg.norm(mol[i]-mol[i-1]))
                                               )
        Xc = np.concatenate((Xc,edge_data[0]),axis=0)
        Yc = np.concatenate((Yc,edge_data[1]),axis=0)
        Zc = np.concatenate((Zc,edge_data[2]),axis=0)

    surcol = np.zeros((Zc-Zc[0]).shape)
    for i in range(len(surcol)):
        surcol[i] = i*np.ones(surcol[i].shape)
    lighting_effects = dict(ambient=0.5, diffuse=0.25, roughness = 0.9, specular=0.65, fresnel=1)
    fig = go.Figure()
    fig.add_trace(go.Surface(x=Xc, y=Yc, z=Zc,
                             showscale=False,
                             surfacecolor=surcol,
                             colorscale='Rainbow',
                             lighting=dict(lighting_effects),
                             lightposition=dict(x=1000,y=1000,z=1000),
                             hoverinfo='none'
                            )
                 )
    colorbar_trace = go.Scatter3d(x=[None],
                                  y=[None],
                                  z=[None],
                                  mode='markers',
                                  marker=dict(
                                      colorscale='Rainbow',
                                      showscale=True,
                                      cmin=surcol[0][0],
                                      cmax=surcol[-1][0],
                                      colorbar=dict(thickness=25,
                                                    tickvals=[surcol[0][0],surcol[-1][0]],
                                                    ticktext=['Start','End'],
                                                    outlinewidth=0)
                                  ),
                                  hoverinfo='none'
                                 )
    fig.add_trace(colorbar_trace)
    fig.update_layout(width=1000,height=1000)
    fig.update_layout(
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
    fig.show()


###--------------------------------------------------------------------------------------------------------####
###                                                                                                        ####
###                                 For Josh - Arron's Additions                                           ####
###                                                                                                        ####
###--------------------------------------------------------------------------------------------------------####

# Post Prediction stuff

# Dealing with carbonara Log file outputs

# DEPRCIATED LOGGING - CHRIS COMPATABLE

def old_log2df(log_file):
    arr = np.loadtxt(log_file, str, skiprows=4)

    columns = ['index', 'step_num', 'total penalty', 'writhe_penalty', 'overlap_penalty', 'contact_penalty', 'time', 'max_fitting_k', 'file_name']
    df_log = pd.DataFrame(arr, columns=columns)

    for c in columns[:-1]:
        df_log[c] = df_log[c].astype(float)


    df_log['time'] = df_log['time'].astype(float) * 1/60 * 10**(-6)


    df_log['fit'] = df_log['total penalty'] -df_log['writhe_penalty'] - df_log['overlap_penalty'] - df_log['contact_penalty']

    return df_log


def old_df2plot(df_log,dict_type, highlight=False):

    fig = go.Figure()


    fig.add_trace(go.Scatter(x=df_log['time'].values,
                             y=df_log[dict_type].values,
                             name=dict_type,
                             mode='lines + markers',
                             line=dict(color="lightblue", width=3) ))

    if highlight:
        fig.add_trace(go.Scatter(x=[df_log['time'].values[highlight]],
                                 y=[df_log[dict_type].values[highlight]],
                                 name='selected',
                                 mode='markers',
                                 marker=dict( symbol="hexagon", color="crimson", size=10
                                            ) ))

    fig.update_layout(template='simple_white')
    # height=600, width=1400
    fig.update_layout(xaxis_title="Time (Minutes)", yaxis_title=dict_type)
    fig.update_layout(margin=dict(l=0, r=0, t=0, b=0))
    fig.update_traces(showlegend=False)
    # fig.show()
    return fig

# NEW LOGGING

def getLogs(MolPath,RunName):
    return glob(MolPath+RunName+'/fitLog*')


def LogFile2df(logFilePath):
    '''
    Reads in Josh's new json format LogFile, into a pandas dataframe.
    Skips the first row as this is just run metadata.
    '''
    lines=[]
    with open(logFilePath) as f:
        for i, line in enumerate(f):
            if i<1:
                continue
            obj = json.loads(line)
            lines.append(obj)
    df = pd.DataFrame(lines)
    return df


def getAcceptableFits(df,cutoff=0.005):
    '''
    Trim Log dataframe to include scatter fits with Chi^2 <= cutoff (usually 0.005).
    '''
    acceptable_df = df[df['ScatterFitFirst']<=cutoff].reset_index()
    return acceptable_df


def viewBestSAXSFit(RunPath,LogFilePath):
    '''
    Plots the current best scattering fit (<= cutoff and largest qMax).
    '''
    df = getAcceptableFits(LogFile2df(LogFilePath))
    print('Best Chi^2 Fit is',df.tail(1)['ScatterFitFirst'].values[0],'for q <=',df.tail(1)['KmaxCurr'].values[0])
    return SAXS_fit_plotter(RunPath+"Saxs.dat",df['ScatterPath'].values[-1],full_q=False)


def viewBestMolChange(RunPath,LogFilePath):
    '''
    Plots the backbone for the current best fit, aligned to the start coords.
    Highlights with the largest change from the start.
    '''
    df = getAcceptableFits(LogFile2df(LogFilePath))
    current_coords = np.genfromtxt(df.tail(1)['MoleculePath'].values[0],skip_footer=1)
    start_coords = np.genfromtxt(RunPath+'coordinates1.dat')
    aligned = rmsd.kabsch_fit(current_coords,start_coords)
    diff_coords = np.array([np.linalg.norm(start_coords[i]-aligned[i]) for i in range(len(start_coords))])
    norm_diff = diff_coords/max(diff_coords)
    fig = go.Figure()
    fig.add_trace(go.Scatter3d(
        x=start_coords[:,0],
        y=start_coords[:,1],
        z=start_coords[:,2],
        name='Start',
        visible='legendonly',
        marker=dict(
            size=1,
            color='red',
        ),
        line=dict(
            color='red',
            width=10
        )
    ))
    fig.add_trace(go.Scatter3d(
        x=aligned[:,0],
        y=aligned[:,1],
        z=aligned[:,2],
        name='Change',
        marker=dict(
            size=1,
            color=norm_diff,
            colorscale='jet',
            colorbar=dict(thickness=20),
        ),
        line=dict(
            color=norm_diff,
            colorscale='jet',
            colorbar=dict(thickness=20),
            width=10
        )
    ))
    fig.add_trace(go.Scatter3d(
        x=aligned[:,0],
        y=aligned[:,1],
        z=aligned[:,2],
        name='Best',
        visible='legendonly',
        marker=dict(
            size=1,
            color='blue',
        ),
        line=dict(
            color='blue',
            width=10
        )
    ))
    fig.update_layout(width=1000,height=1000)
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
    fig.show()


def CollectBestOutputs(MolPath,RunName):
    '''
    Finds the best prediction (chi^2<=0.005 and most qMax) from each run.
    Returns them as a (# of runs, len(CA), 3) tensor.
    '''
    logs = getLogs(MolPath,RunName)
    all_coords_tensor = []
    for log in logs:
        df = getAcceptableFits(LogFile2df(log))
        best_coords = np.genfromtxt(df.tail(1)['MoleculePath'].values[0],skip_footer=1)
        all_coords_tensor.append(best_coords)
    all_coords_tensor = np.array(all_coords_tensor)
    return all_coords_tensor


def CollectBestOutputNames(MolPath,RunName):
    '''
    Collects the file names for the best prediction from each run.
    '''
    logs = getLogs(MolPath,RunName)
    best_output_names = []
    for log in logs:
        df = getAcceptableFits(LogFile2df(log))
        best_output_names.append(df.tail(1)['MoleculePath'].values[0].split('/')[-1][:-8])
    return best_output_names


def AlignBestOutputs(MolPath,RunName):
    '''
    Aligns the best prediction (chi^2<=0.005 and most qMax) from each run
    Returns them as a (# of runs, len(CA), 3) tensor.
    '''
    best_coords_tensor = CollectBestOutputs(MolPath,RunName)
    aligned_coords_tensor = np.zeros_like(best_coords_tensor)
    aligned_coords_tensor[0,:,:] = best_coords_tensor[0,:,:]
    size = best_coords_tensor.shape[0]
    for i in range(1,size):
        aligned_coords_tensor[i,:,:] = rmsd.kabsch_fit(best_coords_tensor[i,:,:],best_coords_tensor[0,:,:])
    return aligned_coords_tensor


def BestOutputRMSD(MolPath,RunName):
    '''
    Computes the pairwise RMSD matrix for the best predictions from each run.
    Returns a (# runs, #runs) array.
    '''
    aligned_best_coords = AlignBestOutputs(MolPath,RunName)
    size = aligned_best_coords.shape[0]
    rmsd_arr = np.zeros((size, size))
    for i in range(size):
        for j in range(i):
            rmsd_arr[i,j] = rmsd_arr[j,i] = rmsd.kabsch_rmsd(aligned_best_coords[i,:,:],aligned_best_coords[j,:,:])
    return rmsd_arr


def PlotBestOutputRMSD(MolPath,RunName):
    '''
    Plots the pairwise RMSD matrix for the best predictions from each run as a heatmap.
    '''
    names = CollectBestOutputNames(MolPath,RunName)
    rmsd_arr = BestOutputRMSD(MolPath,RunName)
    fig = go.Figure(data=go.Heatmap(z=rmsd_arr,x=names,y=names))
    fig.update_layout(template='simple_white')
    fig.update_layout(showlegend = False,
                    width = 750,
                    height = 750,
                    autosize = False
                    )
    fig.show()

# > RMSD 
 
def fit_rms(ref_c,c):

    # move geometric center to the origin
    ref_trans = np.average(ref_c, axis=0)
    ref_c = ref_c - ref_trans
    c_trans = np.average(c, axis=0)
    c = c - c_trans

    # covariance matrix
    C = np.dot(c.T, ref_c)

    # Singular Value Decomposition
    (r1, s, r2) = np.linalg.svd(C)

    # compute sign (remove mirroring)

    if np.linalg.det(C) < 0:
        r2[2,:] *= -1.0
    U = np.dot(r1, r2)

    return (c_trans, U, ref_trans)


def find_rmsd(c1, c2, align_return=False):

    rmsd = 0.0
    c_trans, U, ref_trans = fit_rms(c1, c2)

    new_c2 = np.dot(c2 - c_trans, U) + ref_trans

    rmsd = np.sqrt( np.average( np.sum( ( c1 - new_c2 )**2, axis=1 ) ) )

    if align_return:
        return rmsd, new_c2

    else:
        return rmsd

def align_coords(tensor):

    new_coord_tensor = np.zeros_like(tensor)
    new_coord_tensor[:,:,0] = tensor[:,:,0]
    size = tensor.shape[-1]

    for i in range(1,size):
        rmsd, align_coords = find_rmsd(tensor[:,:,0], tensor[:,:,i], align_return=True)
        new_coord_tensor[:,:,i] = align_coords


    return new_coord_tensor


def coord_tensor_pairwise_rmsd(coord_tensor):

    size = coord_tensor.shape[-1]
    rmsd_arr = np.zeros((size, size))
    # size = 10
    for i in range(size):
        for j in range(i):
            rmsd = find_rmsd(coord_tensor[:,:,i], coord_tensor[:,:,j])
            rmsd_arr[i,j] = rmsd_arr[j,i] = rmsd

    return rmsd_arr


# > older structural comparision

def overlay_coords(tensor):

    fig = go.Figure()

    for i in range(tensor.shape[-1]):
#         print(i)
        fig.add_trace(go.Scatter3d(x=tensor[:,0,i], y=tensor[:,1,i], z=tensor[:,2,i], opacity=0.4,
                               mode='lines',
                               line=dict( width=12,
#                                          color='red',
#                                          color=np.arange(tensor[:,0,i].shape[0]),
                                         colorscale='greys')))


    av_coords = tensor.mean(axis=2)

    fig.add_trace(go.Scatter3d(x=av_coords[:,0], y=av_coords[:,1], z=av_coords[:,2], opacity=0.7,
                               mode='lines',
                               line=dict( width=12,
                                         color='black',
#                                          color=np.arange(tensor[:,0,i].shape[0]),
                                         colorscale='greys')))

    fig['layout']['showlegend'] = False
    fig.update_layout( scene=dict(
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

    fig.update_layout(width=800, height=800)
    fig.update_layout(margin=dict(l=0, r=0, t=0, b=0))
    fig.update_traces(showlegend=False)
    return fig


# Clustering methods + visulisations

#def cluster(rmsd_arr, min_cluster_size=8):
#    clusterer = hdbscan.HDBSCAN(min_cluster_size=min_cluster_size, gen_min_span_tree=True)
#    clusterer.fit(rmsd_arr)
#
#    return clusterer.labels_, clusterer.probabilities_


def visualise_clusters(coord_tensor, labels, best_fit_files):


    test = np.dstack([  coord_tensor[:,:,labels==0], coord_tensor[:,:,labels==1],
                        coord_tensor[:,:,labels==2], coord_tensor[:,:,labels==3],
                        coord_tensor[:,:,labels==-1]])

    rmsd_arr_t = coord_tensor_pairwise_rmsd(test)

    cluster_cumsum = np.cumsum( np.asarray( [sum(labels==0), sum(labels==1), sum(labels==2),

                                         sum(labels==3), sum(labels==-1)] ) )

    cluster_cumsum = np.insert(cluster_cumsum, 0, 0)

    bf_names = []
    for bf in best_fit_files:
        bf_names.append(bf.split('/')[-2]+'/'+bf.split('/')[-1][:-4])

    bf_names = np.asarray(bf_names)
    bf_names_sort = np.concatenate([bf_names[labels==0], bf_names[labels==1], bf_names[labels==2], bf_names[labels==3], bf_names[labels==-1]])

    hovertext = list()
    for yi, yy in enumerate(bf_names_sort):
        hovertext.append(list())
        for xi, xx in enumerate(bf_names_sort):
            hovertext[-1].append('Structure X: {}<br />Structure Y: {}<br /> RMSD: {}'.format(xx, yy, round(rmsd_arr_t[xi, yi],1)))

    fig = go.Figure()
    fig.add_trace(go.Heatmap(z=rmsd_arr_t, colorscale='RdBu_r', text=hovertext, hoverinfo='text'))

    for v in cluster_cumsum:

        fig.add_hline(y=v-0.5, line_width=2, line_color="black")
        fig.add_vline(x=v-0.5, line_width=2, line_color="black")

    fig.update_layout(height=700, width=650)
    return fig, bf_names_sort


def cluster(rmsd_arr):
    '''
    Updated your cluster, think allow_single_cluster is fine for us (if we make good predictions then this could well happen)
    '''
    clusterer = hdbscan.HDBSCAN(metric = 'precomputed',min_cluster_size=3,allow_single_cluster=True)
    clusterer.fit(rmsd_arr)
    return clusterer.labels_, clusterer.probabilities_


# Arron SKMT structure analysis bits

def get_all_edges(curve):
    '''
    As the function name suggestions
    '''
    edges=[]
    for i in range(1,len(curve)):
        edges.append([curve[i-1],curve[i]])
    return edges


def get_sses(ss_file):
    '''
    From a SS fingerprint, returns the SS with the subsection length, eg. ---SSS--- becomes [[-,3],[S,3],[-,3]]]
    '''
    lines = []
    with open(ss_file,'r') as fin:
        for line in fin:
            lines+= [line.split()]
    ss = lines[-1][0]
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
    return sses


def intersect_line_triangle(q1,q2,p1,p2,p3):
    '''
    Checks if an edge q1-q2 intersects a triangle p1-p2-p3
    Used in the SKMT algorithm
    '''
    def signed_tetra_volume(a,b,c,d):
        return np.sign(np.dot(np.cross(b-a,c-a),d-a)/6.0)

    s1 = signed_tetra_volume(q1,p1,p2,p3)
    s2 = signed_tetra_volume(q2,p1,p2,p3)
    if s1 != s2:
        s3 = signed_tetra_volume(q1,q2,p1,p2)
        s4 = signed_tetra_volume(q1,q2,p2,p3)
        s5 = signed_tetra_volume(q1,q2,p3,p1)
        if s3 == s4 and s4 == s5:
            return True
    return False


def skmt(coords,fp_file):
    '''
    Given a CA backbone (coords) and SS fingerprint (fp_file), returns a minimal representation of the
    entanglement of the backbone.
    '''
    ss = get_sses(fp_file)
    splitcurve = []
    index = 0
    for i in ss:
        splitcurve.append(coords[index:index+i[1]])
        index+=i[1]
    newcurve = []
    for i in range(len(splitcurve)):
        for j in range(len(splitcurve[i])):
            newcurve.append(splitcurve[i][j])
    for subsec in range(len(splitcurve)):
        if len(splitcurve[subsec])>2:
            checks = []
            for idx in range(1,len(splitcurve[subsec])-1):
                p1 = 1.25*splitcurve[subsec][0]-splitcurve[subsec][1]
                p2 = splitcurve[subsec][idx]
                p3 = 1.25*splitcurve[subsec][-1]-splitcurve[subsec][-2]
                for edge in get_all_edges(newcurve):
                    q0 = edge[0]
                    q1 = edge[1]
                    checks.append(intersect_line_triangle(q0,q1,p1,p2,p3))
            if not any(checks):
                splitcurve[subsec] = [splitcurve[subsec][0]]
                newcurve = []
                for l in range(len(splitcurve)):
                    for m in range(len(splitcurve[l])):
                        newcurve.append(splitcurve[l][m])
            else:
                idx=2
                while idx<len(splitcurve[subsec]):
                    newcurve = []
                    for i in range(len(splitcurve)):
                        for j in range(len(splitcurve[i])):
                            newcurve.append(splitcurve[i][j])
                    p1 = splitcurve[subsec][idx-2]
                    p2 = splitcurve[subsec][idx-1]
                    p3 = splitcurve[subsec][idx]
                    checks = []
                    for edge in get_all_edges(newcurve):
                        q0 = edge[0]
                        q1 = edge[1]
                        checks.append(intersect_line_triangle(q0,q1,p1,p2,p3))
                    if not any(checks):
                        splitcurve[subsec] = np.delete(splitcurve[subsec],idx-1,axis=0)
                        idx=2
                    else:
                        idx+=1
        else:
            splitcurve[subsec] = [splitcurve[subsec][0]]
            newcurve = []
            for l in range(len(splitcurve)):
                for m in range(len(splitcurve[l])):
                    newcurve.append(splitcurve[l][m])
    newcurve = []
    for i in range(len(splitcurve)):
        for j in range(len(splitcurve[i])):
            newcurve.append(splitcurve[i][j])
    if not np.array_equal(newcurve[-1],coords[-1]):
        newcurve.append(coords[-1])
    return newcurve
    
def write_curve_to_file(curve,outfile_name):
    '''
    As the function name suggests.
    '''
    with open(outfile_name,'w+') as f:
        for i in range(len(curve)-1):
            string = ' '.join(map(str,curve[i]))
            f.write(string)
            f.write('\n')
        f.write(' '.join(map(str,curve[-1])))
        f.close()


def writeBestSKMTCurves(MolPath,RunName):
    '''
    Finds the best predictions from each run, computes their SKMT curve, saves them in a tmp dir.
    '''
    best_coords_tensor = CollectBestOutputs(MolPath,RunName)
    for i in range(best_coords_tensor.shape[0]):
        curve = skmt(best_coords_tensor[i,:,:],MolPath+'fingerPrint1.dat')
        if not os.path.exists(MolPath+RunName+'/tmp'):
            os.makedirs(MolPath+RunName+'/tmp')
        write_curve_to_file(curve,MolPath+RunName+'/tmp/'+str(i+1)+'.xyz')


def extract_nums(text):
    '''
    Chris Logic...
    '''
    for item in text.split(' '):
        try:
            yield float(item)
        except ValueError:
            pass


def toPairs(compareLst):
    '''
    Chris Logic...
    '''
    outLst = []
    for i in range(0,len(compareLst)-2,4):
        outLst.append([intLst(compareLst[i:i+2]),intLst(compareLst[i+2:i+4])])
    outLst.append([compareLst[-2],compareLst[-1]])
    return outLst


def intLst(list_of_floats):
    return [int(item) for item in list_of_floats]


def compare_molecules(fl1_name,fl2_name,cutOff=0.1):
    '''
    Given two xyz coordinate files, computes their writhe fingerprints and finds the largest similar subsections
    See SKMT paper for details...
    '''
    get_data=subprocess.check_output("compareFingerprints "+fl1_name+" "+fl2_name+" " + str(cutOff),shell=True, encoding='utf-8')
    return toPairs(list(extract_nums(get_data)))


def compareBestSKMTCurves(MolPath,RunName):
    '''
    For the best predictions form each run, performs the above writhe comparison
    Returns a (# runs, # runs) matrix with (i,j) the percentage deviation.
    '''
    if not os.path.exists(MolPath+RunName+'/tmp'):
        writeBestSKMTCurves(MolPath,RunName)
    curveFiles = glob(MolPath+RunName+'/tmp/*')
    size = len(curveFiles)
    comp_arr = np.zeros((size, size))
    for i in range(size):
        for j in range(i):
            comp = compare_molecules(curveFiles[i],curveFiles[j])[-1]
            comp_arr[i][j] = comp_arr[j][i] = 1-comp[0]
    return comp_arr


def plotBestSKMTComp(MolPath,RunName):
    '''
    Plots the pairwise difference matrix as a heatmap.
    '''
    names = CollectBestOutputNames(MolPath,RunName)
    comp_arr = compareBestSKMTCurves(MolPath,RunName)
    fig = go.Figure(data=go.Heatmap(z=comp_arr,x=names,y=names))
    fig.update_layout(template='simple_white')
    fig.update_layout(showlegend = False,
                    width = 750,
                    height = 750,
                    autosize = False
                    )
    fig.show()


def overlayCluster(MolPath,RunName,label):
    '''
    Performs HDBScan clustering based on the writhe similarity measure of SKMT curves
    Overlays the backbone coords of cluser "label"
    '''
    comp_arr = compareBestSKMTCurves(MolPath,RunName)
    labels,_ = cluster(comp_arr)
    best_coords_tensor = AlignBestOutputs(MolPath,RunName)
    cluster_coords_tensor = best_coords_tensor[labels==label,:,:]
    fig = go.Figure()
    for i in range(cluster_coords_tensor.shape[0]):
        fig.add_trace(go.Scatter3d(
            x=cluster_coords_tensor[i,:,0],
            y=cluster_coords_tensor[i,:,1],
            z=cluster_coords_tensor[i,:,2],
            opacity=0.4,
            mode='lines',
            line=dict(
                width=12,
                colorscale='greys'
                )
            )
        )
    av_coords = cluster_coords_tensor.mean(axis=0)
    fig.add_trace(go.Scatter3d(
        x=av_coords[:,0],
        y=av_coords[:,1],
        z=av_coords[:,2],
        opacity=0.7,
        mode='lines',
        line=dict(
            width=12,
            color='black',
            colorscale='greys'
            )
        )
    )
    fig['layout']['showlegend'] = False
    fig.update_layout(scene=dict(
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
            showticklabels=False)
            )
        )

    fig.update_layout(width=800, height=800)
    fig.update_layout(margin=dict(l=0, r=0, t=0, b=0))
    fig.update_traces(showlegend=False)
    return fig

def toggle_paired_predictions(script_path):
    with open(script_path, 'r') as f:
        lines = f.readlines()

    with open(script_path, 'w') as f:
        for line in lines:
            if line.strip().startswith("pairedPredictions="):
                line = "pairedPredictions=True\n"
            f.write(line)
