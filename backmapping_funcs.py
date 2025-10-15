import os
import glob
import re
import sys

from tempfile import mkstemp
from os.path import basename

import numpy as np
from modeller import *
from modeller.automodel import *
from modeller.scripts import complete_pdb


from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.Polypeptide import is_aa

import CarbonaraDataTools as cdt
import numpy as np
import json

def renumber_pdb_chains_start_from_1(input_pdb, output_pdb):
    """
    Renumber each chain in a PDB file so residues start from 1.

    Parameters:
        input_pdb (str): Path to the original PDB file
        output_pdb (str): Path to write the renumbered PDB
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", input_pdb)

    for model in structure:
        for chain in model:
            new_resnum = 1
            for res in chain:
                if is_aa(res, standard=True) or res.id[0] == ' ':
                    res.id = (' ', new_resnum, ' ')
                    new_resnum += 1

    io = PDBIO()
    io.set_structure(structure)
    io.save(output_pdb)

def CA2AA(filename, outputname, iterations=1, stout=False):
    # output=None

    # 

    # template for pir file format (shitty format but modeller likes it and I want him to be happy)
    _PIR_TEMPLATE = '\n'.join(['>P1;%s', 'sequence:::::::::', '%s', '*', '', '>P1;model_ca', 'structure:%s:FIRST:@:END:@::::', '*'])

   
    # printing is slowing us down/bogging up the output - preventing output
    old_stdout = sys.stdout
    if stout != True:
        sys.stdout = open(os.devnull, 'w')

    # need a name for a loads of temporary files - naming doesn't matter as they'll be deleted later
    pdb = mkstemp(prefix='.', suffix='.pdb', dir='.', text=True)[1]
    prefix = basename(pdb).rsplit('.', 1)[0]

    # dictionary for mapping residues names to letter representation
    aa_names = {
        'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU',
        'F': 'PHE', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
        'K': 'LYS', 'L': 'LEU', 'M': 'MET', 'N': 'ASN',
        'P': 'PRO', 'Q': 'GLN', 'R': 'ARG', 'S': 'SER',
        'T': 'THR', 'V': 'VAL', 'W': 'TRP', 'Y': 'TYR'
    }
    # The ol' revers-a-roo - didn't want to rewrite ^
    aa_names = {v: k for k, v in aa_names.items()}

    # Initialise the list where all CA info will be appended
    atoms = []

    # Template for finding the CA lines in PDB format
    pattern = re.compile('ATOM.{9}CA .([A-Z]{3}) ([A-Z ])(.{5}).{27}(.{12}).*')

    # reading the pdb + writing new tmp PDB file
    with open(filename, 'r') as f, open(pdb, 'w') as tmp:

        for line in f:

            # stop at end of PDB
            if line.startswith('ENDMDL'):
                break

            else:
                # find matches to 'ATOM CA ...' format (lines in PDB containing CA backbone positions)
                match = re.match(pattern, line)

                # Skipping lines that dont match (return None)
                if match:
                    # append to list of CA atom positions + write to tmp PDB
                    atoms.append(match.groups())
                    tmp.write(line)

        # Quick error check - just in case no CA atoms are present in PDB
        if not len(atoms):
            raise Exception('File %s contains no CA atoms' % filename)

        # # Setup for multiple/broken chains check in PDB
        # chains = [atoms[0][1]]
        # seq = ''
        # rr = int(atoms[0][2]) - 1
        #
        # for a in atoms:
        #     s, c, r = a[:3]
        #
        #     # check for broken chain
        #     if int(r) != int(rr) + 1:
        #         seq += '/'
        #
        #     # move along sequence + update the sequence text representation
        #     rr = r
        #     seq += aa_names[s]
        #
        #     if c not in chains:
        #         chains += c


    atoms_array = np.array(atoms)
    chains_lst = atoms_array[:,1]

    chains = list(np.unique(chains_lst))

    current_chain = chains_lst[0]
    seq = ''
    rr = int(atoms[0][2]) - 1

    for a in atoms:
        s, c, r = a[:3]

        # check for broken chain
        if c != current_chain:
            current_chain = c
            seq += '/'

        # move along sequence + update the sequence text representation
        rr = r
        seq += aa_names[s]

        # if c not in chains:
        #     chains += c

    # temp PIR file
    pir = prefix + '.pir'
    with open(pir, 'w') as f:
        f.write(_PIR_TEMPLATE % (prefix, seq, pdb))

    # Modeller bit - documented standard usage
    env = Environ()
    env.io.atom_files_directory = ['.']
    env.libs.topology.read(file='$(LIB)/top_allh.lib')
    env.libs.parameters.read(file='$(LIB)/par.lib')


    class MyModel(automodel):

        def special_patches(self, aln):
            self.rename_segments(segment_ids=chains)

    # **NEXT*** - adding secondary structure constraints:
    # https://salilab.org/modeller/manual/node28.html


    mdl = MyModel(
        env,
        alnfile=pir,
        knowns='model_ca',
        sequence=prefix,
        assess_methods=assess.DOPE
    )

    mdl.md_level = refine.fast
    mdl.auto_align(matrix_file=prefix + '.mat')
    mdl.starting_model = 1
    mdl.ending_model = int(iterations)
    mdl.final_malign3d = True
    mdl.make()

    # selecting successful models
    models = [m for m in mdl.outputs if m['failure'] is None]

    # Sorting by best to worst model (out of iteration)
    sorted_models = sorted(models, key=lambda d: d['DOPE score'])
    final = sorted_models[0]['name'].rsplit('.', 1)[0] + '_fit.pdb'

    # outputname = filename.split('.')[0] + '_AA.pdb'

    sss = complete_pdb(env, final)
    sss.write(file=outputname, model_format='PDB')

    outfile = sys.stdout

    with open(final) as f:
        a = iter(atoms)
        current = ch = r = t = nl = None
        for line in f:
            if line.startswith('ATOM'):
                res = line[21:27]
                if not current or current != res:
                    current = res
                    ch, r, t = a.__next__()[1:]
                nl = line[:21] + ch + r + line[27:54] + t
                if len(line) > 66:
                    nl += line[66:]
                outfile.write(nl)
            elif line.startswith('TER '):
                outfile.write(line[:22] + nl[22:27] + '\n')
            else:
                outfile.write(line)

    # clean up (theres tonnes of shit modeller has created)
    junk = glob.glob(prefix + '*')
    for j in junk:
        os.remove(j)
        
    sys.stdout = old_stdout



def CA2AA_secondary_slow(filename, outputname, ss_list, iterations=1, stout=False):
    # output=None

    # template for pir file format (shitty format but modeller likes it and I want him to be happy)
    _PIR_TEMPLATE = '\n'.join(['>P1;%s', 'sequence:::::::::', '%s', '*', '', '>P1;model_ca', 'structure:%s:FIRST:@:END:@::::', '*'])

   
    # printing is slowing us down/bogging up the output - preventing output
    old_stdout = sys.stdout
    if stout != True:
        sys.stdout = open(os.devnull, 'w')

    # need a name for a loads of temporary files - naming doesn't matter as they'll be deleted later
    pdb = mkstemp(prefix='.', suffix='.pdb', dir='.', text=True)[1]
    prefix = basename(pdb).rsplit('.', 1)[0]

    # dictionary for mapping residues names to letter representation
    aa_names = {
        'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU',
        'F': 'PHE', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
        'K': 'LYS', 'L': 'LEU', 'M': 'MET', 'N': 'ASN',
        'P': 'PRO', 'Q': 'GLN', 'R': 'ARG', 'S': 'SER',
        'T': 'THR', 'V': 'VAL', 'W': 'TRP', 'Y': 'TYR'
    }
    # The ol' revers-a-roo - didn't want to rewrite ^
    aa_names = {v: k for k, v in aa_names.items()}

    # Initialise the list where all CA info will be appended
    atoms = []

    # Template for finding the CA lines in PDB format
    pattern = re.compile('ATOM.{9}CA .([A-Z]{3}) ([A-Z ])(.{5}).{27}(.{12}).*')

    # reading the pdb + writing new tmp PDB file
    with open(filename, 'r') as f, open(pdb, 'w') as tmp:

        for line in f:

            # stop at end of PDB
            if line.startswith('ENDMDL'):
                break

            else:
                # find matches to 'ATOM CA ...' format (lines in PDB containing CA backbone positions)
                match = re.match(pattern, line)

                # Skipping lines that dont match (return None)
                if match:
                    # append to list of CA atom positions + write to tmp PDB
                    atoms.append(match.groups())
                    tmp.write(line)

        # Quick error check - just in case no CA atoms are present in PDB
        if not len(atoms):
            raise Exception('File %s contains no CA atoms' % filename)

    atoms_array = np.array(atoms)
    chains_lst = atoms_array[:,1]

    chains = list(np.unique(chains_lst))

    current_chain = chains_lst[0]
    seq = ''
    rr = int(atoms[0][2]) - 1

    for a in atoms:
        s, c, r = a[:3]

        # check for broken chain
        if c != current_chain:
            current_chain = c
            seq += '/'

        # move along sequence + update the sequence text representation
        rr = r
        seq += aa_names[s]

        # if c not in chains:
        #     chains += c

    # temp PIR file
    pir = prefix + '.pir'
    with open(pir, 'w') as f:
        f.write(_PIR_TEMPLATE % (prefix, seq, pdb))

    # Modeller bit - documented standard usage
    env = Environ()
    env.io.atom_files_directory = ['.']
    env.libs.topology.read(file='$(LIB)/top_allh.lib')
    env.libs.parameters.read(file='$(LIB)/par.lib')


    # class MyModel(automodel):

    #     def special_patches(self, aln):
    #         self.rename_segments(segment_ids=chains)

    # **NEXT*** - adding secondary structure constraints:
    # https://salilab.org/modeller/manual/node28.html

    class MyModel(automodel):
        def __init__(self, env, alnfile, knowns, sequence, assess_methods, ss_list):
            super().__init__(env, alnfile=alnfile, knowns=knowns, sequence=sequence, assess_methods=assess_methods)
            self.ss_list = ss_list

        def special_patches(self, aln):
            self.rename_segments(segment_ids=chains)

        def special_restraints(self, aln):
            rsr = self.restraints
            at = self.atoms

            if self.ss_list is not None and len(self.ss_list) > 0:
                self.add_secondary_structure_restraints()

        def add_secondary_structure_restraints(self):
            rsr = self.restraints
            current_element = self.ss_list[0]
            start = 1
            for i, ss in enumerate(self.ss_list[1:], start=2):
                if ss != current_element:
                    end = i - 1
                    if current_element == 'H':
                        rsr.add(secondary_structure.alpha(self.residue_range(f'{start}:A', f'{end}:A')))
                    elif current_element == 'S':
                        rsr.add(secondary_structure.strand(self.residue_range(f'{start}:A', f'{end}:A')))
                    start = i
                    current_element = ss
            
            # Add the last element
            end = len(self.ss_list)
            if current_element == 'H':
                rsr.add(secondary_structure.alpha(self.residue_range(f'{start}:A', f'{end}:A')))
            elif current_element == 'S':
                rsr.add(secondary_structure.strand(self.residue_range(f'{start}:A', f'{end}:A')))



    mdl = MyModel(
        env,
        alnfile=pir,
        knowns='model_ca',
        sequence=prefix,
        assess_methods=assess.DOPE,
        ss_list=ss_list
    )

    mdl.md_level = refine.slow
    mdl.auto_align(matrix_file=prefix + '.mat')
    mdl.starting_model = 1
    mdl.ending_model = int(iterations)
    mdl.final_malign3d = True
    mdl.make()

    # selecting successful models
    models = [m for m in mdl.outputs if m['failure'] is None]

    # Sorting by best to worst model (out of iteration)
    sorted_models = sorted(models, key=lambda d: d['DOPE score'])
    final = sorted_models[0]['name'].rsplit('.', 1)[0] + '_fit.pdb'

    # outputname = filename.split('.')[0] + '_AA.pdb'

    sss = complete_pdb(env, final)
    sss.write(file=outputname, model_format='PDB')

    outfile = sys.stdout

    with open(final) as f:
        a = iter(atoms)
        current = ch = r = t = nl = None
        for line in f:
            if line.startswith('ATOM'):
                res = line[21:27]
                if not current or current != res:
                    current = res
                    ch, r, t = a.__next__()[1:]
                nl = line[:21] + ch + r + line[27:54] + t
                if len(line) > 66:
                    nl += line[66:]
                outfile.write(nl)
            elif line.startswith('TER '):
                outfile.write(line[:22] + nl[22:27] + '\n')
            else:
                outfile.write(line)

    # clean up (theres tonnes of shit modeller has created)
    junk = glob.glob(prefix + '*')
    for j in junk:
        os.remove(j)
        
    sys.stdout = old_stdout


def CA2AA_secondary_fast(filename, outputname, ss_list, iterations=1, stout=False):
    # output=None

    # template for pir file format (shitty format but modeller likes it and I want him to be happy)
    _PIR_TEMPLATE = '\n'.join(['>P1;%s', 'sequence:::::::::', '%s', '*', '', '>P1;model_ca', 'structure:%s:FIRST:@:END:@::::', '*'])

   
    # printing is slowing us down/bogging up the output - preventing output
    old_stdout = sys.stdout
    if stout != True:
        sys.stdout = open(os.devnull, 'w')

    # need a name for a loads of temporary files - naming doesn't matter as they'll be deleted later
    pdb = mkstemp(prefix='.', suffix='.pdb', dir='.', text=True)[1]
    prefix = basename(pdb).rsplit('.', 1)[0]

    # dictionary for mapping residues names to letter representation
    aa_names = {
        'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU',
        'F': 'PHE', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
        'K': 'LYS', 'L': 'LEU', 'M': 'MET', 'N': 'ASN',
        'P': 'PRO', 'Q': 'GLN', 'R': 'ARG', 'S': 'SER',
        'T': 'THR', 'V': 'VAL', 'W': 'TRP', 'Y': 'TYR'
    }
    # The ol' revers-a-roo - didn't want to rewrite ^
    aa_names = {v: k for k, v in aa_names.items()}

    # Initialise the list where all CA info will be appended
    atoms = []

    # Template for finding the CA lines in PDB format
    pattern = re.compile('ATOM.{9}CA .([A-Z]{3}) ([A-Z ])(.{5}).{27}(.{12}).*')

    # reading the pdb + writing new tmp PDB file
    with open(filename, 'r') as f, open(pdb, 'w') as tmp:

        for line in f:

            # stop at end of PDB
            if line.startswith('ENDMDL'):
                break

            else:
                # find matches to 'ATOM CA ...' format (lines in PDB containing CA backbone positions)
                match = re.match(pattern, line)

                # Skipping lines that dont match (return None)
                if match:
                    # append to list of CA atom positions + write to tmp PDB
                    atoms.append(match.groups())
                    tmp.write(line)

        # Quick error check - just in case no CA atoms are present in PDB
        if not len(atoms):
            raise Exception('File %s contains no CA atoms' % filename)

    atoms_array = np.array(atoms)
    chains_lst = atoms_array[:,1]

    chains = list(np.unique(chains_lst))

    current_chain = chains_lst[0]
    seq = ''
    rr = int(atoms[0][2]) - 1

    for a in atoms:
        s, c, r = a[:3]

        # check for broken chain
        if c != current_chain:
            current_chain = c
            seq += '/'

        # move along sequence + update the sequence text representation
        rr = r
        seq += aa_names[s]

        # if c not in chains:
        #     chains += c

    # temp PIR file
    pir = prefix + '.pir'
    with open(pir, 'w') as f:
        f.write(_PIR_TEMPLATE % (prefix, seq, pdb))

    # Modeller bit - documented standard usage
    env = Environ()
    env.io.atom_files_directory = ['.']
    env.libs.topology.read(file='$(LIB)/top_allh.lib')
    env.libs.parameters.read(file='$(LIB)/par.lib')


    # class MyModel(automodel):

    #     def special_patches(self, aln):
    #         self.rename_segments(segment_ids=chains)

    # **NEXT*** - adding secondary structure constraints:
    # https://salilab.org/modeller/manual/node28.html

    class MyModel(automodel):
        def __init__(self, env, alnfile, knowns, sequence, assess_methods, ss_list):
            super().__init__(env, alnfile=alnfile, knowns=knowns, sequence=sequence, assess_methods=assess_methods)
            self.ss_list = ss_list

        def special_patches(self, aln):
            self.rename_segments(segment_ids=chains)

        def special_restraints(self, aln):
            rsr = self.restraints
            at = self.atoms

            if self.ss_list is not None and len(self.ss_list) > 0:
                self.add_secondary_structure_restraints()

        def add_secondary_structure_restraints(self):
            rsr = self.restraints
            current_element = self.ss_list[0]
            start = 1
            for i, ss in enumerate(self.ss_list[1:], start=2):
                if ss != current_element:
                    end = i - 1
                    if current_element == 'H':
                        rsr.add(secondary_structure.alpha(self.residue_range(f'{start}:A', f'{end}:A')))
                    elif current_element == 'S':
                        rsr.add(secondary_structure.strand(self.residue_range(f'{start}:A', f'{end}:A')))
                    start = i
                    current_element = ss
            
            # Add the last element
            end = len(self.ss_list)
            if current_element == 'H':
                rsr.add(secondary_structure.alpha(self.residue_range(f'{start}:A', f'{end}:A')))
            elif current_element == 'S':
                rsr.add(secondary_structure.strand(self.residue_range(f'{start}:A', f'{end}:A')))



    mdl = MyModel(
        env,
        alnfile=pir,
        knowns='model_ca',
        sequence=prefix,
        assess_methods=assess.DOPE,
        ss_list=ss_list
    )

    mdl.md_level = refine.fast
    mdl.auto_align(matrix_file=prefix + '.mat')
    mdl.starting_model = 1
    mdl.ending_model = int(iterations)
    mdl.final_malign3d = True
    mdl.make()

    # selecting successful models
    models = [m for m in mdl.outputs if m['failure'] is None]

    # Sorting by best to worst model (out of iteration)
    sorted_models = sorted(models, key=lambda d: d['DOPE score'])
    final = sorted_models[0]['name'].rsplit('.', 1)[0] + '_fit.pdb'

    # outputname = filename.split('.')[0] + '_AA.pdb'

    sss = complete_pdb(env, final)
    sss.write(file=outputname, model_format='PDB')

    outfile = sys.stdout

    with open(final) as f:
        a = iter(atoms)
        current = ch = r = t = nl = None
        for line in f:
            if line.startswith('ATOM'):
                res = line[21:27]
                if not current or current != res:
                    current = res
                    ch, r, t = a.__next__()[1:]
                nl = line[:21] + ch + r + line[27:54] + t
                if len(line) > 66:
                    nl += line[66:]
                outfile.write(nl)
            elif line.startswith('TER '):
                outfile.write(line[:22] + nl[22:27] + '\n')
            else:
                outfile.write(line)

    # clean up (theres tonnes of shit modeller has created)
    junk = glob.glob(prefix + '*')
    for j in junk:
        os.remove(j)
        
    sys.stdout = old_stdout


def CA2AA_secondary_multimer(filename, outputname, ss_list, disulfides=None, iterations=1, stout=False):
    _PIR_TEMPLATE = '\n'.join([
        '>P1;%s',
        'sequence:::::::::',
        '%s',
        '*',
        '',
        '>P1;model_ca',
        'structure:%s:FIRST:@:END:@::::',
        '*'
    ])

    # Suppress output if not requested
    if not stout:
        sys.stdout = open(os.devnull, 'w')

    pdb = mkstemp(prefix='.', suffix='.pdb', dir='.', text=True)[1]
    prefix = basename(pdb).rsplit('.', 1)[0]

    # Reverse map residue names to 1-letter code
    aa_names = {
        'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E',
        'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N',
        'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S',
        'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
    }

    atoms = []
    pattern = re.compile(r'^ATOM.{9}CA\s+([A-Z]{3})\s([A-Z])\s+(\d+).{27}(.{12}).*')

    # Read and write filtered CA-only PDB
    with open(filename, 'r') as f, open(pdb, 'w') as tmp:
        for line in f:
            if line.startswith('ENDMDL'):
                break
            match = re.match(pattern, line)
            if match:
                atoms.append(match.groups())
                tmp.write(line)

    if not atoms:
        raise Exception(f'File {filename} contains no CA atoms')

    atoms_array = np.array(atoms)
    chains_lst = atoms_array[:,1]
    unique_chains = list(np.unique(chains_lst))

    seq = ''
    current_chain = chains_lst[0]

    for i, a in enumerate(atoms):
        resname, chain, resnum = a[:3]
        if chain != current_chain:
            seq += '/'
            current_chain = chain
        seq += aa_names.get(resname.strip(), 'X')

    pir = prefix + '.pir'
    with open(pir, 'w') as f:
        f.write(_PIR_TEMPLATE % (prefix, seq, pdb))

    env = Environ()
    env.io.atom_files_directory = ['.']
    env.libs.topology.read(file='$(LIB)/top_allh.lib')
    env.libs.parameters.read(file='$(LIB)/par.lib')

    class MyModel(automodel):
        def __init__(self, env, alnfile, knowns, sequence, assess_methods, ss_list, disulfides):
            super().__init__(env, alnfile=alnfile, knowns=knowns, sequence=sequence, assess_methods=assess_methods)
            self.ss_list = ss_list
            self.disulfides = disulfides or []

        def special_patches(self, aln):
            # FIX: Extract string names from chains
            seen_chain_ids = sorted({str(res.chain.name) for res in self.residues})
            self.rename_segments(segment_ids=seen_chain_ids)
            
            for res1_str, res2_str in self.disulfides:
                res1_num, chain1 = res1_str.split(':')
                res2_num, chain2 = res2_str.split(':')
                try:
                    res1 = self.residues[f'{int(res1_num)}:{chain1}']
                    res2 = self.residues[f'{int(res2_num)}:{chain2}']
                    self.patch(residue_type='DISU', residues=(res1, res2))
                except KeyError:
                    print(f"Warning: could not find residues {res1_str} or {res2_str} for disulfide bond")
            

        def special_restraints(self, aln):
            if not self.ss_list:
                return
            self.add_secondary_structure_restraints()

        def add_secondary_structure_restraints(self):
            rsr = self.restraints
            i = 1
            start = 1
            current_ss = self.ss_list[0]
            chain = chains_lst[0]

            for i in range(1, len(self.ss_list)):
                if self.ss_list[i] != current_ss or chains_lst[i] != chains_lst[i-1]:
                    end = i
                    if current_ss == 'H':
                        rsr.add(secondary_structure.alpha(self.residue_range(f'{start}:{chain}', f'{end}:{chain}')))
                    elif current_ss == 'S':
                        rsr.add(secondary_structure.strand(self.residue_range(f'{start}:{chain}', f'{end}:{chain}')))
                    start = i + 1
                    current_ss = self.ss_list[i]
                    chain = chains_lst[i]
            # Final stretch
            if current_ss == 'H':
                rsr.add(secondary_structure.alpha(self.residue_range(f'{start}:{chain}', f'{len(self.ss_list)}:{chain}')))
            elif current_ss == 'S':
                rsr.add(secondary_structure.strand(self.residue_range(f'{start}:{chain}', f'{len(self.ss_list)}:{chain}')))

    mdl = MyModel(
        env,
        alnfile=pir,
        knowns='model_ca',
        sequence=prefix,
        assess_methods=assess.DOPE,
        ss_list=ss_list,
        disulfides=disulfides
    )

    mdl.md_level = refine.fast
    mdl.auto_align(matrix_file=prefix + '.mat')
    mdl.starting_model = 1
    mdl.ending_model = int(iterations)
    mdl.final_malign3d = True
    mdl.make()

    models = [m for m in mdl.outputs if m['failure'] is None]
    sorted_models = sorted(models, key=lambda d: d['DOPE score'])
    final = sorted_models[0]['name'].rsplit('.', 1)[0] + '_fit.pdb'

    sss = complete_pdb(env, final)
    sss.write(file=outputname, model_format='PDB')
    renumber_pdb_chains_start_from_1(outputname,outputname)

def backmap_ca_chain(coords_file, fingerprint_file, write_directory, name, ss_constraint=False):

    # write the CA chain into pdb format - note this won't work if non-standard residues are present!
    ca_pdb_output_name = os.path.join(write_directory, name+'_CA.pdb')
    cdt.Carbonara_2_PDB(coords_file, fingerprint_file, ca_pdb_output_name)
    print('Alpha Coordinates pdb written to: ', ca_pdb_output_name)

    aa_pdb_output_name = os.path.join(write_directory, name+'_AA.pdb')

    ss_list = list(np.genfromtxt(fingerprint_file, dtype=str)[2])

    if ss_constraint:
        CA2AA_secondary(ca_pdb_output_name, aa_pdb_output_name, ss_list, iterations=3, stout=False)

    else:
        CA2AA(ca_pdb_output_name, aa_pdb_output_name, iterations=3, stout=False)

    print('All Atomistic pdb written to: ', aa_pdb_output_name)

def backmap_ca_chain_multimer(coords_file, fingerprint_file, write_directory, name,lengths,disulfides=None):

    # write the CA chain into pdb format - note this won't work if non-standard residues are present!
    split_coords_into_chains(coords_file,coords_file, lengths)
    
    ca_pdb_output_name = os.path.join(write_directory, name+'_CA.pdb')
    cdt.Carbonara_2_PDB(coords_file, fingerprint_file, ca_pdb_output_name)
    print('Alpha Coordinates pdb written to: ', ca_pdb_output_name)

    aa_pdb_output_name = os.path.join(write_directory, name+'_AA.pdb')

    ss_list = list(np.genfromtxt(fingerprint_file, dtype=str)[2])

    CA2AA_secondary_multimer(ca_pdb_output_name, aa_pdb_output_name, ss_list, disulfides,iterations=1,stout=False)

    print('All Atomistic pdb written to: ', aa_pdb_output_name)

def read_json_from_file(file_path):
    with open(file_path, 'r') as f:
        log_data = f.read()
    return log_data

def split_coords_into_chains(input_path, output_path, segment_lengths):
    """
    Re-splits a .dat coordinate file with 'End chain ...' lines into new chains
    using the given segment_lengths. Always writes 'End chain' after each block.
    """
    import re

    with open(input_path, 'r') as infile:
        raw_lines = [line.strip() for line in infile if line.strip()]
    
    # Ignore lines that contain 'End chain' (in any form)
    coord_lines = [line for line in raw_lines if not re.search(r'end\s+chain', line, re.IGNORECASE)]

    total_input = len(coord_lines)
    total_expected = sum(segment_lengths)

    print(f"🔍 Found {total_input} coordinates, expecting {total_expected} from segment_lengths")

    if total_input != total_expected:
        raise ValueError(f"Mismatch: {total_input} coords vs {total_expected} expected")

    idx = 0
    with open(output_path, 'w') as outfile:
        for i, L in enumerate(segment_lengths):
            for _ in range(L):
                outfile.write(coord_lines[idx] + '\n')
                idx += 1
            outfile.write(f"End chain {i+1}\n")

    print(f"✅ Wrote {output_path} with {len(segment_lengths)} chains.")


def getFitFiles(directory, threshold="last"):
    # List all files in the directory that contain "fitLog" in their filename
    fitlog_files = [f for f in os.listdir(directory) if 'fitLog' in f]
    fitlog_paths = [os.path.join(directory, f) for f in fitlog_files]
    
    molecule_paths = []
    for fitlog in fitlog_paths:
        log_data = read_json_from_file(fitlog)
        lines = log_data.strip().split('\n')

        if threshold == "last":
            data = json.loads(lines[-1])
            mol_path = data.get("MoleculePath")
            scat_path = data.get("ScatterPath")
            molecule_paths.append([
                _relativize_path(mol_path, directory),
                _relativize_path(scat_path, directory)
            ])
        else:
            for line in lines:
                if not line or line.startswith('{"Run"'):
                    continue
                data = json.loads(line)
                if data.get("ScatterFitFirst", float('inf')) < threshold:
                    mol_path = data.get("MoleculePath")
                    scat_path = data.get("ScatterPath")
                    molecule_paths.append([
                        _relativize_path(mol_path, directory),
                        _relativize_path(scat_path, directory)
                    ])
    return molecule_paths

def _relativize_path(full_path, directory):
    """
    Strips everything before the directory and returns the relative file path from 'directory'.
    If the file is not in 'directory', return the filename joined with directory.
    """
    if not full_path:
        return None
    filename = os.path.basename(full_path)
    return os.path.join(directory, filename)


def generateAllAtomisticFits(directory,run,threshold="last"):
    moleculePaths =getFitFiles(directory+run,threshold)
    [backmap_ca_chain(moleculePaths[i][0], directory+"fingerPrint1.dat", directory+run,moleculePaths[i][0].strip().split('/')[-1].strip().split('xyz.dat')[0], ss_constraint=True) for i in range(len(moleculePaths)) ]

def generateAllAtomisticFitsMultimter(directory,run,lengths,threshold="last",disulfides=None):
    moleculePaths =getFitFiles(directory+run,threshold)
    [backmap_ca_chain_multimer(moleculePaths[i][0], directory+"fingerPrint1.dat", directory+run,moleculePaths[i][0].strip().split('/')[-1].strip().split('xyz.dat')[0],lengths,disulfides) for i in range(len(moleculePaths)) ]

import mdtraj as md
import string

def get_disulfide_distances(pdb_file, disulfide_list):
    traj = md.load(pdb_file)
    topology = traj.topology
    atom_df, _ = topology.to_dataframe()

    # Dynamically assign 'A', 'B', ... to chain indices in order of appearance
    available_letters = list(string.ascii_uppercase)
    chain_id_map = {
        available_letters[i]: chain.index
        for i, chain in enumerate(topology.chains)
    }

    print("User chain label → MDTraj chain index mapping:", chain_id_map)

    # Build SG atom index mapping: (resSeq, internal chain index)
    sg_atom_indices = {}
    for atom in topology.atoms:
        if atom.name == 'SG' and atom.residue.name == 'CYS':
            chain_idx = atom.residue.chain.index
            res_seq = atom.residue.resSeq
            sg_atom_indices[(res_seq, chain_idx)] = atom.index

    # Process each disulfide pair
    results = []
    for res1, res2 in disulfide_list:
        resnum1, chain1 = res1.split(':')
        resnum2, chain2 = res2.split(':')

        if chain1 not in chain_id_map or chain2 not in chain_id_map:
            print(f"Warning: Chain {chain1} or {chain2} not found in chain map.")
            results.append((res1, res2, None))
            continue

        key1 = (int(resnum1), chain_id_map[chain1])
        key2 = (int(resnum2), chain_id_map[chain2])

        idx1 = sg_atom_indices.get(key1)
        idx2 = sg_atom_indices.get(key2)

        if idx1 is None or idx2 is None:
            print(f"Missing SG atom for: {key1 if idx1 is None else ''} {key2 if idx2 is None else ''}")
            results.append((res1, res2, None))
        else:
            dist_nm = md.compute_distances(traj, [[idx1, idx2]])[0][0]
            dist_angstrom = dist_nm * 10
            results.append((res1, res2, dist_angstrom))

    return results

from Bio.PDB import PDBParser

def count_c_alpha_atoms(pdb_filename):
    """
    Counts the number of C-alpha (CA) atoms in a PDB file.

    Parameters:
        pdb_filename (str): Path to the PDB file

    Returns:
        int: Number of CA atoms
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("structure", pdb_filename)
    ca_count = 0

    for model in structure:
        for chain in model:
            for residue in chain:
                if "CA" in residue:
                    ca_count += 1
    return ca_count


