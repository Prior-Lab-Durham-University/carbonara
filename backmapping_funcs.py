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



def CA2AA_secondary(filename, outputname, ss_list, iterations=1, stout=False):
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
