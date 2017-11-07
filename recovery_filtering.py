#!/usr/bin/env python3

"""
All I'm interested in is RMSD to ideal ligand position

Usage:
    recovery_filtering <clean_pdb> <ligand_ID> <matched_pdb_path>

Arguments:
    <clean_pdb>
        Path to cleaned PDB for which binding site is being recovered

    <ligand_ID>
        Three letter ligand ID

    <matched_pdb_path>
        Path to directory containing all matched PDBs
"""
import os
import re

import docopt
import numpy as np
import pandas as pd
from pathos.multiprocessing import ProcessingPool as Pool

from utils import *

class BindingSiteRecoveryFiltering():
    """
    Filter matches for binding site recovery
    """
    # Shitty atom mapping for now {3D2D: params}
    # todo: generate cleaned PDBs with correct atom names...
    REN_atom_mapping = {'C23': 'C18',
                        'O22': 'O1',
                        'C21': 'C15',
                        'C19': 'C14',
                        'O20': 'O2',
                        'C18': 'C13',
                        'C17': 'C12',
                        'C25': 'C17',
                        'C24': 'C16',
                        'C16': 'C11',
                        'C14': 'C10',
                        'N12': 'N1',
                        'C13': 'C1',
                        'C11': 'C2',
                        'C10': 'C3',
                        'C9': 'C4',
                        'C1': 'C9',
                        'C2': 'C8',
                        'C3': 'C7',
                        'O4': 'O3',
                        'C5': 'C6',
                        'O6': 'O4',
                        'C7': 'C19',
                        'C8': 'C5'}

    three_to_one = {'ALA':'A', 'CYS':'C', 'ASP':'D', 'GLU':'E', 'PHE':'F', 'GLY':'G', 'HIS':'H', 'ILE':'I', 'LYS':'K',
                    'LEU':'L', 'MET':'M', 'ASN':'N', 'PRO':'P', 'GLN':'Q', 'ARG':'R', 'SER':'S', 'THR':'T', 'VAL':'V',
                    'TRP':'W', 'TYR':'Y'}

    def __init__(self, cleaned_scaffold, ligand_ID):
        self.cleaned_scaffold = prody.parsePDB(cleaned_scaffold)
        self.ligand_ID = ligand_ID
        self.cleaned_scaffold_ligand = self.cleaned_scaffold.select('resname {}'.format(self.ligand_ID))
        self.ideal_resnum_resname_dict = self._generate_ideal_resnum_resname_mapping()

    def _generate_ideal_resnum_resname_mapping(self):
        """
        Generate mapping between ideal resnum and resname
        :return:
        """
        clean_pridy_hv = self.cleaned_scaffold.getHierView()
        return {residue.getResnum(): self.three_to_one[residue.getResname()] for residue in clean_pridy_hv.iterResidues() if residue.getResname() != self.ligand_ID}

    def calculate_match_ideal_rmsd(self, match_prody):
        """
        Calculate absolute RMSD between ideal and matched ligands
        :return:
        """
        match_prody_ligand = match_prody.select('resname {}'.format(self.ligand_ID))

        # Create two np arrays with coordinates for ideal and matched ligands
        match_prody_atom_list = []

        for atom in self.cleaned_scaffold_ligand:
            current_ideal_atom = atom.getName()
            current_match_atom_name = self.REN_atom_mapping[current_ideal_atom]

            match_prody_atom_list.append(match_prody.select('name {}'.format(current_match_atom_name)).getCoords()[0])

        match_ligand_coords = np.asarray(match_prody_atom_list)
        ideal_ligand_coords = self.cleaned_scaffold_ligand.getCoords()

        return prody.calcRMSD(ideal_ligand_coords, match_ligand_coords)

    def filter_for_correct_residues(self, match_name):
        """
        Determine fraction of match residues that recover the correct identities at scaffold positions
        :return:
        """
        match_residue_block = match_name.split('_')[2]

        match_residue_block_split = re.findall('[A-Z]|\d{1,5}', match_residue_block)
        match_residue_block_tuple = []

        for i in range(int(len(match_residue_block_split) / 2)):
            match_residue_block_tuple.append((match_residue_block_split[i * 2], int(match_residue_block_split[i * 2 + 1])))

        correct_count = 0
        for resname, resnum in match_residue_block_tuple:
            if resname == self.ideal_resnum_resname_dict[resnum]:
                correct_count += 1

        return correct_count / len(match_residue_block_tuple)

if __name__ == '__main__':
    args = docopt.docopt(__doc__)

    filter = BindingSiteRecoveryFiltering(args['<clean_pdb>'], args['<ligand_ID>'])

    def recovery_filtering(pdb):
        match_prody = prody.parsePDB(pdb)
        match_rmsd = filter.calculate_match_ideal_rmsd(match_prody)

        match_name = os.path.basename(os.path.normpath(pdb))
        fraction_correct = filter.filter_for_correct_residues(match_name)

        return {'match_name': match_name,
                'match_rmsd': match_rmsd,
                'fraction_correct': fraction_correct}

    # Multiprocess match evaluation
    process = Pool()
    match_metrics_list_of_dicts = process.map(recovery_filtering, [pdb for pdb in pdb_check(args['<matched_pdb_path>'])])
    process.close()
    process.join()

    # Return passing matcher results
    df = pd.DataFrame(match_metrics_list_of_dicts)
    df.set_index(['match_name'], inplace=True)
    df.to_csv('Recovery_Filter_Results-{}-{}.csv'.format(args['<ligand_ID>'], os.path.basename(os.path.normpath(args['<clean_pdb>']))[:4]))



