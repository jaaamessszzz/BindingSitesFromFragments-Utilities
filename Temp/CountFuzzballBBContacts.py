import io
import os
import sys
import tempfile
import tarfile
import shutil
from pprint import pprint
from multiprocessing import Pool

from docopt import docopt
import pandas as pd
import prody

import pyrosetta
from pyrosetta import rosetta

from BindingSitesFromFragments.motifs import Generate_Constraints

def count_backbone_contacts(fuzzball_path):
    """
    Update design_metrics.csv with number of special_rot rotamers incorporated

    :param fuzzball_path: path to fuzzball.ag.npz
    """
    print(fuzzball_path)

    fuzzball = prody.loadAtoms(fuzzball_path)
    fuzzball_hv = fuzzball.getHierView()
    ligand = fuzzball_hv['X', 1]
    ligand_resname = ligand.getResname()

    constraints = Generate_Constraints(ligand_resname)

    motif_count = 0
    bb_contacts = 0

    for index, motif in enumerate(fuzzball_hv.iterResidues(), start=1):
        if index != 1:
            constraint_dict = constraints.determine_constraint_atoms(motif, ligand)
            if constraint_dict is False: continue
            if constraint_dict['residue']['atom_names'][0] in ['C', 'CA', 'O', 'N']:
                bb_contacts += 1
            motif_count += 1

    return {'motifs': motif_count,
            'bb_contacts': bb_contacts,
            'fuzzball': fuzzball_path}


def fuzzball_count(args):
    """
    Add counts of special_rot rotamers incorporated into designs

    Usage:
      python3 bb_contacts <compound_dir> [options]

    Arguments:
      <compound_dir>                    Path to directory containing BSFF project directories

    Options:
      -c=<cores>, --cores=<cores>       Number of cores to use (Default: 12)
    """
    compound_root = args['<compound_dir>']
    cores = args['--cores'] if args['--cores'] else 12

    fuzzballs_to_parse = list()

    # Gather paths to all Fuzzballs
    for root, dirs, files in os.walk(compound_root):

        # Add designs to dirs_to_parse
        if 'Fuzzballs' in root.split('/') and 'Design' not in root.split('/') and any([file.endswith('.ag.npz') for file in files]):
            if any([a.startswith('_') for a in root.split('/')]):
                continue
            fuzzball_agnpz = [os.path.join(root, file) for file in files if file.endswith('.ag.npz')]
            fuzzballs_to_parse = fuzzballs_to_parse + fuzzball_agnpz

    print(f'Counting backbone contacts for the following fuzzballs:')
    print("\n".join(fuzzballs_to_parse))

    # Assemble pool
    pool = Pool(processes=cores)
    data = pool.map(count_backbone_contacts, fuzzballs_to_parse)
    pool.close()
    pool.join()

    df = pd.DataFrame(data)
    df.to_csv('Fuzzball_bb_contacts.csv')

if __name__ == '__main__':
    argv = sys.argv[1:]
    if argv[0] == 'bb_contacts':
        fuzzball_count(docopt(fuzzball_count.__doc__, argv=argv))