#! /usr/bin/python3
"""
AnnotatedMotifResiduesFromMatch
Pull annotated fuzzball residues for match motif residues

Usage:
    annotations <match_pdb> <fuzzball_dir>

Arguments:
    <match_pdb>             Path to match PDB
    <fuzzball_dir>          Path to Fuzzball directory for match ligand

Options:


"""

import os
import docopt
import prody

args = docopt.docopt(__doc__)
match_pdb_path = args['<match_pdb>']
fuzzball_dir = args['<fuzzball_dir>']

# Parse match name
match_filename, _ = os.path.splitext(os.path.normpath(os.path.basename(match_pdb_path)))
match_filename_split = match_filename.split('-')

# Get correct fuzzball
ccd = match_filename_split[0].split('_')[-2]
ccd_conformer = match_filename_split[0].split('_')[-1]
fuzzball_iteration = match_filename_split[1]
fuzzball_index = match_filename_split[2]

fuzzball_path = os.path.join(fuzzball_dir, f'{ccd}_{ccd_conformer}-{fuzzball_iteration}-{fuzzball_index}.ag.npz')

if not os.path.exists(fuzzball_path):
    raise Exception(f'{fuzzball_path} doesn\'t exist!')

fuzzball_prody = prody.loadAtoms(fuzzball_path)
fuzzball_prody_hv = fuzzball_prody.getHierView()

# Get residues from fuzzball
fuzzball_residue_indicies = [int(a) for a in match_filename_split[-1].split('_')[1:-1]]

# Report annotations
for idx in fuzzball_residue_indicies:
    motif_residue = fuzzball_prody_hv['A', idx]
    print(motif_residue)
    print(motif_residue.getData('contact_source')[0])
