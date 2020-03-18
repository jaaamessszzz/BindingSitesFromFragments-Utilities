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

import pyrosetta
from pyrosetta import rosetta


def fix_sc_filter(args):
    """
    Finds design_metrics.csv in an output directory generated by Design-CompRotSets.py and recalculates shape
    complementarity. This was necessary for designs where the shape complementarity filter failed on the cluster
    due to a missing entry for chlorine atoms in sc_radii.lib.

    Usage:
      fix_sc_filter <design_dir> [options]

    Arguments:
      <design_dir>                      Path to ligand PDB generated by molfile_to_params.py

    Options:
      -c=<cores>, --cores=<cores>       Number of cores to use (Default: 12)
      --compounds=<compounds>           Comma delimited list of compounds to recalculate SC
    """

    design_root = args['<design_dir>']
    cores = args['--cores'] if args['--cores'] else 12
    compounds = list(args['--compounds'].split(',')) if args['--compounds'] else list()

    pyrosetta.init('-ex1 -ex2 -extrachi_cutoff 0 -use_input_sc -mute \
    protocols.simple_filters.ShapeComplementarityFilter core.conformation.Conformation')

    dirs_to_parse = list()
    design_metrics_path_list = list()

    for root, dirs, files in os.walk(design_root):
        _, design_dir = os.path.split(root)

        # Add designs to dirs_to_parse
        if design_dir.startswith('UM_'):

            # Check for specified compounds
            if len(compounds) > 0 and not any([ccd in design_dir.split('_') for ccd in compounds]):
                continue

            dirs_to_parse.append(design_dir)

        if any([parse_dir in root.split('/') for parse_dir in dirs_to_parse]) and 'design_metrics.csv' in files:
            design_metrics_path_list.append(os.path.join(root, 'design_metrics.csv'))

    print(f'Fixing ShapeComplementarity calculations for the following designs:')
    print("\n".join(dirs_to_parse))

    # Assemble pool
    pool = Pool(processes=cores)
    pool.map(fix_sc, design_metrics_path_list)
    pool.close()
    pool.join()


def fix_sc(design_csv_path):

    print(design_csv_path)

    # ShapeComplementarity
    shape_complementarity_filter = rosetta.protocols.simple_filters.ShapeComplementarityFilter()
    shape_complementarity_filter.filtered_sc(0)
    shape_complementarity_filter.filtered_area(0)
    shape_complementarity_filter.jump_id(1)
    shape_complementarity_filter.quick(0)
    shape_complementarity_filter.verbose(0)
    shape_complementarity_filter.write_int_area(0)

    root, design_csv = os.path.split(design_csv_path)
    design_dir, design_subdir = os.path.split(root)
    sc_dict_list = list()

    # Open design metrics as df
    df = pd.read_csv(design_csv_path, index_col=None)
    df = df.drop(['shapecomplementarity'], axis=1)
    source_path, _ = os.path.split(df['path'].iloc[0])

    # Get params
    params_file = [file for file in os.listdir(os.path.join(design_dir, 'Conformers')) if file.endswith('.params')][0]
    params_path = os.path.join(design_dir, 'Conformers', params_file)
    ligand_params = pyrosetta.Vector1([params_path])

    # Get tarball path
    tarball_path = os.path.join(root, f'{design_subdir}.tar.gz')

    # Open tarball
    with tarfile.open(tarball_path, 'r|gz') as tar:

        for filepath in tar:
            # Assumes tarball only contains .pdb files!
            if filepath.isfile():
                tarpdb = tar.extractfile(filepath)
                _, filename = os.path.split(filepath.name)

                # Load params and pose as in 07.03-Ligand-Docking-PyRosetta.ipynb
                complex_pose = pyrosetta.rosetta.core.pose.Pose()
                residue_type_set = complex_pose.conformation().modifiable_residue_type_set_for_conf()
                residue_type_set.read_files_for_base_residue_types(ligand_params)
                complex_pose.conformation().reset_residue_type_set_for_conf(residue_type_set)
                rosetta.core.import_pose.pose_from_pdbstring(complex_pose, tarpdb.read().decode('utf-8'))

                # Get shapecomplementarity
                shape_complementarity_filter.apply(complex_pose)
                shapecomplementarity = shape_complementarity_filter.report_sm(complex_pose)

                sc_dict_list.append({'shapecomplementarity': shapecomplementarity,
                                     'path': os.path.join(source_path, filename)}
                                    )

    sc_df = pd.DataFrame(sc_dict_list)
    merged_df = pd.merge(df, sc_df, how='inner', on='path', sort=True, copy=True)
    merged_df.to_csv(design_csv_path)


def fix_specialrot_count(args):
    """
    Add counts of special_rot rotamers incorporated into designs

    Usage:
      python3 fix_specialrot <design_dir> [options]

    Arguments:
      <design_dir>                      Path to ligand PDB generated by molfile_to_params.py

    Options:
      -c=<cores>, --cores=<cores>       Number of cores to use (Default: 12)
    """
    design_root = args['<design_dir>']
    cores = args['--cores'] if args['--cores'] else 12

    dirs_to_parse = list()
    design_metrics_path_list = list()

    for root, dirs, files in os.walk(design_root):
        _, design_dir = os.path.split(root)

        # Add designs to dirs_to_parse
        if design_dir.startswith('UM_'):
            print(design_dir)
            dirs_to_parse.append(design_dir)

        if any([parse_dir in root.split('/') for parse_dir in dirs_to_parse]) and 'design_metrics.csv' in files:
            design_metrics_path_list.append(os.path.join(root, 'design_metrics.csv'))

    print(f'Fixing special_rot counts for the following designs:')
    print("\n".join(dirs_to_parse))

    # Assemble pool
    pool = Pool(processes=cores)
    pool.map(count_specialrot, design_metrics_path_list)
    pool.close()
    pool.join()

def count_specialrot(design_csv_path):
    """
    Update design_metrics.csv with number of special_rot rotamers incorporated

    :param design_csv_path: path to design_metrics.csv
    """
    print(design_csv_path)

    root, design_csv = os.path.split(design_csv_path)
    design_dir, design_subdir = os.path.split(root)
    specialrot_dict_list = list()

    # Open design metrics as df
    df = pd.read_csv(design_csv_path, index_col=None, usecols=['path',
                                                               'match',
                                                               'bindingstrain',
                                                               'residueie',
                                                               'packstat',
                                                               'heavyburiedunsats',
                                                               'ligand_sasa',
                                                               'hbonds',
                                                               'comprotset',
                                                               'special_rot_weight',
                                                               'holes',
                                                               'shapecomplementarity',])
    source_path, _ = os.path.split(df['path'].iloc[0])
    use_specialrot = df['special_rot_weight'].iloc[0] != 0

    # Get tarball path
    tarball_path = os.path.join(root, f'{design_subdir}.tar.gz')

    # Open tarball
    with tarfile.open(tarball_path, 'r|gz') as tar:

        for filepath in tar:
            # Assumes tarball only contains .pdb files!
            if filepath.isfile():
                tarpdb = tar.extractfile(filepath)
                _, filename = os.path.split(filepath.name)

                for bline in tarpdb:
                    line = bline.decode('utf-8')
                    if line.startswith('weights'):
                        specialrot_weight = float(line.split()[-3])
                    if line.startswith('pose'):
                        if use_specialrot:
                            incorporated_specialrot = float(line.split()[-3]) / specialrot_weight
                        else:
                            incorporated_specialrot = 0

                specialrot_dict_list.append({'incorporated_specialrot': incorporated_specialrot,
                                             'path': os.path.join(source_path, filename)}
                                            )

    specialrot_df = pd.DataFrame(specialrot_dict_list)
    merged_df = pd.merge(df, specialrot_df, how='inner', on='path', sort=True, copy=True)
    merged_df.to_csv(design_csv_path)


if __name__ == '__main__':
    argv = sys.argv[1:]
    if argv[0] == 'fix_sc_filter':
        fix_sc_filter(docopt(fix_sc_filter.__doc__, argv=argv))
    if argv[0] == 'fix_specialrot':
        fix_specialrot_count(docopt(fix_specialrot_count.__doc__, argv=argv))
