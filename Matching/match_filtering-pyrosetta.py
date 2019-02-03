#! /netapp/home/james.lucas/Packages/Python3.6/bin/python
#$ -S /netapp/home/james.lucas/Packages/Python3.6/bin/python
#$ -cwd
#$ -r yes
#$ -l h_rt=240:00:00
#$ -t 1-10025
#$ -j y
#$ -l arch=linux-x64
#$ -l mem_free=15G
#$ -l netapp=1G,scratch=10G

"""
Testing PyRosetta...

Iterate over matches to get pairwise energies for each motif residue with ligand. Also mutate all residues to ALA
except match residues and return largest fa_rep for each motif residue.

Usage:
    pyrosetta_scoring_test match <block_size> <ligand_ID> <params_dir>
    pyrosetta_scoring_test score <input_PDB_dir> <ligand_ID> <params_dir>
    pyrosetta_scoring_test score_matched <ligand_ID> <params_dir>
    pyrosetta_scoring_test consolidate
    pyrosetta_scoring_test move

Arguments:
    match
        Match and score

    score
        Score matches

    consolidate
        Consolidate filter metrics after matching into a single dataframe

    move
        Move quality matches

    <input_PDB_dir>
        Path to directory containing input PDB s

    <params_file>
        Path to params file

    <energy_csv>
        .csv report with energies
"""
import os
import sys
import re
import json
import shutil
import tempfile
import itertools
import subprocess

import docopt
import prody
import pandas as pd
from pandas.io.common import EmptyDataError
import numpy as np

import pyrosetta
from pyrosetta import rosetta

def residue_indicies_from_match_name(match_pdb_name):
    pnc = re.split('_|-|\.', os.path.basename(os.path.normpath(match_pdb_name)))
    motif_residue_ID_list = [a for a in re.split('(\D+)', pnc[2]) if a != '']
    motif_residue_IDs = [int(motif_residue_ID_list[indx + 1]) for indx in range(0, len(motif_residue_ID_list), 2)]
    return motif_residue_IDs


def calculate_match_binding_motif_energy(pose, residue_index_list):
    """
    Calculates the total two-body interaction energy for a given match
    :return:
    """

    ligand_index = pose.size()
    motif_and_ligand_idx = residue_index_list + [ligand_index]
    edges = pose.energies().energy_graph()

    score_list = []
    for tuple in itertools.combinations(motif_and_ligand_idx, 2):

        current_edge = edges.find_energy_edge(tuple[0], tuple[1])
        if current_edge is not None:
            edge_scores = current_edge.fill_energy_map()
            score_list.append(sum([edge_scores[rosetta.core.scoring.fa_rep] * 0.55,
                                   edge_scores[rosetta.core.scoring.fa_atr],
                                   edge_scores[rosetta.core.scoring.hbond_sc],
                                   edge_scores[rosetta.core.scoring.fa_sol],
                                   edge_scores[rosetta.core.scoring.fa_elec]
                                   ]
                                  )
                              )

    return sum(score_list)


def ALAnate_protein(pose, sfxn, residue_index_list):
    """
    Mutates all positions into ALA if not (motif or disulfide)
    :param pose:
    :param reside_index_list:
    :return:
    """
    # --- Mutate all non-motif residues to ALA --- #
    ligand_index = pose.size()

    mutate = rosetta.protocols.simple_moves.MutateResidue()
    mutate.set_res_name('ALA')

    for res in range(1, ligand_index):
        if all([res not in residue_index_list,
                pose.residue(res).name3() not in ['GLY', 'PRO'],
                'disulfide' not in pose.residue(res).name()]):
            mutate.set_target(res)
            mutate.apply(pose)

    sfxn(pose)


def find_binding_motif_clashes(pose, sfxn, residue_index_list):
    """
    Returns the largest fa_rep value for any of the motif residues and ligand where all other residues in the match
    scaffold are mutated to ALA
    :param pose:
    :param residue_index_list:
    :return:
    """
    ligand_index = pose.size()
    motif_and_ligand_idx = residue_index_list + [ligand_index]

    sfxn(pose)  # Redundant
    edges = pose.energies().energy_graph()

    # --- Calculate max fa_rep for each motif residue --- #
    fa_rep_list = []
    fa_sol_list = []

    for motif_res in motif_and_ligand_idx:
        for res in range(1, ligand_index):
            if res not in motif_and_ligand_idx:
                current_edge = edges.find_energy_edge(motif_res, res)
                if current_edge is not None:
                    current_edge.fill_energy_map()
                    fa_rep_list.append(current_edge[rosetta.core.scoring.fa_rep])
                    fa_sol_list.append(current_edge[rosetta.core.scoring.fa_sol])

    return max(fa_rep_list), max(fa_sol_list)


def hbond_satisfaction_filter(pose, sfxn, residue_index_list):
    """
    Iterate over all positions within 10A of the ligand and find rotamers that make Hbonds with the ligand.
    :return:
    """
    ligand_index = pose.size()
    ligand_residue_selector = rosetta.core.select.residue_selector.ChainSelector('X')
    neighborhood_selector = rosetta.core.select.residue_selector.NeighborhoodResidueSelector(ligand_residue_selector, 10, False)

    neighborhood_selector_bool = neighborhood_selector.apply(pose)
    neighborhood_residues_resnums = rosetta.core.select.get_residues_from_subset(neighborhood_selector_bool)
    positions_to_consider = list(set(neighborhood_residues_resnums) - set (residue_index_list))
    print(positions_to_consider)

    # --- Generate backbone-dependent rotamer library for positions in neighborhood_selector --- #

    # Define packertask using neighborhood_selector
    packer_task = rosetta.core.pack.task.TaskFactory.create_packer_task(pose)
    packer_task.restrict_to_residues(neighborhood_selector_bool)

    # Only build rotamers for residues with Hbond donors/acceptors
    restrict_CAAs = rosetta.core.pack.task.operation.RestrictAbsentCanonicalAAS(0, rosetta.utility.vector1_bool(20))
    restrict_CAAs.keep_aas('DEHKNQSTY')
    restrict_CAAs.apply(pose, packer_task)

    # packer_neighbor_graph = create_packer_graph( pose, scfxn, task );
    packer_neighbor_graph = rosetta.core.pack.create_packer_graph(pose, sfxn, packer_task)
    satisfied_atoms = set()

    for residue in positions_to_consider:

        print(f'\nEvaluating position {residue}...\n')

        # Get residue identity
        residue_ID = pose.residue(residue).name3()

        # Built rotamers
        neighborghood_rotamer_set = rosetta.core.pack.rotamer_set.RotamerSetFactory.create_rotamer_set(pose)
        neighborghood_rotamer_set.set_resid(residue)
        neighborghood_rotamer_set.build_rotamers(pose, sfxn, packer_task, packer_neighbor_graph)

        # Apply rotamers to pose and score ligand Hbond
        for rot in range(0, neighborghood_rotamer_set.num_rotamers()):
            rotamer_index = rot + 1
            pose.replace_residue(residue, neighborghood_rotamer_set.rotamer(rotamer_index), False)

            # Get energies
            sfxn(pose)
            edges = pose.energies().energy_graph()
            current_edge = edges.find_energy_edge(residue, ligand_index)

            if current_edge is not None:
                edge_scores = current_edge.fill_energy_map()

                # Evaluate possible clashes with motif residues
                max_fa_rep, max_fa_sol = find_binding_motif_clashes(pose, sfxn, residue_index_list)

                # Get HBond term for residue with ligand
                hbbond_term = edge_scores[rosetta.core.scoring.hbond_sc]

                # Get full score
                gurobi_terms = sum([edge_scores[rosetta.core.scoring.fa_rep] * 0.55,
                                    edge_scores[rosetta.core.scoring.fa_atr],
                                    edge_scores[rosetta.core.scoring.hbond_sc],
                                    edge_scores[rosetta.core.scoring.fa_sol],
                                    edge_scores[rosetta.core.scoring.fa_elec]
                                    ]
                                   )

                if all([hbbond_term <= -0.2, max_fa_rep < 10, gurobi_terms < 10]):
                    contact_distance, ligand_contact_atom = report_contact_atoms(pose, residue)
                    satisfied_atoms.add(ligand_contact_atom)

                    print(f'Position {residue} rotamer {rotamer_index} looks good...')
                    print(contact_distance, ligand_contact_atom)

            # Mutate position back to ALA/GLY/PRO
            mutate = rosetta.protocols.simple_moves.MutateResidue()
            mutate.set_res_name(residue_ID)
            mutate.set_target(residue)
            mutate.apply(pose)

    # Get hbonds from motif residues
    # todo: fix this redundant code block...
    sfxn(pose)
    edges = pose.energies().energy_graph()

    for motif_index in residue_index_list:
        print(f'Motif residue {motif_index}...')
        current_edge = edges.find_energy_edge(residue, motif_index)

        if current_edge is not None:
            edge_scores = current_edge.fill_energy_map()
            hbbond_term = edge_scores[rosetta.core.scoring.hbond_sc]

            if hbbond_term <= -0.2:
                contact_distance, ligand_contact_atom = report_contact_atoms(pose, residue)
                print(f'Hydrogen bond found for motif residue {motif_index}:\tcontact_distance\tligand_contact_atom')
                satisfied_atoms.add(ligand_contact_atom)

    # todo: get backbone hbonds with ligand

    return list(satisfied_atoms)


def report_contact_atoms(pose, resnum):

    # --- Calculate shortest distance between residue and specific atoms on ligand --- #

    # 38E HBond donor/acceptor list
    hbond_atoms = ['O1', 'O2', 'N2']

    ligand_res = pose.residue(pose.size())  # ASSUMPTION: Ligand is last residue in pose
    hbond_residue = pose.residue(resnum)

    contact_distance = 999
    ligand_contact_atom = None

    for lig_atom in hbond_atoms:
        ligand_hbond_atom = ligand_res.atom(lig_atom)
        distances = [np.linalg.norm(ligand_hbond_atom.xyz() - res_atom.xyz()) for res_atom in hbond_residue.atoms()]
        current_contact = min(distances)

        if current_contact < contact_distance:
            contact_distance = current_contact
            ligand_contact_atom = lig_atom

    return contact_distance, ligand_contact_atom


def consolidate_match_stats():
    """
    Consolidates energy states csv files from separate taskid-XXXX directories into a single dataframe
    :return:
    """
    df = pd.DataFrame(columns=['ligand_sasa','match','max_fa_rep','max_fa_sol','satisfied_atoms','total_motif_energy','taskid'])

    for dir in os.scandir(os.getcwd()):
        if dir.name.startswith('taskid-') and dir.is_dir():
            filter_csv = os.path.join(os.getcwd(), dir.name, 'Matches_Filtered-Energy_Stats.csv')
            if os.path.exists(filter_csv):
                try:
                    temp_df = pd.read_csv(filter_csv)

                    taskid = int(dir.name.split('-')[1])
                    temp_df = temp_df.assign(taskid=taskid)
                    df = df.append(temp_df, ignore_index=True)

                except EmptyDataError:
                    continue

    return df


def move_quality_matches(df, src, dst='Matches-Energy_Filtered'):
    """
    Move.
    :return:
    """
    os.makedirs(dst, exist_ok=True)
    set_list = []

    passing_values = {'total_motif_energy': (-5, 'min'),
                      'max_fa_rep': (5, 'min'),
                      'max_fa_sol': (5, 'min'),
                      # 'satisfied_atoms': (2, 'max'),
                      'ligand_sasa': (50, 'min')
                      }

    print('Evaluating matches based on the following criteria:')

    for key, value in passing_values.items():
        print(key)

        if value[1] == 'min':
            passing_rows = df.loc[df[key] <= value[0]].index
        else:
            raise Exception

        set_list.append(set(passing_rows))

    # Get intersection of matches
    final_set = set.intersection(*set_list)
    final_df = df[df.index.isin(list(final_set))]

    # todo: add optional motif residue filter
    final_df['motif_res'] = final_df.apply(lambda x: len(re.split('-|\.', x['match'])[1].split('_')) - 2, axis=1)
    final_df = final_df.loc[final_df['motif_res'] == 5]

    # todo: only return best scoring ~25 matches for each match group
    final_df['unique_match_group'] = final_df.apply(lambda x: '_'.join(x['match'].split('_')[:3] + x['match'].split('_')[4:]), axis=1)

    score_filtered_sub_df = []
    for unique_match, sub_df in final_df.groupby('unique_match_group'):
        score_filtered_sub_df.append(sub_df.nsmallest(25, 'total_motif_energy'))
    final_df = pd.concat(score_filtered_sub_df)

    final_df.to_csv('Matches_Final-Energy_Filtered.csv')

    for index, row in final_df.iterrows():
        match = row['match']
        match_pdb_path = os.path.join(src, 'taskid-{0}'.format(row['taskid']), match)
        shutil.copy(os.path.join(src, match_pdb_path), os.path.join(dst, match))


def score_matches(match_dir, params_file):
    """
    Score things.
    :return:
    """
    pyrosetta.init(options=f"-extra_res_fa {params_file} -ex1 -ex2 -extrachi_cutoff 0")
    sfxn = rosetta.core.scoring.get_score_function()

    list_of_dicts = []
    print(match_dir)
    print([pdb for pdb in os.listdir(match_dir)])

    for match in [pdb for pdb in os.listdir(match_dir) if pdb.endswith('.pdb')]:
        try:
            residue_index_list = residue_indicies_from_match_name(match)

            target_pose = rosetta.core.pose.Pose()
            rosetta.core.import_pose.pose_from_file(target_pose, os.path.join(match_dir, match))
            sfxn(target_pose)

            # Calculate ligand SASA before turning everything to ALA
            ligand_sasa_calc = rosetta.core.scoring.sasa.SasaCalc()
            ligand_sasa_calc.calculate(target_pose)
            ligand_sasa = list(ligand_sasa_calc.get_residue_sasa())[-1]

            # Convert all non-essential positions to ALA
            ALAnate_protein(target_pose, sfxn, residue_index_list)

            # Determine max fa_rep and fa_sol for motif residues
            max_fa_rep, max_fa_sol = find_binding_motif_clashes(target_pose, sfxn, residue_index_list)

            if max_fa_rep >= 20:
                print(f'\nSkipping {match}: high fa_rep in binding motif...\n')
                os.remove(os.path.join(match_dir, match))
                continue

            print(f'\n{match}:\tfa_rep\t{max_fa_rep}\tfa_sol\t{max_fa_sol}\n')

            # If binding motif passes fa_rep and fa_sol filters, determine hbond designability
            # ASSUMPTION : find_binding_motif_clashes() mutated all binding site residues to ALA, or kept GLY/PRO
            satisfied_atoms = hbond_satisfaction_filter(target_pose, sfxn, residue_index_list)

            # Calculate binding motif energy as computed by Gurobi
            total_motif_energy = calculate_match_binding_motif_energy(target_pose, residue_index_list)

            if total_motif_energy > 0:
                print(f'\nSkipping {match}: binding motif energy with gurobi score terms > 0...\n')
                os.remove(os.path.join(match_dir, match))
                continue

            match_info = {'match': match,
                          'total_motif_energy': total_motif_energy,
                          'max_fa_rep': max_fa_rep,
                          'max_fa_sol': max_fa_sol,
                          'satisfied_atoms': ','.join(satisfied_atoms),
                          'ligand_sasa': ligand_sasa
                          }

            list_of_dicts.append(match_info)

        except Exception as e:
            print(f'Something went wrong here...\n{e}\n')

    if list_of_dicts:
        data = pd.DataFrame(list_of_dicts)
        print(data)
        data.to_csv(os.path.join(match_dir, 'Matches_Filtered-Energy_Stats.csv'), index=False)
        return True

    else:
        return False


def match_things(block_size, working_home_dir, working_temp_dir=None):
    """
    Contents from matcher_submission-condensed.py
    :return:
    """
    # SGE_ID and JOD_ID
    sge_task_id = 0
    if "SGE_TASK_ID" in os.environ.keys():
        # sge_task_id - 1 is to allow for easy indexing of lists...
        sge_task_id = int(os.environ["SGE_TASK_ID"]) - 1

    job_id = 0
    if "JOB_ID" in os.environ.keys():
        job_id = int(os.environ["JOB_ID"])

    print('Job id:', job_id)
    print('Task id:', sge_task_id)
    print('Block size:', block_size)

    # --- Start submitting tasks --- #

    matcher_arg_json = json.load(open(os.path.join(working_home_dir, 'matcher_argument_list.json'), 'r'))
    current_arg_block = matcher_arg_json['args'][(sge_task_id * block_size):((sge_task_id + 1) * block_size)]

    # Make directory for each job and switch into it... this is to prevent multiple simultaneous writes to matcher_score.sc
    task_id_dir = 'taskid-{0}'.format(sge_task_id + 1)
    working_temp_taskid_dir = os.path.join(working_temp_dir, task_id_dir)
    os.mkdir(working_temp_taskid_dir)
    os.chdir(working_temp_taskid_dir)

    constraint_json_name = os.path.join(working_home_dir, '{0}-constraint_blocks.json'.format(matcher_arg_json["target_compound"]))

    # --- Prepare directory for task-specific constraint files --- #

    if os.path.exists(constraint_json_name):
        print('{0} found, preparing directory for constraint files...'.format(constraint_json_name))
        temp_cstfiles_dir = os.path.join(working_temp_dir, 'cst_files')
        os.mkdir(temp_cstfiles_dir)
        # os.mkdir(os.path.join(matcher_arg_json['cst_dir'], task_id_dir))

    for block in current_arg_block:

        # For reference...
        # block = [conformer, constraint, current_scaffold, posfile_name, params_name] + [gridlig]

        # --- Generate constraint file from constraint json and dump into matcher_arg_json['cst_dir'] --- #
        cst_file_path = os.path.join(matcher_arg_json['cst_dir'], block[1])

        if os.path.exists(constraint_json_name):
            print('{0} found, generating constraint file...'.format(os.path.basename(constraint_json_name)))
            with open(constraint_json_name, 'r') as cst_json:
                # Get required blocks from constraint file name
                cst_json_json = json.load(cst_json)
                relevant_constraint_blocks = re.split('-|\.', block[1])[1].split('_')[1:]

                # Write generated cst_files to scratch/tmp-#/cst_files
                cst_file_path = os.path.join(temp_cstfiles_dir, block[1])

                with open(cst_file_path, 'w') as cst_file:
                    for cst_block_index in relevant_constraint_blocks:
                        cst_file.write('CST::BEGIN\n')
                        cst_file.write(cst_json_json[block[0]][cst_block_index])
                        cst_file.write('\n')
                        cst_file.write('CST::END\n')

        arg = ['/netapp/home/james.lucas/Packages/Rosetta/main/source/bin/match.linuxgccrelease',
               '-database',
               '/netapp/home/james.lucas/Packages/Rosetta/main/database',
               '-s',
               os.path.join(matcher_arg_json['scaffold_pdb_dir'], block[2]),
               '-match::lig_name',
               matcher_arg_json['target_compound'],
               '-match::scaffold_active_site_residues',
               os.path.join(matcher_arg_json['posfile_dir'], block[3]),
               '-match::geometric_constraint_file',
               cst_file_path,
               '-extra_res_fa',
               os.path.join(matcher_arg_json['params_dir'], block[4]),
               '-output_matches_per_group',
               '1',
               '-match:consolidate_matches',
               '-ex1',
               '-ex2',
               '-extrachi_cutoff',
               '0',
               '-use_input_sc',
               '-euclid_bin_size',
               '1',  # Roland : 1.5
               '-euler_bin_size',
               '10',  # Roland: 15
               '-bump_tolerance',
               '0.5',
               '-out::path',
               working_temp_taskid_dir,
               # matcher_arg_json['output_path'],
               '-match:output_format',
               'PDB',
               '-out:file:scorefile',  # match scores
               '-match_grouper',
               # 'SameSequenceGrouper',  # Two matches belong in the same group if their hits come from the same amino acids at the same scaffold build positions
               'SameRotamerComboGrouper',  # Two matches belong in the same group if their hits come from the same rotamers of the same amino acids at the same scaffold build positions
               '-mute',
               'protocols.idealize']

        if len(block) == 7:
            arg.append('-match::grid_boundary')
            arg.append(os.path.join(matcher_arg_json['gridlig_path'], block[5]))

        print(' '.join(arg))

        rosetta_process = subprocess.Popen(arg, cwd=os.getcwd())
        return_code = rosetta_process.wait()

        print('Task return code:', return_code, '\n')

        # Delete constraint file after it is done being used
        os.remove(cst_file_path)

if __name__ == '__main__':
    args = docopt.docopt(__doc__)
    working_home_dir = os.getcwd()

    # Get params file
    if args['<params_dir>']:
        ligand_ID = args['<ligand_ID>']
        if os.path.exists(os.path.join(args['<params_dir>'], f'{ligand_ID}_conformers.pdb')):
            scoring_params_file = os.path.abspath(os.path.join(args['<params_dir>'], f'{ligand_ID}.params'))
        else:
            scoring_params_file = os.path.abspath(os.path.join(args['<params_dir>'], f'{ligand_ID}_0001.params'))

    if args['match']:
        # --- Move into /scratch directory --- #
        os.chdir('/scratch')
        with tempfile.TemporaryDirectory(dir=os.getcwd()) as working_temp_dir:
            print(f'Temporary files are being written to: {working_temp_dir}')

            # --- Match --- #
            match_things(int(args['<block_size>']), working_home_dir, working_temp_dir=working_temp_dir)

            current_taskid = os.getcwd()  # Current working directory is /scratch/unique_temp_dir/taskid-#
            print(f'Current taskid: {current_taskid}')
            keep_dir = score_matches(current_taskid, scoring_params_file)

            # --- Writing to netapp/home/ directory, deletes taskid-# if none pass fa_rep filter --- #
            # if not keep_dir:
            #     os.chdir('..')
            #     shutil.rmtree(current_taskid)
            #     os.remove(f'{os.path.basename(__file__)}.o{int(os.environ["JOB_ID"])}.{int(os.environ["SGE_TASK_ID"])}')

            # --- Writing to /scratch directory, moves taskid-# to netapp/home/ if any match passes fa_rep filter --- #
            if keep_dir:
                destination_dir = os.path.join(working_home_dir, os.path.basename(current_taskid))
                print(f'Copying data from {current_taskid} to {destination_dir}...')
                shutil.copytree(current_taskid, destination_dir)
            else:
                print('Output was all crap...')

            # todo: clean copied matches to delete matches that fail obvious filter metrics

    if args['score']:
        score_matches(args['<input_PDB_dir>'], scoring_params_file)

    if args['score_matched']:
        os.chdir(os.path.join(os.getcwd(), f'task-{int(os.environ["SGE_TASK_ID"])}'))
        keep_dir = score_matches(os.getcwd(), scoring_params_file)

    if args['consolidate']:
        df = consolidate_match_stats()
        df.to_csv('derp.csv')

    if args['move']:
        df = consolidate_match_stats()
        move_quality_matches(df, os.getcwd())