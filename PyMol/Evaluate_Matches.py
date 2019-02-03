#!/usr/bin/env python3

"""
Creates Pymol sessions to evaluate matches

Usage:
    Evaluate_Matches single (<ideal_motif_dir> | fuzzball <fuzzball_dir>) <Match_PDB> [--monomer]
    Evaluate_Matches batch (<ideal_motif_dir> | fuzzball <fuzzball_dir>) <Match_Dir> [--monomer]

Arguments:
    single
        Create session for a single design

    batch
        Create a new session for each design in an input directory

    <ideal_motif_dir>
        Directory containing ideal binding motif PDBs

    <fuzzball_dir>
        Directory containing fuzzballs

    <Match_PDB>
        Path to match PDB

    <Match_Dir>
        Path to directory containing a bunch of matches

Options:
    --monomer -m
        Monomer matches
"""

import os
import sys
import subprocess
import json
import re

import docopt
import prody


def generate_session(match_pdb_path, ideal_motif_path):
    """
    Generate a pymol session...
    :param match_pdb_path:
    :param ideal_motif_path:
    :return:
    """
    match_pdb_base = os.path.basename(match_pdb_path).split('.')[0]
    binding_motif_base = os.path.basename(ideal_motif_path).split('.')[0]
    compound_id = binding_motif_base[:3]

    cmd_list = ['/Applications/PyMOL.app/Contents/MacOS/PyMOL', '-m', match_pdb_path, ideal_motif_path]

    # --- Initialize --- #
    # Nice settings
    cmd_list += ['-d', 'set ray_trace_mode, 1']
    cmd_list += ['-d', 'set ray_texture, 1']  # Matte 1
    cmd_list += ['-d', 'set ray_shadows, 0']  # Shadows off
    cmd_list += ['-d', 'bg_color white']  # White background
    cmd_list += ['-d', 'set_name {0}, Match'.format(match_pdb_base)]
    cmd_list += ['-d', 'set_name {0}, Binding_Motif'.format(binding_motif_base)]

    # --- Set up comparisons --- #
    # Align
    cmd_list += ['-d', 'align Binding_Motif and resn {0}, Match and resn {0}'.format(compound_id)]

    # Show motif residues
    motif_residues_list_str = [str(a) for a in residue_indicies_from_match_name(match_pdb_base)]
    print(match_pdb_base)
    print(motif_residues_list_str)
    cmd_list += ['-d', 'select Match-motif_residues, Match and resi {0}'.format('+'.join(motif_residues_list_str))]
    cmd_list += ['-d', 'show sticks, Match-motif_residues and (sidechain or name CA)']

    # Show ideal binding motif
    cmd_list += ['-d', 'show sticks, Binding_Motif and (sidechain or name CA)']
    cmd_list += ['-d', 'set stick_radius, 0.1, Binding_Motif']

    # Show CB within 2.5A of ligand
    cmd_list += ['-d', 'select clashing_cb, (Match and resn {0} and not hydrogen) around 2.5 and name CB'.format(compound_id)]
    cmd_list += ['-d', 'show spheres, clashing_cb']
    cmd_list += ['-d', 'set sphere_scale, 0.5']

    # --- Formatting --- #
    # Color everything
    match_prody = prody.parsePDB(match_pdb_path)
    match_chains = list(set(match_prody.getChids()) - set('X'))
    print(match_chains)

    chain_colors = ['marine', 'purple']

    for chain, color in zip(match_chains, chain_colors[:len(match_chains)]):
        cmd_list += ['-d', 'color {0}, Match* and chain {1} and name C*'.format(color, chain)]

    # Color specific selections
    cmd_list += ['-d', 'color yelloworange, Match-motif_residues and name C*']
    cmd_list += ['-d', 'color brightorange, Binding_Motif and name C*']
    cmd_list += ['-d', 'color density, resn {0} and name C*'.format(compound_id)]
    cmd_list += ['-d', 'color yelloworange, Binding_Motif and name C*']
    cmd_list += ['-d', 'color red, clashing_cb']

    # Hydrogens on ligand only
    cmd_list += ['-d', 'hide sticks, hydrogen']
    cmd_list += ['-d', 'show sticks, resn {}'.format(compound_id)]
    cmd_list += ['-d', 'color white, hydrogen']

    # Center
    cmd_list += ['-d', 'center resn {0}'.format(compound_id)]

    # Save and quit
    cmd_list += ['-d', 'save {0}'.format(os.path.join(os.getcwd(), '{0}-pretty.pse'.format(match_pdb_base)))]
    cmd_list += ['-d', 'quit']

    print(cmd_list)
    make_session = subprocess.Popen(cmd_list)
    make_session.wait()


def generate_ideal_motifs(fuzzball_dir, ideal_bs_dir, match_PDB_dir, monomer=False):
    """
    Generate ideal binding motifs for a given
    :param fuzzball_dir:
    :param ideal_bs_dir:
    :return:
    """
    conformer_set = set()
    ideal_motif_set = set()
    fuzzball_dict = {}

    for path, subdirs, files in os.walk(match_PDB_dir):
        # Go through match PDB directory
        for pdb in files:
            if pdb.startswith('UM'):
                # Parse filename to get relevant values
                conformer_name, constraint_resnums = find_conformer_and_constraint_resnums(pdb, monomer=monomer)

                # Load fuzzball if necessary
                if conformer_name not in conformer_set:
                    conformer_set.add(conformer_name)
                    fuzzball_dict[conformer_name] = prody.parsePDB(os.path.join(fuzzball_dir, '{0}-single_pose.pdb'.format(conformer_name)))

                motif_pdb_filename = '{}-{}.pdb'.format(conformer_name, '1_' + '_'.join([str(a) for a in constraint_resnums]))

                # Generate Motif PDB if not already in ideal_motif_set
                if motif_pdb_filename not in ideal_motif_set:
                    ideal_motif_set.add(motif_pdb_filename)

                    current_fuzzball = fuzzball_dict[conformer_name]

                    # Start binding motif with ligand
                    current_binding_motif = current_fuzzball.select('resnum 1')

                    # Add remaining residues
                    for resnum in constraint_resnums:
                        current_binding_motif = current_binding_motif + current_fuzzball.select('resnum {}'.format(resnum))

                    motif_pdb_filename = '{}-{}.pdb'.format(conformer_name, '1_' + '_'.join([str(a) for a in constraint_resnums]))
                    prody.writePDB(os.path.join(ideal_bs_dir, motif_pdb_filename), current_binding_motif)


def find_conformer_and_constraint_resnums(pdb_name, monomer=False):
    """
    Generates name of ideal binding motif from match PDB name
    :param pdb_name:
    :return:
    """
    pdb_split = re.split('_|-|\.', pdb_name)

    ligand_name_index = 5
    conformer_id_index = 6

    if len(pdb_split[4]) == 4:
        ligand_name_index += 1
        conformer_id_index += 1

    ligand_name = pdb_split[ligand_name_index]
    conformer_id = pdb_split[conformer_id_index]
    conformer_name = '{}_{}'.format(ligand_name, conformer_id)

    constraint_resnum_block = re.split('-|\.', pdb_name)[1][:-2]
    constraint_resnums = [int(a) for a in constraint_resnum_block.split('_') if a != ''][1:]

    return conformer_name, constraint_resnums


def pdb_check(dir, base_only=False):
    for file in os.listdir(dir):
        path = os.path.join(dir, file)
        if path.endswith('.pdb'):
            if base_only:
                yield file
            else:
                yield path


def residue_indicies_from_match_name(match_pdb_name):
    pnc = re.split('_|-|\.', os.path.basename(os.path.normpath(match_pdb_name)))
    motif_residue_ID_list = [a for a in re.split('(\D+)', pnc[2]) if a != '']
    motif_residue_IDs = [motif_residue_ID_list[indx + 1] for indx in range(0, len(motif_residue_ID_list), 2)]
    return motif_residue_IDs


if __name__ == '__main__':
    args = docopt.docopt(__doc__)

    # --- Parse Arguments --- #
    match_pdb_path = args['<Match_PDB>']
    match_pdb_dir = args['<Match_Dir>']
    fuzzball_dir = args['<fuzzball_dir>']
    monomer = args['--monomer']

    # --- Find/Generate Ideal Binding Motifs --- #
    # todo: update for Match_PDB option, only works for Match_Dir at the moment...
    if args['fuzzball']:
        ideal_motif_dir = 'Ideal_Binding_Motifs'
        ideal_bs_dir = os.path.join(os.getcwd(), ideal_motif_dir)
        os.makedirs(ideal_bs_dir, exist_ok=True)

        generate_ideal_motifs(fuzzball_dir, ideal_bs_dir, match_pdb_dir, monomer=monomer)

    else:
        ideal_motif_dir = args['<ideal_motif_dir>']

    # --- Generate Pymol Sessions --- #
    if args['single']:
        # Find ideal binding motif
        conformer_name, constraint_resnums = find_conformer_and_constraint_resnums(match_pdb_path, monomer=monomer)
        motif_pdb_filename = '{}-{}.pdb'.format(conformer_name, '1_' + '_'.join([str(a) for a in constraint_resnums]))
        generate_session(os.path.join(match_pdb_dir, match_pdb_path), os.path.join(ideal_motif_dir, motif_pdb_filename))

    if args['batch']:
        for file in os.listdir(match_pdb_dir):
            if file.endswith('.pdb'):
                if not os.path.exists('{}-pretty.pse'.format(file.split('.')[0])):

                    # Find ideal binding motif
                    conformer_name, constraint_resnums = find_conformer_and_constraint_resnums(file, monomer=monomer)
                    motif_pdb_filename = '{}-{}.pdb'.format(conformer_name, '1_' + '_'.join([str(a) for a in constraint_resnums]))
                    generate_session(os.path.join(match_pdb_dir, file), os.path.join(ideal_motif_dir, motif_pdb_filename))