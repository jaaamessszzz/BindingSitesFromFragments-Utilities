#!/usr/bin/env python3

"""
Map special_rot counts onto design positions in a scaffold protein

Usage:
    special_rot <scaffold_dir> <stats_dir>

Arguments:
    <scaffold_dir>      Directory containing scaffolds to map special_rot counts onto
    <stats_dir>         Directory containing special_rot information for scaffolds
"""

import os
import sys
import subprocess
import json
import re
import pandas as pd

import docopt
import prody

def visualize_specialrot(scaffold_path, stats_df):
    """
    derp.
    :return:
    """
    # Consolidate scaffold special_rot stats per position
    scaffold_base = os.path.basename(scaffold_path)[:-4]
    cmd_list = ['pymol', '-m', scaffold_path]

    # --- Initialize --- #
    # Nice settings
    cmd_list += ['-d', 'set ray_trace_mode, 1']
    cmd_list += ['-d', 'set ray_texture, 1'] # Matte 1
    cmd_list += ['-d', 'set ray_shadows, 0'] # Shadows off
    cmd_list += ['-d', 'bg_color white'] # White background
    cmd_list += ['-d', f'set_name {scaffold_base}-clean, Scaffold']

    cmd_list += ['-d', 'show cartoon, all']
    cmd_list += ['-d', 'hide lines, all']
    cmd_list += ['-d', 'color white, all']

    cmd_list += ['-d', 'show sticks, chain X']
    cmd_list += ['-d', 'color skyblue, chain X']

    cmd_list += ['-d', 'set sphere_color, firebrick']
    cmd_list += ['-d', 'set sphere_scale, 0.5']

    cmd_list += ['-d', 'set label_color, black']
    cmd_list += ['-d', 'set label_size, 12']
    cmd_list += ['-d', 'set label_font_id, 13']
    cmd_list += ['-d', 'set label_position, (1.5,1.5,1.5)']

    for position, position_df in stats_df.groupby('Position'):
        cmd_list += ['-d', f'show spheres, resi {position} and name CA']
        total_generated = position_df['RotamersGenerated'].sum()
        total_accepted = position_df['Accepted'].sum()
        cmd_list += ['-d', f'label resi {position} and name CA, "{total_accepted}/{total_generated}"']

    # Set View Better...
    cmd_list += ['-d', 'orient Scaffold']

    # Save and quit
    cmd_list += ['-d', 'save {0}'.format(os.path.join(os.getcwd(), '{0}-annotated.pse'.format(scaffold_base)))]
    cmd_list += ['-d', 'quit']

    print(cmd_list)
    make_session = subprocess.Popen(cmd_list)
    make_session.wait()


if __name__ == '__main__':
    args = docopt.docopt(__doc__)

    scaffold_dir = args['<scaffold_dir>']
    stats_dir = args['<stats_dir>']

    for scaffold in os.listdir(scaffold_dir):
        if scaffold.endswith('.pdb'):
            scaffold_path = os.path.join(scaffold_dir, scaffold)
            scaffold_stats = os.path.join(stats_dir, scaffold[:-4] + '.csv')
            stats_df = pd.read_csv(scaffold_stats)
            visualize_specialrot(scaffold_path, stats_df)
