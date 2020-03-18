#! /usr/bin/python3
"""
Cluster2PyMOL
Converts a cluster into a PyMOL file for generating figures.

Usage:
    Cluster2PyMOL <cluster>

Arguments:
    <cluster>           Path to cluster .ag.npz file

Options:
    None

"""

import os
import subprocess
import tempfile

import docopt
import prody

def generate_session(cluster):
    """
    derp.
    :return:
    """

    basepath, cluster_dir = os.path.split(cluster)
    fragment = os.path.basename(basepath)

    with tempfile.TemporaryDirectory() as tempdir:
        cluster_prody = prody.loadAtoms(cluster)
        residue_list = list()

        for residue in cluster_prody.iterResidues():
            residue_path = os.path.join(tempdir, f"{residue.getData('contact_source')[0]}.pdb")
            residue_list.append(residue_path)
            prody.writePDB(residue_path, residue)

        cmd_list = ['pymol', '-m'] + residue_list

        # --- Initialize --- #
        # Nice settings
        cmd_list += ['-d', 'set ray_trace_mode, 1']
        cmd_list += ['-d', 'set ray_texture, 1'] # Matte 1
        cmd_list += ['-d', 'set ray_shadows, 0'] # Shadows off
        cmd_list += ['-d', 'bg_color white'] # White background

        # Save and quit
        cluster_path, cluster_file = os.path.split(cluster)
        save_path = os.path.join(os.getcwd(), f"{fragment}-{os.path.splitext(cluster_file)[0]}.pse")
        cmd_list += ['-d', f'save {save_path}']
        cmd_list += ['-d', 'quit']

        make_session = subprocess.Popen(cmd_list)
        make_session.wait()

if __name__ == '__main__':
    args = docopt.docopt(__doc__)

    cluster = args['<cluster>']
    generate_session(cluster)
