#!/usr/bin/env python3

"""
Creates Pymol sessions in the style of the Gen9 designs for Tanja to look at

Usage:
    Gen9_Style_Sessions single <Design_json> <Binding_Motif> <Scaffold_PDB> <Design_PDB>
    Gen9_Style_Sessions batch <Design_json> <Binding_Motif> <Scaffold_PDB> <Design_Dir>

Arguments:
    single
        Create session for a single design

    batch
        Create a new session for each design in an input directory

    <Scaffold_PDB>
        PDB to use as template

    <Design_PDB>
        PDB to create session for

    <Design_Dir>
        Directory containing design

    <Binding_Motif>
        Path to ideal binding motif generated with Gurobi

    <Design_json>
        JSON containing all relevant info for designs
"""

import os
import sys
import subprocess
import json
import re

import docopt
import prody

class pymol_sessions():
    def __init__(self, design_json_path, match_pdb_path, binding_motif_path):
        self.design_json = json.load(open(design_json_path, 'r'))
        self.match_pdb_path = match_pdb_path
        self.binding_motif_path = binding_motif_path

    def generate_session(self, design_pdb_path):
        """
        derp.
        :return:
        """
        match_pdb_base = os.path.basename(self.match_pdb_path).split('.')[0]
        design_pdb_base = os.path.basename(design_pdb_path).split('.')[0]
        binding_motif_base = os.path.basename(self.binding_motif_path).split('.')[0]
        compound_id = binding_motif_base[:3]

        cmd_list = ['/Applications/PyMOL.app/Contents/MacOS/PyMOL', '-m', self.match_pdb_path, design_pdb_path, self.binding_motif_path]

        # --- Initialize --- #
        # Nice settings
        cmd_list += ['-d', 'set ray_trace_mode, 1']
        cmd_list += ['-d', 'set ray_texture, 1'] # Matte 1
        cmd_list += ['-d', 'set ray_shadows, 0'] # Shadows off
        cmd_list += ['-d', 'bg_color white'] # White background
        cmd_list += ['-d', 'set_name {0}, Scaffold'.format(match_pdb_base)]
        cmd_list += ['-d', 'set_name {0}, Design'.format(design_pdb_base)]
        cmd_list += ['-d', 'set_name {0}, Binding_Motif'.format(binding_motif_base)]

        # --- Set up comparisons --- #
        # Align
        cmd_list += ['-d', 'align Design, Scaffold']
        cmd_list += ['-d', 'align Binding_Motif and resn {0}, Design and resn {0}'.format(compound_id)]

        # Show designable residues
        designable_residues_list_str = [str(a) for a in self.design_json['design_residue_list']]
        cmd_list += ['-d', 'select Design-designable_residues, Design and resi {0}'.format('+'.join(designable_residues_list_str))]
        cmd_list += ['-d', 'select Scaffold-designable_residues, Scaffold and resi {0}'.format('+'.join(designable_residues_list_str))]
        cmd_list += ['-d', 'show sticks, Design-designable_residues and (sidechain or name CA)']
        cmd_list += ['-d', 'show sticks, Scaffold-designable_residues and (sidechain or name CA)']

        # Show motif residues
        motif_residues_list_str = [str(a[1]) for a in self.design_json['match_residues']]
        cmd_list += ['-d', 'select Design-motif_residues, Design and resi {0}'.format('+'.join(motif_residues_list_str))]
        cmd_list += ['-d', 'show sticks, Design-motif_residues and (sidechain or name CA)']

        cmd_list += ['-d', 'select Scaffold-motif_residues, Scaffold and resi {0}'.format('+'.join(motif_residues_list_str))]
        cmd_list += ['-d', 'show sticks, Scaffold-motif_residues and (sidechain or name CA)']

        # Show ideal binding motif
        cmd_list += ['-d', 'show sticks, Binding_Motif and (sidechain or name CA)']
        cmd_list += ['-d', 'set stick_radius, 0.1, Binding_Motif']

        # --- Show unsats --- #
        unsat_list = self.get_unsats(design_pdb_path)
        for unsat in unsat_list:
            cmd_list += ['-d', 'show sticks, Design and resi {0}'.format(unsat[0])]
            cmd_list += ['-d', 'show spheres, Design and resi {0} and name {1}'.format(unsat[0], unsat[1])]

        all_unsat_selection_string = ['(Design and resi {0} and name {1})'.format(unsat[0], unsat[1]) for unsat in unsat_list]
        cmd_list += ['-d', 'select unsat_atoms, {0}'.format(' or '.join(all_unsat_selection_string))]
        cmd_list += ['-d', 'set sphere_scale, 0.5']

        # --- Set design env objects --- #
        design_env_residues = [str(a) for a in self._design_env(design_pdb_path)]

        cmd_list += ['-d', 'create Design_Env, Design and resi {0}'.format('+'.join(design_env_residues))]
        cmd_list += ['-d', 'show sticks, Design_Env']
        cmd_list += ['-d', 'hide cartoon, Design_Env']
        cmd_list += ['-d', 'disable Design_Env']

        cmd_list += ['-d', 'create Scaffold_Env, Scaffold and resi {0}'.format('+'.join(design_env_residues))]
        cmd_list += ['-d', 'show sticks, Scaffold_Env']
        cmd_list += ['-d', 'hide cartoon, Scaffold_Env']
        cmd_list += ['-d', 'disable Scaffold_Env']

        # --- Design binding site surface --- #
        cmd_list += ['-d', 'select (Design and resn {0}) around 4 and Design'.format(compound_id)]
        cmd_list += ['-d', 'show surface, sele']
        cmd_list += ['-d', 'set transparency, 0.5']

        # --- Formatting --- #
        # Color everything
        cmd_list += ['-d', 'color white, Scaffold* and chain A and name C*']
        cmd_list += ['-d', 'color gray40, Scaffold* and chain B and name C*']
        cmd_list += ['-d', 'color marine, Design* and chain A and name C*']
        cmd_list += ['-d', 'color purple, Design* and chain B and name C*']

        # Color specific selections
        cmd_list += ['-d', 'color yelloworange, Design-motif_residues and name C*']
        cmd_list += ['-d', 'color brightorange, Binding_Motif and name C*']
        cmd_list += ['-d', 'color density, resn {0} and name C*'.format(compound_id)]
        cmd_list += ['-d', 'color yelloworange, Binding_Motif and name C*']

        # Hydrogens on ligand only
        cmd_list += ['-d', 'hide sticks, hydrogen']
        cmd_list += ['-d', 'show sticks, resn {}'.format(compound_id)]
        cmd_list += ['-d', 'color white, hydrogen']

        # Set View
        # my_view = """set_view (\
        # -0.546047747,   -0.557899892,   -0.624951005,\
        # -0.785422504,    0.600424647,    0.150261030,\
        #  0.291410387,    0.572908461,   -0.766049445,\
        #  0.001221966,   -0.000024810,  -81.042785645,\
        # 70.940414429,   21.246437073,    3.970802307,\
        # 64.890853882,   97.411834717,  -20.000000000 )
        # """
        # cmd_list += ['-d', my_view]

        # Set View Better...
        cmd_list += ['-d', 'orient Binding_Motif']

        # Save and quit
        cmd_list += ['-d', 'save {0}'.format(os.path.join(os.getcwd(), '{0}-pretty.pse'.format(design_pdb_base)))]
        cmd_list += ['-d', 'quit']

        print(cmd_list)
        make_session = subprocess.Popen(cmd_list)
        make_session.wait()

    def _design_env(self, design_pdb_path):
        # Get residues within 5A of design/motif residues that can potentially have unsats from designs
        input_PDB_prody = prody.parsePDB(design_pdb_path)

        designable_and_motif_positions = [a[1] for a in self.design_json['match_residues']] + self.design_json['design_residue_list']
        potential_unsat_set = set(designable_and_motif_positions)

        for residue in designable_and_motif_positions:
            ca_resnum_list = [res.getResnum() for res in
                              input_PDB_prody.select('protein within 5 of resnum {}'.format(residue))]
            potential_unsat_set.update(ca_resnum_list)

        return potential_unsat_set

    def get_unsats(self, design_pdb_path):
        """
        Parse design PDBs to get unsats
        :return: unsat_list of (resnum, atom name) tuples
        """
        # Ripped from design_filtering

        potential_unsat_set = self._design_env(design_pdb_path)
        unsat_list = []

        with open(design_pdb_path, 'r') as pdb_file:
            for line in pdb_file:

                # HBUnsats gets written several times... the "real" value is the last one written to file
                if line.startswith('BuriedUnsatHbonds buried_unsat:'):
                    unsat_list = []

                if line.startswith('      Unsatisfied H'):
                    split_line = [a for a in re.split('[ :]', line.strip()) if a]

                    # Only count unsats if they are in potential_unsat_set
                    if int(split_line[6]) in potential_unsat_set:
                        if 'H' in split_line[-1] and split_line[1] == 'Hpol':
                            continue

                        else:
                            unsat_list.append((split_line[6], split_line[8]))

        return unsat_list

if __name__ == '__main__':
    args = docopt.docopt(__doc__)

    match_pdb_path = args['<Scaffold_PDB>']
    design_json_path = args['<Design_json>']
    binding_motif_path = args['<Binding_Motif>']

    do_pymol_things = pymol_sessions(design_json_path, match_pdb_path, binding_motif_path)

    if args['single']:
        do_pymol_things.generate_session(args['<Design_PDB>'])

    if args['batch']:
        design_pdb_dir = args['<Design_Dir>']
        for file in os.listdir(design_pdb_dir):
            if file.endswith('.pdb'):
                print(file)
                do_pymol_things.generate_session(os.path.join(design_pdb_dir, file))