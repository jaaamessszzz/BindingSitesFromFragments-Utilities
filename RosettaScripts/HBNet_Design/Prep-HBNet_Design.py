#! /usr/bin/python3
"""
Prep-HBNet_Design
This script will prepare inputs for Submit-HBNet_Design, namely design positions for HBNet. Settings will be exported
as a json that can be imported by the cluster submission script.

Design Positions - All residues within 10A of ligand except for matched residues
Static Positions - Not {Design Positions}

Usage:
    Prep-HBNet_Design s <matched_pdb_path>
    Prep-HBNet_Design l <matched_pdb_list>
    
Arguments:     
    s
        Generate inputs for a single binding motif match
    
    <matched_pdb_path>
        Path to the matched pdb
        
    l
        Generate inputs for a list of binding motif matches
    
    <matched_pdb_list>
        Path to a text file containing a list of matched pdb paths
        
Options:
    

"""
import docopt
import os
import prody
import re
import json

class Prepare_Designs():
    """
    Class for all your design preparation needs
    """

    def __init__(self, match_pdb_path):
        self.match_pdb_path = match_pdb_path
        self.match_position_list = self.determine_matched_residue_positions()
        self.match_prody = prody.parsePDB(self.match_pdb_path)
        self.compound_id = os.path.basename(os.path.normpath(self.match_pdb_path)).split('_')[5]
        self.design_position_list = self.determine_design_positions()

    def determine_design_positions(self):
        """
        Select all residues with CA within 10A the ligand
        :return: 
        """
        ca_shell_prody = self.match_prody.select('name CA and within 10 of resname {}'.format(self.compound_id))
        ca_shell = [atom.getResnum() for atom in ca_shell_prody]

        return set(ca_shell) - set([a[1] for a in self.match_position_list])

    def determine_matched_residue_positions(self):
        """
        Parse the filename of the match PDB to determine IDs and positions of match residues
        :return: 
        """
        positions_block = os.path.basename(os.path.normpath(self.match_pdb_path)).split('_')[2]
        resnames = [a for a in re.split("[0-9]*", positions_block) if a]
        resnums = [int(a) for a in re.split("[a-zA-Z]*", positions_block) if a]

        return [(a, b) for a, b in zip(resnames, resnums)]

    def dump_json(self):
        """
        Dump all arguments to json
        Since all residues in a pose are default designable and packable, I'm going to use the 
        PreventResiduesFromRepacking Task Operation to set all residues outside of the 10A design shell and the matched
        motif residues to what is essentially NATRO
        :return: 
        """

        static_residue_list = list(set(self.match_prody.getResnums()) - set(self.design_position_list))
        json_dict = {'pdbid': os.path.basename(os.path.normpath(self.match_pdb_path)),
                     'static_residue_list': ','.join([str(a) for a in static_residue_list])
                     }

        print(json_dict)
        # json.dump()


def main():
    args = docopt.docopt(__doc__)

    if args['s']:
        prepare_designs = Prepare_Designs(args['<matched_pdb_path>'])
        design_position_list = prepare_designs.determine_design_positions()
        prepare_designs.dump_json()

if __name__ == '__main__':
    main()