#!/usr/bin/env python3

"""
Fuzzball match iterations are taking way too long and don't reliably produce designable solutions.
Tanja suggested taking three-residue matches and using the rest of the fuzzball to generate design constraints for
the rest of the binding site. This is similar to what Will Hanson (Sagar Lab) suggested at RosettaCon, except Will
suggested fixed AA identities for each designable position while Tanja is suggesting weighted AA preferences.

Tanja suggests generating inverse rotamers for each fuzzball residue and using any hits to bias design preferences.

Need to figure out how to:
1) Generate inverse rotamers
2) Identify inverse rotamer matches
3) Define resfiles or some other way to weight designable residues

Usage:

    fuzzball_design <designable_match_pdb> <corresponding_fuzzball>

Arguments:

    <designable_match_pdb>
        Match PDB that will be used for design

    <corresponding_fuzzball>
        Fuzzball for which the binding motif was matched

Options:


"""

import docopt

def

def generate_inverse_rotamers():
    """
    Generate inverse rotamers
    :return:
    """
    pass

if __name__ == '__main__':
