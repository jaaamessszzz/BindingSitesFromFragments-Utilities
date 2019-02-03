#!/usr/bin/env python3

"""
Generate a library of monomeric/dimeric scaffolds for matching.

Usage:
    generate_scaffold_library (monomer|dimer) (download|<scaffold_dir>) [<trash_scaffolds>]

Arguments:     
    monomer
        generate a scaffold library of monomers
        
    dimer
        generate a scaffold library of dimers

    download
        Download scaffolds

    <scaffold_dir>
        Directory containing raw PDB scaffolds

    <trash_scaffolds>
        Text file containing PDBIDs of scaffolds to skip

Options:
    --wat -w
        wat.
"""
import gzip
import io
import os
import pprint
import shutil
import urllib

import docopt
import numpy as np
import prody
import xmltodict

from Matching.generate_gridlig_posfile_from_PDB import clean_pdb


# Pirated from bsff.fragments
def download_PDBs(monomer=False):
    """
    Download all PDBs with the provided fragment-containing ligand
    :param monomer: monomer if True, otherwise False
    :return: 
    """

    # Generate directory tree for pdbs, gridlig, and posfile
    new_dimer_scaffolds = 'Heterodimer_scaffolds-James'
    pdb_output_dir = os.path.join(new_dimer_scaffolds, 'PDBs')
    gridlig_output_path = os.path.join(new_dimer_scaffolds, 'gridlig')
    posfile_output_path = os.path.join(new_dimer_scaffolds, 'posfile')

    os.makedirs(pdb_output_dir, exist_ok=True)
    os.makedirs(gridlig_output_path, exist_ok=True)
    os.makedirs(posfile_output_path, exist_ok=True)

    for expression_organism in ['ESCHERICHIA COLI', 'saccharomyces cerevisiae'.upper(), 'pichia pastoris'.upper()]:

        # <description>Chain Type: there is a Protein chain but not any DNA or RNA or Hybrid</description>
        # <description>Stoichiometry in biological assembly: Stoichiometry is AB</description>
        # <description># of Disulfide link records: Min is 0 and Max is 4</description>

        # <queryRefinement>
        # <queryRefinementLevel>2</queryRefinementLevel>
        # <conjunctionType>and</conjunctionType>
        # <orgPdbQuery>
        # <queryType>org.pdb.query.simple.CloseContactsQuery</queryType>
        # <min>0</min>
        # <max>4</max>
        # </orgPdbQuery>
        # </queryRefinement>

        # <description>Oligomeric state Search : Min Number of oligomeric state=2 Max Number of oligomeric state=2</description>
        # <description>ExpressionOrganismQuery: entity_src_gen.pdbx_host_org_scientific_name.comparator=contains entity_src_gen.pdbx_host_org_scientific_name.value=ESCHERICHIA COLI </description>
        # <description>Resolution is between 1.0 and 3.0 </description>
        # <description>Sequence Length is between 30 and 400 </description>
        # <description>Representative Structures at 95% Sequence Identity</description>

        print("Searching for proteins expressed in {}".format(expression_organism))

        REST_search_xml = """
        <orgPdbCompositeQuery version="1.0">
        
        <queryRefinement>
        <queryRefinementLevel>0</queryRefinementLevel>
        <orgPdbQuery>
        <queryType>org.pdb.query.simple.ChainTypeQuery</queryType>
        <containsProtein>Y</containsProtein>
        <containsDna>N</containsDna>
        <containsRna>N</containsRna>
        <containsHybrid>N</containsHybrid>
        </orgPdbQuery>
        </queryRefinement>
        
        <queryRefinement>
        <queryRefinementLevel>1</queryRefinementLevel>
        <conjunctionType>and</conjunctionType>
        <orgPdbQuery>
        <queryType>org.pdb.query.simple.StoichiometryQuery</queryType>
        <stoichiometry>AB</stoichiometry>
        </orgPdbQuery>
        </queryRefinement>
        
        <queryRefinement>
        <queryRefinementLevel>3</queryRefinementLevel>
        <conjunctionType>and</conjunctionType>
        <orgPdbQuery>
        <version>head</version>
        <queryType>org.pdb.query.simple.BiolUnitQuery</queryType>
        <oligomeric_statemin>2</oligomeric_statemin>
        <oligomeric_statemax>2</oligomeric_statemax>
        </orgPdbQuery>
        </queryRefinement>
        
        <queryRefinement>
        <queryRefinementLevel>4</queryRefinementLevel>
        <conjunctionType>and</conjunctionType>
        <orgPdbQuery>
        <version>head</version>
        <queryType>org.pdb.query.simple.ExpressionOrganismQuery</queryType>
        <entity_src_gen.pdbx_host_org_scientific_name.comparator>contains</entity_src_gen.pdbx_host_org_scientific_name.comparator>
        <entity_src_gen.pdbx_host_org_scientific_name.value>{0}</entity_src_gen.pdbx_host_org_scientific_name.value>
        </orgPdbQuery>
        </queryRefinement>
        
        <queryRefinement>
        <queryRefinementLevel>5</queryRefinementLevel>
        <conjunctionType>and</conjunctionType>
        <orgPdbQuery>
        <queryType>org.pdb.query.simple.ResolutionQuery</queryType>
        <refine.ls_d_res_high.comparator>between</refine.ls_d_res_high.comparator>
        <refine.ls_d_res_high.max>3.0</refine.ls_d_res_high.max>
        </orgPdbQuery>
        </queryRefinement>
        
        <queryRefinement>
        <queryRefinementLevel>6</queryRefinementLevel>
        <conjunctionType>and</conjunctionType>
        <orgPdbQuery>
        <version>head</version>
        <queryType>org.pdb.query.simple.SequenceLengthQuery</queryType>
        <v_sequence.chainLength.min>25</v_sequence.chainLength.min>
        <v_sequence.chainLength.max>500</v_sequence.chainLength.max>
        </orgPdbQuery>
        </queryRefinement>
        
        <queryRefinement>
        <queryRefinementLevel>7</queryRefinementLevel>
        <conjunctionType>and</conjunctionType>
        <orgPdbQuery>
        <version>head</version>
        <queryType>org.pdb.query.simple.HomologueEntityReductionQuery</queryType>
        <identityCutoff>95</identityCutoff>
        </orgPdbQuery>
        </queryRefinement>
        
        </orgPdbCompositeQuery>
        """.format(expression_organism)

        search_xml = xmltodict.unparse(xmltodict.parse(REST_search_xml), pretty=False)
        search_request = urllib.request.Request('http://www.rcsb.org/pdb/rest/search', data=search_xml.encode())
        search_result = urllib.request.urlopen(search_request, data=None, timeout=300).read()

        result_pdbs = search_result.split()

        pprint.pprint(len(result_pdbs))
        pprint.pprint(result_pdbs)

def generate_scaffolds(scaffold_pdbs, scaffold_output_dir, trash_scaffolds=None, download=False):

    # Keep track of numbers
    total_returned_pdbs = len(scaffold_pdbs)
    acceptable_pdbs = 0
    total_scaffold_units = 0

    # Make separate directories for pdbs/posfiles/gridligs
    pdb_dump_path = os.path.join(scaffold_output_dir, 'PDBs')
    gridlig_dump_path = os.path.join(scaffold_output_dir, 'gridligs')
    posfile_dump_path = os.path.join(scaffold_output_dir, 'posfiles')

    os.makedirs(os.path.join(scaffold_output_dir, 'PDBs'), exist_ok=True)
    os.makedirs(os.path.join(scaffold_output_dir, 'gridligs'), exist_ok=True)
    os.makedirs(os.path.join(scaffold_output_dir, 'posfiles'), exist_ok=True)

    # Download PDB files
    for pdb in scaffold_pdbs:
        print(pdb)
        # pdb_split = pdb.decode("utf-8").split(':')
        pdb_split = pdb.split('.')
        pdb_name = pdb_split[0].upper()
        # pdb_biomol_index = pdb_split[1]

        # --- Check failure conditions --- #

        # Has this PDB already been processed
        already_processed = any([pdb_name in pdb for pdb in os.listdir(pdb_dump_path)])

        if already_processed:
            print('{} already processed!'.format(pdb_name))
            continue

        # Is the scaffold previously identified to be trash?
        if trash_scaffolds:
            with open(trash_scaffolds, 'r') as trashy_trash:
                trash_scaffold_list = [trash.strip().upper() for trash in trashy_trash]

            definitely_trash = pdb_name in trash_scaffold_list

            if definitely_trash:
                print('{} is definitely trash...'.format(pdb_name))
                continue

        # Can prody work with this PDB?
        try:
            pdb_prody = prody.parsePDB(pdb_name)
            pdb_prody_header = prody.parsePDB(pdb_name, header=True, model=0, meta=False)
        except Exception as e:
            print(e)
            continue

        # set of database identifier tuples for a given scaffold...
        processed_heterodimers = set()
        biomol_count = 0

        # Is the scaffold an undesireable classification?
        trash_proteins = ['IMMUNOGLOBULIN', 'IMMUNE', 'MEMBRANE', 'ANTIBODY', 'TOXIN']
        is_this_trash = any([trash_protein in pdb_prody_header['classification'] for trash_protein in trash_proteins])

        if is_this_trash:
            print('{} passed: scaffold is {}!'.format(pdb_name, pdb_prody_header['classification']))
            continue

        # --- Generate Scaffold --- #

        acceptable_pdbs += 1

        for biomolecular_complex_index in pdb_prody_header['biomoltrans']:

            total_scaffold_units += 1

            if len(pdb_prody_header['biomoltrans'][biomolecular_complex_index][0]) == 2:
                chain_1_mapping = pdb_prody_header['biomoltrans'][biomolecular_complex_index][0][0]
                chain_2_mapping = pdb_prody_header['biomoltrans'][biomolecular_complex_index][0][1]
            else:
                continue

            # Get actual chain names for selections
            chain_1_polymer = pdb_prody_header[chain_1_mapping]
            chain_2_polymer = pdb_prody_header[chain_2_mapping]

            chain_1_actual = chain_1_polymer.chid
            chain_2_actual = chain_2_polymer.chid

            # Get DB reference for each chain and check against processed_heterodimers
            # db_tuple = tuple(sorted((chain_1_polymer.dbrefs[0].accession, chain_2_polymer.dbrefs[0].accession)))

            biomol_count += 1
            # processed_heterodimers.add(db_tuple)

            selected_chains_prody = pdb_prody.select('chain {} {}'.format(chain_1_actual, chain_2_actual))
            cleaned_prody, removed_residues = clean_pdb(selected_chains_prody)

            temp_output_stream = io.StringIO()
            prody.writePDBStream(temp_output_stream, cleaned_prody)

            # Get rid of the goddamn remarks
            temp_output_stream.seek(0)
            temp_output_stream.readline()

            temp_final = io.StringIO()
            temp_final.write('REMARK {} Chains {} {}\n'.format(pdb_name, chain_1_actual, chain_2_actual))
            temp_final.write('REMARK Removed_Native_Residues: {}\n'.format(', '.join(removed_residues)))
            temp_final.write('REMARK Chain_{}_NativeSeq: {}\n'.format(chain_1_actual, chain_1_polymer.sequence))
            temp_final.write('REMARK Chain_{}_NativeSeq: {}\n'.format(chain_2_actual, chain_2_polymer.sequence))
            shutil.copyfileobj(temp_output_stream, temp_final)
            temp_final.seek(0)

            with gzip.open(os.path.join(pdb_dump_path, '{}.pdb{}.gz'.format(pdb_name, biomol_count)), mode='wt') as final_out:
                shutil.copyfileobj(temp_final, final_out)

            # Generate posfile
            # all residues with CB within 10A of interface minus glycines and prolines
            interface_residues_A = cleaned_prody.select('(name CB within 10 of chain {}) and not resname PRO GLY and not chain {}'.format(chain_1_actual, chain_1_actual))
            interface_residues_B = cleaned_prody.select('(name CB within 10 of chain {}) and not resname PRO GLY and not chain {}'.format(chain_2_actual, chain_2_actual))
            interface_resnums = sorted(list(interface_residues_A.getResnums()) + list(interface_residues_B.getResnums()))
            print(interface_resnums)

            with open(os.path.join(posfile_dump_path, '{}.pdb{}.pos'.format(pdb_name, biomol_count)), 'w') as posfile:
                posfile.write(' '.join([str(a) for a in interface_resnums]))

            # Generate gridlig... or not
            # Maybe 30A cube centered on interface residues??? Gridlig may be necessary to prevent matcher from crashing...
            all_interface_residues = cleaned_prody.select('resnum {}'.format(' '.join([str(a) for a in interface_resnums])))
            ligand_center = prody.calcCenter(all_interface_residues)

            gridlid_lower_corner = (ligand_center - np.array([15, 15, 15])).tolist()
            with open(os.path.join(gridlig_dump_path, '{}.pdb{}.gridlig'.format(pdb_name, biomol_count)), 'w') as gridlig:
                gridlig.write('NAME: gridlig\n')
                gridlig.write('BASE: {}\n'.format(' '.join([str(coord) for coord in gridlid_lower_corner])))
                gridlig.write('SIZE: 300 300 300\n')
                gridlig.write('LENGTH: 0.1 0.1 0.1')

    # Report
    print('\nCompleted scaffold library!\n')
    print('{0} of {1} PDBs were used to generate {2} unique heterodimer scaffolds'.format(acceptable_pdbs, total_returned_pdbs, total_scaffold_units))

def main():
    args = docopt.docopt(__doc__)

    download_monomers = args['monomer'] if args['monomer'] else args['dimer']
    scaffold_pdbs = download_PDBs(monomer=download_monomers) if args['download'] else os.listdir(args['<scaffold_dir>'])

    # Make scaffold directory
    scaffold_output_dir = os.path.join(os.getcwd(), 'Scaffold_Library_by_James')
    os.makedirs(scaffold_output_dir, exist_ok=True)

    generate_scaffolds(scaffold_pdbs, scaffold_output_dir, args['<trash_scaffolds>'])

if __name__ == '__main__':
    main()
