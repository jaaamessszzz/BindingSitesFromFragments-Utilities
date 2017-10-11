#!/usr/bin/env python3

"""
Generate a library of monomeric/dimeric scaffolds for matching.

Usage:
    generate_scaffold_library (monomer|dimer)

Arguments:     
    monomer
        generate a scaffold library of monomers
        
    dimer
        generate a scaffold library of dimers

Options:
    --wat -w
        wat.
"""
import docopt
import sys
import os
import shutil
import prody
import xmltodict
import urllib
import pprint
from generate_gridlig_posfile_from_PDB import clean_pdb
import io
import gzip

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

    for expression_organism in ['ESCHERICHIA COLI', 'saccharomyces cerevisiae'.upper()]:

        # <description>Chain Type: there is a Protein chain but not any DNA or RNA or Hybrid</description>
        # <description>Stoichiometry in biological assembly: Stoichiometry is AB</description>
        # <description># of Disulfide link records: Min is 0 and Max is 2</description>
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
        <queryRefinementLevel>2</queryRefinementLevel>
        <conjunctionType>and</conjunctionType>
        <orgPdbQuery>
        <queryType>org.pdb.query.simple.CloseContactsQuery</queryType>
        <min>0</min>
        <max>2</max>
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

        if len(result_pdbs) > 0:
            # Download PDB files
            # todo: turn result_pdbs into a set, get heterodimer information from biomoltrans
            for pdb in result_pdbs:
                pdb_split = pdb.decode("utf-8").split(':')
                pdb_name = pdb_split[0].upper()
                pdb_biomol_index = pdb_split[1]

                pdb_prody = prody.parsePDB(pdb_name)
                pdb_prody_header = prody.parsePDB(pdb_name, header=True, model=0, meta=False)

                # set of database identifier tuples for a given scaffold
                processed_heterodimers = set()
                biomol_count = 0

                for biomolecular_complex_index in pdb_prody_header['biomoltrans']:
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
                    db_tuple = tuple(sorted((chain_1_polymer.dbrefs[0].accession, chain_2_polymer.dbrefs[0].accession)))

                    if db_tuple not in processed_heterodimers:
                        biomol_count += 1
                        processed_heterodimers.add(db_tuple)

                        selected_chains_prody = pdb_prody.select('chain {} {}'.format(chain_1_actual, chain_2_actual))
                        cleaned_prody, removed_residues = clean_pdb(selected_chains_prody)

                        temp_output_stream = io.StringIO()
                        prody.writePDBStream(temp_output_stream, cleaned_prody)

                        # Get rid of the goddamn remarks
                        temp_output_stream.seek(0)
                        temp_output_stream.readline()

                        temp_final = io.StringIO()
                        temp_final.write('REMARK {} Chains {} {}\n'.format(pdb_name, chain_1_actual, chain_2_actual))
                        temp_final.write('REMARK Removed_Native_Residues: {}\n'.format(removed_residues))
                        temp_final.write('REMARK Chain_{}_NativeSeq: {}\n'.format(chain_1_actual, chain_1_polymer.sequence))
                        temp_final.write('REMARK Chain_{}_NativeSeq: {}\n'.format(chain_2_actual, chain_2_polymer.sequence))
                        shutil.copyfileobj(temp_output_stream, temp_final)
                        temp_final.seek(0)

                        with gzip.open(os.path.join(pdb_output_dir, '{}.pdb{}.gz'.format(pdb_name, biomol_count)), mode='wt') as final_out:
                            shutil.copyfileobj(temp_final, final_out)

                    else:
                        continue

                    # Generate posfile
                    # all residues with CB within 10A of interface minus glycines and prolines
                    interface_residues_A = cleaned_prody.select('(name CB within 10 of chain {}) and not resname PRO and not chain {}'.format(chain_1_actual, chain_1_actual))
                    interface_residues_B = cleaned_prody.select('(name CB within 10 of chain {}) and not resname PRO and not chain {}'.format(chain_2_actual, chain_2_actual))
                    interface_resnums = sorted(list(interface_residues_A.getResnums()) + list(interface_residues_B.getResnums()))
                    print(interface_resnums)

                    with open(os.path.join(posfile_output_path, '{}.pdb{}.pos'.format(pdb_name, biomol_count)), 'w') as posfile:
                        posfile.write(' '.join([str(a) for a in interface_resnums]))

                    # Generate gridlig... or not
                    # Nah.ge

def main():
    args = docopt.docopt(__doc__)

    if args['monomer']:
        download_PDBs(monomer=True)

    if args['dimer']:
        download_PDBs(monomer=False)

if __name__ == '__main__':
    main()
