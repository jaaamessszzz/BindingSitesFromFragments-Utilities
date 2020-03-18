#!/usr/bin/env python3

"""
Generate hypothetical ligand binding pockets for a given ligand based on experimentally observed protein-ligand
interactions in the PDB.

Usage:  python3 <command> [<args>...]

Supported Supplemental Materials:

    FragmentPlots
    DesignPairPlots
    DesignViolinPlots
    DesignMetricsGrid
    BBSequenceLogos
    ProfileSimilarity
    SequenceRecovery

Options:
    -h --help
"""

import re
import os
import sys
from pprint import pprint

import docopt
import pandas as pd
from jinja2 import Environment, FileSystemLoader

def AllFragmentPlots(args):
    """
    Generate Fragment Plot supplemental materials

    Usage:  python3 FragmentPlots <plot_dir> <fragment_csv> [options]

    Arguments:
      <plot_dir>            Directory containing fragment plots
      <fragment_csv>        Path to Fragment_Inputs.csv

    """

    plot_dir = args['<plot_dir>']
    fragment_inputs = args['<fragment_csv>']

    template_dir = '/home/james/Repositories/BindingSitesFromFragments-Utilities/LatexSuppMat/Templates'
    template_name = 'AllFragmentPlots.tex'

    # http://eosrei.net/articles/2015/11/latex-templates-python-and-jinja2-generate-pdfs
    latex_jinja_env = Environment(
        block_start_string='\BLOCK{',
        block_end_string='}',
        variable_start_string='\VAR{',
        variable_end_string='}',
        comment_start_string='\#{',
        comment_end_string='}',
        line_statement_prefix='%%',
        line_comment_prefix='%#',
        trim_blocks=True,
        autoescape=False,
        loader=FileSystemLoader(template_dir))

    template = latex_jinja_env.get_template(template_name)

    # Assemble dict
    plot_dir_files = os.listdir(plot_dir)
    figure_set = set(fragment.split('-')[0] for fragment in plot_dir_files if fragment.endswith('.png'))
    fragment_sorted = sorted(figure_set, key=lambda x: int(x.split('_')[1]))
    figure_dict_temp = {}
    figure_dict_temp['Fragments'] = {fragment: {} for fragment in fragment_sorted}

    # SMILES
    df = pd.read_csv(fragment_inputs, index_col='Fragment')

    for fragment in fragment_sorted:
        fragment_index = int(fragment.split('_')[1])
        fragment_smiles = df.at[fragment_index, 'SMILES_fragment']
        fragment_recovery = os.path.join(os.path.abspath(plot_dir), f'{fragment}-cluster_recovery-bars.png')
        fragment_occupancy = os.path.join(os.path.abspath(plot_dir), f'{fragment}-cluster_occupancy.png')

        if os.path.exists(os.path.join(plot_dir, fragment_recovery)):
            figure_dict_temp['Fragments'][fragment]['recovery'] = fragment_recovery
        if os.path.exists(os.path.join(plot_dir, fragment_occupancy)):
            figure_dict_temp['Fragments'][fragment]['occupancy'] = fragment_occupancy
        figure_dict_temp['Fragments'][fragment]['SMILES'] = fragment_smiles

    # Lazy enumerate fragments
    figure_dict = {}
    figure_dict['Fragments'] = {}

    for index, fragment in enumerate(fragment_sorted, start=1):
        figure_dict['Fragments'][f'Fragment {index}'] = figure_dict_temp['Fragments'][fragment]

    figure_dict_temp['graphicspath'] = os.path.abspath(plot_dir)

    # Generate tex
    with open(f'SupplementalMaterials-AllFragmentPlots.tex', 'w') as asdf:
        asdf.write(template.render(figure_dict=figure_dict))


def DesignPairPlots(args):
    """
        Generate Fragment Plot supplemental materials

        Usage:  python3 DesignPairPlots <plot_dir> [options]

        Arguments:
          <plot_dir>            Directory containing design pairplots plots

        """

    plot_dir = args['<plot_dir>']

    template_dir = '/home/james/Repositories/BindingSitesFromFragments-Utilities/LatexSuppMat/Templates'
    template_name = 'DesignPairPlots.tex'

    # http://eosrei.net/articles/2015/11/latex-templates-python-and-jinja2-generate-pdfs
    latex_jinja_env = Environment(
        block_start_string='\BLOCK{',
        block_end_string='}',
        variable_start_string='\VAR{',
        variable_end_string='}',
        comment_start_string='\#{',
        comment_end_string='}',
        line_statement_prefix='%%',
        line_comment_prefix='%#',
        trim_blocks=True,
        autoescape=False,
        loader=FileSystemLoader(template_dir))

    template = latex_jinja_env.get_template(template_name)

    # Assemble dict
    plot_dir_files = sorted([design for design in os.listdir(plot_dir) if design.endswith('Subplots.png')])
    print(plot_dir_files)
    figure_set = set(design.split('-')[0] for design in plot_dir_files)

    figure_dict = {}
    for design in figure_set:
        figure_dict[design] = {}
        figure_dict[design]['plots'] = {}

    weight_dict = {'00': 0.0,
                   '15': -1.5,
                   '30': -3.0,
                   '40': -4.0}

    for design in plot_dir_files:
        design_weight_raw = re.split('_|-', design)[-2]
        design_setname = design.split('-')[0]

        figure_dict[design_setname]['plots'][weight_dict[design_weight_raw]] = os.path.join(os.path.abspath(plot_dir), design)

    pprint(figure_dict)
    # Generate tex
    with open(f'SupplementalMaterials-DesignPairPlots.tex', 'w') as asdf:
        asdf.write(template.render(figure_dict=figure_dict))


def DesignViolinPlots(args):
    """
    Generate Fragment Plot supplemental materials

    Usage:  python3 DesignViolinPlots <plot_dir> [options]

    Arguments:
      <plot_dir>            Directory containing design Violin plots
    """
    plot_dir = args['<plot_dir>']

    template_dir = '/home/james/Repositories/BindingSitesFromFragments-Utilities/LatexSuppMat/Templates'
    template_name = 'DesignViolinPlots.tex'

    # http://eosrei.net/articles/2015/11/latex-templates-python-and-jinja2-generate-pdfs
    latex_jinja_env = Environment(
        block_start_string='\BLOCK{',
        block_end_string='}',
        variable_start_string='\VAR{',
        variable_end_string='}',
        comment_start_string='\#{',
        comment_end_string='}',
        line_statement_prefix='%%',
        line_comment_prefix='%#',
        trim_blocks=True,
        autoescape=False,
        loader=FileSystemLoader(template_dir))

    template = latex_jinja_env.get_template(template_name)

    # Assemble dict
    plot_dir_files = sorted([design for design in os.listdir(plot_dir) if design.endswith('.png') and design.startswith('Maintext')])
    figure_dict = {}

    for design in plot_dir_files:
        design_setname = design.split('-')[1]
        figure_dict[design_setname] = os.path.join(os.path.abspath(plot_dir), design)

    pprint(figure_dict)
    # Generate tex
    with open(f'SupplementalMaterials-DesignViolinPlots.tex', 'w') as asdf:
        asdf.write(template.render(figure_dict=figure_dict))


def DesignMetricsGrid(args):
    """
    Generate Fragment Plot supplemental materials

    Usage:  python3 DesignMetricsGrid <plot_dir> [options]

    Arguments:
      <plot_dir>            Directory containing design grid plots
    """
    plot_dir = args['<plot_dir>']

    template_dir = '/home/james/Repositories/BindingSitesFromFragments-Utilities/LatexSuppMat/Templates'
    template_name = 'DesignMetricsGrids.tex'

    # http://eosrei.net/articles/2015/11/latex-templates-python-and-jinja2-generate-pdfs
    latex_jinja_env = Environment(
        block_start_string='\BLOCK{',
        block_end_string='}',
        variable_start_string='\VAR{',
        variable_end_string='}',
        comment_start_string='\#{',
        comment_end_string='}',
        line_statement_prefix='%%',
        line_comment_prefix='%#',
        trim_blocks=True,
        autoescape=False,
        loader=FileSystemLoader(template_dir))

    template = latex_jinja_env.get_template(template_name)

    # Assemble dict
    plot_dir_files = sorted([design for design in os.listdir(plot_dir) if design.endswith('-designmetrics.png') and design.startswith('UM')])
    figure_dict = {}

    for design in plot_dir_files:
        design_setname = design.split('-')[0]
        figure_dict[design_setname] = os.path.join(os.path.abspath(plot_dir), design)

    pprint(figure_dict)
    # Generate tex
    with open(f'SupplementalMaterials-DesignMetricsGrids.tex', 'w') as asdf:
        asdf.write(template.render(figure_dict=figure_dict))


def BBSequenceLogos(args):
    """
    Generate Fragment Plot supplemental materials

    Usage:  python3 BBSequenceLogos <logo_dir> [options]

    Arguments:
      <logo_dir>            Directory containing binding benchmark sequence logos
    """
    logo_dir = args['<logo_dir>']

    template_dir = '/home/james/Repositories/BindingSitesFromFragments-Utilities/LatexSuppMat/Templates'
    template_name = 'BBSequenceLogos.tex'

    # http://eosrei.net/articles/2015/11/latex-templates-python-and-jinja2-generate-pdfs
    latex_jinja_env = Environment(
        block_start_string='\BLOCK{',
        block_end_string='}',
        variable_start_string='\VAR{',
        variable_end_string='}',
        comment_start_string='\#{',
        comment_end_string='}',
        line_statement_prefix='%%',
        line_comment_prefix='%#',
        trim_blocks=True,
        autoescape=False,
        loader=FileSystemLoader(template_dir))

    template = latex_jinja_env.get_template(template_name)

    # Assemble dict
    logo_dir_files = [logo for logo in os.listdir(logo_dir) if logo.endswith('.pdf')]
    logo_set = set(logo.split('_')[0] for logo in logo_dir_files if 'Native' not in logo)

    figure_dict = {logo: {} for logo in logo_set}

    for logo in logo_dir_files:
        complex = logo[:4]
        weight = 'Native' if 'Native' in logo else logo.split('_')[1][:-4]  # lol.

        figure_dict[complex][weight] = os.path.join(os.path.abspath(logo_dir), logo)

    # Add order list
    for complex in figure_dict:
        figure_dict[complex]['order'] = ['0.0', '-0.5', '-1.0', '-1.5', '-2.0', '-2.5', '-3.0', '-3.5', '-4.0',
                                         '-4.5', '-5.0', '-6.0', '-8.0', '-10.0', 'Control', 'Native']

    pprint(figure_dict)
    # Generate tex
    with open(f'SupplementalMaterials-BBSequenceLogos.tex', 'w') as asdf:
        asdf.write(template.render(figure_dict=figure_dict))


def ProfileSimilarity(args):
    """
    Generate ProfileSimilarity supplemental materials

    Usage:  python3 ProfileSimilarity <plot_dir> [options]

    Arguments:
      <plot_dir>            Directory containing binding benchmark ProfileSimilarity plots
    """
    plot_dir = args['<plot_dir>']

    template_dir = '/home/james/Repositories/BindingSitesFromFragments-Utilities/LatexSuppMat/Templates'
    template_name = 'ProfileSimilarityPlots.tex'

    # http://eosrei.net/articles/2015/11/latex-templates-python-and-jinja2-generate-pdfs
    latex_jinja_env = Environment(
        block_start_string='\BLOCK{',
        block_end_string='}',
        variable_start_string='\VAR{',
        variable_end_string='}',
        comment_start_string='\#{',
        comment_end_string='}',
        line_statement_prefix='%%',
        line_comment_prefix='%#',
        trim_blocks=True,
        autoescape=False,
        loader=FileSystemLoader(template_dir))

    template = latex_jinja_env.get_template(template_name)

    # Assemble dict
    plot_dir_files = sorted([design for design in os.listdir(plot_dir) if design.endswith('.png') and design.startswith('ProfileSimilarity-Variable-')])
    figure_dict = {}

    for design in plot_dir_files:
        design_setname = design[:-4].split('-')[-1].upper()
        figure_dict[design_setname] = os.path.join(os.path.abspath(plot_dir), design)

    pprint(figure_dict)
    # Generate tex
    with open(f'SupplementalMaterials-ProfileSimilarity.tex', 'w') as asdf:
        asdf.write(template.render(figure_dict=figure_dict))


def SequenceRecovery(args):
    """
    Generate SequenceRecovery supplemental materials

    Usage:  python3 SequenceRecovery <plot_dir> [options]

    Arguments:
      <plot_dir>            Directory containing binding benchmark SequenceRecovery plots
    """
    plot_dir = args['<plot_dir>']

    template_dir = '/home/james/Repositories/BindingSitesFromFragments-Utilities/LatexSuppMat/Templates'
    template_name = 'SequenceRecoveryPlots.tex'

    # http://eosrei.net/articles/2015/11/latex-templates-python-and-jinja2-generate-pdfs
    latex_jinja_env = Environment(
        block_start_string='\BLOCK{',
        block_end_string='}',
        variable_start_string='\VAR{',
        variable_end_string='}',
        comment_start_string='\#{',
        comment_end_string='}',
        line_statement_prefix='%%',
        line_comment_prefix='%#',
        trim_blocks=True,
        autoescape=False,
        loader=FileSystemLoader(template_dir))

    template = latex_jinja_env.get_template(template_name)

    # Assemble dict
    plot_dir_files = sorted([design for design in os.listdir(plot_dir) if design.endswith('.png') and design.startswith('SeqeuenceRecovery')])
    figure_dict = {}

    for design in plot_dir_files:
        design_setname = design[:-4].split('-')[-1].upper()
        figure_dict[design_setname] = os.path.join(os.path.abspath(plot_dir), design)

    pprint(figure_dict)
    # Generate tex
    with open(f'SupplementalMaterials-SequenceRecovery.tex', 'w') as asdf:
        asdf.write(template.render(figure_dict=figure_dict))


if __name__ == '__main__':
    argv = sys.argv[1:]
    if len(argv) == 0:
        docopt.docopt(__doc__)

    if argv[0] == 'FragmentPlots':
        AllFragmentPlots(docopt.docopt(AllFragmentPlots.__doc__, argv=argv))
    if argv[0] == 'DesignPairPlots':
        DesignPairPlots(docopt.docopt(DesignPairPlots.__doc__, argv=argv))
    if argv[0] == 'DesignViolinPlots':
        DesignViolinPlots(docopt.docopt(DesignViolinPlots.__doc__, argv=argv))
    if argv[0] == 'DesignMetricsGrid':
        DesignMetricsGrid(docopt.docopt(DesignMetricsGrid.__doc__, argv=argv))
    if argv[0] == 'BBSequenceLogos':
        BBSequenceLogos(docopt.docopt(BBSequenceLogos.__doc__, argv=argv))
    if argv[0] == 'ProfileSimilarity':
        ProfileSimilarity(docopt.docopt(ProfileSimilarity.__doc__, argv=argv))
    if argv[0] == 'SequenceRecovery':
        SequenceRecovery(docopt.docopt(SequenceRecovery.__doc__, argv=argv))
