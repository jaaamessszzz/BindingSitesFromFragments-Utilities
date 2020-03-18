import io
import os
import sys
import tempfile
import tarfile
import shutil
from pprint import pprint
from multiprocessing import Pool

from docopt import docopt
import pandas as pd


def comprot_tables(design_stdout):
    """
    Parse standard out to get rotamers found for each scaffold

    :param design_stdout: path to stdout
    """

    stdout_root, stdout_file = os.path.split(design_stdout)
    stdout_root, design_path = os.path.split(stdout_root)

    row_list = list()

    with open(design_stdout, 'r') as stdout:
        position_info = False
        row = {}

        for line in stdout:
            if line.startswith('Position'):
                line_split = line.split()
                row = {'Position': int(line_split[1]),
                       'Restype': line_split[3],
                       'RotamersGenerated': int(line_split[6]),
                       'ContactModes': int(line_split[9])}
                position_info = True
            if 'rotamers accepted' in line and not position_info:
                raise Exception('Position info found without accepted rotamers reported!')
            if 'rotamers accepted' in line and position_info:
                row['Accepted'] = int(line.split()[0])
                row_list.append(row)
                # Reset
                row = {}
                position_info = False

    df = pd.DataFrame(row_list)
    df.to_csv(f'{design_path}.csv')

    # Design = design_path
    # Number of unique design positions
    designable_positions = ','.join([str(a) for a in sorted(df['Position'].unique())])
    unique_positions = len(df['Position'].unique())
    # Total number of rotamers generated
    total_rotamers = df['RotamersGenerated'].sum()
    # Total number of rotamers accepted
    total_accepted = df['Accepted'].sum()
    # Total number of rotamers added to rotamersets
    total_applied = sum([value if value < 50 else 50 for index, value in df['Accepted'].iteritems()])

    df_agg = {'design': design_path,
              'designable_positions': designable_positions,
              'unique_positions': unique_positions,
              'total_rotamers': total_rotamers,
              'total_accepted': total_accepted,
              'total_applied': total_applied
              }

    return df_agg


def complementary_rotamer_tables(args):
    """
    Generate csvs for rotamers found and added to Packer

    Usage:
      python3 comprot_tables <compound_dir> [options]

    Arguments:
      <compound_dir>                    Path to directory containing BSFF project directories

    Options:
      -c=<cores>, --cores=<cores>       Number of cores to use (Default: 12)
    """
    compound_root = args['<compound_dir>']
    cores = args['--cores'] if args['--cores'] else 12

    stdout_to_parse = list()

    # Gather paths to all stdout (just look at jobid == 1)
    for root, dirs, files in os.walk(compound_root):

        # Add designs to dirs_to_parse
        if any([word in root.split('/') for word in ['Design', 'Designs']]) and 'RMSD_0.75' not in root.split('/') and \
                any([file.endswith('.1') and file.startswith('bsff-design.sh.') for file in files]):
            stdout = [os.path.join(root, file) for file in files if (file.endswith('.1') and file.startswith('bsff-design.sh.'))][0]
            stdout_to_parse.append(stdout)

    print(f'Creating complementary rotamer reports for the following:')
    print("\n".join(stdout_to_parse))

    # Assemble pool
    pool = Pool(processes=cores)
    data = pool.map(comprot_tables, stdout_to_parse)
    pool.close()
    pool.join()

    df = pd.DataFrame(data)
    df.to_csv('DesignSpecialRotAgg.csv')


if __name__ == '__main__':
    argv = sys.argv[1:]
    if argv[0] == 'comprot_tables':
        complementary_rotamer_tables(docopt(complementary_rotamer_tables.__doc__, argv=argv))