#! /usr/bin/python3
"""
Quantify the number of design positions where profile similarity improves between design conditions
(i.e. special_rot weights).

Usage:
  QuantifyPSImprovement.py <stats_dir> <condition> [options]

Arguments:
  <stats_dir>             Directory containing design statistics
  <condition>             special_rot weight for 'improved' condition

Options:
    --help                   Halp

"""

import os
import subprocess
import tempfile

import pandas as pd
import docopt
import prody

def quantify_improvement(original_df, improved_df):
    """
    derp.
    :return:
    """
    joined_df = pd.merge(original_df, improved_df, on=['complex', 'position'], how='inner', suffixes=('_orig', '_impr'))
    print(joined_df)
    design_positions = 0
    improved_positions = 0
    for index, row in joined_df.iterrows():
        design_positions += 1
        if row['similarity_impr'] - row['similarity_orig'] > 0.1:
            improved_positions += 1

    print('design_positions:', design_positions)
    print('improved_positions:', improved_positions)


if __name__ == '__main__':
    args = docopt.docopt(__doc__)
    stats_dir = args['<stats_dir>']
    improved_arg = float(args['<condition>'])

    # ASSUMES WE AREN'T LOOKING AT CONTROL!!!!
    improved_str = str(-improved_arg) if improved_arg > 0 else str(improved_arg)

    df = pd.concat([pd.read_csv(os.path.join(stats_dir, stats)) for stats in os.listdir(stats_dir) if stats.endswith('.csv')])

    weight_groups = df.groupby(['weight'])
    original_df = weight_groups.get_group('Control')
    improved_df = weight_groups.get_group(improved_str)

    print(original_df)
    print(improved_df)
    quantify_improvement(original_df, improved_df)