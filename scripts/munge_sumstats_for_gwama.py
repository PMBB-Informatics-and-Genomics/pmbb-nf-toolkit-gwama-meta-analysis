import pandas as pd
from scipy.stats import norm
import numpy as np
import argparse as ap

TEST_ROWS = None

def make_arg_parser():
    parser = ap.ArgumentParser(description=".")
    
    parser.add_argument('-p', '--pheno', required=True, help='Phenotype')
    
    # Add a non-optional list argument for cohorts
    parser.add_argument('-c', '--cohort', required=True, help='Cohort')
    
    # Add an argument for launch_directory
    parser.add_argument('-s', '--sumstats', required=True, help='Path to sumstats file')

    parser.add_argument('--maf', help='MAF filter during munging', type=float)

    # Add an argument for output_directory
    parser.add_argument('-o', '--outDir', default='./', help='Path to output directory. Default: current working directory')

    parser.add_argument('--colnames', required=True, help='File with column name mappings')

    parser.add_argument('--traitType')
    
    return parser


args = make_arg_parser().parse_args()
colnames_file = args.colnames
ss_file = args.sumstats
min_maf = args.maf
pheno = args.pheno
cohort = args.cohort
out_dir = args.outDir
trait_type = args.traitType

output_file = f'{cohort}.{pheno}.munged.txt'

colnames_rows = open(colnames_file).read().splitlines()
col_map = dict(zip([r.split('=')[1] for r in colnames_rows],
                   [r.split('=')[0] for r in colnames_rows]))


chunks = []
for chunk in pd.read_table(ss_file, nrows=TEST_ROWS, chunksize=1E5):
    chunk = chunk.rename(columns=col_map)
    if 'EAF' in chunk.columns:
        chunk = chunk[chunk['EAF'].between(min_maf, 1-min_maf)]
    chunks.append(chunk)


df = pd.concat(chunks)
print(df.columns)
fixed = False

if trait_type == 'quant' and ('BETA' not in df.columns or 'SE' not in df.columns):
    # For QT we need to make sure we have BETA and SE
    raise ValueError('Quantitative trait needs BETA and SE columns')
elif trait_type == 'bin' and 'OR' in df.columns:
    if 'SE' in df.columns:
        df['BETA'] = np.log(df['OR'])
        df['OR_95L'] = np.exp(df['BETA'] - 1.96*df['SE'])
        df['OR_95U'] = np.exp(df['BETA'] + 1.96*df['SE'])
        fixed = True
    elif 'OR_95L' not in df.columns or 'OR_95U' not in df.columns:
        # Handle BT with OR but no Confidence Interval
        raise ValueError('Binary trait needs Odds Ratio with (Confidence Interval or SE)')
    else:
        # We already have what we need
        fixed = True
elif trait_type == 'bin' and 'OR' not in df.columns:
    # For BT we need to make sure we have OR, OR_95L, and OR_95U
    if 'BETA' in df.columns and 'SE' in df.columns:
        df['OR'] = np.exp(df['BETA'])
        df['OR_95L'] = np.exp(df['BETA'] - 1.96*df['SE'])
        df['OR_95U'] = np.exp(df['BETA'] + 1.96*df['SE'])
        fixed = True
    else:
        raise ValueError('Binary trait needs Odds Ratio column plus Confidence Interval OR Beta and SE to convert')
else:
    fixed = True
print(df)

if 'N' not in df.columns:
    if 'N_CASE' in df.columns and 'N_CTRL' in df.columns:
        df['N'] = df['N_CASE'] + df['N_CTRL']

col_order = ['CHR', 'POS', 'MARKERNAME', 'EA', 'NEA', 'OR', 'OR_95L', 'OR_95U', 'BETA', 'SE', 'N', 'EAF']
col_order = [c for c in col_order if c in df.columns]
mandatory = ['CHR', 'POS', 'MARKERNAME', 'EA', 'NEA']
if trait_type == 'bin':
    mandatory.extend(['OR', 'OR_95L', 'OR_95U'])
elif trait_type == 'quant':
    mandatory.extend(['BETA', 'SE'])

all_there = np.all([c in col_order for c in mandatory])

print(df.columns)

if fixed and all_there:
    df[col_order].to_csv(f'{out_dir}/{output_file}', sep='\t', na_rep='NA', index=False)
elif not all_there:
    raise ValueError(f'After attempting to munge, missing columns: {[c for c in mandatory if c not in col_order]}')
else:
    raise ValueError('Unable to properly munge summary statistics')