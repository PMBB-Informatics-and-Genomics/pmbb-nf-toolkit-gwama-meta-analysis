from manhattan_plot import ManhattanPlot
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse as ap
import os


def make_arg_parser():
    parser = ap.ArgumentParser(description=".")

    # Add a non-optional argument for phenotype
    parser.add_argument('-p', '--phenotype', required=True, help='phenotype')
    # Add a non-optional argument for cohort
    parser.add_argument('-c', '--analysis', required=True,
                        help='which meta-analysis group of cohorts')
    # Add an argument for output_directory
    parser.add_argument('-o', '--outDir', default='./',
                        help='Path to output directory. Default: current working directory')
    # Add an argument for summary statistics file
    parser.add_argument('-s', '--sumstats', required=True,
                        help='Path to summary statistics file')
    parser.add_argument('-a', '--annot', required=False, default=None)

    return parser


# parse arguments
args = make_arg_parser().parse_args()
output_dir = args.outDir
output_dir += '/' if output_dir[-1] != '/' else ''
pheno, analysis = args.phenotype, args.analysis
annot_file = args.annot

# pheno_table = args.phenoTable # only used for trait cardinality for plot
sumstats_file = args.sumstats

output_manhattan = f'{output_dir}{analysis}.{pheno}.manhattan.png'
output_qq = f'{output_dir}{analysis}.{pheno}.qq.png'

# Example output file to be used lives in:
# /project/path/to/data/*.gz

# Instantiate manhattan plot object
plot_title = f'GWAMA Meta Manhattan for {analysis}: {pheno.replace("_", " ")}'
mp = ManhattanPlot(sumstats_file, title=plot_title)
mp.load_data()
# clean data, use GWAMA default column names
col_map = {'CHR': '#CHROM', 'variant_id': 'ID', 'p-value': 'P'}
mp.clean_data(col_map=col_map)

if annot_file is not None:
    annot_df = pd.read_csv(annot_file)
    annot_df['ID'] = annot_df['Gene']
    mp.add_annotations(annot_df, extra_cols=['RSID'])

mp.get_thinned_data()

annot_thresh = 1E-5 if np.any(mp.thinned['P'].min() < 1E-5) else np.nanquantile(mp.thinned['P'], 10 / len(mp.thinned))
print(annot_thresh)
print(mp.thinned['P'].describe())

mp.update_plotting_parameters(vertical=True, merge_genes=True, 
                              sig=annot_thresh if not np.any(mp.thinned['P'] < 5E-8) else 5E-8, 
                              sug=annot_thresh, annot_thresh=annot_thresh, ld_block=1E6)
mp.full_plot(save=output_manhattan, extra_cols={'RSID': 'RSID'} if annot_file is not None else {})

# close fig and save
plt.clf()
print(f"Saved Manhattan plot to: {output_manhattan}")

# mp.qq_plot
mp.qq_plot(save=output_qq)
print(f"Saved qq plot to: {output_qq}")
