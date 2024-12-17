import os
import sys
import pandas as pd
import subprocess as sub

# ========> EDIT THE TOP OF THIS SCRIPT

# Decide if you're providing a file name pattern or a table
file_source = 'pattern'
# file_source = 'table'

# ====> Use this section if providing a directory with all files that follow a naming pattern
# Set this variable to the full path containing all of your summary stats
# This should end in a /
sumstats_dir = '/path/to/data/'
# Set this variable to a delimited pattern
# Any variables contributing to the cohort name should be COHORT (dataset, ancestry, sex, batch, etc)
# The phenotype variable should be PHENO
# At the end, indicate the suffix you want to keep
pattern = 'COHORT.IGNORE.PHENO.IGNORE.COHORT.IGNORE.IGNORE.IGNORE.IGNORE.SUFFIX'
delim = '.' # The delimiter for your file name pattern

# ====> Use this section if providing a table
# If your file names don't really follow a pattern you can provide a table
# The table should have at a minimum COHORT_DIR, PHENO, and FILE_PATH columns (in any order)
file_table = '/path/to/table/with/file_info.csv'
table_delim = ','
# This is the suffix you want to use when sym linking
# The suffix should start with a .
new_suffix = '.txt.gz'

# ========> CAUTION DO NOT EDIT PAST THIS LINE


if file_source == 'pattern':
    # Error handling
    if 'COHORT' not in pattern:
        sys.exit('The pattern you provided does not contain any cohort variables')
    if 'PHENO' not in pattern:
        sys.exit('The pattern you provided does not contain any cohort variables')
    
    # Detect the number of variables in the file name
    num_splits = len(pattern.split(delim)) - 1
    # List the files in the directory
    files = pd.Series(os.listdir(sumstats_dir))
    # Split each file name to make a DataFrame
    files_df = files.str.split(delim, n=num_splits, expand=True)
    files_df.columns = pattern.split(delim) # Add columns from the pattern
    files_df['FILE_NAME'] = files # Add file name column
    files_df['FILE_PATH'] = sumstats_dir + files_df['FILE_NAME'] # Add full file path column
    files_df = files_df.dropna(how='any') # Drop any files not matching the pattern (missing values)

    # If there is more than one COHORT value in the pattern,
    # Concatenate the variables in order with underscores
    multiple_cohort_cols = pattern.count('COHORT') > 1
    if multiple_cohort_cols:
        files_df['COHORT_DIR'] = files_df['COHORT'].apply(lambda x: '_'.join(x), axis=1)
    else:
        files_df['COHORT_DIR'] = files_df['COHORT']
elif file_source == 'table':
    files_df = pd.read_table(file_table, sep=table_delim)
    files_df['SUFFIX'] = new_suffix[1:]

    # Error handling
    if 'COHORT_DIR' not in files_df.columns:
        sys.exit('The table you provided does not have a COHORT_DIR variable')
    if 'PHENO' not in files_df.columns:
        sys.exit('The table you provided does not have a PHENO variable')
    if 'FILE_PATH' not in files_df.columns:
        sys.exit('The table you provided does not have a FILE_PATH variable')

print(files_df)

# Construct the new path
files_df['LINK_PATH'] = files_df['COHORT_DIR'] + '/Sumstats/' + files_df['PHENO'] + '.' + files_df['SUFFIX']

# Iterate over the DataFrame and run the symbolic link command
for _, row in files_df.iterrows():
    os.makedirs(f'{row["COHORT_DIR"]}/Sumstats/', exist_ok=True)
    symlink_cmd = f'ln -s {row["FILE_PATH"]} {row["LINK_PATH"]}'
    sub.run(symlink_cmd, shell=True, check=True)