import os
from os import path,listdir 
import pandas as pd
import argparse

# ARGUMENTS

parser = argparse.ArgumentParser(description="PLACEHOLDER")

Mandatory = parser.add_argument_group('Mandatory arguments')

Mandatory.add_argument('-dir', '--directory',
                        dest='folder',
                        action='store',
                        required=True,
                        help="Mandatory Argument:")

Mandatory.add_argument('-typ', '--type',
                        dest='type',
                        action='store',
                        required=True,
                        help="Mandatory Argument:")

Mandatory.add_argument('-cell', '--cell line',
                        dest='cell_line',
                        action='store',
                        required=True,
                        help="Mandatory Argument:")

arguments = parser.parse_args()

# ARGUMENT ERROR LOG

if not arguments.folder:
    raise NameError(
        "ERROR! You need to add a directory.")
else:
    if (os.path.isdir(arguments.folder)):

        path_subfolder_list = [f.path for f in os.scandir(arguments.folder) if f.is_dir()]
    else:
        raise NameError("ERROR! Directory not found.")

# READ THE FILE FOUND IN EACH FOLDER IN THE RESULTS DIRECTORY

file_df = []
cwd = os.getcwd()

for subfolder in path_subfolder_list:

    files = sorted(list(filter(lambda x: x.endswith(".txt"), os.listdir(subfolder))))

    folder_name = path.split(subfolder)[1]

    # results_path = os.path.join(cwd, subfolder)

    for output in files: # Iterate through each file in the folder
        output = os.path.join(subfolder, output)
        results = pd.read_csv(output, sep='\t')
        file_df.append(results)

merged_list = pd.concat(file_df)

# Select the results by type, type being lncrna/onc/sup/neu and their subsets highconf/lowconf/1/2/3. E.g. lncrna_highconf or onc_1.

type_list = merged_list[merged_list['tracks'].str.contains(arguments.type)]

# Consider only genomic elements with significant p-values in both FDR-adjusted count overlap and FDR-adjusted percentage overlap

type_list = type_list[type_list['c_pvalue_adj'] <= 0.05]
type_list = type_list[type_list['p_pvalue_adj'] <= 0.05]

# Select the results by cell line (U87, hela, hek, MCF7) and operation if applicable (DP,DI,DC,DF). E.g. hela_DF.

cell_line_df = type_list[type_list['tracks'].str.contains(arguments.cell_line)]

cell_line_df = cell_line_df.sort_values(by=['scaled_logOR', 'Rank_biserial_cor'], ascending=True)

try:
    cell_line_df['tracks'] = cell_line_df['tracks'].str.split('_', expand=True)[4]
except KeyError:
    cell_line_df['tracks'] = cell_line_df['tracks'].str.split('_', expand=True)[3]

cell_line_df.to_csv(arguments.cell_line +'.txt', sep = '\t', index=None)