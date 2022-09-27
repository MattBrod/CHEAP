import os 
from os import path,listdir 
import numpy as np
import pandas as pd

import argparse

import pybedtools as pbt
import scipy.stats as stats

import statsmodels.stats.multitest as multi
from sklearn.preprocessing import MinMaxScaler

# FUNCTIONS

def element_tests(trck_name, hit_list_name, ab, list_track_names, list_track_results):
    """ This is the function that separates the hit lists further into different types these being lncRNAs and controls (oncogenes, tumour supressors and neutral).
    It iterates through the hit lists to extract information from the merged dataframe and perform the following statistical tests: Fisher Exact test for count 
    overlap, and MannWhitneyU test for percentage of overlap."""
    track_name = hit_list_name + '_' + trck_name
    track_name = track_name.replace('|','_')
    
    hit_by_type = ab.groupby('type')    
    hit_by_type_list = [hit_by_type.get_group(x) for x in hit_by_type.groups]

    for hit_df in hit_by_type_list:

        type_str = hit_df['type'].unique().tolist()
        track_type = type_str[0] + '_' + track_name

        # COUNT OVERLAP

        print("Analysing %s" % track_type)
        
        hit_bed = hit_df[hit_df[hit_list_name]>=1].copy()
        not_hit_bed = hit_df[hit_df[hit_list_name]<1].copy()

        conditions_hit = [
        hit_bed['bp'] > 0,
        hit_bed['bp'] == 0,
        ]

        conditions_no_hit = [
        not_hit_bed['bp'] > 0,
        not_hit_bed['bp'] == 0,
        ]

        choices_hit = ['HO','HNO']
        choices_no_hit = ['NHO','NHNO']

        hit_bed['contig_pos'] = np.select(conditions_hit, choices_hit)
        not_hit_bed['contig_pos'] = np.select(conditions_no_hit, choices_no_hit)

        HO = hit_bed[hit_bed['contig_pos'] == 'HO']['contig_pos'].count()
        NHO = not_hit_bed[not_hit_bed['contig_pos'] == 'NHO']['contig_pos'].count()
        HNO = hit_bed[hit_bed['contig_pos'] == 'HNO']['contig_pos'].count()
        NHNO = not_hit_bed[not_hit_bed['contig_pos'] == 'NHNO']['contig_pos'].count()

        values = [HO, NHO, HNO, NHNO]


        if 0 in values:
            HO = HO + 0.5
            NHO = NHO + 0.5
            HNO = HNO + 0.5
            NHNO = NHNO + 0.5

        contingency_table = [[HO, NHO],[HNO, NHNO]]

        oddsr, count_p = stats.fisher_exact(contingency_table)

        # PERCENTAGE OVERLAP

        hit_bed = hit_df[hit_df[hit_list_name]>=1]
        not_hit_bed = hit_df[hit_df[hit_list_name]<1]

        rbc_flag = False
        try:
            U1, percentage_p = stats.mannwhitneyu(hit_bed['percentage'], not_hit_bed['percentage'])
        except ValueError:
            U1, percentage_p = 0, 1
            rbc_flag = True

        n1, n2 = len(hit_bed['percentage']), len(not_hit_bed['percentage'])

        U2 = n1*n2 - U1

        if rbc_flag == True:
            rbc = 0
        else:
            rbc = 1 - ((2*U2)/(n1*n2))


        track_results = [count_p, oddsr, HO, NHO, HNO, NHNO,  percentage_p, rbc]

        list_track_names.append(track_type)
        list_track_results.append(track_results)

        rbc_flag = False

def output(list_track_r, list_track_n):
    """Output a file containing all the tests results. Columns 2 to 5 are the values used in the contingency table to perform Fisher Exact test. g1 indicates
    that the column is for hit regions and g2 indicates that the column contains the values for non-hit regions. The 'yes' and 'no' found in the column name
    indicate if the region overlaps in the case of 'yes', or not in the case of 'no' with the tested genomic element."""

    final_results = pd.DataFrame(list_track_r, columns=['count_pvalue','OddsRatio', 'g1.yes', 'g1.no', 'g2.yes','g2.no','percentage_pvalue','Rank_biserial_cor'])
    final_results['tracks'] = list_track_n
    final_results = final_results[['tracks','count_pvalue','OddsRatio', 'g1.yes', 'g1.no', 'g2.yes','g2.no', 'percentage_pvalue','Rank_biserial_cor']]
    finalfinal.append(final_results)

def process(track_list, hit_df_list, track_list_name):
    """The function takes the list of merged dataframes containing what genomic elements overlap with the tested exons, the list containing what exons were
    considered hits in each cell line, time point and operation, and the list with the names of the genomic elements so the rows in the resulting dataframe 
    are aproppiately marked. The function goes through all the merged dataframes that originate from the genomic elements BED files found in the same folder.
    The function also iterates through the hit lists."""

    track_name_counter = 0
    for track_df in track_list:
        for df in hit_df_list:

            region_hit = track_df.merge(df) # Merge the genomic elements dataframes with the overlap information with the dataframe containing which exons are hits
            
            hit_list_name = region_hit.columns[-1]

            list_track_names = []
            list_track_results = []

            track_name = track_list_name[track_name_counter]

            # Function that makes the statistical tests to determine if an overlap between genomic regions and lncRNA exons are significant (p-values and size effect measures)
            element_tests(track_name, hit_list_name, region_hit, list_track_names, list_track_results)

            # Output the results as files.
            output(list_track_results, list_track_names)
            
            list_track_names = []
            list_track_results = []

        track_name_counter += 1

# ARGUMENTS

# Region analysis
# python3 hit_pipeline.py -hit OUTPUT_GENES_IBARC_RAW1.csv -ex regions.bed -ct CONTROLS.csv -trc folder

# Region:exon analysis
# python3 hit_pipeline.py -hit OUTPUT_GENES_IBARC_RAW1.csv -ex intersections.bed -ct CONTROLS.csv -trc folder -int REASSIGNMENTS.csv

parser = argparse.ArgumentParser(description="DESCRIPTION PLACEHOLDER")

Mandatory = parser.add_argument_group('Mandatory arguments')
Optional = parser.add_argument_group('Optional arguments')

# This argument requires the raw .CSV output from MAGeCK-IBAR (OUTOUTPUT_GENES_IBARC_RAW.csv).
Mandatory.add_argument('-hit', '--hitlist',
                       dest="hitlist",
                       action="store",
                       required=True,
                       help="MANDATORY ARGUMENT: CSV with the raw file from MAGeCK IBAR")

# This argument requires the BED file containing the genomic coordinates of the tested lncRNA exons (regions.bed). For the by region:exon analysis, the file has
# to be the BED file with the coordinates of the exons and the regions they intersect with (intersections.bed).
Mandatory.add_argument('-ex', '--exons',
                       dest="exons",
                       action="store",
                       required=True,
                       help="MANDATORY ARGUMENT: BED file with genomic coordinates of regions/region:exons")

# This argument requires the CSV file containing the region IDs of the controls used in the CRISPR screening (CONTROLS.csv).
Mandatory.add_argument('-ct', '--controls',
                       dest="controls",
                       action="store",
                       required=True,
                       help="MANDATORY ARGUMENT: CSV file containing which regions in the CRISPR screen are controls and what type of controls they are (oncogenes,tumour supressors and neutral)")

# This argument requires the directory path that contains the BED files of the genomic elements that will be tested for enrichment. In the case of wanting to
# test genomic elements into separate families or groups, they have to be separated into different folders in the directort path provided.
Mandatory.add_argument('-trc', '--tracks',
                       dest="tracks",
                       action="store",
                       required=True,
                       help="MANDATORY ARGUMENT: directory of BED files of genomic elements of interest")

# To perform the by region:exon analysis, this argument needs to be provided with the CSV file containing the specifics of the intersections between regions and
# exons and what type of operation (DC, DF, DP, DI) is being done (REASSIGNMENTS.csv)
Optional.add_argument('-int', '--intersections',
                       dest="intersections",
                       action="store",
                       required=False,
                       help="OPTIONAL ARGUMENT: CSV file containing the intersection of the regions and the experimental exons. It also contains the type of operation of the intersection (DF, DI, DP, DC)")

arguments = parser.parse_args()

# ERROR LOGS

# HITLIST ARGUMENT ERROR
if not arguments.hitlist:
    raise NameError(
        "ERROR! You need to add a hit list file.")
else:
    if (os.path.isfile(arguments.hitlist)):
        if arguments.hitlist.endswith((".csv")):
            arguments.hitlist = os.path.abspath(arguments.hitlist)
        else:
            raise NameError("ERROR! Extension must be .csv")
    else:
        raise NameError("ERROR! Hit list file not found.")

# TRACKS ARGUMENT ERROR
if not arguments.tracks:
    raise NameError(
        "ERROR! You need to add a tracks directory.")
else:
    if (os.path.isdir(arguments.tracks)):

        path_subfolder_list = [f.path for f in os.scandir(arguments.tracks) if f.is_dir()]
    else:
        raise NameError("ERROR! Tracks directory not found.")

# FILE PREPROCESSING

pd.options.mode.chained_assignment = None

finalfinal = [] #list where results will be stored

# Exon file argument
region_bed = pd.read_csv(arguments.exons, sep='\t', header=None)
region_bed.columns = ["chr","start","end","region_exon"]

# Controls file argument
controls = pd.read_csv(arguments.controls)
controls = controls[['region_id','type']]
controls = controls.drop_duplicates().reset_index(drop=True) # File might have duplicated rows

# Merging together REASSIGNMENTS, exons and controls (region:exon analysis)
if arguments.intersections is not None:
    reassignments = pd.read_csv(arguments.intersections)
    reassignments_short = reassignments[['op_type_v2','exon_id','region_id']].copy()
    reassignments_short['region_exon'] = reassignments_short['region_id'] + ';' + reassignments_short['exon_id'] #column to match with intersections file
    reassignments_short = reassignments_short.drop_duplicates().reset_index(drop=True)
    reassignments_coords = region_bed.merge(reassignments_short)
    reassignments_coords['region_exon_op'] = reassignments_coords['region_exon'] + ';' + reassignments_coords['op_type_v2']
    reassignments_coords = reassignments_coords.merge(controls, how='left', on='region_id').drop(columns=['exon_id','region_id'])
    reassignments_coords = reassignments_coords.drop(columns=['region_exon','op_type_v2'])
    reassignments_coords['type'] = reassignments_coords['type'].fillna('lncrna')
    reassignments_coords['type'] = reassignments_coords['type'].replace(['+R','-','+A'], ['sup','neu','onc'])
else: # Merge exons and controls (region analysis)
    region_bed.columns = ["chr","start","end","region_id"]
    region_coords = region_bed.merge(controls, how='left', on='region_id')
    region_coords['type'] = region_coords['type'].fillna('lncrna')
    region_coords['type'] = region_coords['type'].replace(['+R','-','+A'], ['sup','neu','onc'])

# Hit list argument
mageck_output = pd.read_csv(arguments.hitlist)

fdr_thr = 0.25 # FDR threshold to determine what counts as a hit

mageck_output['is_hit'] = mageck_output['neg|fdr'].apply(lambda x:int(x<=fdr_thr))

# Divide hit lists by hit confidence: low confidence (hit in 1 timepoint out of 3) and high confidence (hit in 2 timepoints out of 3)
D = mageck_output.groupby(by=['region_id','cell_line'])['is_hit'].sum().reset_index()
D['lowconf'] = D['is_hit'].apply(lambda x:int(x>0))
D['highconf'] = D['is_hit'].apply(lambda x:int(x>1))

P = D.pivot('region_id','cell_line',['is_hit','lowconf','highconf']).reset_index()
P.columns = ['%s|%s'%x for x in P.columns]
P.rename(columns={'region_id|':'region_id'},inplace=True)

hit_list_df_list = [] #list for storing hit lists divided by cell line, time point and operation (only in region:exon analysis).

# Divide the hit lists by cell line and operation (region:exon analysis).
cell_lines = mageck_output['cell_line'].unique()
num_desired_col = len(P.columns) - len(cell_lines) -1
counter = 0
while counter < num_desired_col:
    second_col = counter + len(cell_lines) + 1
    cell_line_df = P.iloc[:,[0,second_col]]

    if arguments.intersections is not None: # Region:exon analysis
        merged_reassignments = reassignments_short.merge(cell_line_df)

        operations = merged_reassignments['op_type_v2'].unique()

        for ops in operations: # Divide hit lists into the 4 distinct operations (DF, DC, DI, DP)
            op_df = merged_reassignments['op_type_v2'] == ops
            op_df = merged_reassignments[op_df]
            op_df = op_df.add_suffix('_' + ops)
            op_df_short = op_df.iloc[:,-2:]
            op_df_short.rename(columns={ op_df_short.columns[0]: "name" }, inplace = True)
            op_df_short["name"]  = op_df_short["name"] + ';' + ops
            op_df_short = op_df_short.reset_index(drop=True)
            hit_list_df_list.append(op_df_short)
        counter += 1
    else: # Region analysis
        cell_line_df.rename(columns={ cell_line_df.columns[0]: "name" }, inplace = True)
        hit_list_df_list.append(cell_line_df)
        counter += 1

# Divide the hit lists by time point and operation (region:exon analysis).
timepoints = mageck_output['T'].unique()

for number in timepoints:
    number_df = mageck_output.loc[mageck_output['T'] == number].copy()
    number_df['is_hit'] = number_df['neg|fdr'].apply(lambda x:int(x<=fdr_thr))
    number_short_df = number_df[['region_id','cell_line','is_hit']]
    number_pivot = number_short_df.pivot('region_id','cell_line',['is_hit']).reset_index()
    number_pivot.columns = ['%s|%s'%x for x in number_pivot.columns]
    number_pivot.rename(columns={'region_id|':'region_id'},inplace=True)
    number_pivot.columns = number_pivot.columns.str.replace("is_hit", str(number))
    counter = 0
    while counter < (len(number_pivot.columns)-1):
        col = counter + 1
        timepoint_df = number_pivot.iloc[:,[0,col]]

        if arguments.intersections is not None: # Region:exon analysis
            merged_reassignments = reassignments_short.merge(timepoint_df)

            operations = merged_reassignments['op_type_v2'].unique()

            for ops in operations:
                op_df = merged_reassignments['op_type_v2'] == ops
                op_df = merged_reassignments[op_df]
                op_df = op_df.add_suffix('_' + ops)
                op_df_short = op_df.iloc[:,-2:]
                op_df_short.rename(columns={ op_df_short.columns[0]: "name" }, inplace = True)
                op_df_short["name"]  = op_df_short["name"] + ';' + ops
                op_df_short = op_df_short.reset_index(drop=True)
                hit_list_df_list.append(op_df_short)
            counter += 1
        else: # Region analysis
            timepoint_df.rename(columns={ timepoint_df.columns[0]: "name" }, inplace = True)
            hit_list_df_list.append(timepoint_df)
            counter += 1

# Iterate through folders from track directory path to get BED files
for folder in path_subfolder_list:
    files = sorted(list(filter(lambda x: x.endswith(".bed"), os.listdir(folder))))

    folder_name = path.split(folder)[1]

    track_list = [] # List that will contain the track_intersect_simple dataframes from one folder
    track_list_name = [] # List of strings containing the name of the BED file
    for track in files: # Iterate through each BED file in the folder

        track_specific = path.join(folder, track)
        track_BedTool = pd.read_csv(track_specific, sep='\t', header=None)
        track_BedTool = pbt.BedTool.from_dataframe(track_BedTool)
        track_BedTool = track_BedTool.sort()

        # This section intersects the coordinates dataframe (exons+controls or exons+controls+reassignments) with the BED file coordinates.
        # This way we can obtain what regions of the genomic elements overlap with the tested exons, obtaining information such as number of base pairs of
        # overlap and calculating from that the percentage of overlap.
        if arguments.intersections is not None: # Region:exon analysis

            reassignments_bedtool = pbt.BedTool.from_dataframe(reassignments_coords)
            reassignments_bedtool = reassignments_bedtool.sort()

            merged_intersect = reassignments_bedtool.intersect(track_BedTool, sorted = True).sort().merge() # Merge intersections that overlap so they are not counted separately
            track_intersect = reassignments_bedtool.intersect(merged_intersect, wao=True, sorted = True) # 'wao' argument is used to obtain the number of base pairs of overlap between 2 features

            track_intersect = track_intersect.to_dataframe()

            track_intersect[['region_id', 'exon_id', 'operation']] = track_intersect['name'].str.split(';', expand=True)
            track_intersect.columns = ['chrom','start','end','name','type','chrom1','start1','end1','bp','region_id','exon_id','operation']

            track_intersect['bp_region'] = track_intersect['end'] - track_intersect['start']

            track_intersect_simple = track_intersect.groupby(['name']).agg({'chrom': 'first', 'start': 'first', 'end': 'first', 'type': 'first', 'chrom1': 'first', 'start1': 'first', 'end1': 'first', 'bp': 'sum', 'region_id': 'first', 'exon_id': 'first', 'operation': 'first', 'bp_region': 'first'}).reset_index()


            track_intersect_simple = track_intersect_simple.drop_duplicates(subset="name", keep='first', ignore_index=True)

            track_intersect_simple['percentage'] = track_intersect_simple['bp'] / track_intersect_simple['bp_region'] * 100
            

            track_name = track.split('.',)[0]
            
            track_list_name.append(track_name)
            track_list.append(track_intersect_simple)
            pbt.helpers.cleanup() # Due to a bug, if the temp files of pyBEDTools pile up and the library can't read properly the pbt objects
        else: # Region analysis

            reassignments_bedtool = pbt.BedTool.from_dataframe(region_coords)
            reassignments_bedtool = reassignments_bedtool.sort()

            merged_intersect = reassignments_bedtool.intersect(track_BedTool, sorted = True).sort().merge()
            track_intersect = reassignments_bedtool.intersect(merged_intersect, wao=True, sorted = True) # 'wao' argument is used to obtain the number of base pairs of overlap between 2 features

            track_intersect = track_intersect.to_dataframe()

            track_intersect = track_intersect.drop(columns=['strand','thickStart','thickEnd'])

            track_intersect.columns = ['chrom','start','end','name','type','bp']

            track_intersect['bp_region'] = track_intersect['end'] - track_intersect['start']

            track_intersect_simple = track_intersect.groupby(['name']).agg({'chrom': 'first', 'start': 'first', 'end': 'first', 'type': 'first', 'bp': 'sum','bp_region': 'first'}).reset_index()

            track_intersect_simple = track_intersect_simple.drop_duplicates(subset="name", keep='first', ignore_index=True)

            track_intersect_simple['percentage'] = track_intersect_simple['bp'] / track_intersect_simple['bp_region'] * 100

            track_name = track.split('.',)[0]
            
            track_list_name.append(track_name)
            track_list.append(track_intersect_simple)
            pbt.helpers.cleanup()


    # SCRIPT START

    cwd = os.getcwd()

    results_path = os.path.join(cwd, "results")
    if not os.path.exists(results_path):
        os.makedirs(results_path)

    subfolder_results_path = os.path.join(results_path, folder_name)
    if not os.path.exists(subfolder_results_path):
        os.makedirs(subfolder_results_path)

    # The following function will perform the tests with the processed files to obtain the statistic results, all except the scaling of Odds Ratio and Rank
    # Biserial Correlation, along with the FDR correction for the count overlap p-values and percentage of overlap p-values.

    process(track_list, hit_list_df_list, track_list_name)

    merged = pd.concat(finalfinal) # Concatenating together all of the dataframe rows with the results.

    # Scale (log transformation for Odds Ratio and MinMaxScaler for Rank-Biserial Correlation) and FDR correct the results, taking lncRNAs, oncogenes, 
    # tumour supressors and neutral controls separately.

    type_list = (merged['tracks'].str.split('_', expand=True)[0]).unique()

    type_result_dfs = []

    for type in type_list:
        merged_type = merged[merged['tracks'].str.contains(type)].copy()
        
        merged_type["OddsRatio"] = merged_type["OddsRatio"].replace([np.inf, 0], np.nan)

        merged_type = merged_type.dropna()

        if merged_type.empty:
            type_result_dfs.append(merged_type)
        else:
            merged_type["logOR"] = np.log(merged_type["OddsRatio"])
            logOR2D = merged_type['logOR'].values.reshape(-1, 1)
            scaler = MinMaxScaler(feature_range=(-1, 1))
            logOR_scaled = scaler.fit_transform(logOR2D)
            merged_type['scaled_logOR'] = logOR_scaled

            merged_type['c_pvalue_adj'] = multi.fdrcorrection(merged_type['count_pvalue'],alpha=0.05)[1]
            merged_type['p_pvalue_adj'] = multi.fdrcorrection(merged_type['percentage_pvalue'],alpha=0.05)[1]

            type_result_dfs.append(merged_type)

    merged_final = pd.concat(type_result_dfs)

    # FORMATING RESULTS

    cols_float = ['OddsRatio','Rank_biserial_cor', 'logOR', 'scaled_logOR']
    cols_exp = ['count_pvalue','percentage_pvalue', 'c_pvalue_adj', 'p_pvalue_adj']

    for cols in cols_float:
        merged_final[cols] = merged_final[cols].map('{:,.2f}'.format)

    for cols in cols_exp:
        merged_final[cols] = merged_final[cols].map('{:,.2e}'.format)

    final_results_destination = path.join(subfolder_results_path, 'enrichment_results' +'.txt')
    merged_final.to_csv(final_results_destination, sep='\t', index= False)
    finalfinal = []