from collections import Counter
import itertools
import sys
from collections import OrderedDict
import bioframe
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib import cm
from matplotlib.gridspec import GridSpec
from matplotlib.gridspec import GridSpecFromSubplotSpec
import matplotlib.colors as colors
from matplotlib.colors import ListedColormap
import random
import seaborn as sns
import pickle

iteration = sys.argv[1]
alignmentFile = sys.argv[2]
outDataDir = sys.argv[3]

f = open(alignmentFile, 'rb')
overlap_dfs = pickle.load(f)
f.close()

#### This section (below) will change depending on specific samples ####
conditions = [
    't0Mit_R1_T1',
    't2_R1',
    't4DMSO_R1',
    't4ICRF_R1',
    't8DMSO_R1',
    't8ICRF_R1',
    't0Mit_R1_T2',
    't0Mit_R2',
    't2_R2',
    't4DMSO_R2',
    't4ICRF_R2',
    't8DMSO_R2',
    't8ICRF_R2',
    't0Mit_R3',
    't2_R3',
    't4DMSO_R3',
    't4ICRF_R3',
    't8DMSO_R3',
    't8ICRF_R3',   
]

long_names = {
    't0Mit_R1_T1' : 'TI-MC3C-Dpn-t0Mit-4-30',
    't2_R1' : 'TI-MC3C-Dpn-t2-4-30',
    't4DMSO_R1' : 'TI-MC3C-Dpn-t4DMSO-4-30',
    't4ICRF_R1' : 'TI-MC3C-Dpn-t4ICRF-4-30',
    't8DMSO_R1' : 'TI-MC3C-Dpn-t8DMSO-4-30',
    't8ICRF_R1' : 'TI-MC3C-Dpn-t8ICRF-4-30',
    't0Mit_R1_T2' : 'TI-MC3C-Dpn-t0Mit-4-30-T2',
    't0Mit_R2' : 'TI-MC3C-Dpn-t0Mit-4-39',
    't2_R2' : 'TI-MC3C-Dpn-t2-4-39',
    't4DMSO_R2' : 'TI-MC3C-Dpn-t4DMSO-4-39',
    't4ICRF_R2' : 'TI-MC3C-Dpn-t4ICRF-4-39',
    't8DMSO_R2' : 'TI-MC3C-Dpn-t8DMSO-4-39',
    't8ICRF_R2' : 'TI-MC3C-Dpn-t8ICRF-4-39',
    't0Mit_R3' : 'TI-MC3C-Dpn-t0Mit-R3-5-14',
    't2_R3' : 'TI-MC3C-Dpn-t2-R3-5-14',
    't4DMSO_R3' : 'TI-MC3C-Dpn-t4DMSO-R3-5-14',
    't4ICRF_R3' : 'TI-MC3C-Dpn-t4ICRF-R3-5-14',
    't8DMSO_R3' : 'TI-MC3C-Dpn-t8DMSO-R3-5-14',
    't8ICRF_R3' : 'TI-MC3C-Dpn-t8ICRF-R3-5-14',
}

##################################################################

#Shuffle steps within each walk for permutations
#Just shuffle all rows, then resort by query_name - doesn't need to be grouped for the shuffle

permuted_dfs = {}
for cond in conditions:
    df = overlap_dfs[cond]['length_6'].copy()
    df_shuff = df.iloc[np.random.permutation(df.index)].reset_index(drop=True).sort_values(['Query_Name']).reset_index(drop=True)
    df_shuff['Iteration'] = iteration
    permuted_dfs[cond] = df_shuff

    
#Add the annotations needed for splitting walk by first_x read length
annotated_dfs = {}
for cond in conditions:
    df = permuted_dfs[cond].copy() #save more of the columns here
    grouped_df = df.groupby(by = 'Query_Name')
    df['Fragment_Index'] = grouped_df.cumcount()
    
    summary_table = pd.DataFrame()
    summary_table['Fragment_Number'] = grouped_df.size()
    
    df2 = df.merge(summary_table['Fragment_Number'], left_on = 'Query_Name', right_index = True)
    annotated_dfs[cond] = df2  
    
#Annotate walks with largest step size - to use for entanglement analysis - first X walks
stepsize_dfs = {}

for cond in conditions:
    stepsize_dfs[cond] = {}
    df = annotated_dfs[cond].copy()
    df['dist'] = df.mid.diff()
    df['dist'].iloc[np.where(df['Query_Name'] != df['Query_Name'].shift())] = np.nan
    df['dist'].iloc[np.where(df['chrom'] != df['chrom'].shift())] = np.nan
        
    #Adding in whether a step changes chromosomes
    df['Trans_Step'] = df['chrom'] != df['chrom'].shift()
    df['Trans_Step'].iloc[np.where(df['Query_Name'] != df['Query_Name'].shift())] = False 
        
    #adding absolute distance of step as well
    df['Abs_Dist'] = abs(df['dist'].copy())
        
    grouped_df = df.groupby(by = 'Query_Name')
    summary_table1 = pd.DataFrame()

    #size of largest step in walk
    summary_table1['Largest_Step'] = grouped_df['Abs_Dist'].max()
    df = df.join(summary_table1, on = 'Query_Name')
        
    #whether a fragment is part of largest step
    df['Largest_Step_Fragment_Start'] = (df['Abs_Dist'] == df['Largest_Step']).shift(-1)
    df['Largest_Step_Fragment_End'] = (df['Abs_Dist'] == df['Largest_Step'])
        
    #Filter out reads where more than one step is the same size as the largest step - there are a few of these
    summary_table1 = pd.DataFrame()
        
    grouped_df = df.groupby(by = 'Query_Name')
    summary_table1['Num_Largest_Steps'] = grouped_df['Largest_Step_Fragment_Start'].sum().astype(int)
        
    df = df.join(summary_table1, on = 'Query_Name')
    df = df[df['Num_Largest_Steps'] == 1]
        
    #midpoint of start of largest step fragment
    midpoint_start = df.loc[df['Largest_Step_Fragment_Start'] == True, ['Query_Name', 'mid']]
    midpoint_start = midpoint_start.set_index('Query_Name', drop = True)

    #midpoint of end of largest step fragment
    midpoint_end = df.loc[df['Largest_Step_Fragment_End'] == True, ['Query_Name', 'mid']]
    midpoint_end = midpoint_end.set_index('Query_Name', drop = True)

    df = df.join(midpoint_start, on = 'Query_Name', rsuffix = '_Largest_Step_Start')
    df = df.join(midpoint_end, on = 'Query_Name', rsuffix = '_Largest_Step_End')

    #absolute distance from each fragment to start of largest step
    df['Distance_To_Largest_Step_Start'] = abs(df['mid_Largest_Step_Start'].copy() - df['mid'].copy())

    #distance from each fragment to end of largest step
    df['Distance_To_Largest_Step_End'] = abs(df['mid_Largest_Step_End'].copy() - df['mid'].copy())  

    stepsize_dfs[cond] = df
        
#next add annotations to each fragment and summarize walks
overlap_dfs2 = {}
for cond in conditions:
    df = stepsize_dfs[cond].copy()

    #Adding in whether a step changes compartment type, or compartment index
    df['Inter_Comp_Type_Step'] = df['Frag_Comp_Type'] != df['Frag_Comp_Type'].shift()
    df['Inter_Comp_Type_Step'].iloc[np.where(df['Query_Name'] != df['Query_Name'].shift())] = np.nan

    df['Inter_Comp_Index_Step'] = df['Frag_Comp_Index'] != df['Frag_Comp_Index'].shift()
    df['Inter_Comp_Index_Step'].iloc[np.where(df['Query_Name'] != df['Query_Name'].shift())] = np.nan
    overlap_dfs2[cond] = df
        
#next add annotations to each fragment based on how far fragment is from either side of largest step
#next add annotations to each fragment based on how far fragment is from either side of largest step
overlap_dfs3 = {}
for cond in conditions:
    df = overlap_dfs2[cond].copy()

    #is the step within 1/4 of largest step size from start of the largest step?
    df['Near_Largest_Step_Start_Step'] = df['Distance_To_Largest_Step_Start'] < df['Largest_Step']//4

    #is the step within 1/4 of largest step size from end of the largest step?
    df['Near_Largest_Step_End_Step'] = df['Distance_To_Largest_Step_End'] < df['Largest_Step']//4

    #If the fragment is close to one end of the largest step, was the step between the two regions, or within one region?
    df['Inter_Largest_Step_Side_Step'] = df['Near_Largest_Step_Start_Step'] != df['Near_Largest_Step_Start_Step'].shift()
    df['Inter_Largest_Step_Side_Step'].iloc[np.where(df['Query_Name'] != df['Query_Name'].shift())] = np.nan
    df['Inter_Largest_Step_Side_Step'].iloc[np.where((df['Near_Largest_Step_Start_Step'] == False) &
                                                         (df['Near_Largest_Step_End_Step'] == False))] = np.nan

    overlap_dfs3[cond] = df    
    
    
summarized_walks = {}

for cond in conditions:
    summary_table = pd.DataFrame()
    grouped_df = overlap_dfs3[cond].groupby('Query_Name')

    #Number of fragments in walk close to the two sides of the largest step
    summary_table['Near_Largest_Step_Start_Frag_Num'] = grouped_df['Near_Largest_Step_Start_Step'].sum().astype(int)
    summary_table['Near_Largest_Step_End_Frag_Num'] = grouped_df['Near_Largest_Step_End_Step'].sum().astype(int)

    #Number of steps between sides of largest step
    summary_table['Inter_Largest_Step_Side_Step_Num'] = grouped_df['Inter_Largest_Step_Side_Step'].sum().astype(int)

    #Indicate how many fragments are close to one of the two sides
    summary_table['Near_Largest_Step_Either_Side_Frag_Num'] = summary_table['Near_Largest_Step_Start_Frag_Num'] + summary_table['Near_Largest_Step_End_Frag_Num']

    #Sum of distances, cis walks only
    summary_table['Sum_Dists'] = grouped_df['dist'].apply(lambda x: np.sum(np.abs(x)))
    #Set to nan for walks with more than 1 chromosome
    summary_table['Chrom_Number'] = grouped_df['Chrom_Number'].mean().astype(int)
    summary_table.loc[summary_table['Chrom_Number'] > 1, 'Sum_Dists'] = np.nan
    summary_table = summary_table.drop(columns = 'Chrom_Number')

    #Number of inter chromosomal steps
    summary_table['Trans_Steps'] = grouped_df['Trans_Step'].sum().astype(int)

    #Number of inter compartment type steps
    summary_table['Inter_Compartment_Type_Steps'] = grouped_df['Inter_Comp_Type_Step'].sum().astype(int)

    #Number of inter compartment index steps in walk - also includes switches between compartment types
    summary_table['Inter_Compartment_Index_Steps'] = grouped_df['Inter_Comp_Index_Step'].sum().astype(int)

    summarized_walks[cond] = summary_table

#Merge full walks with summarized walks for all walk lengths to save
full_walks_with_summary = {}
for cond in conditions:
    full_walks_with_summary[cond] = overlap_dfs3[cond].merge(summarized_walks[cond], left_on = 'Query_Name', right_on = 'Query_Name')


walks_with_summary_firstx_length_fractions = {}
for cond in conditions:
    df = full_walks_with_summary[cond]
        
    #add max fraction of fragments near to one side of largest step
    df1 = df.groupby(['Query_Name']).agg({
        'Near_Largest_Step_Start_Frag_Num' : 'mean',
        'Near_Largest_Step_End_Frag_Num' : 'mean'}).reset_index()
        
    df1['Max_NearOneLargestStepEnd_FracOfFragments'] = df1[['Near_Largest_Step_Start_Frag_Num', 'Near_Largest_Step_End_Frag_Num']].max(axis = 1)/6
    df_filter = df.merge(df1[['Query_Name', 'Max_NearOneLargestStepEnd_FracOfFragments']], on = 'Query_Name')


    walks_with_summary_firstx_length_fractions[cond] = df_filter
        
#Make filtered summary from permuted_walks_with_summary_filtered
walks_summarized_firstx_length = {}
for cond in conditions:
    walks_summarized_firstx_length[cond] = walks_with_summary_firstx_length_fractions[cond][[
        'Query_Name',
        'Iteration',
        'chrom',
        'Fragment_Number',
        'Chrom_Number',
        'Comp_Type_Number',
        'Comp_Index_Number',
        'Near_Largest_Step_Either_Side_Frag_Num',
        'MaxCoord',
        'MinCoord',
        'Largest_Step',
        'Span',
        'Sum_Dists',
        'Walk_Comp_Type',
        'Trans_Steps',
        'Inter_Compartment_Type_Steps',
        'Inter_Compartment_Index_Steps',
        'Inter_Largest_Step_Side_Step_Num',
        'Max_OneCompIndex_FracOfFragments',
        'Max_OneChrom_FracOfFragments',
        'Max_NearOneLargestStepEnd_FracOfFragments']].drop_duplicates()
        
        
#### This section (below) will change depending on specific samples ####

#pickle dicts to save
f = open(f'{outDataDir}/permutations/220517_MRICRF_R1R2R3_permuted_walks_with_summary_first6fragments.iter{iteration}.pkl', 'wb')
pickle.dump(walks_with_summary_firstx_length_fractions, f)
f.close()

f = open(f'{outDataDir}/permutations/220517_MRICRF_R1R2R3_permuted_walks_summarized_first6fragments.iter{iteration}.pkl', 'wb')
pickle.dump(walks_summarized_firstx_length, f)
f.close()
