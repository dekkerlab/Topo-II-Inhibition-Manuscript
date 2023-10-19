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
import scipy
import pickle
from numpy import diff
import sys

from pandas import read_csv
from sklearn.utils import resample
from sklearn.tree import DecisionTreeClassifier
from sklearn.metrics import accuracy_score
from matplotlib import pyplot

iteration = sys.argv[1]
outDataDir = sys.argv[2]
window_size = float(sys.argv[3])

conditions = [
    't0Mit_R1',
    't2_R1',
    't4DMSO_R1',
    't4ICRF_R1',
    't8DMSO_R1',
    't8ICRF_R1',
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
    't8ICRF_R3'
]

long_names = {
    't0Mit_R1' : 'TI-MC3C-Dpn-t0Mit-4-30',
    't2_R1' : 'TI-MC3C-Dpn-t2-4-30',
    't4DMSO_R1' : 'TI-MC3C-Dpn-t4DMSO-4-30',
    't4ICRF_R1' : 'TI-MC3C-Dpn-t4ICRF-4-30',
    't8DMSO_R1' : 'TI-MC3C-Dpn-t8DMSO-4-30',
    't8ICRF_R1' : 'TI-MC3C-Dpn-t8ICRF-4-30',
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

repdict = {
    't0Mit_R1' : 'R1',
    't2_R1' : 'R1',
    't4DMSO_R1' : 'R1',
    't4ICRF_R1' : 'R1',
    't8DMSO_R1' : 'R1',
    't8ICRF_R1' : 'R1',
    't0Mit_R2' : 'R2',
    't2_R2' : 'R2',
    't4DMSO_R2' : 'R2',
    't4ICRF_R2' : 'R2',
    't8DMSO_R2' : 'R2',
    't8ICRF_R2' : 'R2',
    't0Mit_R3' : 'R3',
    't2_R3' : 'R3',
    't4DMSO_R3' : 'R3',
    't4ICRF_R3' : 'R3',
    't8DMSO_R3' : 'R3',
    't8ICRF_R3' : 'R3',
}

labeldict = {
    't0Mit_R1' : 't0 Prometa',
    't2_R1' : 't2 Ana/Telo',
    't4DMSO_R1' : 't4 G1 DMSO',
    't4ICRF_R1' : 't4 G1 ICRF-193',
    't8DMSO_R1' : 't8 G1 DMSO',
    't8ICRF_R1' : 't8 G1 ICRF-193',
    't0Mit_R2' : 't0 Prometa',
    't2_R2' : 't2 Ana/Telo',
    't4DMSO_R2' : 't4 G1 DMSO',
    't4ICRF_R2' : 't4 G1 ICRF-193',
    't8DMSO_R2' : 't8 G1 DMSO',
    't8ICRF_R2' : 't8 G1 ICRF-193',
    't0Mit_R3' : 't0 Prometa',
    't2_R3' : 't2 Ana/Telo',
    't4DMSO_R3' : 't4 G1 DMSO',
    't4ICRF_R3' : 't4 G1 ICRF-193',
    't8DMSO_R3' : 't8 G1 DMSO',
    't8ICRF_R3' : 't8 G1 ICRF-193',
}

comp_types = ['A', 'B', 'AB']
good_chroms = ['chr4', 'chr14', 'chr17', 'chr18', 'chr20', 'chr21']

f = open(f'{outDataDir}/data/permutations/220517_MRICRF_R1R2R3_permuted_walks_summarized_first6fragments.iter{iteration}.pkl', 'rb')
permuted_walks_summarized = pickle.load(f)
f.close()
    
f = open(f'{outDataDir}/data/permutations/220517_MRICRF_R1R2R3_permuted_walks_with_summary_first6fragments.iter{iteration}.pkl', 'rb')
permuted_walks_with_summary = pickle.load(f)
f.close()

#combining T1 and T2 reads for t0 Mit R1
permuted_walks_with_summary['t0Mit_R1'] = {}
permuted_walks_with_summary['t0Mit_R1'] = permuted_walks_with_summary['t0Mit_R1_T1'].append(
    permuted_walks_with_summary['t0Mit_R1_T2'],
    ignore_index = True) 
    
#combining T1 and T2 reads for t0 Mit R1
permuted_walks_summarized['t0Mit_R1'] = {}
permuted_walks_summarized['t0Mit_R1'] = permuted_walks_summarized['t0Mit_R1_T1'].append(
    permuted_walks_summarized['t0Mit_R1_T2'],
    ignore_index = True) 

#set up windows
start_dist = 0 #start of first window
end_dist = 3e7 #end of last window

window_step = 1e6

Intermingling_Sliding_Window = pd.DataFrame()

#set up windows
for i, start_bp in enumerate(range(int(start_dist), int(end_dist-window_size), int(window_step))):
    end_bp = start_bp + int(window_size)

    #iterate through conditions for each window, all compartments
    for cond in conditions:
        df = permuted_walks_with_summary[cond]
        df['Query_Fragment_Length'] = df['Query_End'] - df['Query_Start']
        grouped_walks = df.groupby('Query_Name')
        walks_min_mapq = grouped_walks.agg({'Mapping_Quality' : 'min'})
        good_walks_mapq = walks_min_mapq[walks_min_mapq['Mapping_Quality'] > 59] #use this to filter for mapq
        walks_frac_map = grouped_walks.agg({'Match_Length' : 'sum',
                                            'Query_Fragment_Length' : 'sum',
                                            'Alignment_Length' : 'sum'
                                           })
        walks_high_frac_map = walks_frac_map[
        (walks_frac_map['Match_Length']/walks_frac_map['Query_Fragment_Length']) > 0.8] #use this to filter for fraction mapped

        df2 = permuted_walks_summarized[cond].copy()
        df_cond = df2[
            (df2['Chrom_Number'] == 1) &
            (df2['Walk_Comp_Type'].isin(['A', 'B', 'AB'])) &
            (df2['Query_Name'].isin(good_walks_mapq.index)) &
            (df2['Query_Name'].isin(walks_high_frac_map.index)) &
            (df2['chrom'].isin(good_chroms)) &
            (df2['Near_Largest_Step_Either_Side_Frag_Num'] == 6) &
            (df2['Largest_Step'] >= start_bp) &
            (df2['Largest_Step'] < end_bp) &
            (df2['Max_NearOneLargestStepEnd_FracOfFragments'] == 5/6)
        ][['Inter_Largest_Step_Side_Step_Num']]
        df_cond['Condition'] = cond
        df_cond['Label'] = labeldict[cond]
        df_cond['Window_Midpoint'] = (start_bp + end_bp)//2
        df_cond['Walk_Comp_Type'] = 'All'
        df_cond['Replicate'] = repdict[cond]

        #add together into one dataframe
        Intermingling_Sliding_Window = Intermingling_Sliding_Window.append(df_cond, ignore_index = True)

        for comp in comp_types:#iterate through comp types
            df_comp = df2[
                (df2['Chrom_Number'] == 1) &
                (df2['Walk_Comp_Type'] == comp) &
                (df2['Query_Name'].isin(good_walks_mapq.index)) &
                (df2['Query_Name'].isin(walks_high_frac_map.index)) &
                (df2['chrom'].isin(good_chroms)) &
                (df2['Near_Largest_Step_Either_Side_Frag_Num'] == 6) &
                (df2['Largest_Step'] >= start_bp) &
                (df2['Largest_Step'] < end_bp) &
                (df2['Max_NearOneLargestStepEnd_FracOfFragments'] == 5/6)
            ][['Inter_Largest_Step_Side_Step_Num']]
            df_comp['Condition'] = cond
            df_comp['Label'] = labeldict[cond]
            df_comp['Window_Midpoint'] = (start_bp + end_bp)//2
            df_comp['Walk_Comp_Type'] = comp
            df_comp['Replicate'] = repdict[cond]
            
            #add together into one dataframe
            Intermingling_Sliding_Window = Intermingling_Sliding_Window.append(df_comp, ignore_index = True)
    
f = open(f'{outDataDir}/data/permutations/220518_MRICRF_R1R2R3_PermutedInterminglingSlidingWindowSweep_window{window_size}.iter{iteration}.pkl', 'wb')
pickle.dump(Intermingling_Sliding_Window, f)
f.close()
    