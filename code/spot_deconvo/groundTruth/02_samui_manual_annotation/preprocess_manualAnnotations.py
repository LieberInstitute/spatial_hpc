#   Standardize and clean manual cell-type labels and re-write CSVs for
#   immediate downstream use
#
#   1. Drop NA labels and use consistent names for cell types
#   2. Ensure exactly 30 labels per cell type per sample
import os
os.chdir('/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/')

import pandas as pd
import numpy as np
import pyhere
from pathlib import Path

################################################################################
#   Variable definitions
################################################################################

df_path = pyhere.here("processed-data", "spot_deconvo", "groundTruth", "02_samui_manual_annotation", "{}" + "_df.csv")
manual_label_path = pyhere.here('processed-data', 'spot_deconvo', 'groundTruth', '02_samui_manual_annotation', "annotations_" + "{}" + "_coords_SCP.csv")
manual_label_path_out = pyhere.here('processed-data', 'spot_deconvo', 'groundTruth', '02_samui_manual_annotation', "annotations_SCP_processed.csv")

#   Define expected cell-type labels (after cleaning), expected columns
#   in the fluorescence-intensity data frame, and expected number of counts
#   for each cell-type label
expected_df_cols = ['id','NeuN', 'TMEM119', 'GFAP', 'OLIG2','area', 'y', 'x']
#expected_label_counts = 30


################################################################################
#   Clean data
################################################################################

#   Read in the list of sample IDs for the ID data
sample_ids = pd.Series(["V12D07-332_D1", "V12D07-335_D1"], dtype=object)
df = pd.DataFrame()

for sample_id in sample_ids:
    #   Determine paths for this sample
    this_df_path = str(df_path).format(sample_id)
    this_manual_label_path = str(manual_label_path).format(sample_id)
       
    #   Read in all required CSVs
    this_df = pd.read_csv(this_df_path, names = expected_df_cols, header=0)
    this_manual_labels = pd.read_csv(this_manual_label_path)
        
    #   Check the columns are all as expected
    #assert all([x == y for x, y in zip(this_df.columns.tolist(), expected_df_cols)])
    
    #   Fix indices (index rows by cell ID)
    this_manual_labels.index = this_manual_labels['id']
    this_df.index = this_df['id']
    
    this_df['label'] = this_manual_labels['label']
    this_df['label_sample'] = this_df['label'] + '_' + sample_id
    df = pd.concat([df, this_df.dropna()])
    
#   Verify we have the correct amount of cells
#    assert(df.shape[0] == len(sample_ids) * expected_num_labels * num_cell_types)
label_counts = df['label'].value_counts()
label_counts        
# neuron       134
# oligo        117
# other         97
# astrocyte     51
# microglia     36
# Name: label, dtype: int64

    
    #   Write a clean copy of both sets of manual labels
    df.to_csv(manual_label_path_out, index = False)
    