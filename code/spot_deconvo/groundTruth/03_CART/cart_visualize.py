import os
os.chdir('/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/')
import pandas as pd
import numpy as np
import copy

import pyhere
from pathlib import Path

import matplotlib.pyplot as plt
import graphviz


#############
## paths
############

mask_path = pyhere.here('processed-data', 'spot_deconvo', 'groundTruth', '01_cellpose', 'final_masks', '{}' + '_DAPI_seg.npy')
plot_dir = pyhere.here('plots','spot_deconvo', 'groundTruth', '03_CART')
img_path = pyhere.here('processed-data', 'spot_deconvo', 'groundTruth', '02_samui_manual_annotation', '{}.tif')
df_path = pyhere.here('processed-data', 'spot_deconvo', 'groundTruth', '03_CART', '{}' + 'cell_metrics.csv')

spaceranger_dirs = pd.read_csv(pyhere.here("code","spot_deconvo","shared_utilities","samples.txt"), sep = '\t', header=None, names = ['SPpath', 'sample_id', 'brain'])
spaceranger_dirs = spaceranger_dirs.iloc[36:].reset_index(drop=True)
sample_ids = spaceranger_dirs.sample_id
sample_id = sample_ids[int(os.environ['SGE_TASK_ID']) - 1]
