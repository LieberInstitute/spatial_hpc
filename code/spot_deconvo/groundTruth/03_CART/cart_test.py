#   Given mean fluorescence intensities for each nucleus, classify cell types
#   present in each spot.

import os
os.chdir('/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/')

import matplotlib.pyplot as plt
import numpy as np
from numpy.random import default_rng
import pandas as pd
import seaborn as sns
import tifffile
from scipy import ndimage
from skimage.measure import regionprops, regionprops_table
import pyhere
from pathlib import Path
import pickle
from scipy.spatial import KDTree
import json

################################################################################
#   Paths
################################################################################
spaceranger_dirs = pd.read_csv(pyhere.here("code","spot_deconvo","shared_utilities","samples.txt"), sep = '\t', header=None, names = ['SPpath', 'sample_id', 'brain'])
spaceranger_dirs = spaceranger_dirs.iloc[36:].reset_index(drop=True)
df_path = pyhere.here('processed-data', 'spot_deconvo', 'groundTruth', '02_samui_manual_annotation', '{}' + '_df.csv')

#   Trained DecisionTreeClassifier path
model_path = pyhere.here('processed-data', 'spot_deconvo', 'groundTruth', '03_CART', 'decision_tree.pkl')

#   Main output: rows are spots and columns are cell types (values are counts)
clusters_path = pyhere.here('processed-data', 'spot_deconvo', 'groundTruth', '03_CART', '{}', 'clusters.csv')

#   Secondary output: rows are cells and columns are metrics/info
cells_path = pyhere.here('processed-data', 'spot_deconvo', 'groundTruth', '03_CART', '{}', 'cell_metrics.csv')

# os.environ['SGE_TASK_ID'] = '1'

sample_id = spaceranger_dirs.sample_id[int(os.environ['SGE_TASK_ID']) - 1]
spot_path = pyhere.here(spaceranger_dirs.SPpath[int(os.environ['SGE_TASK_ID']) - 1],'outs','spatial','tissue_positions.csv')
json_path = pyhere.here(spaceranger_dirs.SPpath[int(os.environ['SGE_TASK_ID']) - 1],'outs','spatial','scalefactors_json.json')
df_path = str(df_path).format(sample_id)

clusters_path = str(clusters_path).format(sample_id)
Path(clusters_path).parents[0].mkdir(parents=True, exist_ok=True)

cells_path = str(cells_path).format(sample_id)
Path(cells_path).parents[0].mkdir(parents=True, exist_ok=True)
