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

#-------------------------------------------------------------------------------
#   Visually verify cell-type calls
#-------------------------------------------------------------------------------
################################################################################
#   Functions
################################################################################

def plot_roi(img, props, indices, vmax: int = 128, pad: int = 25):
    #   Set up the plot
    fig, axs = plt.subplots(nrows=len(indices), ncols=6, figsize=(18, len(indices) * 3))
    axs = axs.flatten()
    
    #   Loop through each nucleus, each of which will be a row in the final plot
    j = 0
    for idx in indices:
        bbox = props[idx]["bbox"]
        roi = props[idx]["image"]
        
        axs[j].imshow(np.pad(roi, (pad, pad), mode="constant", constant_values=0), aspect="equal")
        axs[j].grid(False)
        if j < 6:
            axs[j].set_title("Mask")
        
        j += 1
        
        for i in range(1, 6):
            a = axs[j].imshow(
                img[
                    i-1,
                    max(0, bbox[0] - pad) : min(img.shape[1], bbox[2] + pad),
                    max(0, bbox[1] - pad) : min(img.shape[2], bbox[3] + pad),
                ],
                vmax=vmax,
                aspect="equal",
            )
            plt.colorbar(a, ax = axs[j])
            
            if j < 6:
                axs[j].set_title(names[i])
            
            axs[j].grid(False)
            j += 1
    #
    return fig
