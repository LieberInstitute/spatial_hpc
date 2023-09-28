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
    
#####################################################################################
spaceranger_dirs = pd.read_csv(pyhere.here("code","spot_deconvo","shared_utilities","samples.txt"), sep = '\t', header=None, names = ['SPpath', 'sample_id', 'brain'])
spaceranger_dirs = spaceranger_dirs.iloc[36:].reset_index(drop=True)
sample_id = spaceranger_dirs.sample_id[int(os.environ['SGE_TASK_ID']) - 1]

df_path = str(df_path).format(sample_id)
df = pd.read_csv(df_path)
df.rename({'Unnamed: 0': 'id'}, axis = 1, inplace = True)
df.index = df['id']

mask_path = str(mask_path).format(sample_id)
dat = np.load(mask_path,allow_pickle=True).item()
masks = dat['masks']

props = regionprops(masks)
examples_per_type = 5

rng = default_rng(0)
img_path = str(img_path).format(sample_id)
imgs = tifffile.imread(img_path)
plot_file_type = 'pdf' # 'png'
cell_types = {"NeuN": "neuron","OLIG2": "oligo","TMEM119": "microglia","GFAP": "astrocyte"}

for cell_type in df['cell_type'].unique():
    # Randomly pick 5 distinct rows for this cell type
    indices = rng.choice(df[df['cell_type'] == cell_type].index, examples_per_type, replace=False)
    # Print numeric intensities for these cells
    print(f'Intensities for {examples_per_type} random {cell_type} cells:')
    print(df.loc[indices][list(cell_types.keys())])
    # Plot intensities
    #names = {1: "DAPI", 2: "NeuN", 3: "TMEM119", 4: "GFAP", 5: "OLIG2", 6: "LIP"}
    fig = plot_roi(imgs, props, indices)
    plt.suptitle(f'Cells classified as {cell_type}')
    fig.savefig(os.path.join(plot_dir, f'{cell_type}_{sample_id}.{plot_file_type}'))
    plt.close('all')

