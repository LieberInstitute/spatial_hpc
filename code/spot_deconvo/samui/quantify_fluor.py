#   Given masks from cellpose (from segmenting the DAPI channel of Visium IF
#   images), quantify mean fluorescence in each non-DAPI channel within each
#   nucleus (dilated to include a region around each nucleus), and save a pandas
#   DataFrame with these values for each nucleus.

import os
os.chdir('/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/')
import matplotlib.pyplot as plt
import numpy as np
from numpy.random import default_rng
import pandas as pd
import seaborn as sns
import tifffile
from scipy.spatial import KDTree
from scipy import ndimage
from skimage.measure import regionprops, regionprops_table
import pyhere
from pathlib import Path
import json

################################################################################
#   Variable definitions
################################################################################

#-------------------------------------------------------------------------------
#   Paths
#-------------------------------------------------------------------------------

spaceranger_dirs = pd.read_csv(pyhere.here("code","spot_deconvo","shared_utilities","samples.txt"), sep = '\t', header=None, names = ['SPpath', 'sample_id', 'brain'])
spaceranger_dirs = spaceranger_dirs.iloc[36:].reset_index(drop=True)

img_path = pyhere.here('processed-data', 'Images', 'VistoSeg', 'VSPG', '{}.tif')
mask_path = pyhere.here('processed-data', 'spot_deconvo', 'cellpose', 'final_masks', '{}_DAPI_seg.npy')
spot_path = pyhere.here(spaceranger_dirs.SPpath, 'outs', 'spatial', 'tissue_positions.csv')
scale_path = pyhere.here(spaceranger_dirs.SPpath, 'outs', 'spatial','scalefactors_json.json')

out_df_path = pyhere.here('processed-data', 'spot_deconvo', 'samui')
plot_dir = pyhere.here("plots", "spot_deconvo", "cellpose", "count_cells")
Path(plot_dir).mkdir(parents=True, exist_ok=True)

#-------------------------------------------------------------------------------
#   Dataset-specific variables
#-------------------------------------------------------------------------------

names = {0: "DAPI", 1: "NeuN", 2: "TMEM119", 3: "GFAP", 4: "OLIG2", 5: "LIP"}
plot_file_type = 'png' # 'pdf' is also supported for higher-quality plots


################################################################################
#   Analysis
################################################################################

# os.environ['SGE_TASK_ID'] = '1'

#-------------------------------------------------------------------------------
#   Read in sample info and adjust paths for this particular sample ID
#-------------------------------------------------------------------------------

#   Different sample IDs are used for different files associated with each
#   sample. Determine both forms of the sample ID for this sample and update
#   path variables accordingly

sample_id_img = spaceranger_dirs.sample_id[int(os.environ['SGE_TASK_ID']) - 1]
img_path = str(img_path).format(sample_id_img)
mask_path = str(mask_path).format(sample_id_img)

spot_path = str(spot_path[int(os.environ['SGE_TASK_ID']) - 1])
out_df_path = pyhere.here(out_df_path,str(sample_id_img + '_df.csv'))

Path(out_df_path).parents[0].mkdir(parents=True, exist_ok=True)

#   Path to JSON from spaceranger including spot size for this sample
json_path = scale_path[int(os.environ['SGE_TASK_ID']) - 1]
with open(json_path) as f: 
    json_data = json.load(f)

#-------------------------------------------------------------------------------
#   Read in spot data
#-------------------------------------------------------------------------------

#   Read in all spot data
raw = pd.read_csv(spot_path,header=None,names=["barcode", "included", "row", "col", "x", "y"],)

#   Take only spots that overlap tissue
raw = raw.iloc[raw.included[raw.included == 1].index].reset_index().drop(columns=["included", "index"])

#-------------------------------------------------------------------------------
#   Quantify mean fluorescence for each channel at each nucleus
#-------------------------------------------------------------------------------

#   Load multi-channel image and masks from segmenting DAPI channel
imgs = tifffile.imread(img_path)
dat = np.load(mask_path,allow_pickle=True).item()
masks = dat['masks']

#   Quantify the mean image fluorescence intensity at each nucleus
#   identified by segmenting the DAPI channel. This is done for each
#   (non-lipofuscin, non-DAPI) channel
its = {
    names[i]: regionprops_table(
        masks, intensity_image=imgs[i], properties=["intensity_mean"]
    )["intensity_mean"]
    for i in range(2, 6)
}

#   Create a table containing the centroids and areas of each mask
#   (nucleus), and add this info to the intensities table
general = regionprops_table(masks, properties=["centroid", "area"])
its["area"] = general["area"]
its["x"] = general["centroid-0"]
its["y"] = general["centroid-1"]

df = pd.DataFrame(its)
#-------------------------------------------------------------------------------
#   Exploratory plot: show the distribution of masks over spots
#-------------------------------------------------------------------------------

#   Plot mask spatial distribution vs. spot distribution; there should be
#   quite a bit of overlap
plt.clf()
plt.scatter(raw["x"], raw["y"], 2)
plt.scatter(df["x"], df["y"], 2)
plt.savefig(
    os.path.join(
        plot_dir, f'mask_spot_overlap_{sample_id_img}.{plot_file_type}'
    )
)

#-------------------------------------------------------------------------------
#   Save relevant data
#-------------------------------------------------------------------------------
df.rename(
    {
        'x': 'y',
        'y': 'x',
        'Unnamed: 0': 'id'
    },
    axis = 1, inplace = True
)
