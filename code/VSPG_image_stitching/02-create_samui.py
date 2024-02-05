import os
os.chdir('/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/')

from pyhere import here
from pathlib import Path
import session_info

import numpy as np
import pandas as pd
import json
import sys
from loopy.sample import Sample
import tifffile
from PIL import Image
import re
import matplotlib.pyplot as plt

import scanpy as sc
from rasterio import Affine
from loopy.utils.utils import remove_dupes, Url
import re

sample_info_path = here('processed-data', 'VSPG_image_stitching', 'transformations.csv')
samui_dir = Path(here('processed-data', 'VSPG_image_stitching', this_donor))
samui_dir.mkdir(parents = True, exist_ok = True)
combined_dir = Path(here('processed-data', 'VSPG_image_stitching',f'combined_{this_donor}'))
combined_dir.mkdir(parents = True, exist_ok = True)

img_out_fullres = Path(here(combined_dir, f'combined_{this_donor}.tif'))
tissue_out_path = Path(here(combined_dir, f'tissue_positions_{this_donor}.csv'))
img_out_lowres = Path(here(combined_dir, f'combined_tissue_lowres_{this_donor}.png'))
json_out_path = Path(here(combined_dir, 'scalefactors_json.json'))

#   55-micrometer diameter for Visium spot but 65-micrometer spot diameter used
#   in 'spot_diameter_fullres' calculation for spaceranger JSON. The
#   difference between 55 and 65 does indeed exist and is properly documented,
#   but is likely a bug in the sense the choice was probably unintentional
#   https://kb.10xgenomics.com/hc/en-us/articles/360035487812-What-is-the-size-of-the-spots-on-the-Visium-Gene-Expression-Slide-
#   https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/output/spatial

spg_path = here("processed-data", "VSPG_image_stitching", "spg.h5ad")

SPOT_DIAMETER_M = 55e-6
SPOT_DIAMETER_JSON_M = 65e-6

LOWRES_MAX_SIZE = 1200
BACKGROUND_COLOR = 0

#   Read in sample info and subset to samples of interest
sample_info = pd.read_csv(sample_info_path, index_col = 0)
sample_info = sample_info.loc[
    (sample_info['Brain'] == this_donor) &
    ~sample_info['xml_path'].isna(),# &
    #sample_info['In analysis'],
    :
]