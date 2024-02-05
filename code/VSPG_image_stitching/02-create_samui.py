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
