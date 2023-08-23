from pathlib import Path
import os
os.chdir('/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/')
from pyhere import here
import json

import scanpy as sc

import numpy as np
import pandas as pd
from rasterio import Affine

from loopy.sample import Sample
from loopy.utils.utils import remove_dupes, Url

spot_diameter_m = 55e-6 # 55-micrometer diameter for Visium spot
img_channels = ['DAPI', 'Alexa_488', 'Alexa_555', 'Alexa_594', 'Alexa_647', 'Autofluorescence', 'segmented_DAPI']
default_channels = {'blue': 'DAPI', 'green': 'Alexa_488', 'yellow': 'Alexa_555', 'red': 'Alexa_594', 'magenta': 'Alexa647', 'cyan': 'Autofluorescence', 'white': 'segmented_DAPI'}
#default_gene = 'PPFIA2'

#   Names of continuous features expected to be columns in the observation data (colData) of the AnnData
# spe_cont_features = ['NDAPI', 'CNDAPI', 'PDAPI']
inten_features = ['TMEM119', 'GFAP', 'OLIG2', 'LIP', 'area']
