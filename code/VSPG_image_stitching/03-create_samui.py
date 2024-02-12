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

this_donor = "Br8325"
this_donor = "Br3942"

samui_dir = Path(here('processed-data', 'VSPG_image_stitching', this_donor))
samui_dir.mkdir(parents = True, exist_ok = True)

################################################################################
#   Gather gene-expression data into a DataFrame to later as a feature
################################################################################
spg_path = here("processed-data", "VSPG_image_stitching", "spg.h5ad")
#   Read in AnnData and subset to this_donor
spg = sc.read(spg_path)
#path_groups = spg.obs['path_groups'].cat.categories
spgP = spg[spg.obs['brnum'] == this_donor+'_VSPG', :]
spgP.obs.index = spgP.obs.index.str.replace('_'+this_donor, '')
#   Convert the sparse gene-expression matrix to pandas DataFrame, with the
#   gene symbols as column names
gene_df = pd.DataFrame(
    spgP.X.toarray(),
    index = spgP.obs.index,
    columns = spgP.var['gene_name']
)
#   Some gene symbols are actually duplicated. Just take the first column in
#   any duplicated cases
gene_df = gene_df.loc[: , ~gene_df.columns.duplicated()].copy()

sample_df = spgP.obs[['sample_id', 'domain']].copy()
#sample_df['domain_codes'] = sample_df['domain'].cat.codes
#sample_df['domain_codes'] = sample_df['domain_codes'].astype(int)
#sample_df['domain_codes'] = sample_df['domain_codes'].astype(int)

################################################################################
#   Use the Samui API to create the importable directory for this combined
#   "sample"
################################################################################
img_channels = ['DAPI', 'Alexa_488', 'Alexa_555', 'Alexa_594', 'Alexa_647', 'Autofluorescence']
default_channels = {'blue': 'DAPI', 'green': 'Alexa_488', 'yellow': 'Alexa_555', 'red': 'Alexa_594', 'magenta': 'Alexa647', 'cyan': 'Autofluorescence'}
default_gene = 'SNAP25'

assert default_gene in gene_df.columns, "Default gene not in AnnData"

this_sample = Sample(name = samui_dir.name, path = samui_dir)

tissue_positions_path = Path(here("processed-data", "VSPG_image_stitching", "combined_"+this_donor, f'tissue_positions_{this_donor}.csv'))
tissue_positions = pd.read_csv(tissue_positions_path ,index_col = 0).rename(
                {'pxl_row_in_fullres': 'y', 'pxl_col_in_fullres': 'x'},axis = 1)
tissue_positions.index.name = None
tissue_positions = tissue_positions[['x', 'y']].astype(int)
tissue_positions_filtered = tissue_positions.loc[gene_df.index]
tissue_positions_arranged = tissue_positions.reindex(gene_df.index)
SPOT_DIAMETER_M = 55e-6
m_per_px = 4.971263040387764e-07

this_sample.add_coords(tissue_positions_arranged, name = "coords", mPerPx = m_per_px, size = SPOT_DIAMETER_M)

img_fullres = Path(here("processed-data", "VSPG_image_stitching", "combined_"+this_donor, f'combined_{this_donor}.tif'))

this_sample.add_image(
    tiff = img_fullres,
    channels = img_channels,
    defaultChannels = default_channels,
    scale = m_per_px
)
this_sample.add_chunked_feature(gene_df, name = "Genes", coordName = "coords", dataType = "quantitative")

this_sample.set_default_feature(group = "Genes", feature = default_gene)

this_sample.add_csv_feature(sample_df, name = "Sample Info", coordName = "coords", dataType = "categorical")

this_sample.write()

session_info.show()