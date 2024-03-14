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
from scipy.spatial import KDTree

#this_donor = "Br8325"
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
gene_df.index.name = None

domain_df = spgP.obs[['domain']].copy()
sample_df = sample_df.loc[tissue_positions_filtered.index]
deconvo_df = deconvo_df.loc[tissue_positions_filtered.index]
tissue_positions_arranged = tissue_positions_filtered.reindex(gene_df.index)
################################################################################
#  remove overlapping spots
################################################################################
tissue_positions_path = Path(here("processed-data", "VSPG_image_stitching", "combined_"+this_donor, f'tissue_positions_{this_donor}.csv'))
tissue_positions = pd.read_csv(tissue_positions_path ,index_col = 0).rename({'pxl_row_in_fullres': 'y', 'pxl_col_in_fullres': 'x'},axis = 1)
tissue_positions.index.name = None
tissue_positions = tissue_positions[['x', 'y']].astype(int)
SPOT_DIAMETER_M = 55e-6
m_per_px = 4.971263040387764e-07
 
 # Build KD tree for nearest neighbor search.
kd = KDTree(tissue_positions[["x", "y"]].values)
 # Query the KDTree for pairs of points within the threshold distance
overlapping_pairs = pd.DataFrame([(tissue_positions.index[x], tissue_positions.index[y]) for x, y in kd.query_pairs(150)])

unique_items = set(item.split('-1')[1] for col in overlapping_pairs.columns for item in overlapping_pairs[col])
tissue_positions_f = tissue_positions.loc[~tissue_positions.index.isin(overlapping_pairs[0])]

common_indices = gene_df.index.intersection(tissue_positions_f.index)
tissue_positions_filtered = tissue_positions_f.loc[common_indices]

domain_df = domain_df.loc[tissue_positions_filtered.index]
gene_df = gene_df.loc[tissue_positions_filtered.index]
tissue_positions_arranged = tissue_positions_filtered.reindex(gene_df.index)
 
################################################################################
#   Use the Samui API to create the importable directory for this combined "sample"
################################################################################
img_channels = ['DAPI', 'Alexa_488', 'Alexa_555', 'Alexa_594', 'Alexa_647', 'Autofluorescence']
#default_channels = {'blue': 'DAPI', 'green': 'Alexa_488', 'yellow': 'Alexa_555', 'red': 'Alexa_594', 'magenta': 'Alexa647', 'cyan': 'Autofluorescence'}
default_channels = {'blue': 'DAPI', 'green': 'Alexa_488'}
default_gene = 'SLC17A7'

assert default_gene in gene_df.columns, "Default gene not in AnnData"

notes_md_url = Url('/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/code/VSPG_image_stitching/feature_notes.md')
this_sample = Sample(name = samui_dir.name, path = samui_dir, notesMd = notes_md_url)

this_sample.add_coords(tissue_positions_arranged, name = "coords", mPerPx = m_per_px, size = SPOT_DIAMETER_M)

img_fullres = Path(here("processed-data", "VSPG_image_stitching", "combined_"+this_donor, f'combined_{this_donor}.tif'))

this_sample.add_image(
    tiff = img_fullres,
    channels = img_channels,
    defaultChannels = default_channels,
    scale = m_per_px
)
this_sample.add_csv_feature(sample_df, name = "Capture areas", coordName = "coords", dataType = "categorical")
this_sample.add_csv_feature(domain_df, name = "PRECAST domains", coordName = "coords", dataType = "categorical")
this_sample.add_csv_feature(deconvo_df, name = "broad deconvolution", coordName = "coords", dataType = "quantitative")
this_sample.add_csv_feature(mid_deconvo_df, name = "mid deconvolution", coordName = "coords", dataType = "quantitative")
this_sample.add_csv_feature(fine_deconvo_df, name = "fine deconvolution", coordName = "coords", dataType = "quantitative")

this_sample.add_chunked_feature(gene_df, name = "Genes", coordName = "coords", dataType = "quantitative")

this_sample.set_default_feature(group = "Genes", feature = default_gene)



this_sample.write()

session_info.show()