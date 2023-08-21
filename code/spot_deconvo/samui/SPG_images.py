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
default_gene = 'PPFIA2'

#   Names of continuous features expected to be columns in the observation data (colData) of the AnnData
spe_cont_features = ['NDAPI', 'CNDAPI', 'PDAPI']

spe_path = here('processed-data', 'spot_deconvo', 'shared_utilities', 'spe.h5ad')
IMG_path = here('processed-data', 'spot_deconvo', 'samui', '{}.tif')
spaceranger_dirs = pd.read_csv(here("code","spot_deconvo","shared_utilities","samples.txt"), sep = '\t', header=None, names = ['SPpath', 'sample_id', 'brain'])
spaceranger_dirs =spaceranger_dirs[36:].reset_index(drop=True)
JSON_path = here(spaceranger_dirs.SPpath, 'outs', 'spatial','scalefactors_json.json')
OUT_dir = here('processed-data', 'spot_deconvo', 'samui', '{}')

################################################################################
#   Read in sample info and clean
################################################################################

#   Subset all types of IDs to this sample only
sample_id = spaceranger_dirs.sample_id.iloc[int(os.environ['SGE_TASK_ID']) - 1]

#   Update paths for this sample ID
out_dir = Path(str(OUT_dir).format(sample_id))
json_path = JSON_path.iloc[int(os.environ['SGE_TASK_ID']) - 1]
img_path = Path(str(IMG_path).format(sample_id))

out_dir.mkdir(exist_ok = True)

#   All paths should exist
assert all([x.exists() for x in [out_dir, json_path, img_path]])

################################################################################
#   Read in scale-factors info
################################################################################

#   Read in the spaceranger JSON to calculate meters per pixel for
#   the full-resolution image
with open(json_path, 'r') as f:
    spaceranger_json = json.load(f)

m_per_px = spot_diameter_m / spaceranger_json['spot_diameter_fullres']

################################################################################
#   Gather gene-expression data into a DataFrame to later as a feature
################################################################################

#   Read in AnnData and subset to this sample
spe = sc.read(spe_path)
path_groups = spe.obs['path_groups'].cat.categories
spe = spe[spe.obs['sample_id'] == sample_id_spaceranger, :]
spe.obs.index.name = "barcode"

#   Convert the sparse gene-expression matrix to pandas DataFrame, with the
#   gene symbols as column names
gene_df = pd.DataFrame(
    spe.X.toarray(),
    index = spe.obs.index,
    columns = spe.var['gene_name']
)

#   Some gene symbols are actually duplicated. Just take the first column in
#   any duplicated cases
gene_df = gene_df.loc[: , ~gene_df.columns.duplicated()].copy()

#   Samui seems to break when using > ~ 5,000 genes. Take just the genes where
#   at least 10% of spots have nonzero counts
gene_df = gene_df.loc[:, np.sum(gene_df > 0, axis = 0) > (gene_df.shape[0] * 0.1)].copy()

assert default_gene in gene_df.columns, "Default gene not in AnnData"

print('Using {} genes as features.'.format(gene_df.shape[1]))

################################################################################
#   Split 'path_groups' column into binary columns for each of its values
################################################################################

#   Circumvent a Samui bug (https://github.com/chaichontat/samui/issues/84);
#   turn the categorical column 'path_groups' into several numeric columns with
#   just values of 0 and 1
path_df = pd.DataFrame()
for path_group in path_groups:
    path_df[path_group] = (spe.obs['path_groups'] == path_group).astype(int)

################################################################################
#   Use the Samui API to create the importable directory for this sample
################################################################################

this_sample = Sample(name = sample_id, path = out_dir)

this_sample.add_coords(
    spe.obsm['spatial'].rename(
        columns = {'pxl_col_in_fullres': 'x', 'pxl_row_in_fullres': 'y'}
    ),
    name = "coords", mPerPx = m_per_px, size = spot_diameter_m
)

#   Add the IF image for this sample
this_sample.add_image(
    tiff = img_path, channels = img_channels, scale = m_per_px,
    defaultChannels = default_channels
)

#   Add gene expression results (multiple columns) as a feature
this_sample.add_csv_feature(
    gene_df, name = "Genes", coordName = "coords", dataType = "quantitative"
)

#   Add additional requested observational columns (colData columns)
this_sample.add_csv_feature(
    spe.obs[spe_cont_features], name = "Spot Coverage", coordName = "coords",
    dataType = "quantitative"
)


this_sample.set_default_feature(group = "Genes", feature = default_gene)

this_sample.write()