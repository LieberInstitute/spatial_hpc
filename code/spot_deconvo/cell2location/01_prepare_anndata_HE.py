import sys
import os
os.chdir('/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/')
import scanpy as sc
import anndata
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

import cell2location
from cell2location.utils.filtering import filter_genes

from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42 # enables correct plotting of text
import seaborn as sns

import pyhere
from pathlib import Path
from PIL import Image
import json
import pandas as pd
import pprint
################################################################################
#   Variable definitions
################################################################################

#cell_group = "broad" 
#subtype = "_class"

cell_group = "layer" 
subtype = "_celltype_class1_noHATAGABAAmy."

#sc_path = pyhere.here("processed-data", "spot_deconvo", "shared_utilities","sce_class.h5ad")
sc_path = pyhere.here("processed-data", "spot_deconvo", "shared_utilities","sce_class_noHATAGABA.h5ad")
sp_path = pyhere.here("processed-data", "spot_deconvo", "shared_utilities","spe.h5ad")
#spg_path = pyhere.here("processed-data", "spot_deconvo", "shared_utilities","spg.h5ad")

processed_dir = pyhere.here("processed-data", "spot_deconvo", "cell2location", "HE", cell_group, "newclass")
plot_dir = pyhere.here("plots", "spot_deconvo", "cell2location", "HE", cell_group, "newclass")
Path(plot_dir).mkdir(parents=True, exist_ok=True)
Path(processed_dir).mkdir(parents=True, exist_ok=True)

#   Directory containing hires image and a JSON containing scale factors and spot size for a given sample. 
spaceranger_dirs = pd.read_csv(pyhere.here("code","spot_deconvo","shared_utilities","samples.txt"), sep = '\t', header=None, names = ['SPpath', 'sample_id', 'brain'])
spaceranger_dirs.SPpath = pyhere.here(spaceranger_dirs.SPpath, 'outs', 'spatial')

marker_path = pyhere.here("processed-data", "spot_deconvo", "shared_utilities", "markers_" + cell_group + subtype + ".txt")
#cell_type_var = 'broad.class'
cell_type_var = 'cell.class'

#   Naming conventions used for different columns in the spatial AnnData
sample_id_var = 'sample_id'          # in spatial object only
ensembl_id_var = 'gene_id'           # in both spatial and single-cell objects
gene_symbol_var = 'gene_name'        # in both spatial and single-cell objects
spatial_coords_names = ['pxl_col_in_fullres', 'pxl_row_in_fullres']

plot_file_type = 'pdf'

################################################################################
#   Preprocessing
################################################################################
#adata_ref.obs.iloc[:,27].categories
#adata_ref.obs.iloc[:,27].cat.categories

#  Load AnnDatas
print('Loading AnnDatas...')
adata_vis = sc.read_h5ad(sp_path)
adata_ref = sc.read_h5ad(sc_path)

# update Visium AnnData object to match structure in cell2location tutorial
adata_vis.obs['sample'] = adata_vis.obs[sample_id_var]
# rename genes to ENSEMBL
adata_vis.var['SYMBOL'] = adata_vis.var[gene_symbol_var]
adata_vis.var_names = adata_vis.var[ensembl_id_var]
adata_vis.var_names.name = None

# find mitochondria-encoded (MT) genes
adata_vis.var['MT_gene'] = [gene.startswith('MT-') for gene in adata_vis.var['SYMBOL']]

# remove MT genes for spatial mapping (keeping their counts in the object).
adata_vis.obsm['MT'] = adata_vis[:, adata_vis.var['MT_gene'].values].X.toarray()
adata_vis = adata_vis[:, ~adata_vis.var['MT_gene'].values]

# Spatial AnnData needs unique indices. Rather than using barcode (repeated for every sample), use "key" (barcode + sample ID)
adata_vis.obs_names = adata_vis.obs['key']
adata_vis.obs_names.name = None

# Use ENSEMBL as gene IDs to make sure IDs are unique and correctly matched
adata_ref.var['SYMBOL'] = adata_ref.var[gene_symbol_var]
adata_ref.var.index = adata_ref.var[ensembl_id_var]
adata_ref.var_names = adata_ref.var[ensembl_id_var]
adata_ref.var.index.name = None

#   Subset to marker genes
with open(marker_path, 'r') as f:
    selected = f.read().splitlines()

adata_ref = adata_ref[:, selected].copy()

#-------------------------------------------------------------------------------
#   Attach hi-res images and scaleFactors to spatial AnnData
#-------------------------------------------------------------------------------

adata_vis.uns['spatial'] = {}

for sample_id in adata_vis.obs['sample_id'].cat.categories:
    spaceranger_dir = spaceranger_dirs[spaceranger_dirs.sample_id == sample_id].SPpath.values[0]
    # print(spaceranger_dirs.SPpath.values[0])
    #   Path to JSON from spaceranger including spot size for this sample
    json_path = pyhere.here(str(spaceranger_dir), 'scalefactors_json.json')
    print(json_path)
    
    with open(json_path) as f: 
        json_data = json.load(f)
    
    #   Read in high-res image as numpy array with values in [0, 1] rather than [0, 255], then attach to AnnData object
    img_path = pyhere.here(str(spaceranger_dir),'tissue_hires_image.png')
    img_arr = np.array(Image.open(img_path), dtype = np.float32) / 256
    
    #   Store image and scalefactors in AnnData as squidpy expects
    adata_vis.uns['spatial'][sample_id] = {
        'scalefactors': json_data,
        'images' : { 'hires' : img_arr }
    }

#-------------------------------------------------------------------------------
#   Attach spatialCoords to spatial AnnData
#-------------------------------------------------------------------------------

#   Correct how spatialCoords are stored. Currently, they are a pandas
#   DataFrame, with the columns potentially in the wrong order (depending on the
#   version of SpatialExperiment used in R). We need them as a numpy array.
adata_vis.obsm['spatial'] = np.array(adata_vis.obsm['spatial'][spatial_coords_names])

#-------------------------------------------------------------------------------
#   Replace special characters in some layer groups
#-------------------------------------------------------------------------------

if cell_group == "layer":
    adata_ref.obs[cell_type_var] = pd.Series(
        [x.replace('/', '_') for x in adata_ref.obs[cell_type_var]],
        dtype = 'category', index = adata_ref.obs_names
    )

#-------------------------------------------------------------------------------
#   Save AnnDatas
#-------------------------------------------------------------------------------

if cell_group == 'broad':
    adata_vis.write_h5ad(
        os.path.join(os.path.dirname(processed_dir), 'adata_vis_orig.h5ad')
    )

adata_ref.write_h5ad(
    os.path.join(processed_dir, 'adata_ref_orig.h5ad')
)