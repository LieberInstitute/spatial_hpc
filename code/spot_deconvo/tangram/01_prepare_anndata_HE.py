#   While shared_utilities/02_r2python.R ideally would produce an AnnData totally
#   ready for use with Tangram, in reality there are a few details
#   to adjust.

import os, sys
os.chdir('/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/')

import pyhere
from pathlib import Path

import scanpy as sc
import numpy as np
import pandas as pd
from anndata import AnnData
import skimage
import seaborn as sns
import tangram as tg
from PIL import Image
import json
import pandas as pd

################################################################################
#   Variable definitions
################################################################################

#cell_group = "broad" 
#subtype = "_class"

cell_group = "layer" 
subtype = "_celltype_class1_noHATAGABAAmy"
#-------------------------------------------------------------------------------
#   Paths
#-------------------------------------------------------------------------------

plot_dir = pyhere.here("plots", "spot_deconvo", "tangram", "HE", "2ndRun_newClass", cell_group)
processed_dir = pyhere.here("processed-data", "spot_deconvo", "tangram", "HE", "2ndRun_newClass", cell_group)

#sc_path_in = pyhere.here("processed-data", "spot_deconvo", "shared_utilities", "sce_class.h5ad")
sc_path_in = pyhere.here("processed-data", "spot_deconvo", "shared_utilities", "sce"+ subtype +".h5ad")
#sc_path_in = pyhere.here("processed-data", "spot_deconvo", "shared_utilities", "sce.h5ad")
sp_path_in = pyhere.here("processed-data", "spot_deconvo", "shared_utilities", 'spe.h5ad')
marker_path = pyhere.here("processed-data", "spot_deconvo", "shared_utilities", "markers_" + cell_group + subtype + ".txt")

sc_path_out = pyhere.here(processed_dir, '{}', 'ad_sc.h5ad')
sp_path_out = pyhere.here(processed_dir, '{}', 'ad_sp_orig.h5ad')

#   Directory containing hires image and a JSON containing scale factors and spot size for a given sample.
spaceranger_dirs = pd.read_csv(pyhere.here("code","spot_deconvo","shared_utilities","samples.txt"), sep = '\t', header=None, names = ['SPpath', 'sample_id', 'brain'])
spaceranger_dirs.SPpath = pyhere.here(spaceranger_dirs.SPpath, 'outs', 'spatial')

#-------------------------------------------------------------------------------
#   Dataset-specific variables
#-------------------------------------------------------------------------------

#   Variable name in ad_sp.obs to color by in deconvolution-related plots
cluster_var_plots = 'PRECAST_k18_nnSVG'

#   Variable name in ad_sc.obs representing cell type
cell_type_var = 'cell.class'
#cell_type_var = 'broad.class'
#   Variable name in both ad_sc.var and ad_sp.var containing Ensembl gene ID and
#   variable name in ad_sp.var containing gene symbol
ensembl_id_var = 'gene_id'
gene_symbol_var = 'gene_name'

#   Variable name in ad_sp.obs containing sample ID
sample_id_var = 'sample_id'

#  Genes we want to plot predicted vs. actual expression for
select_genes_names = ['PPFIA2', 'AMPH', 'FNDC1', 'GFRA1', 'KRT17', 'C5orf63', 'GAD2', 'MIF', 'FABP7', 'MAN1A2', 'SFRP2', 'MOBP', 'MAG', 'MTURN', 'PHLDB1', 'ACTA2', 'TTR']

spatial_coords_names = ['pxl_col_in_fullres', 'pxl_row_in_fullres']

################################################################################
#   Preprocessing
################################################################################

#  Load AnnDatas and list of marker genes
print('Loading AnnDatas...')
ad_sp = sc.read_h5ad(sp_path_in)
ad_sp.obs[sample_id_var] = ad_sp.obs[sample_id_var].astype('category')

sample_name = ad_sp.obs[sample_id_var].unique()[int(os.environ['SGE_TASK_ID']) - 1]

#   Different naming conventions are used between sample IDs in ad_sp vs. in
#   file paths for spaceranger files. Grab the spaceranger ID for this sample

spaceranger_dir = spaceranger_dirs[spaceranger_dirs.sample_id == sample_name].SPpath.values[0]

print('Subsetting to just sample {}.'.format(sample_name))
ad_sp = ad_sp[ad_sp.obs[sample_id_var] == sample_name, :]

ad_sc = sc.read_h5ad(sc_path_in)

with open(marker_path, 'r') as f:
    markers = f.read().splitlines()

#  Use Ensembl IDs as gene names for both AnnDatas
ad_sc.var.index = ad_sc.var[ensembl_id_var]
    
#  Note when genes of interest are present in the training set  
select_genes = ad_sp.var[ensembl_id_var][ad_sp.var[gene_symbol_var].isin(select_genes_names)]
assert len(select_genes_names) == len(select_genes)

for i in range(len(select_genes)):
    gene = select_genes[i]
    gene_name = select_genes_names[i]
    
    #  Verify genes of interest were measured in the experiment
    assert gene in ad_sp.var[ensembl_id_var].values
    assert gene in ad_sc.var[ensembl_id_var].values
    
    if gene in markers:
        print('Gene', gene, '(' + gene_name + ') is in the training set.')
    else:
        print('Gene', gene, '(' + gene_name + ') is in the test set.')

tg.pp_adatas(ad_sc, ad_sp, genes=markers)

#   Make sure some variables are categorical, which enables correct coloring of
#   some later plots
ad_sc.obs[cell_type_var] = ad_sc.obs[cell_type_var].astype('category')
ad_sp.obs[cluster_var_plots] = ad_sp.obs[cluster_var_plots].astype('category')

#   Path to JSON from spaceranger including spot size for this sample
json_path = pyhere.here(str(spaceranger_dir), 'scalefactors_json.json')

with open(json_path) as f: 
    json_data = json.load(f)

#   Read in high-res image as numpy array with values in [0, 1] rather than
#   [0, 255], then attach to AnnData object.
img_path = pyhere.here(str(spaceranger_dir),'tissue_hires_image.png')
img_arr = np.array(Image.open(img_path), dtype = np.float32) / 256

#   Store image and scalefactors in AnnData as squidpy expects
ad_sp.uns['spatial'] = {
    sample_name: {
        'scalefactors': json_data,
        'images' : { 'hires' : img_arr }
    }
}

#   Correct how spatialCoords are stored. Currently, they are a pandas
#   DataFrame, with the columns potentially in the wrong order (depending on the
#   version of SpatialExperiment used in R). We need them as a numpy array.
ad_sp.obsm['spatial'] = np.array(
    ad_sp.obsm['spatial'][spatial_coords_names]
)

#-------------------------------------------------------------------------------
#   Save AnnDatas
#-------------------------------------------------------------------------------

#   Ensure output directories exist
Path(os.path.join(plot_dir, sample_name)).mkdir(parents=True, exist_ok=True)
Path(os.path.dirname(str(sc_path_out).format(sample_name))).mkdir(
    parents=True, exist_ok=True
)

#   While the contents of the 'ad_sc' object should be identical regardless of
#   sample, it looks like the order of some variables is random,
#   but alignment and other downstream tasks are dependent on variable
#   ordering. Therefore, while it's a bit wasteful to save many "copies" of
#   'ad_sc' and 'ad_sp' as done here, it simplifies code later by avoiding
#   several order-related complications that would need manual resolution
ad_sc.write_h5ad(str(sc_path_out).format(sample_name))
ad_sp.write_h5ad(str(sp_path_out).format(sample_name))