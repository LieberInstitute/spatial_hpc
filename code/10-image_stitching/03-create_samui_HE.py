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
#this_donor = "Br8325"

this_donor = sys.argv[1:]
this_donor = ''.join(this_donor)

samui_dir = Path(here('processed-data', '10-image_stitching', this_donor))
samui_dir.mkdir(parents = True, exist_ok = True)

sample_info_path = Path(here(
    'processed-data', '10-image_stitching', 'sample_info_clean.csv'
))

sample_info = pd.read_csv(sample_info_path, index_col = 0)
sample_info = sample_info.loc[
    (sample_info['Brain'] == this_donor) &
    ~sample_info['XML file name'].isna(),# &
    #sample_info['In analysis'],
    :
]

################################################################################
#   Gather gene-expression data into a DataFrame to later as a feature
################################################################################
#spe_path = here("processed-data", "spot_deconvo", "shared_utilities", "spe.h5ad")
spe_path = here("processed-data", "10-image_stitching", "spe.h5ad")
#   Read in AnnData and subset to this_donor
spe = sc.read(spe_path)
speP = spe[spe.obs['brnum'] == this_donor, :]


#   Convert the sparse gene-expression matrix to pandas DataFrame, with the
#   gene symbols as column names
gene_df = pd.DataFrame(
    speP.X.toarray(),
    index = speP.obs.index,
    columns = speP.var['gene_name']
)
#   Some gene symbols are actually duplicated. Just take the first column in
#   any duplicated cases
gene_df = gene_df.loc[: , ~gene_df.columns.duplicated()].copy()
gene_df.index.name = None

nmf_df = speP.obs.filter(regex='^(n|f)')
nmf_df = speP.obs.filter(regex='^(?!neuron_cell_body).*nmf.*')
sample_df = speP.obs[['sample_id']].copy()
domain_df = speP.obs[['domain']].copy()

################################################################################
#  deconvo results
################################################################################
deconvo_df = pd.read_csv(Path(here("processed-data", "10-image_stitching", "deconvo.csv")),index_col = 0)
deconvo_df = deconvo_df.set_index('key') 

################################################################################
#  remove overlapping spots
################################################################################
tissue_positions_path = Path(here("processed-data", "10-image_stitching", "tissue_positions_"+this_donor+".csv"))

#tissue_positions_path = Path(here("processed-data", "10-image_stitching", "imageJ", "combined_Br8325", "tissue_positions_"+this_donor+".csv"))
tissue_positions = pd.read_csv(tissue_positions_path ,index_col = 0).rename({'pxl_row_in_fullres': 'y', 'pxl_col_in_fullres': 'x'},axis = 1)
tissue_positions.index.name = None
tissue_positions = tissue_positions[['x', 'y']].astype(int)
SPOT_DIAMETER_M = 55e-6
m_per_px =3.6060572219145977e-07

SPOT_DIAMETER_JSON_M = 65e-6
#json_path = os.path.join(
#    sample_info['spaceranger_dir'].iloc[0], 'scalefactors_json.json'
#)
#with open(json_path, 'r') as f:
#    spaceranger_json = json.load(f)
#m_per_px = SPOT_DIAMETER_JSON_M / spaceranger_json['spot_diameter_fullres']

# json_path = "/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/01_spaceranger/spaceranger-all/V11L05-333_B1/outs/spatial/scalefactors_json.json"
# with open(json_path, 'r') as f:
#      spaceranger_json = json.load(f)
# m_per_px = SPOT_DIAMETER_JSON_M / spaceranger_json['spot_diameter_fullres']

 # Build KD tree for nearest neighbor search.
kd = KDTree(tissue_positions[["x", "y"]].values)
 # Query the KDTree for pairs of points within the threshold distance
overlapping_pairs = pd.DataFrame([(tissue_positions.index[x], tissue_positions.index[y]) for x, y in kd.query_pairs(150)])

unique_items = set(item.split('-1')[1] for col in overlapping_pairs.columns for item in overlapping_pairs[col])
tissue_positions_f = tissue_positions.loc[~tissue_positions.index.isin(overlapping_pairs[0])]

#tissue_positions_f = tissue_positions
common_indices = gene_df.index.intersection(tissue_positions_f.index)
tissue_positions_filtered = tissue_positions_f.loc[common_indices]

nmf_df = nmf_df.loc[tissue_positions_filtered.index]
sample_df = sample_df.loc[tissue_positions_filtered.index]
gene_df = gene_df.loc[tissue_positions_filtered.index]
deconvo_df = deconvo_df.loc[tissue_positions_filtered.index]
domain_df = domain_df.loc[tissue_positions_filtered.index]

tissue_positions_arranged = tissue_positions_filtered.reindex(gene_df.index)
 
################################################################################
#   Use the Samui API to create the importable directory for this combined "sample"
################################################################################
default_gene = 'SLC17A7'
assert default_gene in gene_df.columns, "Default gene not in AnnData"

this_sample = Sample(name = samui_dir.name, path = samui_dir)

this_sample.add_coords(tissue_positions_arranged, name = "coords", mPerPx = m_per_px, size = SPOT_DIAMETER_M)

img_fullres = Path(here("processed-data", "10-image_stitching", "combined_"+this_donor+"_initial.tif"))

#img_fullres = Path(here("processed-data", "10-image_stitching", "imageJ", "combined_Br8325", "combined_"+this_donor+".tif"))
this_sample.add_image(tiff = img_fullres,channels = 'rgb',scale = m_per_px)
this_sample.add_csv_feature(sample_df, name = "Capture areas", coordName = "coords", dataType = "categorical")
this_sample.add_csv_feature(domain_df, name = "Domains", coordName = "coords", dataType = "categorical")
this_sample.add_csv_feature(nmf_df, name = "NMF patterns", coordName = "coords", dataType = "quantitative")
this_sample.add_csv_feature(deconvo_df, name = "Deconvolution", coordName = "coords", dataType = "quantitative")
this_sample.add_chunked_feature(gene_df, name = "Genes", coordName = "coords", dataType = "quantitative")

this_sample.set_default_feature(group = "Genes", feature = default_gene)

this_sample.write()

features_of_interest = [{"feature":"SFRP2","group":"Genes"},
                        {"feature":"CLSTN2","group":"Genes"},
                        {"feature":"APOC1","group":"Genes"},
                        {"feature":"PPFIA2","group":"Genes"},
                        {"feature":"TPM2","group":"Genes"},
                        {"feature":"SLC1A3","group":"Genes"},
                        {"feature":"PRKCG","group":"Genes"},
                        {"feature":"SHTN1","group":"Genes"},
                        {"feature":"nmf79","group":"NMF patterns"},
                        {"feature":"nmf13","group":"NMF patterns"},
                        {"feature":"nmf77","group":"NMF patterns"},
                        {"feature":"fine2broad_RCTD_neuron","group":"Deconvolution"},
                        {"feature":"mid2broad_RCTD_neuron","group":"Deconvolution"}]
with open(here(samui_dir,'sample.json'), 'r') as json_file:
    data = json.load(json_file)

# Replace the "importantFeatures" value
data['overlayParams']['importantFeatures'] = features_of_interest

# Write the modified data back to the JSON file
with open(here(samui_dir,'sample.json'), 'w') as json_file:
    json.dump(data, json_file, indent=4)

session_info.show()