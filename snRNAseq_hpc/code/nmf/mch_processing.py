import pandas as pd

import anndata
import scanpy as sc
from ALLCools.mcds import MCDS
from ALLCools.clustering import log_scale
var_dim_er ='geneslop2k'
chrom_to_remove = ['chrX', 'chrY', 'chrM', 'chrL']
excluded_L1_annot = ['FALSE']
mc_type='CHN'
region_to_subregion = {'CTX': ['MOp','SSp','ACA','AI','RSP','AUD','PTLp','VIS'],
                       'HIP': ['CAa','CAp','DGa','DGp'],
                       'RHP': ['ENT'],
                       'OLF': ['MOB'],
                       'PIR': ['PIRa','PIRp'],
                       'STR': ['STR'],
                       'PAL': ['PAL'],
                       'AMY': ['AMY'],
                       'TH': ['THm','THl','THp'],
                       'HY': ['HY'],
                       'MB': ['SC','MRN','VTA','PAG','IC'],
                       'HB': ['P','MY']
                      }
rs2_subregion_dict = {
    'MOp': ['MOp'],
    'SSp': ['SSp'],
    'ACA': ['ACA'],
    'AI': ['AI'],
    'RSP': ['RSP'],
    'AUD': ['AUD'],
    # Entries from here on follow the same pattern as the first six
    'PTLp': ['PTLp'],
    'VIS': ['VIS'],
    'ENT': ['ENT'],
    'CAa': ['CAa'],
    'CAp': ['CAp'],
    'DGa': ['DGa'],
    'DGp': ['DGp'],
    'PIRa': ['PIRa'],
    'PIRp': ['PIRp'],
    'MOB': ['MOB'],
    'PAL': ['PAL'],
    'STR': ['STR'],
    'AMY': ['AMY'],
    'THl': ['THl'], 
    'THm': ['THm'], 
    'THp': ['THp'],
    'HY': ['HY'],
    'SC': ['SC'],
    'MRN': ['MRN'],
    'VTA': ['VTA'],
    'PAG': ['PAG'],
    'IC': ['IC'],
    'P': ['P'],
    'MY': ['MY']
}

group_name = ['HIP', 'RHP', 'PIR', 'CTX']
selected_ER_slice = []

for region in group_name:
    subregions = region_to_subregion.get(region, [])  # Get the subregions or an empty list if not found
    for subregion in subregions:
        selected_ER_slice.extend(rs2_subregion_dict.get(subregion, []))  # Get the values or an empty list if not found

print(selected_ER_slice)
meta = pd.read_csv('cell_48032_RS2_meta_nooutlier.csv', index_col=0, header=0)
meta

selc = meta.index[meta['Source'].isin(selected_ER_slice) & ~meta['PassTargetFilter'].isin(excluded_L1_annot)]
print(len(selc))

mcds = MCDS.open('CEMBA.epiretro.mcds', var_dim=var_dim_er, use_obs=selc)
mcds

er = mcds.get_adata(mc_type='CHN')
er

log_scale(er, with_mean=True)
er.X = -er.X

er.write_h5ad('rs2_mch_matrix.h5ad')
