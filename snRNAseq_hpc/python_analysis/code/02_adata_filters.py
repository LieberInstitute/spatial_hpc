import scanpy as sc
import anndata as ad

adata = ad.read_h5ad("snRNAseq_hpc/python_analysis/processed-data/adata_raw.h5")

# mitochondrial genes, "MT-" for human, "Mt-" for mouse
adata.var["mt"] = adata.var_names.str.startswith("MT-")
# ribosomal genes
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
# hemoglobin genes
adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")

sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True
)


adata.obs['br_samp'] = [x+" "+y for x, y in zip(adata.obs['brnum'].values, adata.obs['sample'].values)]
sc.pl.violin(adata, "pct_counts_mt", groupby='br_samp', rotation=90)

adata_filt = adata[adata.obs['round']!=3] #118415
sc.pl.violin(adata_filt, "pct_counts_mt", groupby='br_samp', rotation=90)

adata_filt = adata_filt[adata_filt.obs['pct_counts_mt']<20] #118162
sc.pl.violin(adata_filt, "pct_counts_mt", groupby='br_samp', rotation=90)

adata_filt.obs.to_csv("snRNAseq_hpc/python_analysis/processed-data/filtered-matrix_obs_postqc-loose.csv")

adata_pi = adata_filt[adata_filt.obs['sorted']=="PI"]
sc.pl.scatter(adata_pi, "total_counts", "n_genes_by_counts", color="pct_counts_mt")

adata_pinn = adata_filt[adata_filt.obs['sorted']=="PI+NeuN+"]
sc.pl.scatter(adata_pinn, "total_counts", "n_genes_by_counts", color="pct_counts_mt")


import seaborn as sns
import matplotlib.pyplot as plt
sns.distplot(adata_filt[adata_filt.obs['n_genes_by_counts']<5000].obs['n_genes_by_counts'], kde=False, bins=60)
plt.show()

sns.distplot(adata_pi[adata_pi.obs['n_genes_by_counts']<2000].obs['n_genes_by_counts'], kde=False, bins=60)

# universal threshold

checkmin = adata_pi.obs['n_genes_by_counts']<1000
adata_pi.obs[checkmin].groupby('sample', observed=True).size()
sc.pl.violin(adata_pi[checkmin], "n_genes_by_counts", groupby='br_samp', rotation=90)

#sample_specific threshold #https://www.sc-best-practices.org/preprocessing_visualization/quality_control.html

import numpy as np
from scipy.stats import median_abs_deviation
def is_outlier(adata, metric: str, nmads: int):
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier

adata_filt.obs['genes_outlier'] = is_outlier(adata_filt, 'log1p_n_genes_by_counts', 3)
sc.pl.violin(adata_filt[adata_filt.obs['genes_outlier']==True], "n_genes_by_counts", groupby='br_samp', rotation=90)


sns.distplot(adata_pinn[adata_pinn.obs['n_genes_by_counts']<7000].obs['n_genes_by_counts'], kde=False, bins=60)

tmp = adata_pinn[adata_pinn.obs['sample']!="17c-scp"]
sns.distplot(tmp[tmp.obs['n_genes_by_counts']<5000].obs['n_genes_by_counts'], kde=False, bins=60)
checkmin = tmp.obs['n_genes_by_counts']<2000
sc.pl.violin(tmp[checkmin], "n_genes_by_counts", groupby='br_samp', rotation=90)

checkmin = adata_pinn.obs['n_genes_by_counts']<3000
sc.pl.violin(adata_pinn[checkmin], "n_genes_by_counts", groupby='br_samp', rotation=90)
checkmin = adata_pinn.obs['n_genes_by_counts']<1000
adata_pinn.obs[checkmin].groupby('sample', observed=True).size()

adata_pi_filt = adata_pi[adata_pi.obs['n_genes_by_counts']>750]
sc.pl.scatter(adata_pi_filt, "total_counts", "n_genes_by_counts", color="pct_counts_mt")

adata_pinn_filt = adata_pinn[adata_pinn.obs['n_genes_by_counts']>1000]
sc.pl.scatter(adata_pinn_filt, "total_counts", "n_genes_by_counts", color="pct_counts_mt")

# umi counts
sns.distplot(adata_pi_filt[adata_pi_filt.obs['total_counts']<10000].obs['total_counts'], kde=False, bins=60)

tmp = adata_pi_filt[adata_pi_filt.obs['sample']!="16c-scp"]
sns.distplot(tmp[tmp.obs['total_counts']<1000].obs['total_counts'], kde=False, bins=60)
sns.distplot(tmp[tmp.obs['total_counts']<3500].obs['total_counts'], kde=False, bins=60)

checkmin = adata_pi_filt.obs['total_counts']<5000
sc.pl.violin(adata_pi_filt[checkmin], "total_counts", groupby='br_samp', rotation=90)
checkmin = adata_pi_filt.obs['total_counts']<2500
adata_pi_filt.obs[checkmin].groupby('sample', observed=True).size()

sns.distplot(adata_pinn_filt[adata_pinn_filt.obs['total_counts']<10000].obs['total_counts'], kde=False, bins=60)

tmp = adata_pinn_filt[adata_pinn_filt.obs['sample']!="17c-scp"]
sns.distplot(tmp[tmp.obs['total_counts']<10000].obs['total_counts'], kde=False, bins=60)

checkmin = adata_pinn_filt.obs['total_counts']<5000
sc.pl.violin(adata_pinn_filt[checkmin], "total_counts", groupby='br_samp', rotation=90)
checkmin = adata_pinn_filt.obs['total_counts']<3500
adata_pinn_filt.obs[checkmin].groupby('sample', observed=True).size()


adata_pi_filt = adata_pi_filt[adata_pi_filt.obs['total_counts']>2000]
adata_pinn_filt = adata_pinn_filt[adata_pinn_filt.obs['total_counts']>3500]

adata_filt = ad.concat([adata_pi_filt, adata_pinn_filt], label="sample")

adata_filt.obs.to_csv("snRNAseq_hpc/python_analysis/processed-data/filtered-matrix_obs_postqc-strict.csv")
