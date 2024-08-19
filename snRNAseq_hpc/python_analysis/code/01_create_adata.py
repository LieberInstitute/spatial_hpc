import scanpy as sc
import anndata as ad
import pandas as pd
import scglue

sampleNum = [1,2]+list(range(10,28))+[32,33]+list(range(36,40))
sampleNames = [x + "c-scp" for x in map(str, sampleNum)]

r0 = ["snRNAseq_hpc/processed-data/01_align/1c-scp", "snRNAseq_hpc/processed-data/01_align/2c-scp"]

rest = list(range(10,28))+[32,33]+list(range(36,40))
rest = [x + "c-scp" for x in map(str, rest)]
rest_parent = "snRNAseq_hpc/code/01_align/"
rest = [rest_parent + x for x in rest]

sampleList = r0+rest

finalList = dict(zip(sampleNames, sampleList))
testList = dict(zip(sampleNames[0:2], r0))

subdir_raw = "/outs/raw_feature_bc_matrix.h5"
subdir_filt = "/outs/filtered_feature_bc_matrix.h5"

adatas = {}

for sample_id, filename in finalList.items():
	path = filename+subdir_filt
	sample_adata = sc.read_10x_h5(path)
	sample_adata.var_names_make_unique()
	sample_adata.obs['barcode'] = list(sample_adata.obs_names)
	sample_adata.obs['key'] = [sample_id + "-" + x for x in list(sample_adata.obs_names)]
	sample_adata.obs.index = sample_adata.obs['key']
	adatas[sample_id] = sample_adata

adata = ad.concat(adatas, label="sample")
check1 = list(adata.obs_names)
len(check1) != len(set(check1)) #pass if False

# add obs data

mdata = pd.read_csv("snRNAseq_hpc/raw-data/sample_info/snRNAseq_U01_HPC_AllRounds_Master_Spreadsheet_04072023.csv")

mdata = mdata[mdata.columns[0:5]]
mdata = mdata.rename(columns={"Sample #":"sample", "Tissue":"tissue", "Brain":"brnum", "Round":"round", "PI/NeuN":"sorted"})
sampleSub = [x + "c_scp" for x in map(str, sampleNum)]
mdata = mdata.loc[mdata["sample"].isin(sampleSub)]
mdata['sample'] = sampleNames


demo = pd.read_csv("raw-data/sample_info/demographicInfo_Geo.csv")

obs_data = pd.merge(mdata, demo, on="brnum", how="left") 
adata.obs = pd.merge(adata.obs, obs_data, on="sample", how="left")
adata.obs.index = adata.obs['key']

# add var data

scglue.data.get_gene_annotation(adata=adata, gtf="/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-GRCh38-2020-A/genes/genes.gtf", gtf_by="gene_name")

# save

adata.write_h5ad("snRNAseq_hpc/python_analysis/processed-data/adata_raw.h5")
adata.obs.to_csv("snRNAseq_hpc/python_analysis/processed-data/filtered-matrix_obs.csv")
