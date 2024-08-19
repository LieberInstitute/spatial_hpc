library(here)
library(SingleCellExperiment)
library(dplyr)

load(here("snRNAseq_hpc","processed-data","sce","sce_drops_removed.rda"))
obs = read.csv(here("snRNAseq_hpc","python_analysis", "processed-data", "filtered-matrix_obs.csv"))
cmp.df = full_join(as.data.frame(colData(sce)), obs, by=c("Sample"="sample","Barcode"="barcode"), suffix=c("_sce","_adata"))

colSums(is.na(cmp.df))

# nuclei in common = 137617
filter(cmp.df, !is.na(key_sce), !is.na(key_adata)) %>% nrow()
# jacquard coefficient = 0.885
137617/nrow(cmp.df)
