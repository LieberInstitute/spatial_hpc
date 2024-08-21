library(DropletUtils)
library(SingleCellExperiment)
library(dplyr)

sampleNum = c(1:2,10:27,32,33,36:39)
sampleNames = paste0(sampleNum, "c-scp")

sampleList = c("snRNAseq_hpc/processed-data/01_align/1c-scp", "snRNAseq_hpc/processed-data/01_align/2c-scp",
  paste0("snRNAseq_hpc/code/01_align/",sampleNames[3:26]))
subdir = "/outs/filtered_feature_bc_matrix.h5"
paths = paste0(sampleList, subdir)

sce <- read10xCounts(paths, sampleNames, type = "HDF5", col.names = TRUE)
sce$key <- paste(sce$Sample, sce$Barcode, sep="_")

tmp <- read.csv(here("snRNAseq_hpc","raw-data","sample_info",
                       "snRNAseq_U01_HPC_AllRounds_Master_Spreadsheet_04072023.csv"))
tmp = tmp[,1:5]
colnames(tmp) = c("sample","tissue","brnum","round","sorted")
mdata = filter(tmp, sample %in% paste0(sampleNum,"c_scp"))

donorinfo <- read.csv(here("raw-data","sample_info","demographicInfo_Geo.csv"))
mdata = left_join(mdata, donorinfo, by=c("brnum"))
mdata$sample = sampleNames

mdata = left_join(as.data.frame(colData(sce)), mdata, by=c("Sample"="sample"))
identical(sce$key, mdata$key)
colData(sce) <- cbind(colData(sce), mdata[,4:13])

save(sce, file="snRNAseq_hpc/python_analysis/processed-data/sce_filtered-matrix.Rdata")

# limit to adata QC
strict = read.csv("snRNAseq_hpc/python_analysis/processed-data/filtered-matrix_obs_postqc-strict.csv") %>%
  tidyr::separate(br_samp, c("brnum","sample"), sep=" ", remove=FALSE) %>%
  mutate(key=paste(sample, barcode, sep="_"))
length(intersect(sce$key, strict$key))

sce = sce[,sce$key %in% strict$key]
save(sce, file="snRNAseq_hpc/python_analysis/processed-data/sce_postqc-strict.Rdata")
