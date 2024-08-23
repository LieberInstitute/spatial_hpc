library("SingleCellExperiment")
library("scran")
library("scry")

set.seed(800)

load("snRNAseq_hpc/python_analysis/processed-data/sce_postqc-strict.Rdata")

##make sure these are factors
sce$brnum<-factor(sce$brnum)
sce$round<-factor(sce$round)
sce$sorted<-factor(sce$sorted)

sce <- devianceFeatureSelection(sce,
                                assay = "counts", fam = "poisson",
                                sorted = T, batch=sce$brnum)

write.csv(as.data.frame(rowData(sce)), "snRNAseq_hpc/python_analysis/processed-data/rowData_strict_bin-deviance.csv", row.names=F)
save(sce, file="snRNAseq_hpc/python_analysis/processed-data/sce_postqc-strict.Rdata")
