library("SingleCellExperiment")
library("scran")
library("scater")
library("scry")
library("here")
library("sessioninfo")
library("harmony")
library("BiocSingular")


## Load sce
load(here("snRNAseq_hpc","processed-data", "sce", "sce_post_qc.rda"))
dim(sce)

##feature selection using deviance
##let's do a binomial model first
#set.seed(800)
#sce <- devianceFeatureSelection(sce,
#                                assay = "counts", fam = "binomial", sorted = T)
#sce <- save(sce,file = here("snRNAseq_hpc","processed-data", "sce", "sce_deviance.rda"))

# This temp file just used for getting batch-corrected components (drops a variety of entries)
#rownames(sce) <- uniquifyFeatureNames(rowData(sce)$gene_id,
#                                      rowData(sce)$gene_name)
#}
#pdf(here("snRNAseq_hpc","plots", "build_sce", "binomial_deviance.pdf"))
#plot(sort(rowData(sce)$binomial_deviance, decreasing = T),
#     type = "l", xlab = "ranked genes",
#     ylim=c(0,500000),
#     ylab = "binomial deviance", main = "Feature Selection with Deviance"
#)
#abline(v = 2000, lty = 2, col = "red")
#abline(v = 3000, lty = 2, col = "blue")
#abline(v = 4000, lty = 2, col = "green")
#abline(v = 5000, lty = 2, col = "black")
#dev.off()
#
hdg<-rownames(counts(sce))[1:5000]
set.seed(913)
message("running nullResiduals - ", Sys.time())
res<-counts(sce)[rownames(counts(sce)) %in% hdg,]
resids <- nullResiduals(res,
                     fam = "binomial", 
                     type = "deviance"
)

set.seed(915)
message("running PCA - ", Sys.time())
reducedDim(sce,"PCA") <- runPCA(resids, 
                          ncomponents = 100,
                          BSPARAM = BiocSingular::IrlbaParam()
)

## Why multi match as well?
# sce <- multiBatchNorm(sce, batch=sce$sample_short)

#### Reduce Dims####

#message("running TSNE - ", Sys.time())
#sce <- runTSNE(sce, dimred = "PCA")
set.seed(916)
message("running UMAP - ", Sys.time())
sce <- runUMAP(sce, dimred = "PCA")

pdf(here("snRNAseq_hpc","plots","UMAP_uncorrected_by_sample.pdf"))
plotUMAP(sce,colour_by='Sample')
dev.off()

message("Saving Data - ", Sys.time())
save(sce, file = here("snRNAseq_hpc","processed-data", "sce", "sce.rda"))
# sgejobs::job_single('normalize_step1_glm', create_shell = TRUE, queue= 'bluejay', memory = '75G', command = "Rscript normalize_step1_glm.R")

# Run harmony
message("running Harmony - ", Sys.time())
sce <- RunHarmony(sce, group.by.vars = "Sample", verbose = TRUE)


#### TSNE & UMAP ####
#set.seed(602)
#message("running TSNE - ", Sys.time())
#sce <- runTSNE(sce, dimred = "HARMONY")

message("running UMAP - ", Sys.time())
sce <- runUMAP(sce, dimred = "HARMONY")

pdf(here("snRNAseq_hpc","plots","UMAP_corrected_by_sample.pdf"))
plotUMAP(sce,colour_by='Sample')
dev.off()

message("Done UMAP - Saving data...", Sys.time())

save(sce, file = here("snRNAseq_hpc","processed-data", "sce", "sce_harmony.rda"))


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
