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
sce$brnum<-factor(tab$Brain[match(sce$Sample,tab[,1])])
sce$round<-factor(tab$round[match(sce$Sample,tab[,1])])
sce$sort<-factor(tab$PI.NeuN[match(sce$Sample,tab[,1])])
##rm problematic 17c-scp nuclei
assay(sce,'binomial_deviance_residuals')<-NULL
sce$discard<-ifelse(sce$Sample %in% "17c-scp" & sce$detected < 4000,T,F)
sce<-sce[,sce$discard==F]

##feature selection using deviance
##let's do a poisson model
set.seed(800)
sce <- devianceFeatureSelection(sce,
                                assay = "counts", fam = "poisson",
                                sorted = T,batch=sce$brnum)
#sce <- save(sce,file = here("snRNAseq_hpc","processed-data", "sce", "sce_deviance.rda"))

# This temp file just used for getting batch-corrected components (drops a variety of entries)
#rownames(sce) <- uniquifyFeatureNames(rowData(sce)$gene_id,
#                                      rowData(sce)$gene_name)
#}
#pdf(here("snRNAseq_hpc","plots", "build_sce", "binomial_deviance.pdf"))
plot(sort(rowData(sce)$poisson_deviance, decreasing = T),
     type = "l", xlab = "ranked genes",
     ylim=c(0,1000000),
     ylab = "poisson deviance", main = "Feature Selection with Deviance"
)
abline(v = 2000, lty = 2, col = "red")
abline(v = 2500, lty = 2, col = "pink")
abline(v = 3000, lty = 2, col = "blue")
abline(v = 4000, lty = 2, col = "green")
abline(v = 5000, lty = 2, col = "black")
#dev.off()
#

hdg<-rownames(counts(sce))[1:3000]
res<-sce[rownames(counts(sce)) %in% hdg,]
set.seed(913)
# Initialize an empty list to store the results.
res_list <- list()

# Loop over the subsets defined by brnum.
for (i in seq_along(splitit(res$brnum))) {
    # Subset the sce object.
    res_temp <- res[, splitit(res$brnum)[[i]]]

    message("running nullResiduals - ", Sys.time())

    # Apply the nullResiduals function.
    res_temp <- nullResiduals(res_temp,
                              fam = "poisson",
                              type = "pearson"
    )

    # Add the result to the list.
    res_list[[i]] <- res_temp
}

# Combine the results back into a single sce object.
res_combined <- do.call(cbind, res_list)
idx <- match(colnames(res), colnames(res_combined))

# Use the index to order the columns of 'res_combined'
res_combined <- res_combined[,idx]

# Check that the column names are the same as in the original res object.
all(colnames(res) == colnames(res_combined))

res<-res_combined

##Running GLMPCA
set.seed(101)
res<-GLMPCA(res, 100, assay="counts",minibatch='memoized')
fit<-metadata(sce2)$glmpca
# # Initialize an empty list to store results for each batch
# residuals_list <- list()
#
# # Specify the fam and type variables explicitly
# fam = "binomial"
# type = "deviance"
#
# # Run the nullResiduals() function for each batch
# for(b in levels(res$brnum)){
#     idx <- which(res$brnum == b)
#     temp_res <- nullResiduals(res[, idx], fam = fam, type = type, assay='counts')
#     residuals_list[[b]] <- temp_res



set.seed(915)
message("running PCA - ", Sys.time())
res <- scater::runPCA(res,
                      exprs_values='poisson_pearson_residuals',
                      ncomponents=100,
                      BSPARAM=IrlbaParam())


## Why multi match as well?
# sce <- multiBatchNorm(sce, batch=sce$sample_short)

#### Reduce Dims####

#message("running TSNE - ", Sys.time())
#sce <- runTSNE(sce, dimred = "PCA")
reducedDim(sce,'PCA')<-reducedDim(res,'PCA')
set.seed(12)
message("running UMAP - ", Sys.time())
sce <- runUMAP(sce, dimred = "PCA",n_dimred=100)
reducedDim(sce,'UMAP2')<-umap

pdf(here("snRNAseq_hpc","plots","UMAP_uncorrected_by_sample.pdf"))
plotUMAP(sce,colour_by='Sample')
dev.off()

message("Saving Data - ", Sys.time())
save(sce, file = here("snRNAseq_hpc","processed-data", "sce", "sce.rda"))
# sgejobs::job_single('normalize_step1_glm', create_shell = TRUE, queue= 'bluejay', memory = '75G', command = "Rscript normalize_step1_glm.R")


# Run harmony
message("running MNN - ", Sys.time())
set.seed(1788)
 #sce <- harmony::RunHarmony(sce, group.by.vars = "round",
  #                          verbose = TRUE,max.iter.harmony=30)

mnn<-batchelor::reducedMNN(reducedDim(sce,'PCA'),batch=sce$round,k=50)


#### TSNE & UMAP ####
#set.seed(602)
#message("running TSNE - ", Sys.time())
#sce <- runTSNE(sce, dimred = "HARMONY")
reducedDim(sce,'MNN')<-mnn$corrected
set.seed(125)
message("running UMAP - ", Sys.time())
sce <- runUMAP(sce, dimred = "MNN")

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
