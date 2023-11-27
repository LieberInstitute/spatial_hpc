library(RcppML)
library(pheatmap)
library(SingleCellExperiment)
load(file=here::here('snRNAseq_hpc','processed-data','sce','sce_final.rda')
###Run NMF
x<-RcppML::nmf(assay(sce,'logcounts'),
    k=120,
    tol = 1e-06,
    maxit = 1000,
    verbose = T,
    L1 = 0.1,
    seed = 1135,
    mask_zeros = FALSE,
    diag = TRUE,
    nonneg = TRUE
)

loads<-t(x$h)
colData(sce)<-cbind(colData(sce),loads)

data<-as.data.frame(sce$fine.type)
colnames(data)<-'cellType'
onehot_cellType <-  dcast(data = data, rownames(data) ~ cellType, length)
rownames(onehot_cellType)<-onehot_cellType[,1]
onehot_cellType[,1]<-as.numeric(onehot_cellType[,1])
onehot_cellType<-onehot_cellType[order(onehot_cellType[,1],decreasing=F),]
onehot_cellType[,1]<-NULL
#pats<-colData(spe)[,c(109:183)]
heat<-cor(onehot_cellType,as.data.frame(loads))
heat<-heat[,-c(1,2,3,15,24)]
pheatmap(heat)

cvnmf<-cross_validate_nmf(
    logcounts(sce),
    ranks=c(5,10,50,100,125,150,200),
    n_replicates = 2,
    tol = 1e-03,
    maxit = 100,
    verbose = 3,
    L1 = 0.1,
    L2 = 0,
    threads = 0,
    test_density = 0.2
)

save(sce,file=here::here(SCE_PATH_HERE))
save(x,file=here::here(NMF_PATH_HERE))



## Step 1: Initialize an empty matrix to store the ratios
#ratios <- matrix(0, nrow = nrow(pats), ncol = ncol(pats), dimnames = dimnames(pats))
#
## Step 2: Iterate over each pattern
#for (pattern in colnames(pats)) {
#    # Find the indices of the columns for all patterns except the current pattern
#    other_patterns <- colnames(pats)[colnames(pats) != pattern]
#
#    # Find the maximum loading for each gene across all other patterns
#    max_loadings <- apply(pats[, other_patterns], 1, max)
#
#    # Calculate the ratios for the current pattern
#    pattern_loadings <- pats[, pattern]
#    ratios[, pattern] <- pattern_loadings / max_loadings
#}

## Step 4: Find the top 100 genes for each pattern based on the ratios
#top_genes <- lapply(colnames(ratios), function(pattern) {
#    # Sort the ratios for the current pattern in descending order
#    sorted_ratios <- sort(ratios[, pattern], decreasing = TRUE)
#
#    # Get the top 100 genes for the current pattern
#    top_100_genes <- names(sorted_ratios)[1:100]
#
#    # Return the top 100 genes
#    top_100_genes
#})
#
## Print the top 100 genes for each pattern
#for (i in seq_along(top_genes)) {
#    pattern <- colnames(pats)[i]
#    cat("Top 100 genes for", pattern, ":", "\n")
#    print(top_genes[[i]])
#    cat("\n")
#}
#
#loads<-colData(spe)[,c(109:183)]
#precast
#
#
#cell_types <- spe$cellType
## Calculate median expression for each gene within each cell type
#median_expression_list <- lapply(unique(cell_types), function(ct) {
#ct_indexes <- which(cell_types == ct)
# ct_expression <- gene_expression[, ct_indexes]
#rowMedians(ct_expression, na.rm = TRUE) })
#
#    > # Combine all medians into a matrix
#    > median_expression <- do.call(cbind, median_expression_list)
#> colnames(median_expression) <- unique(cell_types)
#>
#    > # Calculate pairwise correlation between cell types
#    > correlation_matrix <- cor(median_expression)
#>
#    > # Generate a heatmap
#    > pheatmap(correlation_matrix)
#
#
## install.packages(c("Matrix", "DelayedMatrixStats", "pheatmap"))  # Uncomment to install the packages if they're not installed
#library(Matrix)
#library(DelayedMatrixStats)
#library(pheatmap)
#
## Assuming that the gene expression data is in an assay named "counts" for sce2
#gene_expression_sce2 <- assays(sce2)$counts
#
## Get cell types
#cell_types_sce2 <- colData(sce2)$cellType
#
## Calculate median expression for each gene within each cell type for sce2
#median_expression_list_sce2 <- lapply(unique(cell_types_sce2), function(ct) {
#    ct_indexes <- which(cell_types_sce2 == ct)
#    ct_expression <- gene_expression_sce2[, ct_indexes]
#    rowMedians(ct_expression, na.rm = TRUE)
#})
#
## Similarly, do this for spe object
#gene_expression_spe <- assays(spe)$counts
#cell_types_spe <- colData(spe)$cellType
#
#median_expression_list_spe <- lapply(unique(cell_types_spe), function(ct) {
#    ct_indexes <- which(cell_types_spe == ct)
#    ct_expression <- gene_expression_spe[, ct_indexes]
#    rowMedians(ct_expression, na.rm = TRUE)
#})
#
## Combine all medians into a matrix for sce2 and spe
#median_expression_sce2 <- do.call(cbind, median_expression_list_sce2)
#colnames(median_expression_sce2) <- unique(cell_types_sce2)
#
#median_expression_spe <- do.call(cbind, median_expression_list_spe)
#colnames(median_expression_spe) <- unique(cell_types_spe)
#
## Calculate pairwise correlation between cell types
#correlation_matrix <- cor(cbind(median_expression_sce2, median_expression_spe))
#
## Generate a heatmap
#pheatmap(correlation_matrix)
