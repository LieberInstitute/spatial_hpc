###############################
# spatial_HPC project
# Plot DE analysis results
# Anthony Ramnauth, Dec 21 2022
###############################

setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")

suppressPackageStartupMessages({
    library(SpatialExperiment)
    library(here)
    library(SingleCellExperiment)
    library(scater)
    library(scran)
    library(RColorBrewer)
    library(Polychrome)
    library(viridis)
    library(dplyr)
    library(ComplexHeatmap)
    library(sessioninfo)
})

## load spe data
load(file = here::here("processed-data", "08_pseudobulk", "PRECAST",
    "spe_pseudo_captureArea_wo_9-15-NA_Fncells50.Rdata"))

# Load modeling results
modeling_results <- readRDS(file = here::here("processed-data", "08_pseudobulk", "PRECAST",
    "modeling_results.rds"))

################################################################
# Make a .csv list of enriched genes for each PRECAST Cluster
################################################################

# Make a data frame summary
clust_1_enriched_summary <- data.frame(
    gene_id = modeling_results$enrichment$ensembl,
    gene_name = modeling_results$enrichment$gene,
    t_stat = modeling_results$enrichment$t_stat_1,
    p_val = modeling_results$enrichment$p_value_1,
    FDR = modeling_results$enrichment$fdr_1
)

clust_1_enriched_summary <- clust_1_enriched_summary %>%
    filter(FDR < 0.05) %>%
    filter(p_val < 0.05)

clust_1_enriched_summary <- arrange(clust_1_enriched_summary, desc(t_stat))

# directory to save lists
dir_outputs <- here("processed-data", "08_pseudobulk", "PRECAST")
fn_out_1 <- file.path(dir_outputs, "Clust_1_enriched_results")

# Export summary as .csv file
write.csv(clust_1_enriched_summary,fn_out_1, row.names = FALSE)

################################################
# Plot top enriched genes per PRECAST Cluster
################################################

# Get mean expression
mat <- assays(spe_pseudo)$logcounts

# filter
gIndex <- rowMeans(mat) > 0.2 # find the genes for which the mean expression is greater than 0.2
mat_filter <- mat[gIndex, ] # subset matrix on just those genes.  want to remove lowly expressed genes.

# Extract the p-values
pvals <- modeling_results$enrichment[, 14:26]
rownames(pvals) <- rownames(mat_filter)

# Extract the t-statistics
t_stat <- modeling_results$enrichment[, 1:13]
rownames(t_stat) <- rownames(mat_filter)

# Extract the FDRs
fdrs <- modeling_results$enrichment[, 27:39]
rownames(fdrs) <- rownames(mat_filter)

### pick top 10 genes per cluster:sample
cluster_specific_indices <- mapply(
    function(t, p, f) {
        oo <- order(t, decreasing = TRUE)[1:10]
    },
    as.data.frame(t_stat),
    as.data.frame(pvals),
    as.data.frame(fdrs)
)
cluster_ind <- unique(as.numeric(cluster_specific_indices))

# Add logcounts from indexed from top genes
exprs_heatmap <- assays(spe_pseudo)[[2]][cluster_ind, ]
rownames(exprs_heatmap) <- rowData(spe_pseudo)$gene_name[cluster_ind]
colnames(exprs_heatmap) <- paste("logcount", 1:281, sep = "")

# Configure column order to match age groups per BayesSpace cluster
Bayes_age_order <- c(
6, 8, 1, 2, 7, 3, 4, 5,
14, 16, 9, 10, 15, 11, 12, 13,
22, 24, 17, 18, 23, 19, 20, 21,
30, 32, 25, 26, 31, 27, 28, 29,
38, 40, 33, 34, 39, 35, 36, 37,
46, 48, 41, 42, 47, 43, 44, 45,
54, 56, 49, 50, 55, 51, 52, 53
    )

# convert to z-scores
scale_rows = function(x){
    m = apply(x, 1, mean, na.rm = T)
    s = apply(x, 1, sd, na.rm = T)
    return((x - m) / s)
}

exprs_heatmap <- scale_rows(exprs_heatmap)

# Add colors for PRECAST clusters
cols <- Polychrome::palette36.colors(13)
names(cols) <- sort(unique(spe_pseudo$PRECAST_k15))

# Plot heatmap of logcounts for clusters and samples
pdf(file = here::here("plots", "08_pseudobulk", "PRECAST", "PRECAST_enrichment_heatmap_all.pdf"),
    width = 12, height = 10)
Heatmap(exprs_heatmap,
    name = "z-score",
    top_annotation = HeatmapAnnotation(PRECAST_cluster = spe_pseudo$PRECAST_k15,
    col = list(PRECAST_cluster = c("1" = "black", "2" = "yellow", "3" = "purple", "4" = "orange",
            "5" = "cyan", "6" = "red", "7" = "tan", "8" = "navyblue", "10" = "green",
        "11" = "pink", "12" = "blue", "13" = "grey", "14" = "brown"))),
    column_title = "Top 10 transcripts from Differential Enrichment of PRECAST Clusters",
    show_column_names = FALSE,
    cluster_columns = FALSE,
    column_split = spe_pseudo$PRECAST_k15,
    row_split = 13,
    row_title = NULL,
    row_names_gp = gpar(fontsize = 7)
    )
dev.off()


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
