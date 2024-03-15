#####################################################
# spatial_HPC & spatial_lifespan_DG projects
# NMF pattern projection onto lifespan dataset
# Anthony Ramnauth, Ded 0792023
#####################################################

suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(SpatialExperiment)
    library(spatialLIBD)
    library(here)
    library(scater)
    library(scran)
    library(CoGAPS)
    library(RcppML)
    library(projectR)
    library(reshape2)
    library(ComplexHeatmap)
    library(circlize)
    library(dplyr)
})

# load NMF pattern
load("nmf_final.rda")

# Load SPE
spe <- readRDS(here::here("processed-data", "Cell_Type_Deconvolution", "harmony_spe_cell2loc.rds"))

# ## Set gene names as row names for nmf patterns
rownames(spe) <- rowData(spe)$gene_name

# First add dominant cell types
cell_props <- as.data.frame(colData(spe)[, c(44:68)])
for (col in 1:ncol(cell_props)){
    colnames(cell_props)[col] <-  sub("meanscell_abundance_w_sf_", "", colnames(cell_props)[col])
}

# Find dominant celltype for each spot
p <- cell_props %>%
 rowwise() %>%
 mutate(row_max = names(.)[which.max(c_across(everything()))])

colData(spe)$dominant_cell_types <- p$row_max

# Project Erik's human snRNAseq NMF patterns onto DG lifespan Visium data

set.seed(1029)
i<-intersect(rownames(spe),rownames(x@w))
loadings<-x@w
loadings<-loadings[rownames(loadings) %in% i,]
spe2<-spe[rownames(spe) %in% i,]
loadings<-loadings[match(rownames(spe2),rownames(loadings)),]
proj<-project(logcounts(spe2),loadings,L1=0)

###scale projected weights
h <- (t(proj) / x@d)

###add to colData
colData(spe)<-cbind(colData(spe),h)

# Check overlap with dominant cell types

data<-as.data.frame(spe$dominant_cell_types)
colnames(data)<-'dominant_cell_types'
onehot_cellType <-  dcast(data = data, rownames(data) ~ dominant_cell_types, length)
rownames(onehot_cellType)<-onehot_cellType[,1]
onehot_cellType[,1]<-as.numeric(onehot_cellType[,1])
onehot_cellType<-onehot_cellType[order(onehot_cellType[,1],decreasing=F),]
onehot_cellType[,1]<-NULL

correlation_celltype_projection <- as.matrix(cor(onehot_cellType, h))

colnames(correlation_celltype_projection) <- colnames(x@w)

# Find NAs
which(is.na(correlation_celltype_projection), arr.ind=TRUE)

#         row col
#Astro_1     1   2
#Astro_2     2   2
#CA1_N       3   2
#CA2_N       4   2
#CA3_N       5   2
#COP         6   2
#Endoth      7   2
#GC          8   2
#InN_LAMP5   9   2
#InN_LHX6   10   2
#InN_MEIS2  11   2
#InN_NR2F2  12   2
#InN_PV     13   2
#InN_SST    14   2
#InN_VIP    15   2
#Macro      16   2
#Microglia  17   2
#Mossy      18   2
#Myeloid    19   2
#Oligo      20   2
#OPC        21   2
#Pericyte   22   2
#SMC        23   2
#T_Cell     24   2
#VLMC       25   2
#Astro_1     1   3
#Astro_2     2   3
#CA1_N       3   3
#CA2_N       4   3
#CA3_N       5   3
#COP         6   3
#Endoth      7   3
#GC          8   3
#InN_LAMP5   9   3
#InN_LHX6   10   3
#InN_MEIS2  11   3
#InN_NR2F2  12   3
#InN_PV     13   3
#InN_SST    14   3
#InN_VIP    15   3
#Macro      16   3
#Microglia  17   3
#Mossy      18   3
#Myeloid    19   3
#Oligo      20   3
#OPC        21   3
#Pericyte   22   3
#SMC        23   3
#T_Cell     24   3
#VLMC       25   3

# All from column 2 & 3 (NMF2 & NMF3) so remove those columns

correlation_celltype_projection2 <- correlation_celltype_projection[,-2]
correlation_celltype_projection3 <- correlation_celltype_projection2[,-2]

# Take the transpose
correlation_celltype_projection3 <- t(correlation_celltype_projection3)

# Reorganize columns to make the other plots from mouse and macaque cell types

list <- c("GC", "Mossy", "CA3_N", "CA2_N", "CA1_N", "InN_SST", "InN_PV", "InN_VIP",
    "InN_LAMP5", "InN_LHX6", "InN_MEIS2", "InN_NR2F2", "Astro_1", "Astro_2",
    "Oligo", "COP", "OPC", "Microglia", "Macro", "T_Cell", "Myeloid", "Pericyte",
    "SMC", "VLMC", "Endoth")

correlation_celltype_projection3 <- correlation_celltype_projection3[,list]

col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

pdf(file = here::here("plots","NMF_snRNA_to_lifespanDG_Visium_heatmap.pdf"),
    width = 8.5, height = 9)

Heatmap(correlation_celltype_projection3,
    name = "corr",
    col = col_fun,
    cluster_columns = FALSE,
    row_names_gp = gpar(fontsize = 6))

dev.off()

# Print spot plots of interesting NMF patterns onto lifespan data

# order spe observations according to age
spe <- spe[, order(spe$age)]

NMF_patterns <- colnames(colData(spe)[70:169])


# Plot neurogenesis marker genes on tissue
for (i in NMF_patterns) {
    vis_grid_gene(
        spe = spe,
        geneid = i,
        pdf = here::here("plots", "NMF_proj_lifespan", paste0(gsub("; ", "_", i), ".pdf")),
        minCount = 0,
        viridis = FALSE,
        point_size = 1.5
    )
}

saveRDS(spe,
    file = here::here("spe_lifespan_DG_NMF_projections.rds"))

