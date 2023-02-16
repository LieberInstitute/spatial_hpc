###############################
# spatial_HPC project
# Plot Volcano plots from DE
# Anthony Ramnauth, Jan 05 2022
###############################

setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")

suppressPackageStartupMessages({
    library('SingleCellExperiment')
    library('SpatialExperiment')
    library('here')
    library('jaffelab')
    library('scater')
    library('scran')
    library('readxl')
    library('Polychrome')
    library('cluster')
    library('limma')
    library('sessioninfo')
    library('ggplot2')
    library('ggrepel')
    library('EnhancedVolcano')
})

## load DE data
load(file = here::here("processed-data", "08_pseudobulk", "PRECAST",
    "DE_list_captureArea_adjBrnum.Rdata"))

# load spe_pseudo object
load(file = here::here("processed-data", "08_pseudobulk", "PRECAST",
    "spe_pseudo_captureArea_wo_9-15-NA_Fncells50.Rdata"))

# Set up data frame for PRECAST clusters that will be used for volcano plots

############
## cluster 1
############

res1 = eb0_list[[1]]

# extract p-values
p_vals1 <- res1$p.value[, 2]

# calculate adjusted p-values (FDRs)
fdrs1 <- p.adjust(p_vals1, method = "fdr")

# store corresponding gene names for convenience later
stopifnot(length(fdrs1) == length(rowData(spe_pseudo)$gene_id))
stopifnot(all(names(fdrs1) == rowData(spe_pseudo)$gene_id))

fdrs_gene_ids <- rowData(spe_pseudo)$gene_id
fdrs_gene_names <- rowData(spe_pseudo)$gene_name

names(fdrs1) <- fdrs_gene_names

# extract log fold change
logfc1 <- res1$coefficients[,2]

stopifnot(length(fdrs1) == length(logfc1))
stopifnot(all(fdrs_gene_ids == names(logfc1)))

names(logfc1) <- names(fdrs1)

# identify significant genes (low FDR and high logFC)
thresh_fdr <- 0.1
thresh_logfc <- log2(2)
sig1 <- (fdrs1 < thresh_fdr) & (abs(logfc1) > thresh_logfc)

# number of significant genes
table(sig1)
#FALSE  TRUE
# 11304  1056

df1 <- data.frame(
  gene_name = names(fdrs1),
  logFC = logfc1,
  FDR = fdrs1,
  sig = sig1
)

############
## cluster 2
############

res2 = eb0_list[[2]]

# extract p-values
p_vals2 <- res2$p.value[, 2]

# calculate adjusted p-values (FDRs)
fdrs2 <- p.adjust(p_vals2, method = "fdr")

# store corresponding gene names for convenience later
stopifnot(length(fdrs2) == length(rowData(spe_pseudo)$gene_id))
stopifnot(all(names(fdrs2) == rowData(spe_pseudo)$gene_id))

names(fdrs2) <- fdrs_gene_names

# extract log fold change
logfc2 <- res2$coefficients[,2]

stopifnot(length(fdrs2) == length(logfc2))
stopifnot(all(fdrs_gene_ids == names(logfc2)))

names(logfc2) <- names(fdrs2)

# identify significant genes (low FDR and high logFC)
sig2 <- (fdrs2 < thresh_fdr) & (abs(logfc2) > thresh_logfc)

# number of significant genes
table(sig2)
#FALSE  TRUE
# 11743    617

df2 <- data.frame(
  gene_name = names(fdrs2),
  logFC = logfc2,
  FDR = fdrs2,
  sig = sig2
)

############
## cluster 3
############

res3 = eb0_list[[3]]

# extract p-values
p_vals3 <- res3$p.value[, 2]

# calculate adjusted p-values (FDRs)
fdrs3 <- p.adjust(p_vals3, method = "fdr")

# store corresponding gene names for convenience later
stopifnot(length(fdrs3) == length(rowData(spe_pseudo)$gene_id))
stopifnot(all(names(fdrs3) == rowData(spe_pseudo)$gene_id))

names(fdrs3) <- fdrs_gene_names

# extract log fold change
logfc3 <- res3$coefficients[,2]

stopifnot(length(fdrs3) == length(logfc3))
stopifnot(all(fdrs_gene_ids == names(logfc3)))

names(logfc3) <- names(fdrs3)

# identify significant genes (low FDR and high logFC)
sig3 <- (fdrs3 < thresh_fdr) & (abs(logfc3) > thresh_logfc)

# number of significant genes
table(sig3)
#FALSE  TRUE
# 12341    19

df3 <- data.frame(
  gene_name = names(fdrs3),
  logFC = logfc3,
  FDR = fdrs3,
  sig = sig3
)

############
## cluster 4
############

res4 = eb0_list[[4]]

# extract p-values
p_vals4 <- res4$p.value[, 2]

# calculate adjusted p-values (FDRs)
fdrs4 <- p.adjust(p_vals4, method = "fdr")

# store corresponding gene names for convenience later
stopifnot(length(fdrs4) == length(rowData(spe_pseudo)$gene_id))
stopifnot(all(names(fdrs4) == rowData(spe_pseudo)$gene_id))

names(fdrs4) <- fdrs_gene_names

# extract log fold change
logfc4 <- res4$coefficients[,2]

stopifnot(length(fdrs4) == length(logfc4))
stopifnot(all(fdrs_gene_ids == names(logfc4)))

names(logfc4) <- names(fdrs4)

# identify significant genes (low FDR and high logFC)
sig4 <- (fdrs4 < thresh_fdr) & (abs(logfc4) > thresh_logfc)

# number of significant genes
table(sig4)
#FALSE  TRUE
# 11419   941

df4 <- data.frame(
  gene_name = names(fdrs4),
  logFC = logfc4,
  FDR = fdrs4,
  sig = sig4
)

############
## cluster 5
############

res5 = eb0_list[[5]]

# extract p-values
p_vals5 <- res5$p.value[, 2]

# calculate adjusted p-values (FDRs)
fdrs5 <- p.adjust(p_vals5, method = "fdr")

# store corresponding gene names for convenience later
stopifnot(length(fdrs5) == length(rowData(spe_pseudo)$gene_id))
stopifnot(all(names(fdrs5) == rowData(spe_pseudo)$gene_id))

names(fdrs5) <- fdrs_gene_names

# extract log fold change
logfc5 <- res5$coefficients[,2]

stopifnot(length(fdrs5) == length(logfc5))
stopifnot(all(fdrs_gene_ids == names(logfc5)))

names(logfc5) <- names(fdrs5)

# identify significant genes (low FDR and high logFC)
sig5 <- (fdrs5 < thresh_fdr) & (abs(logfc5) > thresh_logfc)

# number of significant genes
table(sig5)
#FALSE  TRUE
# 10845  1515

df5 <- data.frame(
  gene_name = names(fdrs5),
  logFC = logfc5,
  FDR = fdrs5,
  sig = sig5
)

############
## cluster 6
############

res6 = eb0_list[[6]]

# extract p-values
p_vals6 <- res6$p.value[, 2]

# calculate adjusted p-values (FDRs)
fdrs6 <- p.adjust(p_vals6, method = "fdr")

# store corresponding gene names for convenience later
stopifnot(length(fdrs6) == length(rowData(spe_pseudo)$gene_id))
stopifnot(all(names(fdrs6) == rowData(spe_pseudo)$gene_id))

names(fdrs6) <- fdrs_gene_names

# extract log fold change
logfc6 <- res6$coefficients[,2]

stopifnot(length(fdrs6) == length(logfc6))
stopifnot(all(fdrs_gene_ids == names(logfc6)))

names(logfc6) <- names(fdrs6)

# identify significant genes (low FDR and high logFC)
sig6 <- (fdrs6 < thresh_fdr) & (abs(logfc6) > thresh_logfc)

# number of significant genes
table(sig6)
#FALSE  TRUE
# 10400  1960

df6 <- data.frame(
  gene_name = names(fdrs6),
  logFC = logfc6,
  FDR = fdrs6,
  sig = sig6
)

############
## cluster 7
############

res7 = eb0_list[[7]]

# extract p-values
p_vals7 <- res7$p.value[, 2]

# calculate adjusted p-values (FDRs)
fdrs7 <- p.adjust(p_vals7, method = "fdr")

# store corresponding gene names for convenience later
stopifnot(length(fdrs7) == length(rowData(spe_pseudo)$gene_id))
stopifnot(all(names(fdrs7) == rowData(spe_pseudo)$gene_id))

names(fdrs7) <- fdrs_gene_names

# extract log fold change
logfc7 <- res7$coefficients[,2]

stopifnot(length(fdrs7) == length(logfc7))
stopifnot(all(fdrs_gene_ids == names(logfc7)))

names(logfc7) <- names(fdrs7)

# identify significant genes (low FDR and high logFC)
sig7 <- (fdrs7 < thresh_fdr) & (abs(logfc7) > thresh_logfc)

# number of significant genes
table(sig7)
#FALSE  TRUE
# 10059  2301

df7 <- data.frame(
  gene_name = names(fdrs7),
  logFC = logfc7,
  FDR = fdrs7,
  sig = sig7
)

############
## cluster 8
############

res8 = eb0_list[[8]]

# extract p-values
p_vals8 <- res8$p.value[, 2]

# calculate adjusted p-values (FDRs)
fdrs8 <- p.adjust(p_vals8, method = "fdr")

# store corresponding gene names for convenience later
stopifnot(length(fdrs8) == length(rowData(spe_pseudo)$gene_id))
stopifnot(all(names(fdrs8) == rowData(spe_pseudo)$gene_id))

names(fdrs8) <- fdrs_gene_names

# extract log fold change
logfc8 <- res8$coefficients[,2]

stopifnot(length(fdrs8) == length(logfc8))
stopifnot(all(fdrs_gene_ids == names(logfc8)))

names(logfc8) <- names(fdrs8)

# identify significant genes (low FDR and high logFC)
sig8 <- (fdrs8 < thresh_fdr) & (abs(logfc8) > thresh_logfc)

# number of significant genes
table(sig8)
#FALSE  TRUE
# 10000  2360

df8 <- data.frame(
  gene_name = names(fdrs8),
  logFC = logfc8,
  FDR = fdrs8,
  sig = sig8
)

#############
## cluster 10
#############

res10 = eb0_list[[9]]

# extract p-values
p_vals10 <- res10$p.value[, 2]

# calculate adjusted p-values (FDRs)
fdrs10 <- p.adjust(p_vals10, method = "fdr")

# store corresponding gene names for convenience later
stopifnot(length(fdrs10) == length(rowData(spe_pseudo)$gene_id))
stopifnot(all(names(fdrs10) == rowData(spe_pseudo)$gene_id))

names(fdrs10) <- fdrs_gene_names

# extract log fold change
logfc10 <- res10$coefficients[,2]

stopifnot(length(fdrs10) == length(logfc10))
stopifnot(all(fdrs_gene_ids == names(logfc10)))

names(logfc10) <- names(fdrs10)

# identify significant genes (low FDR and high logFC)
sig10 <- (fdrs10 < thresh_fdr) & (abs(logfc10) > thresh_logfc)

# number of significant genes
table(sig10)
#FALSE  TRUE
# 11718   642

df10 <- data.frame(
  gene_name = names(fdrs10),
  logFC = logfc10,
  FDR = fdrs10,
  sig = sig10
)

#############
## cluster 11
#############

res11 = eb0_list[[10]]

# extract p-values
p_vals11 <- res11$p.value[, 2]

# calculate adjusted p-values (FDRs)
fdrs11 <- p.adjust(p_vals11, method = "fdr")

# store corresponding gene names for convenience later
stopifnot(length(fdrs11) == length(rowData(spe_pseudo)$gene_id))
stopifnot(all(names(fdrs11) == rowData(spe_pseudo)$gene_id))

names(fdrs11) <- fdrs_gene_names

# extract log fold change
logfc11 <- res11$coefficients[,2]

stopifnot(length(fdrs11) == length(logfc11))
stopifnot(all(fdrs_gene_ids == names(logfc11)))

names(logfc11) <- names(fdrs11)

# identify significant genes (low FDR and high logFC)
sig11 <- (fdrs11 < thresh_fdr) & (abs(logfc11) > thresh_logfc)

# number of significant genes
table(sig11)
#FALSE  TRUE
# 9337  3023

df11 <- data.frame(
  gene_name = names(fdrs11),
  logFC = logfc11,
  FDR = fdrs11,
  sig = sig11
)

#############
## cluster 12
#############

res12 = eb0_list[[11]]

# extract p-values
p_vals12 <- res12$p.value[, 2]

# calculate adjusted p-values (FDRs)
fdrs12 <- p.adjust(p_vals12, method = "fdr")

# store corresponding gene names for convenience later
stopifnot(length(fdrs12) == length(rowData(spe_pseudo)$gene_id))
stopifnot(all(names(fdrs12) == rowData(spe_pseudo)$gene_id))

names(fdrs12) <- fdrs_gene_names

# extract log fold change
logfc12 <- res12$coefficients[,2]

stopifnot(length(fdrs12) == length(logfc12))
stopifnot(all(fdrs_gene_ids == names(logfc12)))

names(logfc12) <- names(fdrs12)

# identify significant genes (low FDR and high logFC)
sig12 <- (fdrs12 < thresh_fdr) & (abs(logfc12) > thresh_logfc)

# number of significant genes
table(sig12)
#FALSE  TRUE
# 9863  2497

df12 <- data.frame(
  gene_name = names(fdrs12),
  logFC = logfc12,
  FDR = fdrs12,
  sig = sig12
)

#############
## cluster 13
#############

res13 = eb0_list[[12]]

# extract p-values
p_vals13 <- res13$p.value[, 2]

# calculate adjusted p-values (FDRs)
fdrs13 <- p.adjust(p_vals13, method = "fdr")

# store corresponding gene names for convenience later
stopifnot(length(fdrs13) == length(rowData(spe_pseudo)$gene_id))
stopifnot(all(names(fdrs13) == rowData(spe_pseudo)$gene_id))

names(fdrs13) <- fdrs_gene_names

# extract log fold change
logfc13 <- res13$coefficients[,2]

stopifnot(length(fdrs13) == length(logfc13))
stopifnot(all(fdrs_gene_ids == names(logfc13)))

names(logfc13) <- names(fdrs13)

# identify significant genes (low FDR and high logFC)
sig13 <- (fdrs13 < thresh_fdr) & (abs(logfc13) > thresh_logfc)

# number of significant genes
table(sig13)
#FALSE  TRUE
# 7205  5155

df13 <- data.frame(
  gene_name = names(fdrs13),
  logFC = logfc13,
  FDR = fdrs13,
  sig = sig13
)

#############
## cluster 14
#############

res14 = eb0_list[[13]]

# extract p-values
p_vals14 <- res14$p.value[, 2]

# calculate adjusted p-values (FDRs)
fdrs14 <- p.adjust(p_vals14, method = "fdr")

# store corresponding gene names for convenience later
stopifnot(length(fdrs14) == length(rowData(spe_pseudo)$gene_id))
stopifnot(all(names(fdrs14) == rowData(spe_pseudo)$gene_id))

names(fdrs14) <- fdrs_gene_names

# extract log fold change
logfc14 <- res14$coefficients[,2]

stopifnot(length(fdrs14) == length(logfc14))
stopifnot(all(fdrs_gene_ids == names(logfc14)))

names(logfc14) <- names(fdrs14)

# identify significant genes (low FDR and high logFC)
sig14 <- (fdrs14 < thresh_fdr) & (abs(logfc14) > thresh_logfc)

# number of significant genes
table(sig14)
#FALSE  TRUE
# 11754   606

df14 <- data.frame(
  gene_name = names(fdrs14),
  logFC = logfc14,
  FDR = fdrs14,
  sig = sig14
)

# Volcano plots for all 6 PRECAST clusters

pdf(file = here::here("plots", "08_pseudobulk","PRECAST",
    "PRECAST_pseudoBulk_DE_volcano.pdf"),
    width = 8.5, height = 8)

EnhancedVolcano(df1,
    lab = df1$gene,
    x = 'logFC',
    y = 'FDR',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 FDR",
    legendLabels = c('Not sig.','Log (base 2) FC','FDR',
      'FDR & Log (base 2) FC'),
    title = "PRECAST_HPC",
    subtitle = "Cluster_1 vs. all_others"
    )

EnhancedVolcano(df2,
    lab = df2$gene,
    x = 'logFC',
    y = 'FDR',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 FDR",
    legendLabels = c('Not sig.','Log (base 2) FC','FDR',
      'FDR & Log (base 2) FC'),
    title = "PRECAST_HPC",
    subtitle = "Cluster_2 vs. all_others"
    )

EnhancedVolcano(df3,
    lab = df3$gene,
    x = 'logFC',
    y = 'FDR',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 FDR",
    legendLabels = c('Not sig.','Log (base 2) FC','FDR',
      'FDR & Log (base 2) FC'),
    title = "PRECAST_HPC",
    subtitle = "Cluster_3 vs. all_others"
    )

EnhancedVolcano(df4,
    lab = df4$gene,
    x = 'logFC',
    y = 'FDR',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 FDR",
    legendLabels = c('Not sig.','Log (base 2) FC','FDR',
      'FDR & Log (base 2) FC'),
    title = "PRECAST_HPC",
    subtitle = "Cluster_4 vs. all_others"
    )

EnhancedVolcano(df5,
    lab = df5$gene,
    x = 'logFC',
    y = 'FDR',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 FDR",
    legendLabels = c('Not sig.','Log (base 2) FC','FDR',
      'FDR & Log (base 2) FC'),
    title = "PRECAST_HPC",
    subtitle = "Cluster_5 vs. all_others"
    )

EnhancedVolcano(df6,
    lab = df6$gene,
    x = 'logFC',
    y = 'FDR',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 FDR",
    legendLabels = c('Not sig.','Log (base 2) FC','FDR',
      'FDR & Log (base 2) FC'),
    title = "PRECAST_HPC",
    subtitle = "Cluster_6 vs. all_others"
    )

EnhancedVolcano(df7,
    lab = df7$gene,
    x = 'logFC',
    y = 'FDR',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 FDR",
    legendLabels = c('Not sig.','Log (base 2) FC','FDR',
      'FDR & Log (base 2) FC'),
    title = "PRECAST_HPC",
    subtitle = "Cluster_7 vs. all_others"
    )

EnhancedVolcano(df8,
    lab = df8$gene,
    x = 'logFC',
    y = 'FDR',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 FDR",
    legendLabels = c('Not sig.','Log (base 2) FC','FDR',
      'FDR & Log (base 2) FC'),
    title = "PRECAST_HPC",
    subtitle = "Cluster_8 vs. all_others"
    )

EnhancedVolcano(df10,
    lab = df10$gene,
    x = 'logFC',
    y = 'FDR',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 FDR",
    legendLabels = c('Not sig.','Log (base 2) FC','FDR',
      'FDR & Log (base 2) FC'),
    title = "PRECAST_HPC",
    subtitle = "Cluster_10 vs. all_others"
    )

EnhancedVolcano(df11,
    lab = df11$gene,
    x = 'logFC',
    y = 'FDR',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 FDR",
    legendLabels = c('Not sig.','Log (base 2) FC','FDR',
      'FDR & Log (base 2) FC'),
    title = "PRECAST_HPC",
    subtitle = "Cluster_11 vs. all_others"
    )

EnhancedVolcano(df12,
    lab = df12$gene,
    x = 'logFC',
    y = 'FDR',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 FDR",
    legendLabels = c('Not sig.','Log (base 2) FC','FDR',
      'FDR & Log (base 2) FC'),
    title = "PRECAST_HPC",
    subtitle = "Cluster_12 vs. all_others"
    )

EnhancedVolcano(df13,
    lab = df13$gene,
    x = 'logFC',
    y = 'FDR',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 FDR",
    legendLabels = c('Not sig.','Log (base 2) FC','FDR',
      'FDR & Log (base 2) FC'),
    title = "PRECAST_HPC",
    subtitle = "Cluster_13 vs. all_others"
    )

EnhancedVolcano(df14,
    lab = df14$gene,
    x = 'logFC',
    y = 'FDR',
    FCcutoff = 1,
    pCutoff = 0.049,
    ylab = "-log10 FDR",
    legendLabels = c('Not sig.','Log (base 2) FC','FDR',
      'FDR & Log (base 2) FC'),
    title = "PRECAST_HPC",
    subtitle = "Cluster_14 vs. all_others"
    )

dev.off()

pdf(file = here::here("plots", "08_pseudobulk","PRECAST",
    "PRECAST_pseudoBulk_DE_cluster_3_& 13_volcano.pdf"),
    width = 8.5, height = 8)

EnhancedVolcano(df3,
    lab = df3$gene,
    x = 'logFC',
    y = 'FDR',
    FCcutoff = 0.01,
    pCutoff = 0.049,
    ylab = "-log10 FDR",
    legendLabels = c('Not sig.','Log (base 2) FC','FDR',
      'FDR & Log (base 2) FC'),
    title = "PRECAST_HPC",
    subtitle = "Cluster_3 vs. all_others"
    )

EnhancedVolcano(df13,
    lab = df13$gene,
    x = 'logFC',
    y = 'FDR',
    FCcutoff = 0.01,
    pCutoff = 0.049,
    ylab = "-log10 FDR",
    legendLabels = c('Not sig.','Log (base 2) FC','FDR',
      'FDR & Log (base 2) FC'),
    title = "PRECAST_HPC",
    subtitle = "Cluster_13 vs. all_others"
    )

dev.off()

