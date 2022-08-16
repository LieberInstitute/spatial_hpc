setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")

library('SingleCellExperiment')
library('here')
library('jaffelab')
library('scater')
library('scran')
library('pheatmap')
library('readxl')
library('Polychrome')
library('cluster')
library('limma')
library('sessioninfo')
library('ggplot2')
library('ggrepel')

## load spe data
load(file = here::here("processed-data", "08_pseudobulk", "mbkmeans", "spe_pseudo_captureArea_mbkmeans17.Rdata"))

# boxplots of spots per cluster
pdf(file = here::here("plots", "08_pseudobulk", "mbkmeans", "ncells_per_captureArea_mbkmeans17.pdf"))
boxplot(ncells ~ spe_pseudo$mbkmeans, data = colData(spe_pseudo))
dev.off()


## Extract the data
mat <- assays(spe_pseudo)$logcounts

# make mat_formula
# var_oi = paste0("mbkmeans_harmony_",k)
var_oi <- "mbkmeans"
covars <- c("age", "sex")
mat_formula <- eval(str2expression(paste("~", "0", "+", var_oi, "+", paste(covars, collapse = " + "))))

# make sure everything is  a factor
colData(spe_pseudo)[[var_oi]] <- as.factor(colData(spe_pseudo)[[var_oi]])
colData(spe_pseudo)$age <- as.numeric(colData(spe_pseudo)$age)
colData(spe_pseudo)$sex <- as.factor(colData(spe_pseudo)$sex)
colData(spe_pseudo)$brnum <- as.factor(colData(spe_pseudo)$brnum)

## Compute correlation
## Adapted from https://github.com/LieberInstitute/Visium_IF_AD/blob/7973fcebb7c4b17cc3e23be2c31ac324d1cc099b/code/10_spatial_registration/01_spatial_registration.R#L134-L150
mod <- model.matrix(mat_formula, data = colData(spe_pseudo))
message(Sys.time(), " running duplicateCorrelation()")
corfit <- duplicateCorrelation(mat, mod, block = spe_pseudo$sample_id)
message("Detected correlation: ", corfit$consensus.correlation)

######### ENRICHMENT t-stats ######################
## Adapted from https://github.com/LieberInstitute/HumanPilot/blob/7049cd42925e00b187c0866f93409196dbcdd526/Analysis/Layer_Guesses/layer_specificity.R#L1423-L1443
# cluster_idx <- splitit(spe_pseudo$mbkmeans_harmony_9) #split by layers not path_grups
cluster_idx <- splitit(colData(spe_pseudo)[, var_oi])

message(Sys.time(), " running the enrichment model")
eb0_list <- lapply(cluster_idx, function(x) {
  res <- rep(0, ncol(spe_pseudo))
  res[x] <- 1
  res_formula <- paste("~", "res", "+", paste(covars, collapse = " + "))
  m <- with(
    colData(spe_pseudo),
    model.matrix(eval(str2expression(res_formula)))
  )
  eBayes(
    lmFit(
      mat,
      design = m,
      block = spe_pseudo$sample_id,
      correlation = corfit$consensus.correlation
    )
  )
})

res = eb0_list[[15]]
# extract p-values
p_vals <- res$p.value[, 2]

# calculate adjusted p-values (FDRs)
fdrs <- p.adjust(p_vals, method = "fdr")

# store corresponding gene names for convenience later
stopifnot(length(fdrs) == length(rowData(spe_pseudo)$gene_id))
stopifnot(all(names(fdrs) == rowData(spe_pseudo)$gene_id))

fdrs_gene_ids <- rowData(spe_pseudo)$gene_id
fdrs_gene_names <- rowData(spe_pseudo)$gene_name

names(fdrs) <- fdrs_gene_names

# ----------------------------------------------
# volcano plot: standard significance thresholds
# ----------------------------------------------

logfc <- res$coefficients[,2]

stopifnot(length(fdrs) == length(logfc))
stopifnot(all(fdrs_gene_ids == names(logfc)))

names(logfc) <- names(fdrs)

# identify significant genes (low FDR and high logFC)
thresh_fdr <- 0.05
thresh_logfc <- log2(2)
sig <- (fdrs < thresh_fdr) & (abs(logfc) > thresh_logfc)

# number of significant genes
table(sig)

df <- data.frame(
  gene = names(fdrs), 
  FDR = fdrs, 
  logFC = logfc, 
  sig = sig
)

pal <- c("black", "red")


# volcano plot without labels
pdf(file = here::here("plots", "08_pseudobulk","mbkmeans", "pseudobulk_captureArea_DE_volcano_DGfrommbkmeans17.pdf"), width = 4.5, height = 4)
ggplot(df, aes(x = logFC, y = -log10(FDR), color = sig)) + 
  geom_point(size = 0.1) + 
  geom_point(data = df[df$sig, ], size = 0.5) + 
  scale_color_manual(values = pal, guide = "none") + 
  geom_hline(yintercept = -log10(0.05), lty = "dashed", color = "royalblue") + 
  geom_vline(xintercept = -log2(2), lty = "dashed", color = "royalblue") + 
  geom_vline(xintercept = log2(2), lty = "dashed", color = "royalblue") + 
  ggtitle("DG - cluster 15 vs. all other clusters") + 
  theme_bw() + 
  theme(plot.title = element_text(face = "bold"), 
        panel.grid.minor = element_blank())

# -------------------------------------------------------
# volcano plot: highly associated significance thresholds
# -------------------------------------------------------

# summarize results with volcano plot

# highly associated significance thresholds

# identify significant genes (low FDR and high logFC)
thresh_fdr <- 1e-13
thresh_logfc <- log2(3)
highlyassoc <- (fdrs < thresh_fdr) & (abs(logfc) > thresh_logfc)

# number of significant genes
table(highlyassoc)


df <- data.frame(
  gene = names(fdrs), 
  FDR = fdrs, 
  logFC = logfc, 
  highlyassoc = highlyassoc
)

pal <- c("black", "red")


# volcano plot without labels
ggplot(df, aes(x = logFC, y = -log10(FDR), color = highlyassoc)) + 
  geom_point(size = 0.1) + 
  geom_point(data = df[df$highlyassoc, ], size = 0.5) + 
  scale_color_manual(values = pal, guide = "none") + 
  geom_hline(yintercept = -log10(1e-13), lty = "dashed", color = "royalblue") + 
  geom_vline(xintercept = -log2(3), lty = "dashed", color = "royalblue") + 
  geom_vline(xintercept = log2(3), lty = "dashed", color = "royalblue") + 
  ggtitle("DG - cluster 15 vs. all other clusters") + 
  theme_bw() + 
  theme(plot.title = element_text(face = "bold"), 
        panel.grid.minor = element_blank())

# volcano plot with labels
set.seed(123)
ggplot(df, aes(x = logFC, y = -log10(FDR), color = highlyassoc, label = gene)) + 
  geom_point(size = 0.1) + 
  geom_point(data = df[df$highlyassoc, ], size = 0.5) + 
  geom_text_repel(data = df[df$highlyassoc, ], 
                  size = 1.5, nudge_y = 0.1, 
                  force = 0.1, force_pull = 0.1, min.segment.length = 0.1, 
                  max.overlaps = 20) + 
  scale_color_manual(values = pal, guide = "none") + 
  geom_hline(yintercept = -log10(1e-13), lty = "dashed", color = "royalblue") + 
  geom_vline(xintercept = -log2(3), lty = "dashed", color = "royalblue") + 
  geom_vline(xintercept = log2(3), lty = "dashed", color = "royalblue") + 
  ggtitle("DG - cluster 15 vs. all other clusters") + 
  theme_bw() + 
  theme(plot.title = element_text(face = "bold"), 
        panel.grid.minor = element_blank())

dev.off()

