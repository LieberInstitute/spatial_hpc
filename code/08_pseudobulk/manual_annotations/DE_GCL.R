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
library('dplyr')

#load(file = here::here("processed-data", "08_pseudobulk", "manual_annotations", "spe_pseudo_captureArea_wo_CP-THAL-CTX_Fncells50_Fol.Rdata"))
#load(file = here::here("processed-data", "08_pseudobulk", "manual_annotations", "DE_eb0_list_captureArea.Rdata"))

load(file = here::here("processed-data", "08_pseudobulk", "manual_annotations", "spe_pseudo_brain_wo_CP-THAL-CTX_Fncells200.Rdata"))
load(file = here::here("processed-data", "08_pseudobulk", "manual_annotations",  "DE_eb0_list_brain.Rdata"))

res = eb0_list$'GCL'
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
  log2FC = logfc, 
  sig = sig
)

pal <- c("black", "red")

# volcano plot without labels
pdf(file = here::here("plots", "08_pseudobulk","manual_annotations", "GCL", "captureArea_volcano.pdf"), width = 4.5, height = 4)
#pdf(file = here::here("plots", "08_pseudobulk","manual_annotations", "GCL", "brain_volcano.pdf"), width = 4.5, height = 4)
ggplot(df, aes(x = log2FC, y = -log10(FDR), color = sig)) + 
  geom_point(size = 0.1) + 
  geom_point(data = df[df$sig, ], size = 0.5) + 
  scale_color_manual(values = pal, guide = "none") + 
  geom_hline(yintercept = -log10(0.05), lty = "dashed", color = "royalblue") + 
  geom_vline(xintercept = -log2(2), lty = "dashed", color = "royalblue") + 
  geom_vline(xintercept = log2(2), lty = "dashed", color = "royalblue") + 
  ggtitle("GCL vs. all other manually annotated regions") + 
  theme_bw() + 
  theme(plot.title = element_text(face = "bold"), 
        panel.grid.minor = element_blank())

# -------------------------------------------------------
# volcano plot: highly associated significance thresholds
# -------------------------------------------------------

# summarize results with volcano plot

# highly associated significance thresholds

# identify significant genes (low FDR and high logFC)
thresh_fdr <- 1e-10 # brain
thresh_logfc <- log2(3) # brain

#thresh_fdr <- 1e-18 # captureArea
#thresh_logfc <- log2(3) # captureArea

highlyassoc <- (fdrs < thresh_fdr) & (abs(logfc) > thresh_logfc)

# number of significant genes
table(highlyassoc)

df <- data.frame(
  gene = names(fdrs), 
  FDR = fdrs, 
  log2FC = logfc, 
  highlyassoc = highlyassoc
)

pal <- c("black", "red")

temp = df[which(df$highlyassoc == TRUE),]
UP = temp[which(temp$log2FC>0),]
UP = UP %>% arrange(FDR)

DOWN = temp[which(temp$log2FC<0),]
DOWN = DOWN %>% arrange(FDR)

# volcano plot without labels
ggplot(df, aes(x = log2FC, y = -log10(FDR), color = highlyassoc)) + 
  geom_point(size = 0.1) + 
  geom_point(data = df[df$highlyassoc, ], size = 0.5) + 
  scale_color_manual(values = pal, guide = "none") + 
  geom_hline(yintercept = -log10(1e-18), lty = "dashed", color = "royalblue") + 
  geom_vline(xintercept = -log2(3), lty = "dashed", color = "royalblue") + 
  geom_vline(xintercept = log2(3), lty = "dashed", color = "royalblue") + 
  ggtitle("GCL vs. all other manually annotated regions") + 
  theme_bw() + 
  theme(plot.title = element_text(face = "bold"), 
        panel.grid.minor = element_blank())

# volcano plot with labels
set.seed(123)
ggplot(df, aes(x = log2FC, y = -log10(FDR), color = highlyassoc, label = gene)) + 
  geom_point(size = 0.1) + 
  geom_point(data = df[df$highlyassoc, ], size = 0.5) + 
  geom_text_repel(data = df[df$highlyassoc, ], 
                  size = 1.5, nudge_y = 0.1, 
                  force = 0.1, force_pull = 0.1, min.segment.length = 0.1, 
                  max.overlaps = 20) + 
  scale_color_manual(values = pal, guide = "none") + 
  geom_hline(yintercept = -log10(1e-18), lty = "dashed", color = "royalblue") + 
  geom_vline(xintercept = -log2(3), lty = "dashed", color = "royalblue") + 
  geom_vline(xintercept = log2(3), lty = "dashed", color = "royalblue") + 
  ggtitle("GCL vs. all other manually annotated regions") + 
  theme_bw() + 
  theme(plot.title = element_text(face = "bold"), 
        panel.grid.minor = element_blank())

dev.off()


#save(UP,DOWN, file = here::here("processed-data", "08_pseudobulk", "manual_annotations", "DEgenes_brain_GCL.Rdata"))
save(UP,DOWN, file = here::here("processed-data", "08_pseudobulk", "manual_annotations", "DEgenes_captureArea_GCL.Rdata"))


