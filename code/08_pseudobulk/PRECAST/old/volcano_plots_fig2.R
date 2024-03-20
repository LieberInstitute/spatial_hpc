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


thresh_fdr <- 0.05
thresh_logfc <- log2(2)
fdrs_gene_ids <- rowData(spe_pseudo)$gene_id
fdrs_gene_names <- rowData(spe_pseudo)$gene_name

# Set up data frame for PRECAST clusters that will be used for volcano plots
for (i in names(eb0_list)) {
    # Extract results
    res = eb0_list[[i]]

    # Extract p-values
    p_vals <- res$p.value[, 2]

    # Calculate adjusted p-values (FDRs)
    fdrs <- p.adjust(p_vals, method = "fdr")

    # Check and rename fdrs
    stopifnot(length(fdrs) == length(fdrs_gene_ids))
    stopifnot(all(names(fdrs) == fdrs_gene_ids))
    names(fdrs) <- fdrs_gene_names

    # Extract log fold change
    logfc <- res$coefficients[,2]

    # Check and rename logfc
    stopifnot(length(fdrs) == length(logfc))
    stopifnot(all(fdrs_gene_ids == names(logfc)))
    names(logfc) <- names(fdrs)

    # Identify significant genes (low FDR and high logFC)
    sig <- (fdrs < thresh_fdr) & (abs(logfc) > thresh_logfc)

    # Number of significant genes
    print(paste("Cluster", i))
    print(table(sig))

    # Create data frame
    df_list[[i]] <- data.frame(
        gene_name = names(fdrs),
        logFC = logfc,
        FDR = fdrs,
        sig = sig
    )
}

pdf(file = here::here("plots", "08_pseudobulk","PRECAST",
                      "PRECAST_pseudoBulk_DE_volcano_PRECASTk17_nnSVG_rep.pdf"),
    width = 8.5, height = 8)

for(i in seq_along(df_list)) {
    print(EnhancedVolcano(df_list[[i]],
                    lab = df_list[[i]]$gene_name,
                    x = 'logFC',
                    y = 'FDR',
                    FCcutoff = 1,
                    pCutoff = 0.049,
                    ylab = "-log10 FDR",
                    legendLabels = c('Not sig.','Log (base 2) FC','FDR',
                                     'FDR & Log (base 2) FC'),
                    title = "PRECAST_HPC",
                    subtitle = paste0("Cluster_", names(df_list[i]), " vs. all_others")
    ))
pdf(file=here::here('plots','figures','figure_2','neuropil_volcano.pdf'),h=5,w=5)
EnhancedVolcano(reg2,
                     lab = reg2$gene,
                    selectLab = c('MT-ND4','ETNPPL','SLC1A3','MT-ND5','DUSP5','ATP1B2'),
                     x = 'logFC_Neuropil',
                     y = 'fdr_Neuropil',
                     FCcutoff = 1,
                     pCutoff = 0.01,
                     ylab = "-log10 FDR",
                     legendLabels = c('n.s.','log2FC > 1','FDR < 0.01',
                                      'DEGs'),
                     title = 'Neuropil',
                     subtitle = NULL,
                    caption=NULL,
                    xlab='log2 fold change',
                    drawConnectors = T,
                    pointSize=1
    )+theme( # Modify theme to remove grid lines, ticks, tick labels and set text color to black
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(color = "black"),
        axis.ticks = element_blank(),
        axis.text = element_blank()
    )

}

dev.off()

pdf(file=here::here('plots','figures','figure_2','neuron_volcano.pdf'),h=5,w=5)
EnhancedVolcano(reg2,
                lab = reg2$gene,
                selectLab = c('PREPL','CLSTN3','ENO2','KLC1','ATPV1B2',
                              'RPL13A','RPS6','IFITM3','RPS3','RPS9'),
                max.overlaps=100,
                x = 'logFC_Neuron',
                y = 'fdr_Neuron',
                FCcutoff = 1,
                pCutoff = 0.01,
                ylab = "-log10 FDR",
                legendLabels = c('n.s.','log2FC > 1','FDR < 0.01',
                                 'DEGs'),
                title = 'Neuron',
                subtitle = NULL,
                caption=NULL,
                xlab='log2 fold change',
                drawConnectors = T,
                pointSize=1
)+theme( # Modify theme to remove grid lines, ticks, tick labels and set text color to black
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(color = "black"),
    axis.ticks = element_blank(),
    axis.text = element_blank()
)
dev.off()

pdf(file=here::here('plots','figures','figure_2','wm_volcano.pdf'),h=5,w=5)
EnhancedVolcano(reg2,
                lab = reg2$gene,
               selectLab = c('ABCA2','SHTN1','TMEM165','DBNDD2',
                              'PTMS','DKK3'),
                max.overlaps=100,
                x = 'logFC_WM',
                y = 'fdr_WM',
                FCcutoff = 1,
                pCutoff = 0.01,
                ylab = "-log10 FDR",
                legendLabels = c('n.s.','log2FC > 1','FDR < 0.01',
                                 'DEGs'),
                title = 'WM',
                subtitle = NULL,
                caption=NULL,
                xlab='log2 fold change',
                drawConnectors = T,
                pointSize=1
)+theme( # Modify theme to remove grid lines, ticks, tick labels and set text color to black
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(color = "black"),
    axis.ticks = element_blank(),
    axis.text = element_blank()
)
dev.off()


pdf(file=here::here('plots','figures','figure_2','vasc_volcano2.pdf'),h=5,w=5)
EnhancedVolcano(reg2,
                lab = reg2$gene,
                selectLab = c('TPM2','TAGLN','FLNA','MYH11',
                              'KIF5C','ANK3'),
                max.overlaps=100,
                x = 'logFC_Vascular',
                y = 'fdr_Vascular',
                FCcutoff = 1,
                pCutoff = 0.01,
                ylab = "-log10 FDR",
                legendLabels = c('n.s.','log2FC > 1','FDR < 0.01',
                                 'DEGs'),
                title = 'Vascular',
                subtitle = NULL,
                caption=NULL,
                xlab='log2 fold change',
                drawConnectors = T,
                pointSize=1
)+theme( # Modify theme to remove grid lines, ticks, tick labels and set text color to black
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(color = "black"),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    legend.position='right'
)
dev.off()
