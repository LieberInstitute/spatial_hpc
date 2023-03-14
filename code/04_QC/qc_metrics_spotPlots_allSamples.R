setwd('/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/')
suppressPackageStartupMessages(library("spatialLIBD"))
suppressPackageStartupMessages(library("ggspavis"))
suppressPackageStartupMessages(library("gridExtra"))
suppressPackageStartupMessages(library("here"))

load(here("processed-data", "04_QC", "spe_QC_allSamples.Rdata"), verbose = TRUE)
spe$discard_auto_br <- spe$low_sum_br | spe$low_detected_br
spe$discard_auto_id <- spe$low_sum_id | spe$low_detected_id

## QC plot of tissue spots discarded by brain with low sum and features
pdf(here("plots", "04_QC", "QC_discard_lowFeaturesCount_by_brain_allSamples.pdf"), width = 21, height = 10)
samples <- unique(colData(spe)[, c("sample_id", "brnum")])
rownames(samples) <- NULL

for (i in 1:36) {
    p <- vis_clus(
        spe = spe,
        sampleid = samples$sample_id[i],
        clustervar = "discard_auto_br",
        colors = c("FALSE" = "grey90", "TRUE" = "red"),
        point_size = 2,
        ... = paste0("_", samples$brnum[i])
    )

    p1 <- plotVisium(spe[, which(spe$sample_id == samples$sample_id[i])], spots = FALSE)

    grid.arrange(p, p1, nrow = 1)
}
dev.off()

pdf(here("plots", "04_QC", "QC_discard_lowFeaturesCount_by_brain_grid_allSamples.pdf"), width = 21, height = 20)
brains = unique(spe$brnum)

for (i in seq_along(brains)){
    speb <- spe[, which(spe$brnum == brains[i])]
    samples <- unique(speb$sample_id)
    samples
    
    if (length(samples) == 1){
        p1 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "discard_auto_br", colors = c("FALSE" = "grey90", "TRUE" = "red"), point_size = 2,... = paste0("_",brains[i]) )
    } else if (length(samples) == 2){
        p1 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "discard_auto_br", colors = c("FALSE" = "grey90", "TRUE" = "red"), point_size = 2,... = paste0("_",brains[i]) )
        p2 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "discard_auto_br", colors = c("FALSE" = "grey90", "TRUE" = "red"), point_size = 2,... = paste0("_",brains[i]) )
        grid.arrange(p1, p2, nrow = 2)
    } else if (length(samples) == 3){
        p1 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "discard_auto_br", colors = c("FALSE" = "grey90", "TRUE" = "red"), point_size = 2,... = paste0("_",brains[i]) )
        p2 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "discard_auto_br", colors = c("FALSE" = "grey90", "TRUE" = "red"), point_size = 2,... = paste0("_",brains[i]) )
        p3 <- vis_clus(spe = speb, sampleid = samples[3], clustervar = "discard_auto_br", colors = c("FALSE" = "grey90", "TRUE" = "red"), point_size = 2,... = paste0("_",brains[i]) )
        grid.arrange(p1, p2, p3, nrow = 2)
    } else if (length(samples) == 4){
        p1 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "discard_auto_br", colors = c("FALSE" = "grey90", "TRUE" = "red"), point_size = 2,... = paste0("_",brains[i]) )
        p2 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "discard_auto_br", colors = c("FALSE" = "grey90", "TRUE" = "red"), point_size = 2,... = paste0("_",brains[i]) )
        p3 <- vis_clus(spe = speb, sampleid = samples[3], clustervar = "discard_auto_br", colors = c("FALSE" = "grey90", "TRUE" = "red"), point_size = 2,... = paste0("_",brains[i]) )
        p4 <- vis_clus(spe = speb, sampleid = samples[4], clustervar = "discard_auto_br", colors = c("FALSE" = "grey90", "TRUE" = "red"), point_size = 2,... = paste0("_",brains[i]) )
        grid.arrange(p1, p2, p3, p4, nrow = 2)
    } else if (length(samples) == 5){
        p1 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "discard_auto_br", colors = c("FALSE" = "grey90", "TRUE" = "red"), point_size = 2,... = paste0("_",brains[i]) )
        p2 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "discard_auto_br", colors = c("FALSE" = "grey90", "TRUE" = "red"), point_size = 2,... = paste0("_",brains[i]) )
        p3 <- vis_clus(spe = speb, sampleid = samples[3], clustervar = "discard_auto_br", colors = c("FALSE" = "grey90", "TRUE" = "red"), point_size = 2,... = paste0("_",brains[i]) )
        p4 <- vis_clus(spe = speb, sampleid = samples[4], clustervar = "discard_auto_br", colors = c("FALSE" = "grey90", "TRUE" = "red"), point_size = 2,... = paste0("_",brains[i]) )
        p5 <- vis_clus(spe = speb, sampleid = samples[5], clustervar = "discard_auto_br", colors = c("FALSE" = "grey90", "TRUE" = "red"), point_size = 2,... = paste0("_",brains[i]) )
        grid.arrange(p1, p2, p3, p4, p5, nrow = 2)}
}

dev.off()

## QC plot of tissue spots discarded by bcapture area with low sum and features
pdf(here("plots", "04_QC", "QC_discard_lowFeaturesCount_by_capture_area_allSamples.pdf"), width = 21, height = 10)
samples <- unique(colData(spe)[, c("sample_id", "brnum")])
rownames(samples) <- NULL

for (i in 1:36) {
    p <- vis_clus(
        spe = spe,
        sampleid = samples$sample_id[i],
        clustervar = "discard_auto_id",
        colors = c("FALSE" = "grey90", "TRUE" = "red"),
        point_size = 2,
        ... = paste0("_", samples$brnum[i])
    )

    p1 <- plotVisium(spe[, which(spe$sample_id == samples$sample_id[i])], spots = FALSE)

    grid.arrange(p, p1, nrow = 1)
}

dev.off()

pdf(here("plots", "04_QC", "QC_discard_lowFeaturesCount_by_capture_area_grid_allSamples.pdf"), width = 21, height = 20)
brains = unique(spe$brnum)
for (i in seq_along(brains)){
    speb <- spe[, which(spe$brnum == brains[i])]
    samples <- unique(speb$sample_id)
    samples
    
    if (length(samples) == 1){
        p1 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "discard_auto_id", colors = c("FALSE" = "grey90", "TRUE" = "red"), point_size = 2,... = paste0("_",brains[i]) )
    } else if (length(samples) == 2){
        p1 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "discard_auto_id", colors = c("FALSE" = "grey90", "TRUE" = "red"), point_size = 2,... = paste0("_",brains[i]) )
        p2 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "discard_auto_id", colors = c("FALSE" = "grey90", "TRUE" = "red"), point_size = 2,... = paste0("_",brains[i]) )
        grid.arrange(p1, p2, nrow = 2)
    } else if (length(samples) == 3){
        p1 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "discard_auto_id", colors = c("FALSE" = "grey90", "TRUE" = "red"), point_size = 2,... = paste0("_",brains[i]) )
        p2 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "discard_auto_id", colors = c("FALSE" = "grey90", "TRUE" = "red"), point_size = 2,... = paste0("_",brains[i]) )
        p3 <- vis_clus(spe = speb, sampleid = samples[3], clustervar = "discard_auto_id", colors = c("FALSE" = "grey90", "TRUE" = "red"), point_size = 2,... = paste0("_",brains[i]) )
        grid.arrange(p1, p2, p3, nrow = 2)
    } else if (length(samples) == 4){
        p1 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "discard_auto_id", colors = c("FALSE" = "grey90", "TRUE" = "red"), point_size = 2,... = paste0("_",brains[i]) )
        p2 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "discard_auto_id", colors = c("FALSE" = "grey90", "TRUE" = "red"), point_size = 2,... = paste0("_",brains[i]) )
        p3 <- vis_clus(spe = speb, sampleid = samples[3], clustervar = "discard_auto_id", colors = c("FALSE" = "grey90", "TRUE" = "red"), point_size = 2,... = paste0("_",brains[i]) )
        p4 <- vis_clus(spe = speb, sampleid = samples[4], clustervar = "discard_auto_id", colors = c("FALSE" = "grey90", "TRUE" = "red"), point_size = 2,... = paste0("_",brains[i]) )
        grid.arrange(p1, p2, p3, p4, nrow = 2)
    } else if (length(samples) == 5){
        p1 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "discard_auto_id", colors = c("FALSE" = "grey90", "TRUE" = "red"), point_size = 2,... = paste0("_",brains[i]) )
        p2 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "discard_auto_id", colors = c("FALSE" = "grey90", "TRUE" = "red"), point_size = 2,... = paste0("_",brains[i]) )
        p3 <- vis_clus(spe = speb, sampleid = samples[3], clustervar = "discard_auto_id", colors = c("FALSE" = "grey90", "TRUE" = "red"), point_size = 2,... = paste0("_",brains[i]) )
        p4 <- vis_clus(spe = speb, sampleid = samples[4], clustervar = "discard_auto_id", colors = c("FALSE" = "grey90", "TRUE" = "red"), point_size = 2,... = paste0("_",brains[i]) )
        p5 <- vis_clus(spe = speb, sampleid = samples[5], clustervar = "discard_auto_id", colors = c("FALSE" = "grey90", "TRUE" = "red"), point_size = 2,... = paste0("_",brains[i]) )
        grid.arrange(p1, p2, p3, p4, p5, nrow = 2)}
}

dev.off()


## QC plot of tissue spots discarded by brain with high mito spots
pdf(here("plots", "04_QC", "QC_discard_highMito_by_brain_allSamples.pdf"), width = 21, height = 20)
for (i in seq_along(brains)){
    speb <- spe[, which(spe$brnum == brains[i])]
    samples <- unique(speb$sample_id)
    samples
    
    if (length(samples) == 1){
        p1 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "high_mito_br", colors = c("FALSE" = "grey90", "TRUE" = "red"), point_size = 2,... = paste0("_",brains[i]) )
    } else if (length(samples) == 2){
        p1 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "high_mito_br", colors = c("FALSE" = "grey90", "TRUE" = "red"), point_size = 2,... = paste0("_",brains[i]) )
        p2 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "high_mito_br", colors = c("FALSE" = "grey90", "TRUE" = "red"), point_size = 2,... = paste0("_",brains[i]) )
        grid.arrange(p1, p2, nrow = 2)
    } else if (length(samples) == 3){
        p1 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "high_mito_br", colors = c("FALSE" = "grey90", "TRUE" = "red"), point_size = 2,... = paste0("_",brains[i]) )
        p2 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "high_mito_br", colors = c("FALSE" = "grey90", "TRUE" = "red"), point_size = 2,... = paste0("_",brains[i]) )
        p3 <- vis_clus(spe = speb, sampleid = samples[3], clustervar = "high_mito_br", colors = c("FALSE" = "grey90", "TRUE" = "red"), point_size = 2,... = paste0("_",brains[i]) )
        grid.arrange(p1, p2, p3, nrow = 2)
    } else if (length(samples) == 4){
        p1 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "high_mito_br", colors = c("FALSE" = "grey90", "TRUE" = "red"), point_size = 2,... = paste0("_",brains[i]) )
        p2 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "high_mito_br", colors = c("FALSE" = "grey90", "TRUE" = "red"), point_size = 2,... = paste0("_",brains[i]) )
        p3 <- vis_clus(spe = speb, sampleid = samples[3], clustervar = "high_mito_br", colors = c("FALSE" = "grey90", "TRUE" = "red"), point_size = 2,... = paste0("_",brains[i]) )
        p4 <- vis_clus(spe = speb, sampleid = samples[4], clustervar = "high_mito_br", colors = c("FALSE" = "grey90", "TRUE" = "red"), point_size = 2,... = paste0("_",brains[i]) )
        grid.arrange(p1, p2, p3, p4, nrow = 2)
    } else if (length(samples) == 5){
        p1 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "high_mito_br", colors = c("FALSE" = "grey90", "TRUE" = "red"), point_size = 2,... = paste0("_",brains[i]) )
        p2 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "high_mito_br", colors = c("FALSE" = "grey90", "TRUE" = "red"), point_size = 2,... = paste0("_",brains[i]) )
        p3 <- vis_clus(spe = speb, sampleid = samples[3], clustervar = "high_mito_br", colors = c("FALSE" = "grey90", "TRUE" = "red"), point_size = 2,... = paste0("_",brains[i]) )
        p4 <- vis_clus(spe = speb, sampleid = samples[4], clustervar = "high_mito_br", colors = c("FALSE" = "grey90", "TRUE" = "red"), point_size = 2,... = paste0("_",brains[i]) )
        p5 <- vis_clus(spe = speb, sampleid = samples[5], clustervar = "high_mito_br", colors = c("FALSE" = "grey90", "TRUE" = "red"), point_size = 2,... = paste0("_",brains[i]) )
        grid.arrange(p1, p2, p3, p4, p5, nrow = 2)}
}

dev.off()

## QC plot of tissue spots discarded by capture with high mito spots
pdf(here("plots", "04_QC", "QC_discard_highMito_by_capture_area_allSamples.pdf"), width = 21, height = 20)
for (i in seq_along(brains)){
    speb <- spe[, which(spe$brnum == brains[i])]
    samples <- unique(speb$sample_id)
    samples
    
    if (length(samples) == 1){
        p1 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "high_mito_id", colors = c("FALSE" = "grey90", "TRUE" = "red"), point_size = 2,... = paste0("_",brains[i]) )
    } else if (length(samples) == 2){
        p1 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "high_mito_id", colors = c("FALSE" = "grey90", "TRUE" = "red"), point_size = 2,... = paste0("_",brains[i]) )
        p2 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "high_mito_id", colors = c("FALSE" = "grey90", "TRUE" = "red"), point_size = 2,... = paste0("_",brains[i]) )
        grid.arrange(p1, p2, nrow = 2)
    } else if (length(samples) == 3){
        p1 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "high_mito_id", colors = c("FALSE" = "grey90", "TRUE" = "red"), point_size = 2,... = paste0("_",brains[i]) )
        p2 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "high_mito_id", colors = c("FALSE" = "grey90", "TRUE" = "red"), point_size = 2,... = paste0("_",brains[i]) )
        p3 <- vis_clus(spe = speb, sampleid = samples[3], clustervar = "high_mito_id", colors = c("FALSE" = "grey90", "TRUE" = "red"), point_size = 2,... = paste0("_",brains[i]) )
        grid.arrange(p1, p2, p3, nrow = 2)
    } else if (length(samples) == 4){
        p1 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "high_mito_id", colors = c("FALSE" = "grey90", "TRUE" = "red"), point_size = 2,... = paste0("_",brains[i]) )
        p2 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "high_mito_id", colors = c("FALSE" = "grey90", "TRUE" = "red"), point_size = 2,... = paste0("_",brains[i]) )
        p3 <- vis_clus(spe = speb, sampleid = samples[3], clustervar = "high_mito_id", colors = c("FALSE" = "grey90", "TRUE" = "red"), point_size = 2,... = paste0("_",brains[i]) )
        p4 <- vis_clus(spe = speb, sampleid = samples[4], clustervar = "high_mito_id", colors = c("FALSE" = "grey90", "TRUE" = "red"), point_size = 2,... = paste0("_",brains[i]) )
        grid.arrange(p1, p2, p3, p4, nrow = 2)
    } else if (length(samples) == 5){
        p1 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "high_mito_id", colors = c("FALSE" = "grey90", "TRUE" = "red"), point_size = 2,... = paste0("_",brains[i]) )
        p2 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "high_mito_id", colors = c("FALSE" = "grey90", "TRUE" = "red"), point_size = 2,... = paste0("_",brains[i]) )
        p3 <- vis_clus(spe = speb, sampleid = samples[3], clustervar = "high_mito_id", colors = c("FALSE" = "grey90", "TRUE" = "red"), point_size = 2,... = paste0("_",brains[i]) )
        p4 <- vis_clus(spe = speb, sampleid = samples[4], clustervar = "high_mito_id", colors = c("FALSE" = "grey90", "TRUE" = "red"), point_size = 2,... = paste0("_",brains[i]) )
        p5 <- vis_clus(spe = speb, sampleid = samples[5], clustervar = "high_mito_id", colors = c("FALSE" = "grey90", "TRUE" = "red"), point_size = 2,... = paste0("_",brains[i]) )
        grid.arrange(p1, p2, p3, p4, p5, nrow = 2)}
}

dev.off()

# to look at all samples in single page

# vis_grid_clus(
#     spe = spe,
#     clustervar = "high_mito_br",
#     pdf = here::here("plots", "04_QC", "QC_discard_mito_by_brain.pdf"),
#     spatial = FALSE,
#     sort_clust = FALSE,
#     colors = c("FALSE" = "grey90", "TRUE" = "red"),
#     point_size = 2
# )
# 
# vis_grid_clus(
#     spe = spe,
#     clustervar = "high_mito_id",
#     pdf = here::here("plots", "04_QC", "QC_discard_mito_by_capture_area.pdf"),
#     spatial = FALSE,
#     sort_clust = FALSE,
#     colors = c("FALSE" = "grey90", "TRUE" = "red"),
#     point_size = 2
# )
