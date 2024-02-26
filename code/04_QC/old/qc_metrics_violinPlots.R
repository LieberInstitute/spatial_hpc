setwd('/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/')
suppressPackageStartupMessages(library("SpatialExperiment"))
suppressPackageStartupMessages(library("here"))
suppressPackageStartupMessages(library("sessioninfo"))
suppressPackageStartupMessages(library("scater"))
suppressPackageStartupMessages(library("gridExtra"))

load(here("processed-data", "04_QC", "spe_QC.Rdata"), verbose = TRUE)

#### QC plots brain####
pdf(here("plots", "04_QC", "QC_outliers_brain.pdf"), width = 21)
## Mito rate
plotColData(spe, x = "brnum", y = "subsets_Mito_percent", colour_by = "high_mito_br") +
    ggtitle("Mito Precent by brain") +
    facet_wrap(~ spe$brnum, scales = "free_x", nrow = 1)

## low sum
plotColData(spe, x = "brnum", y = "sum", colour_by = "low_sum_br") +
    scale_y_log10() +
    ggtitle("Total count by brain") +
    facet_wrap(~ spe$brnum, scales = "free_x", nrow = 1)

## low detected
plotColData(spe, x = "brnum", y = "detected", colour_by = "low_detected_br") +
    scale_y_log10() +
    ggtitle("Detected features by brain") +
    facet_wrap(~ spe$brnum, scales = "free_x", nrow = 1)

# detected features Vs Mito rate
plotColData(spe,
    x = "detected", y = "subsets_Mito_percent",
    colour_by = "discard_auto_br", point_size = 2.5, point_alpha = 0.5) +
    ggtitle("Detected features Vs Mito rate by brain")

brains = unique(spe$brnum)

for (i in seq_along(brains)){
    speb <- spe[, which(spe$brnum == brains[i])]
    samples <- unique(speb$sample_id)
    samples
    
    if (length(samples) == 1){
        p1 <- plotColData(speb[, which(speb$sample_id == samples[1])], 
                          x = "detected", y = "subsets_Mito_percent",
                          colour_by = "discard_auto_br", point_size = 2.5, point_alpha = 0.5
                         ) + ggtitle(paste0(samples[1], "_", brains[i]))
    } else if (length(samples) == 2){
        p1 <- plotColData(speb[, which(speb$sample_id == samples[1])],
                          x = "detected", y = "subsets_Mito_percent",
                          colour_by = "discard_auto_br", point_size = 2.5, point_alpha = 0.5
                          ) + ggtitle(paste0(samples[1], "_", brains[i]))
        p2 <- plotColData(speb[, which(speb$sample_id == samples[2])],
                          x = "detected", y = "subsets_Mito_percent",
                          colour_by = "discard_auto_br", point_size = 2.5, point_alpha = 0.5
                          ) + ggtitle(paste0(samples[2], "_", brains[i]))
        grid.arrange(p1, p2, nrow = 2)
    } else if (length(samples) == 3){
        p1 <- plotColData(speb[, which(speb$sample_id == samples[1])],
                          x = "detected", y = "subsets_Mito_percent",
                          colour_by = "discard_auto_br", point_size = 2.5, point_alpha = 0.5
                          ) + ggtitle(paste0(samples[1], "_", brains[i]))
        p2 <- plotColData(speb[, which(speb$sample_id == samples[2])],
                          x = "detected", y = "subsets_Mito_percent",
                          colour_by = "discard_auto_br", point_size = 2.5, point_alpha = 0.5
                          ) + ggtitle(paste0(samples[2], "_", brains[i]))
        p3 <- plotColData(speb[, which(speb$sample_id == samples[3])],
                          x = "detected", y = "subsets_Mito_percent",
                          colour_by = "discard_auto_br", point_size = 2.5, point_alpha = 0.5
                          ) + ggtitle(paste0(samples[3], "_", brains[i]))
        grid.arrange(p1, p2, p3, nrow = 2)
    } else if (length(samples) == 4){
        p1 <- plotColData(speb[, which(speb$sample_id == samples[1])],
                          x = "detected", y = "subsets_Mito_percent",
                          colour_by = "discard_auto_br", point_size = 2.5, point_alpha = 0.5
                          ) + ggtitle(paste0(samples[1], "_", brains[i]))
        p2 <- plotColData(speb[, which(speb$sample_id == samples[2])],
                          x = "detected", y = "subsets_Mito_percent",
                          colour_by = "discard_auto_br", point_size = 2.5, point_alpha = 0.5
                          ) + ggtitle(paste0(samples[2], "_", brains[i]))
        p3 <- plotColData(speb[, which(speb$sample_id == samples[3])],
                          x = "detected", y = "subsets_Mito_percent",
                          colour_by = "discard_auto_br", point_size = 2.5, point_alpha = 0.5
                          ) + ggtitle(paste0(samples[3], "_", brains[i]))
        p4 <- plotColData(speb[, which(speb$sample_id == samples[4])],
                          x = "detected", y = "subsets_Mito_percent",
                          colour_by = "discard_auto_br", point_size = 2.5, point_alpha = 0.5
                          ) + ggtitle(paste0(samples[4], "_", brains[i]))
        grid.arrange(p1, p2, p3, p4, nrow = 2)
    } else if (length(samples) == 5){
        p1 <- plotColData(speb[, which(speb$sample_id == samples[1])],
                          x = "detected", y = "subsets_Mito_percent",
                          colour_by = "discard_auto_br", point_size = 2.5, point_alpha = 0.5
                          ) + ggtitle(paste0(samples[1], "_", brains[i]))
        p2 <- plotColData(speb[, which(speb$sample_id == samples[2])],
                          x = "detected", y = "subsets_Mito_percent",
                          colour_by = "discard_auto_br", point_size = 2.5, point_alpha = 0.5
                          ) + ggtitle(paste0(samples[2], "_", brains[i]))
        p3 <- plotColData(speb[, which(speb$sample_id == samples[3])],
                          x = "detected", y = "subsets_Mito_percent",
                          colour_by = "discard_auto_br", point_size = 2.5, point_alpha = 0.5
                          ) + ggtitle(paste0(samples[3], "_", brains[i]))
        p4 <- plotColData(speb[, which(speb$sample_id == samples[4])],
                          x = "detected", y = "subsets_Mito_percent",
                          colour_by = "discard_auto_br", point_size = 2.5, point_alpha = 0.5
                          ) + ggtitle(paste0(samples[4], "_", brains[i]))
        p5 <- plotColData(speb[, which(speb$sample_id == samples[5])],
                          x = "detected", y = "subsets_Mito_percent",
                          colour_by = "discard_auto_br", point_size = 2.5, point_alpha = 0.5
                          ) + ggtitle(paste0(samples[5], "_", brains[i]))
        
        grid.arrange(p1, p2, p3, p4, p5, nrow = 2)}
}

# Detected features vs total count
plotColData(spe,
    x = "detected", y = "sum",
    colour_by = "discard_auto_br", point_size = 2.5, point_alpha = 0.5) +
    ggtitle("Detected features Vs total count by brain")

for (i in seq_along(brains)){
    speb <- spe[, which(spe$brnum == brains[i])]
    samples <- unique(speb$sample_id)
    samples
    
    if (length(samples) == 1){
        p1 <- plotColData(speb[, which(speb$sample_id == samples[1])], 
                          x = "detected", y = "sum",
                          colour_by = "discard_auto_br", point_size = 2.5, point_alpha = 0.5
                          ) + ggtitle(paste0(samples[1], "_", brains[i]))
    } else if (length(samples) == 2){
        p1 <- plotColData(speb[, which(speb$sample_id == samples[1])],
                          x = "detected", y = "sum",
                          colour_by = "discard_auto_br", point_size = 2.5, point_alpha = 0.5
                          ) + ggtitle(paste0(samples[1], "_", brains[i]))
        p2 <- plotColData(speb[, which(speb$sample_id == samples[2])],
                          x = "detected", y = "sum",
                          colour_by = "discard_auto_br", point_size = 2.5, point_alpha = 0.5
                          ) + ggtitle(paste0(samples[2], "_", brains[i]))
        grid.arrange(p1, p2, nrow = 2)
    } else if (length(samples) == 3){
        p1 <- plotColData(speb[, which(speb$sample_id == samples[1])],
                          x = "detected", y = "sum",
                          colour_by = "discard_auto_br", point_size = 2.5, point_alpha = 0.5
                          ) + ggtitle(paste0(samples[1], "_", brains[i]))
        p2 <- plotColData(speb[, which(speb$sample_id == samples[2])],
                          x = "detected", y = "sum",
                          colour_by = "discard_auto_br", point_size = 2.5, point_alpha = 0.5
                          ) + ggtitle(paste0(samples[2], "_", brains[i]))
        p3 <- plotColData(speb[, which(speb$sample_id == samples[3])],
                          x = "detected", y = "sum",
                          colour_by = "discard_auto_br", point_size = 2.5, point_alpha = 0.5
                          ) + ggtitle(paste0(samples[3], "_", brains[i]))
        grid.arrange(p1, p2, p3, nrow = 2)
    } else if (length(samples) == 4){
        p1 <- plotColData(speb[, which(speb$sample_id == samples[1])],
                          x = "detected", y = "sum",
                          colour_by = "discard_auto_br", point_size = 2.5, point_alpha = 0.5
                          ) + ggtitle(paste0(samples[1], "_", brains[i]))
        p2 <- plotColData(speb[, which(speb$sample_id == samples[2])],
                          x = "detected", y = "sum",
                          colour_by = "discard_auto_br", point_size = 2.5, point_alpha = 0.5
                          ) + ggtitle(paste0(samples[2], "_", brains[i]))
        p3 <- plotColData(speb[, which(speb$sample_id == samples[3])],
                          x = "detected", y = "sum",
                          colour_by = "discard_auto_br", point_size = 2.5, point_alpha = 0.5
                          ) + ggtitle(paste0(samples[3], "_", brains[i]))
        p4 <- plotColData(speb[, which(speb$sample_id == samples[4])],
                          x = "detected", y = "sum",
                          colour_by = "discard_auto_br", point_size = 2.5, point_alpha = 0.5
                          ) + ggtitle(paste0(samples[4], "_", brains[i]))
        grid.arrange(p1, p2, p3, p4, nrow = 2)
    } else if (length(samples) == 5){
        p1 <- plotColData(speb[, which(speb$sample_id == samples[1])],
                          x = "detected", y = "sum",
                          colour_by = "discard_auto_br", point_size = 2.5, point_alpha = 0.5
                          ) + ggtitle(paste0(samples[1], "_", brains[i]))
        p2 <- plotColData(speb[, which(speb$sample_id == samples[2])],
                          x = "detected", y = "sum",
                          colour_by = "discard_auto_br", point_size = 2.5, point_alpha = 0.5
                          ) + ggtitle(paste0(samples[2], "_", brains[i]))
        p3 <- plotColData(speb[, which(speb$sample_id == samples[3])],
                          x = "detected", y = "sum",
                          colour_by = "discard_auto_br", point_size = 2.5, point_alpha = 0.5
                          ) + ggtitle(paste0(samples[3], "_", brains[i]))
        p4 <- plotColData(speb[, which(speb$sample_id == samples[4])],
                          x = "detected", y = "sum",
                          colour_by = "discard_auto_br", point_size = 2.5, point_alpha = 0.5
                          ) + ggtitle(paste0(samples[4], "_", brains[i]))
        p5 <- plotColData(speb[, which(speb$sample_id == samples[5])],
                          x = "detected", y = "sum",
                          colour_by = "discard_auto_br", point_size = 2.5, point_alpha = 0.5
                          ) + ggtitle(paste0(samples[5], "_", brains[i]))
        grid.arrange(p1, p2, p3, p4, p5, nrow = 2)}
}

dev.off()


#### QC plots sample_id ####
pdf(here("plots", "04_QC", "QC_outliers_capture_area.pdf"), width = 21)
## Mito rate
plotColData(spe, x = "sample_id", y = "subsets_Mito_percent", colour_by = "high_mito_id") +
    ggtitle("Mito Precent by capture area") +
    facet_wrap(~ spe$brnum, scales = "free_x", nrow = 1) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

## low sum
plotColData(spe, x = "sample_id", y = "sum", colour_by = "low_sum_id") +
    scale_y_log10() +
    ggtitle("Total count by capture area") +
    facet_wrap(~ spe$brnum, scales = "free_x", nrow = 1) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

## low detected
plotColData(spe, x = "sample_id", y = "detected", colour_by = "low_detected_id") +
    scale_y_log10() +
    ggtitle("Detected features by capture area") +
    facet_wrap(~ spe$brnum, scales = "free_x", nrow = 1) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# detected features Vs Mito rate
plotColData(spe,
    x = "detected", y = "subsets_Mito_percent",
    colour_by = "discard_auto_id", point_size = 2.5, point_alpha = 0.5) +
    ggtitle("Detected features Vs Mito rate by Capture area")

for (i in seq_along(brains)){
    speb <- spe[, which(spe$brnum == brains[i])]
    samples <- unique(speb$sample_id)
    samples
    
    if (length(samples) == 1){
        p1 <- plotColData(speb[, which(speb$sample_id == samples[1])], 
                          x = "detected", y = "subsets_Mito_percent",
                          colour_by = "discard_auto_id", point_size = 2.5, point_alpha = 0.5
        ) + ggtitle(paste0(samples[1], "_", brains[i]))
    } else if (length(samples) == 2){
        p1 <- plotColData(speb[, which(speb$sample_id == samples[1])],
                          x = "detected", y = "subsets_Mito_percent",
                          colour_by = "discard_auto_id", point_size = 2.5, point_alpha = 0.5
        ) + ggtitle(paste0(samples[1], "_", brains[i]))
        p2 <- plotColData(speb[, which(speb$sample_id == samples[2])],
                          x = "detected", y = "subsets_Mito_percent",
                          colour_by = "discard_auto_id", point_size = 2.5, point_alpha = 0.5
        ) + ggtitle(paste0(samples[2], "_", brains[i]))
        grid.arrange(p1, p2, nrow = 2)
    } else if (length(samples) == 3){
        p1 <- plotColData(speb[, which(speb$sample_id == samples[1])],
                          x = "detected", y = "subsets_Mito_percent",
                          colour_by = "discard_auto_id", point_size = 2.5, point_alpha = 0.5
        ) + ggtitle(paste0(samples[1], "_", brains[i]))
        p2 <- plotColData(speb[, which(speb$sample_id == samples[2])],
                          x = "detected", y = "subsets_Mito_percent",
                          colour_by = "discard_auto_id", point_size = 2.5, point_alpha = 0.5
        ) + ggtitle(paste0(samples[2], "_", brains[i]))
        p3 <- plotColData(speb[, which(speb$sample_id == samples[3])],
                          x = "detected", y = "subsets_Mito_percent",
                          colour_by = "discard_auto_id", point_size = 2.5, point_alpha = 0.5
        ) + ggtitle(paste0(samples[3], "_", brains[i]))
        grid.arrange(p1, p2, p3, nrow = 2)
    } else if (length(samples) == 4){
        p1 <- plotColData(speb[, which(speb$sample_id == samples[1])],
                          x = "detected", y = "subsets_Mito_percent",
                          colour_by = "discard_auto_id", point_size = 2.5, point_alpha = 0.5
        ) + ggtitle(paste0(samples[1], "_", brains[i]))
        p2 <- plotColData(speb[, which(speb$sample_id == samples[2])],
                          x = "detected", y = "subsets_Mito_percent",
                          colour_by = "discard_auto_id", point_size = 2.5, point_alpha = 0.5
        ) + ggtitle(paste0(samples[2], "_", brains[i]))
        p3 <- plotColData(speb[, which(speb$sample_id == samples[3])],
                          x = "detected", y = "subsets_Mito_percent",
                          colour_by = "discard_auto_id", point_size = 2.5, point_alpha = 0.5
        ) + ggtitle(paste0(samples[3], "_", brains[i]))
        p4 <- plotColData(speb[, which(speb$sample_id == samples[4])],
                          x = "detected", y = "subsets_Mito_percent",
                          colour_by = "discard_auto_id", point_size = 2.5, point_alpha = 0.5
        ) + ggtitle(paste0(samples[4], "_", brains[i]))
        grid.arrange(p1, p2, p3, p4, nrow = 2)
    } else if (length(samples) == 5){
        p1 <- plotColData(speb[, which(speb$sample_id == samples[1])],
                          x = "detected", y = "subsets_Mito_percent",
                          colour_by = "discard_auto_id", point_size = 2.5, point_alpha = 0.5
        ) + ggtitle(paste0(samples[1], "_", brains[i]))
        p2 <- plotColData(speb[, which(speb$sample_id == samples[2])],
                          x = "detected", y = "subsets_Mito_percent",
                          colour_by = "discard_auto_id", point_size = 2.5, point_alpha = 0.5
        ) + ggtitle(paste0(samples[2], "_", brains[i]))
        p3 <- plotColData(speb[, which(speb$sample_id == samples[3])],
                          x = "detected", y = "subsets_Mito_percent",
                          colour_by = "discard_auto_id", point_size = 2.5, point_alpha = 0.5
        ) + ggtitle(paste0(samples[3], "_", brains[i]))
        p4 <- plotColData(speb[, which(speb$sample_id == samples[4])],
                          x = "detected", y = "subsets_Mito_percent",
                          colour_by = "discard_auto_id", point_size = 2.5, point_alpha = 0.5
        ) + ggtitle(paste0(samples[4], "_", brains[i]))
        p5 <- plotColData(speb[, which(speb$sample_id == samples[5])],
                          x = "detected", y = "subsets_Mito_percent",
                          colour_by = "discard_auto_id", point_size = 2.5, point_alpha = 0.5
        ) + ggtitle(paste0(samples[5], "_", brains[i]))
        
        grid.arrange(p1, p2, p3, p4, p5, nrow = 2)}
}

# Detected features vs total count

plotColData(spe,
    x = "detected", y = "sum",
    colour_by = "discard_auto_id", point_size = 2.5, point_alpha = 0.5) +
    ggtitle("Detected features Vs total count by Capture area")

for (i in seq_along(brains)){
    speb <- spe[, which(spe$brnum == brains[i])]
    samples <- unique(speb$sample_id)
    samples
    
    if (length(samples) == 1){
        p1 <- plotColData(speb[, which(speb$sample_id == samples[1])], 
                          x = "detected", y = "sum",
                          colour_by = "discard_auto_id", point_size = 2.5, point_alpha = 0.5
        ) + ggtitle(paste0(samples[1], "_", brains[i]))
    } else if (length(samples) == 2){
        p1 <- plotColData(speb[, which(speb$sample_id == samples[1])],
                          x = "detected", y = "sum",
                          colour_by = "discard_auto_id", point_size = 2.5, point_alpha = 0.5
        ) + ggtitle(paste0(samples[1], "_", brains[i]))
        p2 <- plotColData(speb[, which(speb$sample_id == samples[2])],
                          x = "detected", y = "sum",
                          colour_by = "discard_auto_id", point_size = 2.5, point_alpha = 0.5
        ) + ggtitle(paste0(samples[2], "_", brains[i]))
        grid.arrange(p1, p2, nrow = 2)
    } else if (length(samples) == 3){
        p1 <- plotColData(speb[, which(speb$sample_id == samples[1])],
                          x = "detected", y = "sum",
                          colour_by = "discard_auto_id", point_size = 2.5, point_alpha = 0.5
        ) + ggtitle(paste0(samples[1], "_", brains[i]))
        p2 <- plotColData(speb[, which(speb$sample_id == samples[2])],
                          x = "detected", y = "sum",
                          colour_by = "discard_auto_id", point_size = 2.5, point_alpha = 0.5
        ) + ggtitle(paste0(samples[2], "_", brains[i]))
        p3 <- plotColData(speb[, which(speb$sample_id == samples[3])],
                          x = "detected", y = "sum",
                          colour_by = "discard_auto_id", point_size = 2.5, point_alpha = 0.5
        ) + ggtitle(paste0(samples[3], "_", brains[i]))
        grid.arrange(p1, p2, p3, nrow = 2)
    } else if (length(samples) == 4){
        p1 <- plotColData(speb[, which(speb$sample_id == samples[1])],
                          x = "detected", y = "sum",
                          colour_by = "discard_auto_id", point_size = 2.5, point_alpha = 0.5
        ) + ggtitle(paste0(samples[1], "_", brains[i]))
        p2 <- plotColData(speb[, which(speb$sample_id == samples[2])],
                          x = "detected", y = "sum",
                          colour_by = "discard_auto_id", point_size = 2.5, point_alpha = 0.5
        ) + ggtitle(paste0(samples[2], "_", brains[i]))
        p3 <- plotColData(speb[, which(speb$sample_id == samples[3])],
                          x = "detected", y = "sum",
                          colour_by = "discard_auto_id", point_size = 2.5, point_alpha = 0.5
        ) + ggtitle(paste0(samples[3], "_", brains[i]))
        p4 <- plotColData(speb[, which(speb$sample_id == samples[4])],
                          x = "detected", y = "sum",
                          colour_by = "discard_auto_id", point_size = 2.5, point_alpha = 0.5
        ) + ggtitle(paste0(samples[4], "_", brains[i]))
        grid.arrange(p1, p2, p3, p4, nrow = 2)
    } else if (length(samples) == 5){
        p1 <- plotColData(speb[, which(speb$sample_id == samples[1])],
                          x = "detected", y = "sum",
                          colour_by = "discard_auto_id", point_size = 2.5, point_alpha = 0.5
        ) + ggtitle(paste0(samples[1], "_", brains[i]))
        p2 <- plotColData(speb[, which(speb$sample_id == samples[2])],
                          x = "detected", y = "sum",
                          colour_by = "discard_auto_id", point_size = 2.5, point_alpha = 0.5
        ) + ggtitle(paste0(samples[2], "_", brains[i]))
        p3 <- plotColData(speb[, which(speb$sample_id == samples[3])],
                          x = "detected", y = "sum",
                          colour_by = "discard_auto_id", point_size = 2.5, point_alpha = 0.5
        ) + ggtitle(paste0(samples[3], "_", brains[i]))
        p4 <- plotColData(speb[, which(speb$sample_id == samples[4])],
                          x = "detected", y = "sum",
                          colour_by = "discard_auto_id", point_size = 2.5, point_alpha = 0.5
        ) + ggtitle(paste0(samples[4], "_", brains[i]))
        p5 <- plotColData(speb[, which(speb$sample_id == samples[5])],
                          x = "detected", y = "sum",
                          colour_by = "discard_auto_id", point_size = 2.5, point_alpha = 0.5
        ) + ggtitle(paste0(samples[5], "_", brains[i]))
        grid.arrange(p1, p2, p3, p4, p5, nrow = 2)}
}


dev.off()
