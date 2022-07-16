setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")

suppressPackageStartupMessages({
    library(SpatialExperiment)
    library(here)
    library(ggplot2)
    library(ggnewscale)
    library(spatialLIBD)
    library(sessioninfo)
})

# Create directory for QC plots
dir_plots <- here::here("plots", "marker_genes")
dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)

load(file = here::here("processed-data", "05_Batch_correction", "spe_harmony.Rdata"))

# Find marker genes
human_markers <-
    c(
        "SNAP25",
        "SLC17A7",
        "GAD1",
        "GAD2",
        "MBP",
        "MOBP",
        "RELN",
        "AQP4",
        "CCK",
        "HPCAL1",
        "PROX1",
        "NECAB1",
        "MPPED1",
        "SLC17A6",
        "TNNT2",
        "CALB1",
        "GABRQ",
        "THOC3",
        "KCNJ4",
        "SFRP2",
        "APLNR",
        "CD44",
        "TMEM155",
        "SYT13",
        "SLC25A22",
        "SLIT1",
        "SYT4",
        "APOC1",
        "WIF1",
        "MTRNR2L10"
    )

# Locate the marker genes
human_markers_search <- rowData(spe)$gene_search[match(human_markers, rowData(spe)$gene_name)]

# Plot marker genes on tissue
for (i in human_markers_search) {
    vis_grid_gene(
        spe = spe,
        geneid = i,
        pdf = here::here("plots", "marker_genes", paste0(gsub("; ", "_", i), ".pdf")),
        assayname = "logcounts",
        minCount = 0,
        viridis = TRUE,
        alpha = 0.5,
        sample_order = unique(spe$sample_id),
        point_size = 1
    )
}
