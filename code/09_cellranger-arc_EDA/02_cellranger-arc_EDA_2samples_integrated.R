########################################################################
## 1st EDA with scCellRanger-ARC
## Authors. CSC + HT / LIEBER
## Date. March 3rd, 2023
## Reference Link. 
## Use this first: https://stuartlab.org/signac/articles/pbmc_multiomic.html 
## https://stuartlab.org/signac/articles/pbmc_vignette.html
########################################################################

## Reproducibility information
library("sessioninfo")
print('Reproducibility information:')
# Last modification
Sys.time()
#"2023-04-04 12:42:26 EDT"
proc.time()
options(width = 120)
session_info()

# > print('Reproducibility information:')
# [1] "Reproducibility information:"
# > # Last modification
#     > Sys.time()
# [1] "2023-03-27 10:08:00 EDT"
# > #"2023-03-22 14:51:04 EDT"
#     > proc.time()
# user   system  elapsed 
# 188.633   11.724 2061.430 
# > options(width = 120)
# > session_info()
# ─ Session info ──────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.2.3 Patched (2023-03-25 r84072)
# os       CentOS Linux 7 (Core)
# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_US.UTF-8
# tz       US/Eastern
# date     2023-03-27
# pandoc   2.19.2 @ /jhpce/shared/jhpce/core/conda/miniconda3-4.11.0/envs/svnR-4.2.x/bin/pandoc

# load libraries
library(Signac)
# ctype    en_US.UTF-8
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(tidyverse)
library(here)
library(GenomeInfoDb)
library(dplyr)
library(patchwork)
library(scCustomize)          # split by group the VPlots - Seurat complement
set.seed(1234)

# 1) LOAD TWO RNA & ATAC DATA COMBINED (hippocampus samples 42_1 and 42_4)
sc_hippo <- Read10X_h5("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/rafael_rotation/cellranger_rerun/42_1/outs/filtered_feature_bc_matrix.h5")
# raw matrix
#sc_hippo_raw <- Read10X_h5("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/rafael_rotation/cellranger_rerun/42_1/outs/raw_feature_bc_matrix.h5")

sc_hippo2 <- Read10X_h5("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/rafael_rotation/cellranger_rerun/42_4/outs/filtered_feature_bc_matrix.h5")
class(sc_hippo)
# Read10X_h5() read count matrix from 10X CellRanger hdf5 file. This can be used to read both scATAC-seq and scRNA-seq matrices.

######## extract RNA and ATAC data, plus Gene Annotation for hg38 #######
str(sc_hippo)   # read modalities
show(sc_hippo)
# sample 42_1
rna_counts <- sc_hippo$`Gene Expression`
atac_counts <- sc_hippo$Peaks
# sample 42_1 raw data
#rna_counts_raw <- sc_hippo_raw$`Gene Expression`
#atac_counts_raw <- sc_hippo_raw$Peaks
head(rna_counts, n=5)
#head(rna_counts_raw, n=5)
# compare size 
#all.equal.raw(rna_counts, rna_counts_raw). #  Mean relative difference: 271.9585
# sample 42_4
rna_counts2 <- sc_hippo2$`Gene Expression`
atac_counts2 <- sc_hippo2$Peaks

# Create .RData object to load at any time: count matrices sample 42_1 ans 42_4
#save(rna_counts, atac_counts, file = "/fastscratch/myscratch/csoto/counts_sample_42_1/arc_atac_42_1.RData")
#load("/fastscratch/myscratch/csoto/counts_sample_42_1/arc_atac_42_1.RData")
#save(rna_counts2, atac_counts2, file = "/fastscratch/myscratch/csoto/counts_sample_42_4/arc_atac_42_4.RData")
#load("/fastscratch/myscratch/csoto/counts_sample_42_1/arc_atac_42_4.RData")

# In this MTX  each row represents a peak, predicted to represent a region of open chromatin.
fragpath <- "/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/rafael_rotation/cellranger_rerun/42_1/outs/atac_fragments.tsv.gz"
fragpath2 <- "/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/rafael_rotation/cellranger_rerun/42_4/outs/atac_fragments.tsv.gz"
# fragments.tsv is the full list of all unique fragments across all single cells. It contains all fragments associated with each single cell, as opposed to only fragments that map to peaks.

####### Add meta-data to calculate stats to add blacklist ratio and fraction of reads in peaks  ######### 
# For anyone having this issue using cellranger-atac-2.0.0 - meta data is now labeled as singlecell.csv
metadata_42_1 <- read.csv(
    file = "/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/rafael_rotation/cellranger_rerun/42_1/outs/per_barcode_metrics.csv",
    header = TRUE,
    row.names = 1
)
meta_tmp = c('atac_peak_region_fragments','atac_fragments')
meta1 = metadata_42_1[meta_tmp]
metadata_42_4 <- read.csv(
    file = "/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/rafael_rotation/cellranger_rerun/42_4/outs/per_barcode_metrics.csv",
    header = TRUE,
    row.names = 1
)
meta2 = metadata_42_4[meta_tmp]
head(meta2, n=3)

######## Create a Seurat object containing the RNA ########
# Initialize the Seurat object with filtered_feature_bc_matrix (non-normalized data) & meta.data with per barcode metrics attached
hippo <- CreateSeuratObject(
  counts = rna_counts,
  assay = "RNA",
  project = "hippo-42_1",
  meta.data = meta1
)
class(hippo)
# check metadata
head(hippo@meta.data)
hippo2 <- CreateSeuratObject(
    counts = rna_counts2,
    assay = "RNA",
    project = "hippo-42_4",
    meta.data = meta2
)
head(hippo2@meta.data)
show(hippo)
show(hippo2)

# Create .RData object with the Seurat rna-seq assays for samples 42_1 and 42_4
#save(rna_counts, atac_counts, hippo, file = "/fastscratch/myscratch/csoto/counts_sample_42_1/rna_assay_42_1.RData")
load("/fastscratch/myscratch/csoto/counts_sample_42_1/rna_assay_42_1.RData")
#save(rna_counts2, atac_counts2, hippo2, file = "/fastscratch/myscratch/csoto/counts_sample_42_4/rna_assay_42_4.RData")
load("/fastscratch/myscratch/csoto/counts_sample_42_4/rna_assay_42_4.RData")

# Merge two seurat objects
hippo.combined <- merge(hippo, y = hippo2, add.cell.ids = c("Sample_42_1", "Sample_42_4"), project = "HIPPO")
class(hippo.combined)
show(hippo.combined)
head(colnames(hippo.combined))
table(hippo.combined$orig.ident)
head(hippo.combined, n=3)
tail(hippo.combined, n=3)

# Create seurat object with the samples 42_1 and 42_4 merged (only rna counts mtx)
#save(hippo.combined, file = "/fastscratch/myscratch/csoto/merged_samples_42/rna_merged_assay_42.RData")
load("/fastscratch/myscratch/csoto/merged_samples_42/rna_merged_assay_42.RData")

######## QA metrics for the 'Gene Expression' assay  ########
hippo.combined[["percent.mt"]] <- PercentageFeatureSet(hippo.combined, pattern = "^MT-")
head(hippo.combined@meta.data)                 # Access cell-level meta-data / head(hippo[[]])
tail(hippo.combined[["percent.mt"]][])         # Access feature-level meta-data

# Visualize QC metrics as a violin plot
head(hippo2)
VlnPlot(object = hippo.combined, features = c('nCount_RNA','nFeature_RNA','percent.mt'), 
        group.by = "orig.ident")
head(hippo.combined)

# Scatter plot to visualize feature-feature relationships, across the set of single cells.
#  sample() gets rndom samples and Permutations
# hippo.combined$sample_id <- sample(c("Sample_42_1", "Sample_42_4"), size = ncol(hippo.combined), replace = TRUE)
Split_FeatureScatter(seurat_object = hippo.combined, feature1 = "nCount_RNA", feature2="nFeature_RNA", 
                     split.by = "orig.ident")
#Split_FeatureScatter(seurat_object = hippo.combined, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", 
#                    split.by = "sample_id")
Split_FeatureScatter(seurat_object = hippo.combined, feature1 = "nCount_RNA", feature2 = "percent.mt", 
                     split.by = 'orig.ident')

######## ATAC assay: peaks in standard chromosomes were used for analysis of both samples ########
# StringToGRanges(), Convert a genomic coordinate string to a GRanges object
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))      
grange.counts2 <- StringToGRanges(rownames(atac_counts2), sep = c(":", "-"))
length(grange.counts)
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
grange.use2 <- seqnames(grange.counts2) %in% standardChromosomes(grange.counts2)
class(grange.use)       #S4Vector object / booleam
# length(grange.use)
# slotNames(grange.use)
# sub_ranges = window(grange.use,92220,92279)

length(atac_counts)
head(atac_counts)
atac_counts <- atac_counts[as.vector(grange.use), ]
atac_counts2 <- atac_counts2[as.vector(grange.use2), ]
rownames(atac_counts)

########  Create ATAC assay and add it to the Seurat object ######## 
# Get gene annotations for hg38 and extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
show(annotations)   
#View(head(annotations,n=5))
#names(genomeStyles('Homo_sapiens'))
seqlevelsStyle(annotations) <- "UCSC"
length(extractSeqlevelsByGroup(species = 'Homo_sapiens', style = 'UCSC', group = 'all')) #group: sex, auto, circular
genome(annotations) <- "hg38"
show(annotations)

chrom_assay <- CreateChromatinAssay(counts = atac_counts,
                                    sep = c(":", "-"), fragments = fragpath, 
                                    annotation = annotations)
# Non-filtered assays by QCs - ATAC assay for sample 42_1 
show(chrom_assay)
hippo[["ATAC"]] <- chrom_assay
head(hippo@meta.data)
chrom_assay2 <- CreateChromatinAssay(counts = atac_counts2,
                                    sep = c(":", "-"), fragments = fragpath2, 
                                    annotation =annotations)
# Non-filtered assays by QCs - ATAC assay for sample 42_4
hippo2[["ATAC"]] <- chrom_assay2
head(hippo2@meta.data)

# Create .RData object with the rna-seq and atac-seq assays for samples 42_1 and 42_4
#save(chrom_assay, file = "/fastscratch/myscratch/csoto/counts_sample_42_1/chrom_assay_42_1.RData")
load("/fastscratch/myscratch/csoto/counts_sample_42_1/chrom_assay_42_1.RData")
#save(chrom_assay2, file = "/fastscratch/myscratch/csoto/counts_sample_42_4/chrom_assay_42_4.RData")
load("/fastscratch/myscratch/csoto/counts_sample_42_4/chrom_assay_42_4.RData")

#save(hippo, file = "/fastscratch/myscratch/csoto/counts_sample_42_1/rna_atac_assay_42_1.RData")
load("/fastscratch/myscratch/csoto/counts_sample_42_1/rna_atac_assay_42_1.RData")
#save(hippo2, file = "/fastscratch/myscratch/csoto/counts_sample_42_4/rna_atac_assay_42_4.RData")
load("/fastscratch/myscratch/csoto/counts_sample_42_4/rna_atac_assay_42_1.RData")

# Add the gene information to the object  ----  FUNCTIONS DUPLICATES ?? --- HEDIA 
show(annotations)
# Annotations of the object are set
Annotation(hippo[["ATAC"]]) <- annotations
Annotation(hippo2[["ATAC"]]) <- annotations
hippo2[["ATAC"]]
#Cells(hippo)
head(hippo)

######## Quality control to ATAC in sample 42_1 ########
DefaultAssay(hippo) <- "ATAC"
DefaultAssay(hippo2) <- "ATAC"

# Insert size distribution plot, called fragment length histogram 
# Still working on this chunk
FragmentHistogram(object = hippo)
FragmentHistogram(object = hippo2)
#FragmentHistogram(object = chrom_assay)
#FragmentHistogram(object = chrom_assay2)
#################################################################

# Signac QA measure to calculate the strength of the nucleosome signal per cell
hippo <- NucleosomeSignal(hippo)
head(hippo, n=5)

# Signac QA measure to compute the transcription start site (TSS) enrichment score for each cell, as defined by ENCODE.
hippo <- TSSEnrichment(hippo, fast = FALSE)    # If fast = False, it computes the TSS enrichment scores, storing the base-resolution matrix of integration counts at each site.
head(hippo,n=5)
# add blacklist ratio and fraction of reads in peaks
hippo$blacklist_fraction <- FractionCountsInRegion(
    object = hippo,
    assay = 'ATAC',
    regions = blacklist_hg38
)
# add blacklist ratio and fraction of reads in peaks
hippo$pct_reads_in_peaks <- hippo$atac_peak_region_fragments / hippo$atac_fragments * 100
hippo$blacklist_ratio <- hippo$blacklist_fraction / hippo$atac_peak_region_fragments
VlnPlot(hippo, 
        features = c("pct_reads_in_peaks","blacklist_ratio"), ncol = 2)
# Classified the TSS enrichment scores in two groups
hippo$high.tss <- ifelse(hippo$TSS.enrichment > 2, 'High', 'Low')
head(hippo, n =3)
TSSPlot(hippo, group.by = 'high.tss') + NoLegend()

# Group by cells with high or low nucleosomal signal strength. You can see that cells that are outliers for the mononucleosomal / nucleosome-free ratio
hippo$nucleosome_group <- ifelse(hippo$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = hippo, group.by = 'nucleosome_group')
head(hippo,n=5)
# Volcano plot for quality control 
# pdf(here("plots", "09_cellranger-arc_EDA", "Volcanoplot_QC_beforefiltering.pdf"))
# Visualize QC metrics as a violin plot / raw data 
VlnPlot(hippo, features = c("nCount_ATAC", "nFeature_ATAC", "nucleosome_signal", "TSS.enrichment"), ncol = 4)
# VlnPlot(object = hippo, 
#   features = c("nCount_ATAC", "nFeature_ATAC", "nucleosome_signal", "TSS.enrichment"),
#   ncol = 4,
#   log = TRUE,
#   pt.size = 0
# )
# dev.off()
 
######## Filter out low quality cells in the RNA-Seq and ATAC assays ########
# Filter cells that get a standard profile based on the QC results 
hippo2 <- subset(hippo, 
                 subset = percent.mt < 5 #subset = nFeature_RNA > 200 & nFeature_RNA < 7500 
                 & nucleosome_signal < 4
                 & TSS.enrichment > 1)
head(hippo2,n=5)
VlnPlot(hippo2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
VlnPlot(hippo2, features = c("nCount_ATAC", "nFeature_ATAC", "nucleosome_signal", "TSS.enrichment"), ncol = 4)
plot1 <- FeatureScatter(hippo2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(hippo2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Volcano plot for quality control 
# pdf(here("plots", "09_cellranger-arc_EDA", "Volcanoplot_QC_afterfiltering.pdf"))
# VlnPlot(
#   object = hippo,
#   features = c("nCount_RNA", "nCount_ATAC", "percent.mt"),
#   ncol = 3,
#   log = TRUE,
#   pt.size = 0
# )
# dev.off()

# no need to perform ATAC/RNA integration and label transfer since the two measurements are made in the same cell
# Dimensionality reduction
# RNA analysis
DefaultAssay(hippo) <- "RNA"

hippo <- SCTransform(hippo, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:30, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
# Warning: The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
# To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'

# ATAC analysis
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(hippo) <- "ATAC"
hippo <- FindTopFeatures(hippo, min.cutoff = "q0")
hippo <- RunTFIDF(hippo)
hippo <- RunSVD(hippo)
hippo <- RunUMAP(hippo, reduction = 'lsi', dims = 2:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

#Non-linear dimension reduction and clustering
# WNN graph, representing a weighted combination of RNA and ATAC-seq modalities
hippo <- FindMultiModalNeighbors(hippo, reduction.list = list("pca", "lsi"), dims.list = list(1:30, 2:30))
hippo <- RunUMAP(hippo, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
hippo <- FindClusters(hippo, graph.name = "wsnn", algorithm = 3, resolution = 0.5,verbose = TRUE)

# visualization
p1 <- DimPlot(hippo, reduction = "umap.rna", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(hippo, reduction = "umap.atac", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(hippo, reduction = "wnn.umap", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))

# findmarkers
#Map over the cluster levels and run FindMarkers for each cluster
DefaultAssay(hippo) <- "RNA"
cluster.markers <- map(levels(hippo$seurat_clusters[[1]]), ~ FindMarkers(hippo, ident.1 = .x, min.pct = 0.25))

#Add cluster number to each marker using map and mutate
cluster.markers <- map2(cluster.markers, seq_along(cluster.markers), ~ mutate(.x, cluster = .y - 1, gene = rownames(.x)))

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
DefaultAssay(hippo) <- "RNA"

cluster.markers.all <- FindAllMarkers(hippo, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 1)
cluster.markers.all %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

# top 5 
cluster.markers.all %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top5

#Linking peaks to genes
# DefaultAssay(hippo) <- "ATAC"
# 
# # first compute the GC content for each peak
# hippo <- RegionStats(hippo, genome = BSgenome.Hsapiens.UCSC.hg38,assay = "ATAC")

# link peaks to genes taking an example two genes
# hippo <- LinkPeaks(
#   object = hippo,
#   peak.assay = "ATAC",
#   expression.assay = "SCT",
#   genes.use = c("PRAMEF2", "PLCH2")
# )
# 
# p4 <- CoveragePlot(
#   object = hippo,
#   region = "PEX10",
#   features = "PEX10",
#   expression.assay = "SCT",
#   extend.upstream = 500,
#   extend.downstream = 500
# )
# 
# p5 <- CoveragePlot(
#   object = hippo,
#   region = "PLCH2",
#   features = "PLCH2",
#   expression.assay = "SCT",
#   extend.upstream = 500,
#   extend.downstream = 500
# )
# 
# patchwork::wrap_plots(p4, p5, ncol = 1)

# link peaks to differentially expressed features (cluster biomarkers)

pdf(file=here("plots", "09_cellranger-arc_EDA", "CoveragePlot_Top5markers_percluster.pdf"))

walk(seq_along(top5$gene), ~ {
  message(paste0("Processing gene ",.x, " ", top5$gene[.x]))
  
  tryCatch({
    DefaultAssay(hippo) <- "ATAC"
    hippo <- RegionStats(hippo, genome = BSgenome.Hsapiens.UCSC.hg38, assay = "ATAC")
    hippo <- LinkPeaks(object = hippo, peak.assay = "ATAC", expression.assay = "SCT", genes.use = top5$gene[.x])
    message(paste0("Generating coverage plot for gene ", top5$gene[.x]))
    
p <- CoveragePlot(
      object = hippo,
      region = top5$gene[.x],
      features = top5$gene[.x],
      expression.assay = "SCT",
      extend.upstream = 500,
      extend.downstream = 500
    )

print(p)
    
    message(paste0("Gene ", .x , " ", top5$gene[.x], " completed successfully"))
    message(paste0(length(top5$gene)-.x, " ", "genes still need to be processed"))
  }, error = function(e) {
    message(paste0("Error occurred while processing gene ", top5$gene[.x], ": ", e$message))
  })
  
})

dev.off()


# Find differentially accessible peaks between clusters
# change back to working with peaks instead of gene activities
DefaultAssay(hippo) <- "ATAC"
clusters <- unique(levels(hippo$seurat_clusters[[1]]))
pairwise <- combn(clusters, 2)

# all pairwise comparisons
# a logistic regression framework to determine differentially accessible regions
da_peaks <- map(1:ncol(pairwise), ~{
  markers <- FindMarkers(hippo, ident.1 = pairwise[1, .x], ident.2 = pairwise[2, .x],test.use = 'LR') 
  comparisons <- pairwise[, .x]
  markers$comparison <- paste(comparisons[1], comparisons[2], sep = '_')
  markers
})

da_peaks.df <- do.call(rbind, da_peaks)

da_peaks.df

# Find all markers
cluster.peaks.all <- FindAllMarkers(hippo, only.pos = TRUE, min.pct = 0.1)
cluster.markers.all %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)

# top 5 
cluster.peaks.all %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) -> top5