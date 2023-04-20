########################################################################
## 1st EDA with scCellRanger-ARC
## Authors. CSC + HT / LIEBER
## Date. March 3rd, 2023
## Reference Link. 
## Use this first: https://stuartlab.org/signac/articles/pbmc_multiomic.html 
## https://stuartlab.org/signac/articles/pbmc_vignette.html
########################################################################

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
sc_hippo2 <- Read10X_h5("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/rafael_rotation/cellranger_rerun/42_4/outs/filtered_feature_bc_matrix.h5")

# For QC. The raw matrix (unfiltered) is required to assess the total UMI_counts of a cell barcode. UMI elbow plots require the unfiltered count mtx to project the inflection point.
sc_hippo_raw <- Read10X_h5("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/rafael_rotation/cellranger_rerun/42_1/outs/raw_feature_bc_matrix.h5")
sc_hippo_raw2 <- Read10X_h5("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/rafael_rotation/cellranger_rerun/42_4/outs/raw_feature_bc_matrix.h5")

class(sc_hippo)
# Read10X_h5() read count matrix from 10X CellRanger hdf5 file. This can be used to read both scATAC-seq and scRNA-seq matrices.

######## extract RNA and ATAC data, plus Gene Annotation for hg38 #######
str(sc_hippo)   # read modalities
show(sc_hippo)
###############.  Filtered count mtx  ###############
# sample 42_1
rna_counts <- sc_hippo$`Gene Expression`
atac_counts <- sc_hippo$Peaks
head(rna_counts, n=5)
# sample 42_4
rna_counts2 <- sc_hippo2$`Gene Expression`
atac_counts2 <- sc_hippo2$Peaks
###############.  UnFiltered count mtx  #############
# sample 42_1 and 42_4 raw data
rna_counts_raw <- sc_hippo_raw$`Gene Expression`
atac_counts_raw <- sc_hippo_raw$Peaks
rna_counts2_raw <- sc_hippo_raw2$`Gene Expression`
atac_counts2_raw <- sc_hippo_raw2$Peaks
# compare size 
#all.equal.raw(rna_counts, rna_counts_raw). #  Mean relative difference: 271.9585

# Create .RData object to load at any time: count matrices sample 42_1 ans 42_4
#save(rna_counts, atac_counts, file = "/fastscratch/myscratch/csoto/counts_sample_42_1/arc_atac_42_1.RData")
#load("/fastscratch/myscratch/csoto/counts_sample_42_1/arc_atac_42_1.RData")
#save(rna_counts2, atac_counts2, file = "/fastscratch/myscratch/csoto/counts_sample_42_4/arc_atac_42_4.RData")
#load("/fastscratch/myscratch/csoto/counts_sample_42_1/arc_atac_42_4.RData")

############## Create inflection point plots based on UMI read number ##########
# Use the raw mtx instead of the filtered cellranger output
#class(rna_counts2_raw)
# All genes x all barcodes
counts <- Matrix::colSums(rna_counts1_raw) # calculate total UMI read number for each cell barcode
class(counts)
countdf <- as.data.frame(counts) %>%
    as_tibble(rownames = "barcode") %>%
    filter(counts >= 1) %>% # throw out cell barcodes with x or less UMI, Cells<500 is the std low cut off threshold 
    arrange(desc(counts)) %>% # arrange by descending order
    mutate(rank = 1:n()) # rank
head(countdf) # barcodes now ranked by UMI counts
tail(countdf)
ggplot(countdf, aes(x = rank, y = counts)) +
    geom_point() +
    labs(x = "barcodes (Unfiltered-log10)", y = "UMI_counts ") +
    theme_classic() +
    scale_x_log10() +
    scale_y_log10() +
    geom_hline(yintercept = 12000, col = 'darkgreen') +     #SAMPLE 42_1
    geom_hline(yintercept = 2364, col = 'dodgerblue') +    #SAMPLE 42_1  
#    geom_hline(yintercept = 9000, col = 'darkgreen') +     #SAMPLE 42_4
#    geom_hline(yintercept = 615, col = 'dodgerblue') +    #SAMPLE 42_4       
    ggtitle("Sample 42_1: Mean Genes Per Cell") 
legend("bottomleft", legend = c("Knee (12000)", "Inflection (2364)"), col = c("darkgreen", "dodgerblue"), lty = 2, cex = 0.8) #SAMPLE 42_1  
#legend("bottomleft", legend = c("Knee (9000)", "Inflection (615)"), col = c("darkgreen", "dodgerblue"), lty = 2, cex = 0.8) #SAMPLE 42_4  

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
#save(rna_counts, atac_counts, hippo, file = "/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/09_cellranger-arc_EDA/rna_assay_42_1.RData")
load("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/09_cellranger-arc_EDA/rna_assay_42_1.RData")
#save(rna_counts2, atac_counts2, hippo2, file = "/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/09_cellranger-arc_EDA/rna_assay_42_4.RData")
load("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/09_cellranger-arc_EDA/rna_assay_42_4.RData")

# Add number of genes per UMI for each cell to metadata
hippo$log10GenesPerUMI <- log10(hippo$nFeature_RNA) / log10(hippo$nCount_RNA)
hippo2$log10GenesPerUMI <- log10(hippo2$nFeature_RNA) / log10(hippo2$nCount_RNA)

# Merge the two seurat objects to compare one each other
hippo.combined <- merge(hippo, y = hippo2, add.cell.ids = c("Sample_42_1", "Sample_42_4"), project = "HIPPO")
# class(hippo.combined)
# show(hippo.combined)
# head(colnames(hippo.combined))
# table(hippo.combined$orig.ident)
# head(hippo.combined, n=3)
# tail(hippo.combined, n=3)
hippo.combined[["percent.mt"]] <- PercentageFeatureSet(hippo.combined, pattern = "^MT-")
head(hippo.combined@meta.data)                 # Access cell-level meta-data / head(hippo[[]])
tail(hippo.combined[["percent.mt"]][])         # Access feature-level meta-data

######## QA metrics for the 'Gene Expression' assay  ########

# Visualize number of cells per sample
genes_per_cell <- as.data.frame(hippo.combined[[]])
head(genes_per_cell,n=5)
genes_per_cell %>% 
    ggplot(aes(x=orig.ident, fill=orig.ident)) + 
    geom_bar(alpha = 0.7) +
    theme_classic() +
#    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    theme(plot.title = element_text(hjust=0.5)) +   #, face="bold"
    xlab("Sample") +
    ggtitle("Number of Cells Per Sample")

# Visualize the number UMIs/transcripts per cell
genes_per_cell %>% 
    ggplot(aes(color=orig.ident, x=nCount_RNA, fill=orig.ident)) + 
    geom_density(alpha = 0.2) + 
    scale_x_log10() + 
    theme_classic() +
    theme(plot.title = element_text(hjust=0.5)) +   #, face="bold"
    ylab("Log10 Cell density") + 
    xlab("UMI/transcripts") +
    ggtitle("Number UMI (transcripts) Per Cell") +
    geom_vline(xintercept = 500)     # should generally be above 500

# # Visualize the distribution of genes detected per cell via histogram
# genes_per_cell %>% 
#     ggplot(aes(color=orig.ident, x=nFeature_RNA, fill= orig.ident)) + 
#     geom_density(alpha = 0.2) + 
#     theme_classic() +
#     scale_x_log10() + 
#     geom_vline(xintercept = 300)

# Visualize the distribution of genes detected per cell via boxplot
genes_per_cell %>% 
    ggplot(aes(x=orig.ident, y=log10(nFeature_RNA), fill=orig.ident)) + 
    geom_boxplot(alpha = 0.7) + 
    theme_classic() +
    theme(axis.text.x = element_text(vjust = 1, hjust=1)) +
    theme(plot.title = element_text(hjust=0.5)) +
    ylab("Log10(Number Genes)") + 
    xlab("") +
    ggtitle("Number Cells vs Number Genes")

# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
colnames(genes_per_cell)
genes_per_cell %>% 
    ggplot(aes(x=nCount_RNA, y=nFeature_RNA, color=percent.mt)) + 
    geom_point() + 
    scale_colour_gradient(low = "gray90", high = "black") +
    stat_smooth(method=lm) +
    scale_x_log10() + 
    scale_y_log10() + 
    theme_classic() +
    geom_vline(xintercept = 500) +
    geom_hline(yintercept = 250) +
    facet_wrap(~orig.ident) + 
    ylab("Log10(Number Genes)") + 
    xlab("Number Cells") +
    ggtitle("Cells vs Genes By Fraction of ^MT reads")
############## Barcode inflection with Seurat functions to compare knee and inflection points ################ 

# CalculateBarcodeInflections(), calculates an adaptive inflection point ("knee") of the barcode distribution for each sample group. This is useful for determining a threshold for removing low-quality samples.
# Calculated inflection points & test several thresholds derived from the barcode-rank distribution for the nCount_RNA 
col_to_use = "nCount_RNA"
hippo_rank <- CalculateBarcodeInflections(
    hippo.combined,
    barcode.column = col_to_use,
    group.column = "orig.ident",
    threshold.low = 6,
    #    threshold.high = 300000
)
hippo_rank@tools$CalculateBarcodeInflections$inflection_points   # Get the inflection points in the two samples   
# Lables for nCounts
s42_1 = paste('S 42_1 Inflection point:', as.character(hippo_rank@tools$CalculateBarcodeInflections$inflection_points$nCount_RNA[1]),
              ' Rank:', as.character(hippo_rank@tools$CalculateBarcodeInflections$inflection_points$rank[1]))
s42_4 = paste('S 42_4 Inflection point:', as.character(hippo_rank@tools$CalculateBarcodeInflections$inflection_points$nCount_RNA[2]),
              ' Rank:', as.character(hippo_rank@tools$CalculateBarcodeInflections$inflection_points$rank[2]))
############# Here jump to plot the BarcodeInflectionsPlot() below

# Calculated inflection points & test several thresholds derived from the barcode-rank distribution for the nFeature 
col_to_use = "nFeature_RNA"
hippo_rank <- CalculateBarcodeInflections(
    hippo.combined,
    barcode.column = col_to_use,
    group.column = "orig.ident",
#    threshold.low = 6,
#    threshold.high = 300000
)
# Labels for nFeature_RNA
s42_1 = paste('S 42_1 Inflection point:', as.character(hippo_rank@tools$CalculateBarcodeInflections$inflection_points$nFeature_RNA[1]),
              ' Rank:', as.character(hippo_rank@tools$CalculateBarcodeInflections$inflection_points$rank[1]))
s42_4 = paste('S 42_4 Inflection point:', as.character(hippo_rank@tools$CalculateBarcodeInflections$inflection_points$nFeature_RNA[2]),
              ' Rank:', as.character(hippo_rank@tools$CalculateBarcodeInflections$inflection_points$rank[2]))
############# Here jump to plot the BarcodeInflectionsPlot() below

# Calculated inflection points & test several thresholds derived from the barcode-rank distribution for the ^MT% Content
col_to_use = "percent.mt"
hippo_rank <- CalculateBarcodeInflections(
    hippo.combined,
    barcode.column = col_to_use,
    group.column = "orig.ident",
    threshold.low = 6,
    #    threshold.high = 300000
)
# Labels for percent.mt
s42_1 = paste('S 42_1 Inflection point:', as.character(round(hippo_rank@tools$CalculateBarcodeInflections$inflection_points$percent.mt[1],2)), '%',
              ' Rank:', as.character(hippo_rank@tools$CalculateBarcodeInflections$inflection_points$rank[1]))
s42_4 = paste('S 42_4 Inflection point:', as.character(round(hippo_rank@tools$CalculateBarcodeInflections$inflection_points$percent.mt[2],2)), '%',
              ' Rank:', as.character(hippo_rank@tools$CalculateBarcodeInflections$inflection_points$rank[2]))

hippo_rank@tools$CalculateBarcodeInflections$inflection_points   # Get the inflection points in the two samples   

# Plot for any calculation distribution and trace the inflection points for any feature described before
BarcodeInflectionsPlot(hippo_rank) 
legend("top", legend = c(s42_1, s42_4), col = c("red", "green"), cex = 0.7) #SAMPLE 42_1  

# Create seurat object with the samples 42_1 and 42_4 merged (only rna counts mtx)
#save(hippo.combined, file = "/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/09_cellranger-arc_EDA/rna_assay_merged_42.RData")
load("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/09_cellranger-arc_EDA/rna_assay_merged_42.RData")

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
#save(chrom_assay, file ="/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/09_cellranger-arc_EDA/chrom_assay_42_1.RData" )
load("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/09_cellranger-arc_EDA/chrom_assay_42_1.RData")
#save(chrom_assay2, file = "/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/09_cellranger-arc_EDA/chrom_assay_42_4.RData")
load("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/09_cellranger-arc_EDA/chrom_assay_42_4.RData")
# composed objects (rna+chrome)
show(hippo)
#save(hippo, file = "/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/09_cellranger-arc_EDA/rna_chrom_assay_42_1.RData")
load("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/09_cellranger-arc_EDA/rna_chrom_assay_42_1.RData")
#save(hippo2, file = "/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/09_cellranger-arc_EDA/rna_chrom_assay_42_4.RData")
load("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/09_cellranger-arc_EDA/rna_chrom_assay_42_4.RData")

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
