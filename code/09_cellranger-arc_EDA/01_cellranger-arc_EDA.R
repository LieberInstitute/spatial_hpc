########################################################################
## 1st EDA with scCellRanger-ARC
## Authors. CSC + HT / LIEBER
## Date. March 3rd, 2023
## Reference Link. 
## Use this first: https://stuartlab.org/signac/articles/pbmc_multiomic.html 
## https://stuartlab.org/signac/articles/pbmc_vignette.html
########################################################################

# Last modification: CSC
Sys.time()
#"2023-03-14 12:02:09 EDT" 
#"2023-03-15 16:00:57 EDT"

# load libraries
library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(tidyverse)
library(here)
library(GenomeInfoDb)
library(dplyr)
library(patchwork)

set.seed(1234)

# 1) LOAD RNA & ATAC DATA
sc_hippo <- Read10X_h5("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/rafael_rotation/cellranger_rerun/42_1/outs/filtered_feature_bc_matrix.h5")
# Read10X_h5() read count matrix from 10X CellRanger hdf5 file. This can be used to read both scATAC-seq and scRNA-seq matrices.
# In this MTX  each row represents a peak, predicted to represent a region of open chromatin.
fragpath <- "/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/rafael_rotation/cellranger_rerun/42_1/outs/atac_fragments.tsv.gz"
# fragments.tsv is the full list of all unique fragments across all single cells.
# it contains all fragments associated with each single cell, as opposed to only fragments that map to peaks.

# Get gene annotations for hg38 and extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
# change to UCSC style since the data was mapped to hg38
show(annotations)   # GRanger(), By default the show method displays 5 head and 5 tail elements
View(head(annotations,n=5))
#names(genomeStyles('Homo_sapiens'))
seqlevelsStyle(annotations) <- "UCSC"
#extractSeqlevelsByGroup(species="Homo_sapiens", style="UCSC", group="sex")
#extractSeqlevelsByGroup(species="Homo_sapiens", style="UCSC", group="auto")
genome(annotations) <- "hg38"
show(annotations)

######## extract RNA and ATAC data, plus Gene Annotation for hg38 #######
str(sc_hippo)   # read modalities
head(sc_hippo)
rna_counts <- sc_hippo$`Gene Expression`
atac_counts <- sc_hippo$Peaks
View(head(rna_counts,n=5))

######## Create a Seurat object containing the RNA ########
# Initialize the Seurat object with the raw (non-normalized data). 
hippo <- CreateSeuratObject(
  counts = rna_counts,
  assay = "RNA",
  project = "hippo-42_1"    # , min.cells = 3, min.features = 200)
)
class(hippo)
# check metadata
head(hippo@meta.data)

######## QA Standard metrics  ########

Idents(hippo) <- "hippo-42_1"
hippo[["percent.mt"]] <- PercentageFeatureSet(hippo, pattern = "^MT-")
head(hippo@meta.data)                 # Access cell-level meta-data / head(hippo[[]])
head(hippo[["percent.mt"]][])         # Access feature-level meta-data

# Visualize QC metrics as a violin plot
VlnPlot(hippo, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# Scatter plot to visualize feature-feature relationships, across the set of single cells.
plot1 <- FeatureScatter(hippo, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(hippo, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

######## ATAC assay: peaks in standard chromosomes were used for analysis ########
class(atac_counts)
show(atac_counts)
grange.counts <- StringToGRanges(rownames(atac_counts), 
                                 sep = c(":", "-"))      # StringToGRanges(), Convert a genomic coordinate string to a GRanges object
length(grange.counts)
View(head(grange.counts,n=10))
grange.use <- seqnames(grange.counts) %in% 
    standardChromosomes(grange.counts)
class(grange.use)       #S4Vector object / boolean
length(grange.use)
#View(!grange.use)
#View(seqnames(grange.counts)[!grange.use])
atac_counts <- atac_counts[as.vector(grange.use), ]
View(head(atac_counts,n=10))

# Annotation custom: something is not working yet:
## > CoveragePlot(
##   +   object = hippo,
##   +   region = "PEX10",
##   +   features = "PEX10",
##   +   expression.assay = "SCT",
##   +   extend.upstream = 500,
##   +   extend.downstream = 500
##   + )
## [W::hts_idx_load3] The index file is older than the data file: /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/rafael_rotation/cellranger_rerun/42_1/outs/atac_fragments.tsv.gz.tbi
## Error in annotation[annotation$type == "body", ] : 
##   incorrect number of dimensions
# annotations <- rtracklayer::import("/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/genes/genes.gtf.gz")
## annotations <- annotations[annotations$type == "gene"]
# annotations <- keepStandardChromosomes(annotations, pruning.mode = "coarse")
# genome(annotations) <- "hg38"
# annotations$gene_biotype <- tolower(annotations$gene_type)

## Takes forever to run!
# annotations$tx_id <- annotations$transcript_id
# annotation_txdb <- GenomicFeatures::makeTxDbFromGRanges(annotations)
# 
# annotation_txdb <- GenomicFeatures::makeTxDbFromGFF("/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/genes/genes.gtf.gz")
# annotation_chunk <- lapply(seq_along(annotations), function(x) { biovizBase::crunch(obj = annotation_txdb, which = annotations[x]) })
# annotation_chunk_all <- do.call(c, annotation_chunk)
# annotation_chunk_all$gene_name <- NA
# annotation_chunk_all$gene_biotype <- NA
# m <- match(annotation_chunk_all$gene_id, annotations$gene_id)
# annotation_chunk_all$gene_name[!is.na(m)] <- annotations$gene_id[m[!is.na(m)]]
# annotation_chunk_all$gene_biotype[!is.na(m)] <- annotations$gene_id[m[!is.na(m)]]

# faq recommendation
# gtf <- rtracklayer::import('/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/genes/genes.gtf.gz')
# gene.coords <- gtf[gtf$type == 'gene']
# seqlevelsStyle(gene.coords) <- 'UCSC'
# gene.coords <- keepStandardChromosomes(gene.coords, pruning.mode = 'coarse')
# annotations <- gtf

########  Create ATAC assay and add it to the object ######## 
chrom_assay <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotations
)
chrom_assay
class(chrom_assay)
# Unfiltered
hippo[["ATAC"]] <- chrom_assay
# Filtered
# hippo2[["ATAC"]] <- chrom_assay
View(head(hippo))

# Add the gene information to the object
show(annotations)
# Annotations of the object are set 
Annotation(hippo[["ATAC"]]) <- annotations
head(hippo)
# Annotation(hippo2[["ATAC"]]) <- annotations

######## Quality control to ATAC ########
DefaultAssay(hippo) <- "ATAC"

# Signac QA measure to calculate the strength of the nucleosome signal per cell.
hippo <- NucleosomeSignal(hippo)
# Signac QA measure to compute the transcription start site (TSS) enrichment score for each cell, as defined by ENCODE.
hippo <- TSSEnrichment(hippo, fast = FALSE)    # If fast = False, it computes the TSS enrichment scores, storing the base-resolution matrix of integration counts at each site.
head(hippo,n=5)
# # add blacklist ratio and fraction of reads in peaks ---------->need to add 'peak_region_fragments'
# hippo$pct_reads_in_peaks <- hippo$peak_region_fragments / hippo$passed_filters * 100    
# hippo$blacklist_ratio <- hippo$blacklist_region_fragments / hippo$peak_region_fragments
# hippo$high.tss <- ifelse(hippo$TSS.enrichment > 2, 'High', 'Low')
# TSSPlot(pbmc, group.by = 'high.tss') + NoLegend()
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
# Filter cells that have unique feature counts over 2,500 or less than 200, & cells that have >5% mitochondrial counts
hippo2 <- subset(hippo, 
                 subset = nFeature_RNA > 200 & nFeature_RNA < 7500 
                 & percent.mt < 5
                 & nucleosome_signal < 4
                 & TSS.enrichment > 2)
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