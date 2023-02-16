############################################################
# spatial_HPC project
# PCA 2vs1 for capture area of PRECAST clusters & QC metrics
# Anthony Ramnauth, Feb 13 2022
############################################################

setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")

suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(spatialLIBD)
  library(here)
  library(scuttle)
  library(scater)
  library(scran)
  library(dplyr)
  library(PCAtools)
  library(schex)
  library(viridis)
  library(sessioninfo)
  library(gridExtra)
})

# Load SPE
load(file = here::here("processed-data", "06_clustering", "PRECAST", "spe_modify_PRECAST_k15.Rdata"))

# Maddy left the spots that failed QC, discard the failed spots QC'ed by sample_id
spe <- spe[, colData(spe)$discard_auto_id == FALSE]
dim(spe)

# Remove N/A cluster from PRECAST labeling
spe = spe[, which(spe$PRECAST_k15 != "<NA>")]

# Feature selection
dec <- modelGeneVar(spe, block = spe$sample_id)
top_hvgs <- getTopHVGs(dec, prop = 0.1)
length(top_hvgs)
head(top_hvgs)
top.hvgs.fdr5 <- getTopHVGs(dec, fdr.threshold = 0.05)
length(top.hvgs.fdr5)
head(top_hvgs)
top.hvgs.fdr1 <- getTopHVGs(dec, fdr.threshold = 0.01)
length(top.hvgs.fdr1)
head(top_hvgs)

# Dimensionality reduction
set.seed(12345)
spe <- runPCA(spe, subset_row = top_hvgs, ncomponents = 50)

pdf(file = here::here("plots","06_clustering", "PRECAST", "PC2vsPC1_wCP.pdf"), width = 14, height = 14)

plotPCA(spe, colour_by = "PRECAST_k15", ncomponents = 2, point_size = 1) +
    theme(text = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.line = element_line(size=2))

plotPCA(spe, colour_by = "sum_gene", ncomponents = 2, point_size = 1) +
    theme(text = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.line = element_line(size=2))

plotPCA(spe, colour_by = "sum_umi", ncomponents = 2, point_size = 1) +
    theme(text = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.line = element_line(size=2))

plotPCA(spe, colour_by = "subsets_Mito_percent", ncomponents = 2, point_size = 1) +
    theme(text = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.line = element_line(size=2))

dev.off()

# Remove CP clusters
spe = spe[, which(spe$PRECAST_k15 != "9")]
spe = spe[, which(spe$PRECAST_k15 != "15")]

# Feature selection again without CP
dec <- modelGeneVar(spe, block = spe$sample_id)
top_hvgs <- getTopHVGs(dec, prop = 0.1)
length(top_hvgs)
head(top_hvgs)
top.hvgs.fdr5 <- getTopHVGs(dec, fdr.threshold = 0.05)
length(top.hvgs.fdr5)
head(top_hvgs)
top.hvgs.fdr1 <- getTopHVGs(dec, fdr.threshold = 0.01)
length(top.hvgs.fdr1)
head(top_hvgs)

# Dimensionality reduction again without CP clusters
set.seed(12345)
spe <- runPCA(spe, subset_row = top_hvgs, ncomponents = 50)

# Create new column in colData for either cluster 8, 13, or other to map onto
# PC 2

df <- cbind.data.frame(
  colData(spe),
  spatialCoords(spe),
  reducedDim(spe, "PCA")
)

df <- df %>%
    mutate(blah = case_when(
        PRECAST_k15 == 8 ~ "Clust_8",
        PRECAST_k15 == 13 ~ "Clust_13",
        PRECAST_k15 == 6 ~ "Clust_6"
    ))

#Replace NA values with 0
df$blah[is.na(df$blah)] <- 0

spe$blah <- as.factor(df$blah)

# Use schex to circumvent overplotting of spots

hex <- make_hexbin(spe, nbins = 100,
                   dimension_reduction = "PCA", use_dims=c(1,2))

colors <- as.vector(Polychrome::palette36.colors(13))

pal <- c("#5A5156", "#FE00FA", "#16FF32", "#3283FE")

label_df <- make_hexbin_label(hex, col="PRECAST_k15")

label_blah <- make_hexbin_label(hex, col="blah")

pdf(file = here::here("plots","06_clustering", "PRECAST", "PC2vsPC1_woCP.pdf"), width = 14, height = 14)

plotPCA(spe, colour_by = "PRECAST_k15", ncomponents = 2, point_size = 1) +
    theme(text = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.line = element_line(size=2))

plot_hexbin_meta(hex, col = "PRECAST_k15", action = "majority",
                xlab = "PC1", ylab = "PC2", color = colors) +
    ggrepel::geom_label_repel(data = label_df, aes(x=x, y=y, label = label),
    colour="black",  label.size = NA, fill = NA) +
    theme(text = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.line = element_line(size=2))

plotPCA(spe, colour_by = "sum_gene", ncomponents = 2, point_size = 1) +
    theme(text = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.line = element_line(size=2))

plotPCA(spe, colour_by = "sum_umi", ncomponents = 2, point_size = 1) +
    theme(text = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.line = element_line(size=2))

plotPCA(spe, colour_by = "subsets_Mito_percent", ncomponents = 2, point_size = 1) +
    theme(text = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.line = element_line(size=2)) +
    guides(color = guide_legend(title = "Mito %"))

plot_hexbin_meta(hex, col = "blah", action = "majority",
                xlab = "PC1", ylab = "PC2", color = pal) +
    theme(text = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.line = element_line(size=2))

dev.off()


## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
