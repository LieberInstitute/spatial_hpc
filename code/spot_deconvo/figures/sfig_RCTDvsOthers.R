setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")
library("here")
library("jaffelab")
library("SpatialExperiment")
library("sessioninfo")
library("tidyverse")
library("reshape")

grp = "broad"
#grp = "layer"

Dr = here("processed-data", "spot_deconvo", "RCTD", "3rdRun_newClass_deconvoMarkers", grp)
csv_files = data.frame(files = list.files(Dr, pattern = "clusters.csv" , recursive = TRUE, full.names = TRUE))
counts_list <- lapply(csv_files$files, function(file_path) 
{data <- read.csv(file_path, row.names = NULL)  
return(data)})
RCTD <- do.call(rbind, counts_list)
RCTD$sample <-sapply(strsplit(RCTD$X,"-1_"), `[`, 2)
RCTD$sample <-sapply(strsplit(RCTD$sample,"_Br"), `[`, 1)

Dr = here("processed-data", "spot_deconvo", "RCTD", "2ndRun_newClass_RCTDmarkers", grp)
csv_files = data.frame(files = list.files(Dr, pattern = "clusters.csv" , recursive = TRUE, full.names = TRUE))
csv_files = csv_files[33:40,]
counts_list <- lapply(csv_files, function(file_path) 
{data <- read.csv(file_path, row.names = NULL)  
return(data)})
RCTD_rctd <- do.call(rbind, counts_list)
RCTD_rctd$sample <-sapply(strsplit(RCTD_rctd$X,"-1_"), `[`, 2)
RCTD_rctd$sample <-sapply(strsplit(RCTD_rctd$sample,"_Br"), `[`, 1)

Dr = here("processed-data", "spot_deconvo", "cell2location", "IF", "2ndRun_newClass", grp)
csv_files = data.frame(files = list.files(Dr, pattern = "clusters.csv" , recursive = TRUE, full.names = TRUE))
counts_list <- lapply(csv_files$files, function(file_path) 
{data <- read.csv(file_path, row.names = NULL)  
return(data)})
cell2location <- do.call(rbind, counts_list)

Dr = here("processed-data", "spot_deconvo", "tangram", "IF", "2ndRun_newClass", grp)
csv_files = data.frame(files = list.files(Dr, pattern = "clusters.csv" , recursive = TRUE, full.names = TRUE))
counts_list <- lapply(csv_files$files, function(file_path) 
{data <- read.csv(file_path, row.names = NULL)  
return(data)})
tangram <- do.call(rbind, counts_list)
tangram$sample <-sapply(strsplit(tangram$key,"-1_"), `[`, 2)
tangram$sample <-sapply(strsplit(tangram$sample,"_Br"), `[`, 1)


library(ggplot2)

# Create data frames for each table
RCTD_counts <- table(RCTD$sample)
RCTD_counts <- table(RCTD_rctd$sample)
cell2location_counts <- table(cell2location$sample)
tangram_counts <- table(tangram$sample)

# Combine counts into a single data frame
counts_df <- data.frame(
  Sample = names(RCTD_counts),
  RCTD = as.numeric(RCTD_counts),
  cell2location = as.numeric(cell2location_counts),
  Tangram = as.numeric(tangram_counts)
)

# Melt the data frame for plotting
counts_melted <- reshape2::melt(counts_df, id.vars = "Sample", variable.name = "Dataset")

png(here("plots","spot_deconvo","figures", "RCTDbarplot_deconvo.png"), width = 1150, height = 450, units = "px") 
#png(here("plots","spot_deconvo","figures", "RCTDbarplot_rctd.png"), width = 1150, height = 450, units = "px") 

p = ggplot(counts_melted, aes(x = Sample, y = value, fill = Dataset)) +
  geom_bar(stat = "identity", position = "dodge", color = "black") +  # Add border color
  theme_bw() +
  labs(title = "User defined markers", x = "Capture array", y = "Number of spots", fill = "Tool") +
#  labs(title = "RCTD markers", x = "Capture array", y = "Number of spots", fill = "Tool") +
  scale_fill_manual(values = c('black', 'grey40', 'white')) +
  theme(
    text = element_text(size = 30, colour = "black"),
    axis.text = element_text(size = 15, colour = "black"),
    panel.grid.minor = element_blank(), 
    panel.grid.major = element_blank(),
    legend.text = element_text(size = 24, colour = "black"),
    plot.title = element_text(hjust = 0.5)
  )
print(p)
dev.off()
