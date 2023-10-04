setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")
suppressPackageStartupMessages({
  library("here")
  library("SpatialExperiment")
  library("rtracklayer")
  library("lobstr")
  library("sessioninfo")
  suppressPackageStartupMessages(library("SingleCellExperiment"))
  suppressPackageStartupMessages(library("basilisk"))
  suppressPackageStartupMessages(library("zellkonverter"))
})

spg_in <- here("processed-data", "spot_deconvo", "shared_utilities" , "spg.rds")
out <- here("processed-data", "spot_deconvo", "tangram")

spg = readRDS(spg_in)

#### gather Counts results
Dr = here("processed-data", "spot_deconvo", "groundTruth", "03_CART")
csv_files = data.frame(files = list.files(Dr, pattern = "clusters.csv" , recursive = TRUE, full.names = TRUE))
csv_files$sample_id = sapply(strsplit(csv_files$files,"/"), `[`, 11)
counts_list <- list()
for (file_path in csv_files$files) {
  data <- read.csv(file_path, row.names = NULL)  # Read the CSV file
  counts_list[[file_path]] <- data  # Store the data frame in the list with the file path as the key
}

temp_df <- do.call(rbind, counts_list)
rownames(temp_df) <- NULL
colData(spg)$key = gsub("_Br3942", "", colData(spg)$key)
colData(spg)$key = gsub("_Br8325", "", colData(spg)$key)
counts_match <- match(colData(spg)$key, temp_df$key)
counts <-temp_df[counts_match,]
colData(spg)$count = counts$n_cells
  
source(here("code", "spot_deconvo", "shared_utilities","write_anndata.R"))

colData(spg)$dateImg = as.character(colData(spg)$dateImg)
write_anndata(spg, paste0(out,"spg.h5ad"))
