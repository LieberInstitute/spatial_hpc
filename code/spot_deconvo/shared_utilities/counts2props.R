setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")
library("here")
library("sessioninfo")
library("SpatialExperiment")
library('ggplot2')
library('dplyr')

spe = readRDS(here("processed-data", "spot_deconvo", "shared_utilities", "spe.rds"))
dat = as.data.frame(colData(spe)) %>% select("key", "cluster_collapsed")

## cell2location
tool =  "cell2location"
Dr = here("processed-data", "spot_deconvo",tool, "HE", "2ndRun_newClass", "layer")
csv_files = data.frame(files = list.files(Dr, pattern = "clusters.csv" , recursive = TRUE, full.names = TRUE))
csv_files$sample_id = sapply(strsplit(csv_files$files,"/"), `[`, 13)
counts_list <- list()
for (file_path in csv_files$files) {
  data <- read.csv(file_path, row.names = NULL)  # Read the CSV file
  counts_list[[file_path]] <- data  # Store the data frame in the list with the file path as the key
}

temp_df <- do.call(rbind, counts_list)
rownames(temp_df) <- NULL
celltypes = names(temp_df)[3:length(names(temp_df))]
temp_df$count = rowSums(temp_df[,celltypes])

counts = temp_df %>% select("key")
for (celltype in celltypes) {
  counts[,celltype] = temp_df[,celltype]/temp_df$count
}

counts$tool = tool
counts_match <- match(dat$key, counts$key)
counts_info <-counts[counts_match,]
counts_info <- subset(counts_info, select = -key)
cell2location <- cbind(dat, counts_info)

## tangram
tool =  "tangram"
Dr = here("processed-data", "spot_deconvo",tool, "HE", "2ndRun_newClass", "layer")
Dr = here("processed-data", "spot_deconvo",tool, "HE", "2ndRun_newClass", "layer")
csv_files = data.frame(files = list.files(Dr, pattern = "clusters.csv" , recursive = TRUE, full.names = TRUE))

csv_files$sample_id = sapply(strsplit(csv_files$files,"/"), `[`, 13)
counts_list <- list()
for (file_path in csv_files$files) {
  data <- read.csv(file_path, row.names = NULL)  # Read the CSV file
  counts_list[[file_path]] <- data  # Store the data frame in the list with the file path as the key
}

temp_df <- do.call(rbind, counts_list)
rownames(temp_df) <- NULL
temp_df$key = substr(temp_df$key, 1, nchar(temp_df$key) - 14)
celltypes = names(temp_df)[3:length(names(temp_df))]
temp_df$count = rowSums(temp_df[,celltypes])

counts = temp_df %>% select("key")

for (celltype in celltypes) {
  counts[,celltype] = temp_df[,celltype]/temp_df$count
}

which(is.na(counts), arr.ind=TRUE)
counts <- replace(counts, is.na(counts), 0)
counts$key <- gsub("Br2720", "V12F14-051", counts$key)
which(is.na(counts), arr.ind=TRUE)

counts_match <- match(dat$key, counts$key)
counts_info <-counts[counts_match,]
which(is.na(counts_info), arr.ind=TRUE)

counts_info <- subset(counts_info, select = -key)
tangram <- cbind(dat, counts_info)
tangram$tool = "tangram"

column_order <- colnames(cell2location)
tangram <- tangram[, column_order]

dat = rbind(cell2location, tangram)

## for box plots

ggplot(dat, aes(x = dat$cluster_collapsed, y = dat$GC)) + geom_boxplot(aes(fill=tool))

