setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")
library("here")
library("jaffelab")
library("SpatialExperiment")
library("sessioninfo")
library("tidyverse")
library("reshape")

Dr = here("processed-data", "spot_deconvo", "RCTD", "2ndRun_newClass_RCTDmarkers", "broad")
csv_files = data.frame(files = list.files(Dr, pattern = "clusters.csv" , recursive = TRUE, full.names = TRUE))
csv_files = data.frame(files = csv_files[c(1:32,41:44),])
counts_list <- lapply(csv_files$files, function(file_path) 
    {data <- read.csv(file_path, row.names = NULL)  
    return(data)})
temp_df <- do.call(rbind, counts_list)
rownames(temp_df) <- NULL
temp_df$X <- gsub("Br2720", "V12F14-051", temp_df$X)
# Initialize and populate mid2broad_RCTD
mid2broad_RCTD <- data.frame(
    key = sapply(strsplit(temp_df$X, "_Br"), `[`, 1), # Extracting the 'key' column
    neuron = temp_df$ExcN + temp_df$InhN,             # Summing Excitatory and Inhibitory neurons
    other = temp_df$CSF + temp_df$OPC + temp_df$Vascular, # Summing other cell types
    microglia = temp_df$Micro_Macro_T,               # Microglia counts
    astrocyte = temp_df$Astro,                       # Astrocytes
    oligo = temp_df$Oligo                            # Oligodendrocytes
)
which(is.na(mid2broad_RCTD), arr.ind=TRUE)

colnames(mid2broad_RCTD)[2:6] <- paste('mid2broad_RCTD', colnames(mid2broad_RCTD)[2:6], sep = "_")

Dr = here("processed-data", "spot_deconvo", "RCTD", "2ndRun_newClass_RCTDmarkers", "layer")
csv_files = data.frame(files = list.files(Dr, pattern = "clusters.csv" , recursive = TRUE, full.names = TRUE))
csv_files = data.frame(files = csv_files[c(1:32,41:44),])
counts_list <- lapply(csv_files$files, function(file_path) 
    {data <- read.csv(file_path, row.names = NULL)  
    return(data)})
temp_df <- do.call(rbind, counts_list)
rownames(temp_df) <- NULL
temp_df$X <- gsub("Br2720", "V12F14-051", temp_df$X)

inhb = c("GABA.LAMP5","GABA.MGE","GABA.CGE")
exct = c("Thal","L5","CA2.4","GC","L2_3.PrS.PaS","CA1_ProS","L6_6b","Sub.2","Sub.1","L2_3.Prs.Ent","Cajal")
csf = c("Choroid","Ependy")

# Initialize and populate mid2broad_RCTD
fine2broad_RCTD <- data.frame(
    key = sapply(strsplit(temp_df$X, "_Br"), `[`, 1), # Extracting the 'key' column
    neuron = rowSums(temp_df[,c(inhb,exct)]),             # Summing Excitatory and Inhibitory neurons
    other = rowSums(temp_df[,c(csf)]), # Summing other cell types
    microglia = temp_df$Micro_Macro_T,               # Microglia counts
    astrocyte = temp_df$Astro,                       # Astrocytes
    oligo = temp_df$Oligo                            # Oligodendrocytes
)

colnames(fine2broad_RCTD)[2:6] <- paste('fine2broad_RCTD', colnames(fine2broad_RCTD)[2:6], sep = "_")
# Get the intersection of keys
common_keys <- intersect(mid2broad_RCTD$key, fine2broad_RCTD$key)

# Filter both dataframes to retain only rows with common keys
mid2broad_RCTD <- mid2broad_RCTD[mid2broad_RCTD$key %in% common_keys, ]
fine2broad_RCTD <- fine2broad_RCTD[fine2broad_RCTD$key %in% common_keys, ]

# Ensure the dataframes are aligned by key
mid2broad_RCTD <- mid2broad_RCTD[order(mid2broad_RCTD$key), ]
fine2broad_RCTD <- fine2broad_RCTD[order(fine2broad_RCTD$key), ]

# Combine the dataframes
df <- cbind(mid2broad_RCTD, fine2broad_RCTD[, 2:6])
which(is.na(df), arr.ind=TRUE)
write.csv(df, file = here("processed-data","VSPG_image_stitching", "deconvoHE.csv"))