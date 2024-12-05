setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")
library("here")
library("sessioninfo")
library("SpatialExperiment")
library('ggplot2')
library('dplyr')
library('gridExtra')

spg_in <- here("processed-data", "spot_deconvo", "shared_utilities" , "spg.rds")
spg = readRDS(spg_in)
dat = as.data.frame(colData(spg)) %>% select("key", "array_row", "array_col", "sample_id")
dat$key = sapply(strsplit(dat$key,"_Br"), `[`, 1)
celltypes = c("oligo","other","neuron", "microglia", "astrocyte")

## CART
Dr = here("processed-data", "spot_deconvo", "groundTruth", "03_CART")
csv_files = data.frame(files = list.files(Dr, pattern = "clusters.csv" , recursive = TRUE, full.names = TRUE))
csv_files$sample_id = sapply(strsplit(csv_files$files,"/"), `[`, 11)
counts_list <- lapply(seq_along(csv_files$files), function(i) 
{data <- read.csv(csv_files$files[i], row.names = NULL) 
data$sample_id <- csv_files$sample_id[i]
return(data)})
temp_df <- do.call(rbind, counts_list)
rownames(temp_df) <- NULL

temp_df$count = temp_df$n_cells
temp_df[, celltypes] <- lapply(temp_df[, celltypes], function(col) col/temp_df$count)
temp_df$tool = "CART"

rmv = setdiff(temp_df$key, dat$key)
temp_df = temp_df[!(temp_df$key %in% rmv),]
counts_match <- match(dat$key, temp_df$key)
temp_df <-temp_df[counts_match,]
temp_df[which(temp_df$count == 0),celltypes] = 0
which(is.na(temp_df), arr.ind=TRUE)
CART <- cbind(dat, temp_df)
CART <- CART[,c("sample_id", "tool", "key", celltypes)]
which(is.na(CART), arr.ind=TRUE)

group = "broad"
celltypes1 = c("Oligo","Micro_Macro_T","InhN","ExcN","Astro","Vascular","OPC","CSF")

## cell2location
tool =  "cell2location"
Dr = here("processed-data", "spot_deconvo",tool, "IF", "2ndRun_newClass", group)
csv_files = data.frame(files = list.files(Dr, pattern = "clusters.csv" , recursive = TRUE, full.names = TRUE))
csv_files$sample_id = sapply(strsplit(csv_files$files,"/"), `[`, 13)
counts_list <- lapply(seq_along(csv_files$files), function(i) 
{data <- read.csv(csv_files$files[i], row.names = NULL) 
data$sample_id <- csv_files$sample_id[i]
return(data)})
temp_df <- do.call(rbind, counts_list)
rownames(temp_df) <- NULL

temp_df$count = rowSums(temp_df[,celltypes1])
temp_df$neuron = temp_df$ExcN + temp_df$InhN
temp_df$other = temp_df$CSF + temp_df$OPC + temp_df$Vascular
temp_df$microglia = temp_df$Micro_Macro_T
temp_df$astrocyte = temp_df$Astro
temp_df$oligo = temp_df$Oligo

temp_df[, celltypes] <- lapply(temp_df[, celltypes], function(col) col/temp_df$count)
temp_df$tool = tool
temp_df$key = sapply(strsplit(temp_df$key,"_Br"), `[`, 1)
counts_match <- match(dat$key, temp_df$key)
temp_df <-temp_df[counts_match,]
which(is.na(temp_df), arr.ind=TRUE)

cell2location <- cbind(dat, temp_df)
cell2location <- cell2location[,c("sample_id", "tool", "key", celltypes)]
which(is.na(cell2location), arr.ind=TRUE)

## tangram
tool =  "tangram"
Dr = here("processed-data", "spot_deconvo",tool, "IF", "2ndRun_newClass", group)
csv_files = data.frame(files = list.files(Dr, pattern = "clusters.csv" , recursive = TRUE, full.names = TRUE))

csv_files$sample_id = sapply(strsplit(csv_files$files,"/"), `[`, 13)
counts_list <- lapply(seq_along(csv_files$files), function(i) 
{data <- read.csv(csv_files$files[i], row.names = NULL) 
data$sample_id <- csv_files$sample_id[i]
return(data)})
temp_df <- do.call(rbind, counts_list)
rownames(temp_df) <- NULL
temp_df$key = substr(temp_df$key, 1, nchar(temp_df$key) - 21)

temp_df$count = rowSums(temp_df[,celltypes1])
temp_df$neuron = temp_df$ExcN + temp_df$InhN
temp_df$other = temp_df$CSF + temp_df$OPC + temp_df$Vascular
temp_df$microglia = temp_df$Micro_Macro_T
temp_df$astrocyte = temp_df$Astro
temp_df$oligo = temp_df$Oligo

temp_df[, celltypes] <- lapply(temp_df[, celltypes], function(col) col/temp_df$count)
which(is.na(temp_df), arr.ind=TRUE)

temp_df[which(temp_df$count == 0),celltypes] = 0
which(is.na(temp_df), arr.ind=TRUE)

temp_df$tool = tool
counts_match <- match(dat$key, temp_df$key)
temp_df <-temp_df[counts_match,]
which(is.na(temp_df), arr.ind=TRUE)

tangram <- cbind(dat, temp_df)
tangram <- tangram[,c("sample_id", "tool", "key", celltypes)]
which(is.na(tangram), arr.ind=TRUE)

## RCTD
tool =  "RCTD"
Dr = here("processed-data", "spot_deconvo",tool,"2ndRun_newClass_RCTDmarkers", group)
csv_files = data.frame(files = list.files(Dr, pattern = "clusters.csv" , recursive = TRUE, full.names = TRUE))
csv_files = data.frame(files = csv_files[c(33:40),])
csv_files$sample_id = sapply(strsplit(csv_files$files,"/"), `[`, 12)
counts_list <- lapply(seq_along(csv_files$files), function(i) 
{data <- read.csv(csv_files$files[i], row.names = NULL) 
data$sample_id <- csv_files$sample_id[i]
return(data)})
temp_df <- do.call(rbind, counts_list)
rownames(temp_df) <- NULL

temp_df$key = sapply(strsplit(temp_df$X,"_Br"), `[`, 1)
which(is.na(temp_df), arr.ind=TRUE)

temp_df$neuron = temp_df$ExcN + temp_df$InhN
temp_df$other = temp_df$CSF + temp_df$OPC + temp_df$Vascular
temp_df$microglia = temp_df$Micro_Macro_T
temp_df$astrocyte = temp_df$Astro
temp_df$oligo = temp_df$Oligo

temp_df$tool = tool
counts_match <- match(dat$key, temp_df$key)
temp_df <-temp_df[counts_match,]
which(is.na(temp_df), arr.ind=TRUE)

RCTD <- cbind(dat, temp_df)
RCTD <- RCTD[,c("sample_id", "tool", "key", celltypes)]
which(is.na(RCTD), arr.ind=TRUE)

dat1 = rbind(CART, cell2location, tangram, RCTD)
which(is.na(dat1), arr.ind=TRUE)

write.csv(dat1, file = here("processed-data","VSPG_image_stitching", paste0(group,"2broad_deconvoVSPG.csv")))

############################################################################################################
group = "layer"
inhb = c("GABA.LAMP5","GABA.MGE","GABA.CGE")
exct = c("Thal","L5","CA2.4","GC","L2_3.PrS.PaS","CA1_ProS","L6_6b","Sub.2","Sub.1","L2_3.Prs.Ent","Cajal")
csf = c("Choroid","Ependy")
celltypes1 = c(inhb, exct, csf, "Micro_Macro_T", "Astro", "Oligo")


## cell2location
tool =  "cell2location"
Dr = here("processed-data", "spot_deconvo",tool, "IF", "2ndRun_newClass", group)
csv_files = data.frame(files = list.files(Dr, pattern = "clusters.csv" , recursive = TRUE, full.names = TRUE))
csv_files$sample_id = sapply(strsplit(csv_files$files,"/"), `[`, 13)
counts_list <- lapply(seq_along(csv_files$files), function(i) 
{data <- read.csv(csv_files$files[i], row.names = NULL) 
data$sample_id <- csv_files$sample_id[i]
return(data)})
temp_df <- do.call(rbind, counts_list)
rownames(temp_df) <- NULL

temp_df$count = rowSums(temp_df[,celltypes1])
temp_df$neuron = rowSums(temp_df[,c(inhb,exct)])
temp_df$other = rowSums(temp_df[,c(csf)])
temp_df$microglia = temp_df$Micro_Macro_T
temp_df$astrocyte = temp_df$Astro
temp_df$oligo = temp_df$Oligo

temp_df[, celltypes] <- lapply(temp_df[, celltypes], function(col) col/temp_df$count)
temp_df$tool = tool
temp_df$key = sapply(strsplit(temp_df$key,"_Br"), `[`, 1)
counts_match <- match(dat$key, temp_df$key)
temp_df <-temp_df[counts_match,]
which(is.na(temp_df), arr.ind=TRUE)

cell2location <- cbind(dat, temp_df)
cell2location <- cell2location[,c("sample_id", "tool", "key", celltypes)]
which(is.na(cell2location), arr.ind=TRUE)

## tangram
tool =  "tangram"
Dr = here("processed-data", "spot_deconvo",tool, "IF", "2ndRun_newClass", group)
csv_files = data.frame(files = list.files(Dr, pattern = "clusters.csv" , recursive = TRUE, full.names = TRUE))

csv_files$sample_id = sapply(strsplit(csv_files$files,"/"), `[`, 13)
counts_list <- lapply(seq_along(csv_files$files), function(i) 
{data <- read.csv(csv_files$files[i], row.names = NULL) 
data$sample_id <- csv_files$sample_id[i]
return(data)})
temp_df <- do.call(rbind, counts_list)
rownames(temp_df) <- NULL
temp_df$key = substr(temp_df$key, 1, nchar(temp_df$key) - 21)

temp_df$count = rowSums(temp_df[,celltypes1])
temp_df$neuron = rowSums(temp_df[,c(inhb,exct)])
temp_df$other = rowSums(temp_df[,c(csf)])
temp_df$microglia = temp_df$Micro_Macro_T
temp_df$astrocyte = temp_df$Astro
temp_df$oligo = temp_df$Oligo

temp_df[, celltypes] <- lapply(temp_df[, celltypes], function(col) col/temp_df$count)
which(is.na(temp_df), arr.ind=TRUE)

temp_df[which(temp_df$count == 0),celltypes] = 0
which(is.na(temp_df), arr.ind=TRUE)

temp_df$tool = tool
counts_match <- match(dat$key, temp_df$key)
temp_df <-temp_df[counts_match,]
which(is.na(temp_df), arr.ind=TRUE)

tangram <- cbind(dat, temp_df)
tangram <- tangram[,c("sample_id", "tool", "key", celltypes)]
which(is.na(tangram), arr.ind=TRUE)

## RCTD
tool =  "RCTD"
Dr = here("processed-data", "spot_deconvo",tool,"2ndRun_newClass_RCTDmarkers", group)
csv_files = data.frame(files = list.files(Dr, pattern = "clusters.csv" , recursive = TRUE, full.names = TRUE))
csv_files = data.frame(files = csv_files[c(33:40),])
csv_files$sample_id = sapply(strsplit(csv_files$files,"/"), `[`, 12)
counts_list <- lapply(seq_along(csv_files$files), function(i) 
{data <- read.csv(csv_files$files[i], row.names = NULL) 
data$sample_id <- csv_files$sample_id[i]
return(data)})
temp_df <- do.call(rbind, counts_list)
rownames(temp_df) <- NULL

temp_df$key = sapply(strsplit(temp_df$X,"_Br"), `[`, 1)
which(is.na(temp_df), arr.ind=TRUE)

temp_df$neuron = rowSums(temp_df[,c(inhb,exct)])
temp_df$other = rowSums(temp_df[,c(csf)])
temp_df$microglia = temp_df$Micro_Macro_T
temp_df$astrocyte = temp_df$Astro
temp_df$oligo = temp_df$Oligo

temp_df$tool = tool
counts_match <- match(dat$key, temp_df$key)
temp_df <-temp_df[counts_match,]
which(is.na(temp_df), arr.ind=TRUE)

RCTD <- cbind(dat, temp_df)
RCTD <- RCTD[,c("sample_id", "tool", "key", celltypes)]
which(is.na(RCTD), arr.ind=TRUE)

dat1 = rbind(CART, cell2location, tangram, RCTD)
which(is.na(dat1), arr.ind=TRUE)

write.csv(dat1, file = here("processed-data","VSPG_image_stitching", paste0(group,"2broad_deconvoVSPG.csv")))