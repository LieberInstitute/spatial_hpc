setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")
suppressPackageStartupMessages(library("SingleCellExperiment"))
suppressPackageStartupMessages(library("DeconvoBuddies"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("sessioninfo"))
suppressPackageStartupMessages(library("here"))
suppressPackageStartupMessages(library("spatialLIBD"))

Dr <- here("processed-data","spot_deconvo","shared_utilities")
marker_stats_layer = readRDS(here(Dr,paste0("marker_stats_layer_celltype_class1_noHATAGABAAmy.rds")))
markers_layer <- as.data.frame(marker_stats_layer |>filter(rank_ratio <= 25,ratio > 1))
markers_layer$cellTypeResolution = "Fine"

marker_stats_broad = readRDS(here(Dr,paste0("marker_stats_broad_class.rds")))
markers_broad <- as.data.frame(marker_stats_broad |>filter(rank_ratio <= 25,ratio > 1))
markers_broad$cellTypeResolution = "Mid"

markers = rbind(markers_broad, markers_layer)
write.csv(markers, file = here("plots", "spot_deconvo", "markergenes_stats.csv"))

##########
library("jaffelab")
library("SpatialExperiment")
library("sessioninfo")
library("tidyverse")
library("reshape")

#grp = "broad"
grp = "layer"
tools <- c("CART", "tangram", "cell2location", "RCTD")
df_list <- list()
for (tool in tools){
  Dr = here("processed-data", "spot_deconvo", tool, "IF", "2ndRun_newClass", grp)
  if (tool == "CART") {Dr = here("processed-data", "spot_deconvo", "groundTruth", "03_CART")}
  if (tool == "RCTD") {Dr = here("processed-data", "spot_deconvo", tool, "2ndRun_newClass_RCTDmarkers", grp)}
  csv_files = data.frame(files = list.files(Dr, pattern = "clusters.csv" , recursive = TRUE, full.names = TRUE))
  if (tool == "RCTD") {csv_files = data.frame(files = csv_files[33:40,])}
  counts_list <- lapply(csv_files$files, function(file_path) 
  {data <- read.csv(file_path, row.names = NULL)  
  return(data)})
  temp_df <- do.call(rbind, counts_list)
  rownames(temp_df) <- NULL
  if (tool == "RCTD"){temp_df$key = temp_df$X}
  temp_df$key = sapply(strsplit(temp_df$key,"_Br"), `[`, 1)
  if (tool == "CART"){
    temp_df$count = temp_df$n_cells
  }else if (grp == "broad"){ 
    celltypes = c("Oligo", "Micro_Macro_T", "InhN", "ExcN", "Astro", "Vascular", "OPC", "CSF")
    temp_df$count = rowSums(temp_df[,celltypes])
    temp_df$neuron = temp_df$ExcN + temp_df$InhN
    temp_df$other = temp_df$CSF + temp_df$OPC + temp_df$Vascular
    temp_df$microglia = temp_df$Micro_Macro_T
    temp_df$astrocyte = temp_df$Astro
    temp_df$oligo = temp_df$Oligo
  }else if (grp == "layer"){
    inhb = c("GABA.LAMP5","GABA.MGE","GABA.CGE")
    exct = c("Thal","L5","CA2.4","GC","L2_3.PrS.PaS","CA1_ProS","L6_6b","Sub.2","Sub.1","L2_3.Prs.Ent","Cajal")
    csf = c("Choroid","Ependy")
    celltypes = c(inhb, exct, csf, "Micro_Macro_T", "Astro", "Oligo")
    temp_df$count = rowSums(temp_df[,celltypes])
    temp_df$neuron = rowSums(temp_df[,c(inhb,exct)])
    temp_df$other = rowSums(temp_df[,c(csf)])
    temp_df$microglia = temp_df$Micro_Macro_T
    temp_df$astrocyte = temp_df$Astro
    temp_df$oligo = temp_df$Oligo
  }
  celltypes = c('oligo', 'other', 'neuron', 'microglia', 'astrocyte')
  if (tool == "RCTD"){temp_df = temp_df}else{
    temp_df[, celltypes] <- lapply(temp_df[, celltypes], function(col) col / temp_df$count)}
  temp_df$sample <-sapply(strsplit(temp_df$key,"-1_"), `[`, 2)
  temp_df[which(temp_df$count == 0),celltypes] = 0
  temp_df = temp_df[,c("key","sample",celltypes)]
  temp_df$tool = tool
  df_list[[length(df_list) + 1]] <- temp_df
}
final_df <- do.call(rbind, df_list)
which(is.na(final_df), arr.ind=TRUE)

df = melt(final_df, id.variables = c("key", "sample_id", "tool"))

CART = df[which(df$tool == "CART"),]
RCTD = df[which(df$tool == "RCTD"),]
tangram = df[which(df$tool == "tangram"),]
cell2location = df[which(df$tool == "cell2location"),]
rm = setdiff(CART$key, RCTD$key)
CART = CART[!(CART$key %in% rm), ]

sorted_keys <- CART$key
RCTD <- RCTD[RCTD$key %in% sorted_keys, ]
RCTD$actual = CART$value
cell2location <- cell2location[cell2location$key %in% sorted_keys, ]
cell2location$actual = CART$value
tangram <- tangram[tangram$key %in% sorted_keys, ]
tangram$actual = CART$value

# CART = final_df[which(final_df$tool == "CART"),]
# RCTD = final_df[which(final_df$tool == "RCTD"),]
# tangram = final_df[which(final_df$tool == "tangram"),]
# cell2location = final_df[which(final_df$tool == "cell2location"),]
# rm = setdiff(CART$key, RCTD$key)
# CART = CART[!(CART$key %in% rm), ]
# 
# sorted_keys <- CART$key
# RCTD <- RCTD[RCTD$key %in% sorted_keys, ]
# cell2location <- cell2location[cell2location$key %in% sorted_keys, ]
# tangram <- tangram[tangram$key %in% sorted_keys, ]

df1 = rbind(RCTD, cell2location, tangram)

metrics_df = df1 |> filter(variable != "other")|>
  group_by(tool, variable, sample) |>
  summarize(
    corr = cor(value, actual),
    rmse = mean((value - actual)**2)**0.5
  ) 

metrics_df_layer = as.data.frame(metrics_df)
metrics_df_layer$cellTypeResolution = "Fine"

# metrics_df_broad = as.data.frame(metrics_df)
# metrics_df_broad$cellTypeResolution = "Mid"

metrics = rbind(metrics_df_broad, metrics_df_layer)
write.csv(metrics, file = here("plots", "spot_deconvo", "CARTvsTool_stats.csv"))


