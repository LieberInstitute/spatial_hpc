setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")
library("here")
library("jaffelab")
library("SpatialExperiment")
library("sessioninfo")
library("tidyverse")
library("reshape")

spe = readRDS(here("processed-data", "spot_deconvo", "shared_utilities", "spe.rds"))
dat = as.data.frame(colData(spe)) %>% select("key", "cluster_collapsed")
count = as.data.frame(colData(spe)) %>% select("key", "count", "cluster_collapsed")

#grp = "broad"
grp = "layer"
tools <- c("tangram", "cell2location", "RCTD")
df_list <- list()
for (tool in tools){
  Dr = here("processed-data", "spot_deconvo", tool, "HE", "2ndRun_newClass", grp)
  if (tool == "RCTD") {Dr = here("processed-data", "spot_deconvo", tool, "2ndRun_newClass_RCTDmarkers", grp)}
  csv_files = data.frame(files = list.files(Dr, pattern = "clusters.csv" , recursive = TRUE, full.names = TRUE))
  if (tool == "RCTD") {csv_files = data.frame(files = csv_files[-c(33:40),])}
  counts_list <- lapply(csv_files$files, function(file_path) 
  {data <- read.csv(file_path, row.names = NULL)  
  return(data)})
  temp_df <- do.call(rbind, counts_list)
  rownames(temp_df) <- NULL
  if (tool == "RCTD"){temp_df$key = temp_df$X}
  temp_df$key = sapply(strsplit(temp_df$key,"_Br"), `[`, 1)
  if (grp == "broad"){ 
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
  if (tool == "RCTD"){
    temp_df <- temp_df[temp_df$key %in% count$key, ]
    cnt = count$count[count$key %in% temp_df$key]
    temp_df[, celltypes] <- lapply(temp_df[, celltypes], function(col) col*cnt)
    temp_df$count = cnt}
  if (tool == "tangram") {
  temp_df$sample <-sapply(strsplit(temp_df$key,"-1_"), `[`, 2)
  temp_df$sample <-sapply(strsplit(temp_df$sample,"_V"), `[`, 1)
  temp_df$key <- sub("(_V1.*)","","AAACAAGTATCTCCCA-1_V10B01-085_A1_V10B01-085_A1")
  temp_df$key = paste0(temp_df$key,"_",temp_df$sample)
  }
  temp_df$sample <-sapply(strsplit(temp_df$key,"-1_"), `[`, 2)
  temp_df[which(temp_df$count == 0),celltypes] = 0
  temp_df = temp_df[,c("key","sample",celltypes)]
  temp_df$tool = tool
  
  df_list[[length(df_list) + 1]] <- temp_df
}
final_df <- do.call(rbind, df_list)
which(is.na(final_df), arr.ind=TRUE)

write.csv(final_df, file = here("processed-data","10-image_stitching", paste0(grp,"2broad_deconvoHE.csv")))