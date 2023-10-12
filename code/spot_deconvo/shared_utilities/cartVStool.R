setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")
library("here")
library("jaffelab")
library("SpatialExperiment")
library("sessioninfo")
library("tidyverse")
library("reshape")

tools <- c("CART", "tangram", "cell2location", "RCTD")
df_list <- list()
  for (tool in tools){
    Dr = here("processed-data", "spot_deconvo", tool, "IF", "2ndRun_newClass", "broad")
    if (tool == "CART") {Dr = here("processed-data", "spot_deconvo", "groundTruth", "03_CART")}
    if (tool == "RCTD") {Dr = here("processed-data", "spot_deconvo", tool, "2ndRun_newClass_RCTDmarkers", "broad")}
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
    }else{
      celltypes = c("Oligo", "Micro_Macro_T", "InhN", "ExcN", "Astro", "Vascular", "OPC", "CSF")
      temp_df$count = rowSums(temp_df[,celltypes])
      temp_df$neuron = temp_df$ExcN + temp_df$InhN
      temp_df$other = temp_df$CSF + temp_df$OPC + temp_df$Vascular
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

CART = final_df[which(final_df$tool == "CART"),]
RCTD = final_df[which(final_df$tool == "RCTD"),]
tangram = final_df[which(final_df$tool == "tangram"),]
cell2location = final_df[which(final_df$tool == "cell2location"),]
rm = setdiff(CART$key, RCTD$key)
CART = CART[!(CART$key %in% rm), ]

sorted_keys <- CART$key
RCTD <- RCTD[RCTD$key %in% sorted_keys, ]
cell2location <- cell2location[cell2location$key %in% sorted_keys, ]
tangram <- tangram[tangram$key %in% sorted_keys, ]

RCTD_list <- split(RCTD[celltypes], RCTD$sample)
CART_list <- split(CART[celltypes], CART$sample)
tangram_list <- split(tangram[celltypes], tangram$sample)
cell2location_list <- split(cell2location[celltypes], cell2location$sample)

# Define a function to calculate RMSE
calculate_rmse <- function(observed, predicted) {signif(sqrt(mean((observed - predicted)^2)),3)}

# Calculate correlations and RMSE for each sample
results <- lapply(names(RCTD_list), function(sample_id) {
  corr_RCTD_CART <- round(cor(RCTD_list[[sample_id]], CART_list[[sample_id]]),2)
  corr_tangram_CART <- round(cor(tangram_list[[sample_id]], CART_list[[sample_id]]),2)
  corr_cell2location_CART <- round(cor(cell2location_list[[sample_id]], CART_list[[sample_id]]),2)
  
  rmse_RCTD_CART <- sapply(celltypes, function(col) {calculate_rmse(RCTD_list[[sample_id]][[col]], CART_list[[sample_id]][[col]])})
  rmse_tangram_CART <- sapply(celltypes, function(col) {calculate_rmse(tangram_list[[sample_id]][[col]], CART_list[[sample_id]][[col]])})
  rmse_cell2location_CART <- sapply(celltypes, function(col) {calculate_rmse(cell2location_list[[sample_id]][[col]], CART_list[[sample_id]][[col]])})
  #sapply(celltypes, function(col) {calculate_rmse(RCTD[[col]], CART[[col]])})
  
  data.frame(
    sample_id = sample_id,
    celltypes = celltypes,
    Corr_RCTD_CART = diag(corr_RCTD_CART),
    Corr_tangram_CART = diag(corr_tangram_CART),
    Corr_cell2location_CART = diag(corr_cell2location_CART),
    RMSE_RCTD_CART = rmse_RCTD_CART,
    RMSE_tangram_CART = rmse_tangram_CART,
    RMSE_cell2location_CART = rmse_cell2location_CART
  )
})

# Combine the results into a single data frame
results <- do.call(rbind, results)
results = results[which(results$celltypes != "other"),]

RCTD = results %>% select(sample_id = sample_id, celltypes = celltypes,Corr = Corr_RCTD_CART,RMSE = RMSE_RCTD_CART)
RCTD$tool = "RCTD"
cell2location = results %>% select(sample_id = sample_id, celltypes = celltypes,Corr = Corr_cell2location_CART, RMSE = RMSE_cell2location_CART)
cell2location$tool = "cell2location"
tangram = results %>% select(sample_id = sample_id, celltypes = celltypes,Corr = Corr_tangram_CART, RMSE = RMSE_tangram_CART)
tangram$tool = "tangram"

df = rbind(RCTD, cell2location, tangram)
# Load the ggplot2 library
library(ggplot2)
library(gridExtra)

p = ggplot(df, aes(x = Corr, y = RMSE, color = sample_id, shape = celltypes))+geom_point(size = 3)+facet_wrap(~tool)

rctd_corr = mean(RCTD$Corr)
rctd_rmse = mean(RCTD$RMSE)
tan_corr = mean(tangram$Corr)
tan_rmse = mean(tangram$RMSE)
cell_corr = mean(cell2location$Corr)
cell_rmse = mean(cell2location$RMSE)

#gridplot = grid.arrange(grobs = plot_list, nrow = length(celltypes))
ggsave(here("plots","spot_deconvo","shared_utilities","cartVStool.pdf"), plot = p, width = 12, height = 6)

