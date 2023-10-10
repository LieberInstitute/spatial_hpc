setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")
library("here")
library("jaffelab")
library("SpatialExperiment")
library("sessioninfo")
library("tidyverse")
library("reshape")

group <- c("broad", "layer")
tools <- c("tangram", "cell2location", "RCTD")

spg_in <- here("processed-data", "spot_deconvo", "shared_utilities", "spg.rds")
spg = readRDS(spg_in)
# unqiue(spg$brnum)
# Br3942_VSPG Br8325_VSPG
temp_df$key = sapply(strsplit(temp_df$key,"_Br"), `[`, 1)
celltypes = c('oligo', 'other', 'neuron', 'microglia', 'astrocyte')

#### gather CART results
tool =  "CART"
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
#which(is.na(temp_Df), arr.ind=TRUE)
counts = temp_df %>% select("key")
for (celltype in celltypes) {
  counts[,celltype] = temp_df[,celltype]/temp_df$n_cells
}
counts[which(temp_df$n_cells == 0),-c(1)] = 0

counts_match <- match(colData(spg)$key, counts$key)
CART <-counts[counts_match,]

CART$sample_id <-sapply(strsplit(CART$key,"-1_"), `[`, 2)
which(is.na(CART), arr.ind=TRUE)
#counts$tool = tool
#CART <- counts |> melt(id.vars = "key", variable.name = "celltype",value.name = "count")

# colnames(counts)[colnames(counts) %in% celltypes] <- paste0(tool, "_", colnames(counts)[colnames(counts) %in% celltypes])
# counts_match <- match(colData(spg)$key, counts$key)
# counts_info <-counts[counts_match,]
# counts_info <- subset(counts_info, select = -key)
# colData(spg) <- cbind(colData(spg), counts_info)

df_list <- list()
for (grp in group){
  
  for (tool in tools){
    Dr = here("processed-data", "spot_deconvo", tool, "IF", "2ndRun_newClass", grp)
    csv_files = data.frame(files = list.files(Dr, pattern = "clusters.csv" , recursive = TRUE, full.names = TRUE))
    counts_list <- lapply(csv_files$files, function(file_path) 
    {data <- read.csv(file_path, row.names = NULL)  
    return(data)})
    temp_df <- do.call(rbind, counts_list)
    rownames(temp_df) <- NULL
    temp_df$key = sapply(strsplit(temp_df$key,"_Br"), `[`, 1)
    
    if (grp == "broad"){
      celltypes = c("Oligo", "Micro_Macro_T", "InhN", "ExcN", "Astro", "Vascular", "OPC", "CSF")
      temp_df$count = rowSums(temp_df[,celltypes])
      temp_df$neuron = temp_df$ExcN + temp_df$InhN
      temp_df$other = temp_df$CSF + temp_df$OPC + temp_df$Vascular
      temp_df$microglia = temp_df$Micro_Macro_T
      temp_df$astrocyte = temp_df$Astro
      temp_df$oligo = temp_df$Oligo
    }else{
      neurons = c("GABA.LAMP5","GABA.MGE","GABA.CGE","Thal","L5","CA2.4","GC","L2_3.PrS.PaS",
                  "CA1_ProS","L6_6b","Sub.2","Sub.1","L2_3.Prs.Ent","Cajal")
      others= c("Choroid","Ependy","Vascular","OPC")
      celltypes = c(neurons, others, "Micro_Macro_T", "Astro", "Oligo")
      temp_df$count = rowSums(temp_df[,celltypes])
      temp_df$neuron = rowSums(temp_df[,neurons])
      temp_df$other = rowSums(temp_df[,others])
      temp_df$microglia = temp_df$Micro_Macro_T
      temp_df$astrocyte = temp_df$Astro
      temp_df$oligo = temp_df$Oligo
    }
    
    celltypes = c('oligo', 'other', 'neuron', 'microglia', 'astrocyte')
    #  temp_df[, celltypes] <- lapply(temp_df[, celltypes], function(col) col / temp_df$count)
    temp_df$sample <-sapply(strsplit(temp_df$key,"-1_"), `[`, 2)
    temp_df$group = grp
    temp_df[which(temp_df$count == 0),celltypes] = 0
    temp_df = temp_df[,c("key","sample",celltypes,"group")]
    temp_df$tool = tool
    df_list[[length(df_list) + 1]] <- temp_df
  }
}
final_df <- do.call(rbind, df_list)
which(is.na(final_df), arr.ind=TRUE)

## RCTD ####
tool =  "RCTD"
Dr = here("processed-data", "spot_deconvo", tool, "2ndRun_newClass_RCTDmarkers", "broad")
csv_files = data.frame(files = list.files(Dr, pattern = "clusters.csv" , recursive = TRUE, full.names = TRUE))
csv_files = csv_files[33:40,]
counts_list <- list()
for (file_path in csv_files) {
  data <- read.csv(file_path, row.names = NULL)  # Read the CSV file
  counts_list[[file_path]] <- data  # Store the data frame in the list with the file path as the key
}

temp_df <- do.call(rbind, counts_list)
rownames(temp_df) <- NULL
colnames(temp_df)[colnames(temp_df) == "X"] = "key"
temp_df$key = gsub("_Br3942", "", temp_df$key)
temp_df$key = gsub("_Br8325", "", temp_df$key)
temp_df$neuron = temp_df$ExcN + temp_df$InhN
temp_df$other = temp_df$CSF + temp_df$OPC + temp_df$Vascular
temp_df$microglia = temp_df$Micro_Macro_T
temp_df$astrocyte = temp_df$Astro
temp_df$oligo = temp_df$Oligo
temp_df <- temp_df[, colnames(temp_df) %in% c(celltypes,'key')]
which(is.na(temp_df), arr.ind=TRUE)

counts_match <- match(colData(spg)$key, temp_df$key)
RCTD <-temp_df[counts_match,]
which(is.na(RCTD), arr.ind=TRUE)
RCTD = RCTD[,c("key",celltypes)]

RCTD$sample_id <-sapply(strsplit(RCTD$key,"-1_"), `[`, 2)
which(is.na(RCTD), arr.ind=TRUE)
#counts$tool = tool
#RCTD <- counts |> melt(id.vars = "key", variable.name = "celltype",value.name = "count")

# colnames(temp_df)[colnames(temp_df) %in% celltypes] <- paste0(tool, "_", colnames(temp_df)[colnames(temp_df) %in% celltypes])
# counts_match <- match(colData(spg)$key, temp_df$key)
# counts_info <-temp_df[counts_match,]
# counts_info <- subset(counts_info, select = -key)
# colData(spg) <- cbind(colData(spg), counts_info)

#columns_to_compare <- c("oligo", "other", "neuron", "microglia", "astrocyte")

RCTD_list <- split(RCTD[celltypes], RCTD$sample_id)
CART_list <- split(CART[celltypes], CART$sample_id)
tangram_list <- split(tangram[celltypes], tangram$sample_id)
cell2location_list <- split(cell2location[celltypes], cell2location$sample_id)

# Define a function to calculate RMSE
calculate_rmse <- function(observed, predicted) {signif(sqrt(mean((observed - predicted)^2)),3)}

# Calculate correlations and RMSE for each sample
results <- lapply(names(RCTD_list), function(sample_id) {
  corr_RCTD_CART <- cor(RCTD_list[[sample_id]], CART_list[[sample_id]])
  corr_tangram_CART <- cor(tangram_list[[sample_id]], CART_list[[sample_id]])
  corr_cell2location_CART <- cor(cell2location_list[[sample_id]], CART_list[[sample_id]])
  
  rmse_RCTD_CART <- sapply(celltypes, function(col) {calculate_rmse(RCTD_list[[sample_id]][[col]], CART_list[[sample_id]][[col]])})
  rmse_tangram_CART <- sapply(celltypes, function(col) {calculate_rmse(tangram_list[[sample_id]][[col]], CART_list[[sample_id]][[col]])})
  rmse_cell2location_CART <- sapply(celltypes, function(col) {calculate_rmse(cell2location_list[[sample_id]][[col]], CART_list[[sample_id]][[col]])})
  #sapply(celltypes, function(col) {calculate_rmse(RCTD[[col]], CART[[col]])})
  
  data.frame(
    Sample_ID = sample_id,
    celltypes = celltypes,
    Corr_RCTD_CART = round(diag(corr_RCTD_CART),2),
    Corr_tangram_CART = round(diag(corr_tangram_CART),2),
    Corr_cell2location_CART = round(diag(corr_cell2location_CART),2),
    RMSE_RCTD_CART = rmse_RCTD_CART,
    RMSE_tangram_CART = rmse_tangram_CART,
    RMSE_cell2location_CART = rmse_cell2location_CART
  )
})

# Combine the results into a single data frame
results_df <- do.call(rbind, results)

# Load the ggplot2 library
library(ggplot2)
library(gridExtra)

results_df = results_df[which(results_df$celltypes != "other"),]
plot_list = list()
# Create a scatterplot of Correlation vs. RMSE with sample_id as shape
plot_list[[1]] = ggplot(results_df, aes(x = Corr_RCTD_CART, y = RMSE_RCTD_CART, color = Sample_ID, shape = celltypes))+geom_point(size = 3) +
  labs(x = "Correlation (RCTD vs. CART)",y = "RMSE (RCTD vs. CART)", title = paste0("corr = ",mean(results_df$Corr_RCTD_CART), " & RMSE = ",  mean(results_df$RMSE_RCTD_CART))) 
plot_list[[2]] = ggplot(results_df, aes(x = Corr_tangram_CART, y = RMSE_tangram_CART, color = Sample_ID, shape = celltypes))+geom_point(size = 3) +
  labs(x = "Correlation (Tangram vs. CART)",y = "RMSE (Tangram vs. CART)", title = paste0("corr = ",mean(results_df$Corr_tangram_CART), " & RMSE = ",  mean(results_df$RMSE_tangram_CART)))
plot_list[[3]] = ggplot(results_df, aes(x = Corr_cell2location_CART, y = RMSE_cell2location_CART, color = Sample_ID, shape = celltypes))+geom_point(size = 3) +
  labs(x = "Correlation (Cell2location vs. CART)",y = "RMSE (Cell2location vs. CART)", title = paste0("corr = ",mean(results_df$Corr_cell2location_CART), " & RMSE = ",  mean(results_df$RMSE_cell2location_CART)))

#gridplot = grid.arrange(grobs = plot_list, nrow = length(celltypes))
ggsave(here("plots","spot_deconvo","shared_utilities","cartVStool.pdf"), plot = marrangeGrob(plot_list, nrow=1, ncol=1),  width = 10, height = 8)

