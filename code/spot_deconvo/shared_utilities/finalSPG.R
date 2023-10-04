setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")
library("here")
library("jaffelab")
library("SpatialExperiment")
library("sessioninfo")
library("tidyverse")

cell_groups <- c("broad", "layer")
deconvo_tools <- c("tangram", "cell2location", "RCTD")

spg_in <- here("processed-data", "spot_deconvo", "shared_utilities", "spg.rds")
spg_out <- here("processed-data", "spot_deconvo", "shared_utilities", "spg_final.rds")

spg = readRDS(spg_in)
# unqiue(spg$brnum)
# Br3942_VSPG Br8325_VSPG
colData(spg)$key = gsub("_Br3942", "", colData(spg)$key)
colData(spg)$key = gsub("_Br8325", "", colData(spg)$key)
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
counts = temp_df %>% select("key")
for (celltype in celltypes) {
  counts[,celltype] = temp_df[,celltype]/temp_df$n_cells
}
counts$tool = tool
CART_df_long <- counts |> melt(id.vars = "key", variable.name = "celltype",value.name = "count")

# colnames(counts)[colnames(counts) %in% celltypes] <- paste0(tool, "_", colnames(counts)[colnames(counts) %in% celltypes])
# counts_match <- match(colData(spg)$key, counts$key)
# counts_info <-counts[counts_match,]
# counts_info <- subset(counts_info, select = -key)
# colData(spg) <- cbind(colData(spg), counts_info)

## cell2location
tool =  "cell2location"
Dr = here("processed-data", "spot_deconvo", tool, "IF", "2ndRun_newClass", "broad")
csv_files = data.frame(files = list.files(Dr, pattern = "clusters.csv" , recursive = TRUE, full.names = TRUE))
csv_files$sample_id = sapply(strsplit(csv_files$files,"/"), `[`, 13)
counts_list <- list()
for (file_path in csv_files$files) {
  data <- read.csv(file_path, row.names = NULL)  # Read the CSV file
  counts_list[[file_path]] <- data  # Store the data frame in the list with the file path as the key
}

temp_df <- do.call(rbind, counts_list)
rownames(temp_df) <- NULL
temp_df$key = gsub("_Br3942", "", temp_df$key)
temp_df$key = gsub("_Br8325", "", temp_df$key)
temp_df$neuron = temp_df$ExcN + temp_df$InhN
temp_df$other = temp_df$CSF + temp_df$OPC + temp_df$Vascular
temp_df$microglia = temp_df$Micro_Macro_T
temp_df$astrocyte = temp_df$Astro
temp_df$oligo = temp_df$Oligo
temp_df <- temp_df[, colnames(temp_df) %in% c(celltypes,'key')]
temp_df$count = rowSums(temp_df[,celltypes])

counts = temp_df %>% select("key")
for (celltype in celltypes) {
  counts[,celltype] = temp_df[,celltype]/temp_df$count
}

counts$tool = tool
cell2location <- counts |> melt(id.vars = c("key","tool"), variable.name = "celltype",value.name = "count")

# colnames(counts)[colnames(counts) %in% celltypes] <- paste0(tool, "_", colnames(counts)[colnames(counts) %in% celltypes])
# counts_match <- match(colData(spg)$key, counts$key)
# counts_info <-counts[counts_match,]
# counts_info <- subset(counts_info, select = -key)
# colData(spg) <- cbind(colData(spg), counts_info)

## RCTD
tool =  "RCTD"
Dr = here("processed-data", "spot_deconvo", tool, "2ndRun_newClass_RCTDmarkers", "broad")
csv_files = data.frame(files = list.files(Dr, pattern = "clusters.csv" , recursive = TRUE, full.names = TRUE))
csv_files = csv_files[30:35,]
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

counts$tool = tool
RCTD_df_long <- counts |> melt(id.vars = "key", variable.name = "celltype",value.name = "count")

# colnames(temp_df)[colnames(temp_df) %in% celltypes] <- paste0(tool, "_", colnames(temp_df)[colnames(temp_df) %in% celltypes])
# counts_match <- match(colData(spg)$key, temp_df$key)
# counts_info <-temp_df[counts_match,]
# counts_info <- subset(counts_info, select = -key)
# colData(spg) <- cbind(colData(spg), counts_info)
df = rbind(CART_df,cell2location_df,RCTD_df)
df = subset(colData(spg), select = c(1,15,88:102))
df <- na.omit(df)
which(is.na(df), arr.ind=TRUE)
temp = list()
for (sample in unique(df$sample_id)) {
  data = data.frame(celltype = c("oligo", "neuron", "microglia", "astrocyte"),  RMSE_RCTD = rep(NA, 4), corr_RCTD = rep(NA, 4), RMSE_cell2location = rep(NA, 4),corr_cell2location = rep(NA, 4), sample = rep(sample,4))
  dat = df[which(df$sample_id == sample),]
  data[1,"RMSE_RCTD"] = sqrt(mean((dat$CART_oligo - dat$RCTD_oligo)^2))
  data[1,"corr_RCTD"] = round(cor(dat$RCTD_oligo,dat$CART_oligo),2)
  data[4,"RMSE_RCTD"] = sqrt(mean((dat$CART_astrocyte - dat$RCTD_astrocyte)^2))
  data[4,"corr_RCTD"] = cor(dat$CART_astrocyte,dat$RCTD_astrocyte)
  data[2,"RMSE_RCTD"] = sqrt(mean((dat$CART_neuron - dat$RCTD_neuron)^2))
  data[2,"corr_RCTD"] = cor(dat$CART_neuron,dat$RCTD_neuron)
  data[3,"RMSE_RCTD"] = sqrt(mean((dat$CART_oligo - dat$RCTD_oligo)^2))
  data[3,"corr_RCTD"] = cor(dat$CART_oligo,dat$RCTD_oligo)
  
  data[1,"RMSE_cell2location"] = sqrt(mean((dat$CART_oligo - dat$cell2location_oligo)^2))
  data[1,"corr_cell2location"] = cor(dat$CART_oligo, dat$cell2location_oligo)
  data[4,"RMSE_cell2location"] = sqrt(mean((dat$CART_astrocyte - dat$cell2location_astrocyte)^2))
  data[4,"corr_cell2location"] = cor(dat$CART_astrocyte,dat$cell2location_astrocyte)
  data[2,"RMSE_cell2location"] = sqrt(mean((dat$CART_neuron - dat$cell2location_neuron)^2))
  data[2,"corr_cell2location"] = cor(dat$CART_neuron, dat$cell2location_neuron)
  data[3,"RMSE_cell2location"] = sqrt(mean((dat$CART_microglia - dat$cell2location_microglia)^2))
  data[3,"corr_cell2location"] = cor(dat$CART_microglia, dat$cell2location_microglia)
  temp[[sample]] = data
}

temp_df <- do.call(rbind, temp)


metrics_df <- df |>
  group_by(tool) |>
  summarize(
    corr = round(cor(CART, actual), 2),
    rmse = signif(mean((observed - actual)**2)**0.5, 3)
  ) |>
  ungroup()


ggplot(temp_df,aes( x = RMSE_RCTD, y = corr_RCTD, color = celltype, shape = sample))+geom_point()


across_spots <- function(count_df, plot_name, x_angle = 0) {
  #   Compute metrics for each deconvolution tool: correlation between
  #   observed and actual values as well as RMSE
  metrics_df <- count_df |>
    group_by(deconvo_tool) |>
    summarize(
      corr = round(cor(observed, actual), 2),
      rmse = signif(mean((observed - actual)**2)**0.5, 3)
    ) |>
    ungroup()
  
  #   Add KL divergence from observed section-wide cell-type proportions to
  #   ground-truth proportions, averaged across sections
  metrics_df <- count_df |>
    #   Ensure counts or proportions are normalized to add to 1 across
    #   cell types
    group_by(sample_id, deconvo_tool) |>
    mutate(
      observed = observed / sum(observed),
      actual = actual / sum(actual),
    ) |>
    #   Compute each term in the sum for KL divergence
    group_by(sample_id, deconvo_tool, cell_type) |>
    summarize(kl_piece = actual * log(actual / observed)) |>
    #   Add all terms to form the sum for each sample
    group_by(sample_id, deconvo_tool) |>
    summarize(kl = sum(kl_piece)) |>
    #   Take the mean across samples to form one value per tool
    group_by(deconvo_tool) |>
    summarize(kl = round(mean(kl), 2)) |>
    left_join(metrics_df)
  
  #   Improve labels for plotting
  metrics_df$corr <- paste("Cor =", metrics_df$corr)
  metrics_df$rmse <- paste("RMSE =", metrics_df$rmse)
  metrics_df$kl <- paste("KL Div. =", metrics_df$kl)
  
  p <- ggplot(count_df) +
    geom_point(
      aes(x = observed, y = actual, shape = sample_id, color = cell_type)
    ) +
    facet_wrap(~deconvo_tool) +
    geom_abline(
      intercept = 0, slope = 1, linetype = "dashed", color = "red"
    ) +
    #   Correlation label
    geom_text(
      data = metrics_df,
      mapping = aes(
        x = max(count_df$observed),
        y = 0.05 * max(count_df$actual),
        label = corr
      ),
      hjust = 1, vjust = 0
    ) +
    #   RMSE label
    geom_text(
      data = metrics_df,
      mapping = aes(
        x = max(count_df$observed),
        y = 0.15 * max(count_df$actual),
        label = rmse
      ),
      hjust = 1, vjust = 0
    ) +
    #   KL divergence label
    geom_text(
      data = metrics_df,
      mapping = aes(
        x = max(count_df$observed),
        y = 0.25 * max(count_df$actual),
        label = kl
      ),
      hjust = 1, vjust = 0
    ) +
    scale_color_manual(values = cell_type_labels) +
    scale_x_continuous(
      limits = c(0, max(count_df$observed) * 1.05),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      limits = c(0, max(count_df$actual) * 1.05),
      expand = c(0, 0)
    ) +
    labs(
      x = "Software-Estimated", y = "CART-Calculated",
      color = "Cell Type", shape = "Sample ID"
    ) +
    theme_bw(base_size = 15) +
    theme(axis.text.x = element_text(angle = x_angle))
  
  pdf(file.path(plot_dir, plot_name), height = 4, width = 10)
  print(p)
  dev.off()
}


observed_df_long <- observed_df |>
  melt(
    id.vars = added_colnames, variable.name = "cell_type",
    value.name = "count"
  ) |>
  pivot_wider(
    names_from = obs_type, values_from = count,
  ) |>
  #   Use an ordered factor for plotting cell types
  mutate(cell_type = factor(cell_type, levels = cell_types, order = TRUE))
