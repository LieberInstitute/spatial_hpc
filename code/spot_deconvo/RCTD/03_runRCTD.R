setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")

suppressPackageStartupMessages(library("sessioninfo"))
suppressPackageStartupMessages(library("here"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("spacexr"))
suppressPackageStartupMessages(library("spatialLIBD"))
library('gridExtra')
library('viridis')

#  Paths
cell_group = "broad"
subtype = "_class"
cell_type_var = 'broad.class'
#cell_types = unique(colData(sce)$broad.class)

# cell_group = "layer"
# subtype = "_celltype_class1_noHATAGABAAmy"
# cell_type_var = 'cell.class'
#cell_types = colData(sce)$cell.class
Ncol = 4

Dr <- here("processed-data","spot_deconvo","RCTD","2ndRun_newClass_RCTDmarkers", cell_group)
plots = here("plots","spot_deconvo","RCTD","2ndRun_newClass_RCTDmarkers", cell_group)
#   Load objects
spaceranger_dirs = read.csv(file.path(here::here("code","VistoSeg","code","samples.txt")), header = FALSE, sep = '\t', stringsAsFactors = FALSE, col.names = c('SPpath','sample_id','brain'))
#spaceranger_dirs = spaceranger_dirs[1:36,]
spaceranger_dirs = spaceranger_dirs[37:44,]
sample_ids = spaceranger_dirs$sample_id

sample_id = sample_ids[as.numeric(Sys.getenv("SGE_TASK_ID"))]
myRCTD = readRDS(here(Dr,sample_id,"myRCTD.rds"))
celltypes = levels(myRCTD@reference@cell_types)

# results <- myRCTD@results
# norm_weights = normalize_weights(results$weights)
# cell_type_names <- myRCTD@cell_type_info$info[[2]] #list of cell type names
# spatialRNA <- myRCTD@spatialRNA
# resultsdir <- 'RCTD_Plots' ## you may change this to a more accessible directory on your computer.
# dir.create(resultsdir)
# 
# plot_all_cell_types(results$results_df, spatialRNA@coords, cell_type_names, resultsdir)
# 

myRCTD1 <- run.RCTD(myRCTD, doublet_mode = 'multi')
saveRDS(myRCTD1, here(Dr,sample_id,paste0(sample_id,"_multi.rds")))

## multi mode ## 
print('Examining multi mode all weights results')

results = myRCTD1@results
weights = lapply(results, function(x) x$all_weights)
weights_df <- data.frame(do.call(rbind, weights))
colnames(weights_df)[colnames(weights_df) == "CA2.4"] = "CA2-4"
norm_weights <- normalize_weights(weights_df)
coords = myRCTD1@spatialRNA@coords

plot_list <- lapply(celltypes, function(i) {
  print(i)  
  ggplot(coords, aes(x = coords$x, y=coords$y, color = weights_df[,i])) + labs(title = i, x="", y="") + 
    geom_point(size = 0.5)+scale_color_gradientn(colours = viridis(10, option = "magma"), limits = c(0,1)) +
    scale_y_reverse()+ theme(legend.key.width = unit(0.1, "cm"))+labs(color = "")
})

gridplot = grid.arrange(grobs = plot_list, ncol = Ncol)
ggsave(here(plots,sample_id,"multi_allWeight.pdf"), plot = gridplot, width = 18, height = 8)
weights_df = as.data.frame(weights_df)
rownames(weights_df) = rownames(myRCTD1@spatialRNA@coords)
write.csv(weights_df, here(Dr, sample_id, "clusters_allWeights.csv"))

print('Examining multi mode sub weights results')
results = myRCTD1@results
weights = lapply(results, function(x) x$sub_weights)

df = weights[[1]]
missing_columns <- setdiff(celltypes, names(df))
if (length(missing_columns) > 0) {for (col_name in missing_columns) {df[[col_name]] <- 0 }}
df <- as.data.frame(t(df))
if (length(missing_columns) == length(celltypes)) {df = select(df, -c(setdiff(colnames(df), celltypes)))}

for (i in 2:length(weights)) {
  # Extract the list
  current_list <- weights[[i]]
  missing_columns <- setdiff(celltypes, names(current_list))
  if (length(missing_columns) > 0) {for (col_name in missing_columns) {current_list[[col_name]] <- 0 }}
  current_list = as.data.frame(t(current_list))
  if (length(missing_columns) == length(celltypes)) {current_list = select(current_list, -c(setdiff(colnames(current_list), celltypes)))}
  current_list = current_list[,colnames(df)]
  # Combine the current list with the combined dataframe by column names
  df <- rbind(df, current_list)
}


plot_list <- lapply(celltypes, function(i) {
  ggplot(coords, aes(x = coords$x, y=coords$y , color = df[,i])) + labs(title = i, x="", y="") + 
    geom_point(size = 0.5)+scale_color_gradientn(colours = viridis(10, option = "magma"), limits = c(0,1)) +
    scale_y_reverse()+ theme(legend.key.width = unit(0.1, "cm"))+labs(color = "")
})

gridplot = grid.arrange(grobs = plot_list, ncol = Ncol)
ggsave(here(plots,sample_id,"multi_subWeight.pdf"), plot = gridplot, width = 18, height = 8)
rownames(df) = rownames(myRCTD1@spatialRNA@coords)
write.csv(df, here(Dr, sample_id, "clusters_subWeights.csv"))

## full mode
myRCTD2 <- run.RCTD(myRCTD, doublet_mode = 'full')
saveRDS(myRCTD2, here(Dr,sample_id,paste0(sample_id,"_full.rds")))

barcodes <- colnames(myRCTD2@spatialRNA@counts)
weights <- myRCTD2@results$weights
norm_weights <- normalize_weights(weights)

plot_list <- lapply(celltypes, function(i) {
  ggplot(coords, aes(x = coords$x, y=coords$y , color = norm_weights[,i])) + labs(title = i, x="", y="") + 
    geom_point(size = 0.5)+scale_color_gradientn(colours = viridis(10, option = "magma"), limits = c(0,1)) +
    scale_y_reverse()+ theme(legend.key.width = unit(0.1, "cm"))+labs(color = "")
})

gridplot = grid.arrange(grobs = plot_list, ncol = Ncol)
ggsave(here(plots,sample_id,"full_weight.pdf"), plot = gridplot, width = 18, height = 8)

write.csv(as.data.frame(as.matrix(weights)), here(Dr, sample_id, "clusters.csv"))

