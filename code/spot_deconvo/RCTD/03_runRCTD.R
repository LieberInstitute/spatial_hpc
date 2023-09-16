setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")

suppressPackageStartupMessages(library("sessioninfo"))
suppressPackageStartupMessages(library("here"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("spacexr"))

#  Paths
Dr <- here("processed-data","spot_deconvo","RCTD","2ndRun_newClass")
cell_group = "broad"
subtype = "_class"
cell_type_var = 'broad.class'
cell_types = colData(sce)$broad.class

# cell_group = "layer"
# subtype = "_celltype_class1_noHATAGABAAmy"
# cell_type_var = 'cell.class'
# cell_types = colData(sce)$cell.class

out = here("processed-data","spot_deconvo","RCTD","2ndRun_newClass")

#   Load objects
sce = readRDS(here(Dr,"sce_class.rds"))
# sce = readRDS(here(Dr,"sce_class_noHATAGABA.rds"))
spe = readRDS(here(Dr,"spe.rds"))

myRCTD <- run.RCTD(myRCTD, doublet_mode = 'multi', MAX_MULTI_TYPES = 5)

results <- myRCTD@results
norm_weights = normalize_weights(results$weights) 
cell_type_names <- myRCTD@cell_type_info$info[[2]] #list of cell type names
spatialRNA <- myRCTD@spatialRNA
resultsdir <- 'RCTD_Plots' ## you may change this to a more accessible directory on your computer.
dir.create(resultsdir)

plot_all_cell_types(results$results_df, spatialRNA@coords, cell_type_names, resultsdir) 

