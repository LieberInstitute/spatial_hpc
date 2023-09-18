setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")

suppressPackageStartupMessages(library("sessioninfo"))
suppressPackageStartupMessages(library("here"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("spacexr"))
suppressPackageStartupMessages(library("spatialLIBD"))
library(gridExtra)

#  Paths
cell_group = "broad"
subtype = "_class"
cell_type_var = 'broad.class'
#cell_types = unique(colData(sce)$broad.class)

# cell_group = "layer"
# subtype = "_celltype_class1_noHATAGABAAmy"
# cell_type_var = 'cell.class'
# cell_types = colData(sce)$cell.class
Ncol = 4

Dr <- here("processed-data","spot_deconvo","RCTD","2ndRun_newClass",cell_group)
plots = here("plots","spot_deconvo","RCTD","2ndRun_newClass",cell_group)
#   Load objects
spaceranger_dirs = read.csv(file.path(here::here("code","VistoSeg","code","samples.txt")), header = FALSE, sep = '\t', stringsAsFactors = FALSE, col.names = c('SPpath','sample_id','brain'))
spaceranger_dirs = spaceranger_dirs[1:36,]
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

