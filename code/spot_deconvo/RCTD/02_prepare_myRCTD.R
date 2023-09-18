setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")

suppressPackageStartupMessages(library("SingleCellExperiment"))
suppressPackageStartupMessages(library("SpatialExperiment"))
suppressPackageStartupMessages(library("sessioninfo"))
suppressPackageStartupMessages(library("here"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("spacexr"))

#  Paths
Dr <- here("processed-data","spot_deconvo","shared_utilities")
processed_out = here("processed-data","spot_deconvo","RCTD","2ndRun_newClass")
plots_out = here("plots","spot_deconvo","RCTD","2ndRun_newClass")
cell_group = "broad"

spaceranger_dirs = read.csv(file.path(here::here("code","VistoSeg","code","samples.txt")), header = FALSE, sep = '\t', stringsAsFactors = FALSE, col.names = c('SPpath','sample_id','brain'))
spaceranger_dirs = spaceranger_dirs[1:36,]
sample_ids = spaceranger_dirs$sample_id

sample_id = sample_ids[as.numeric(Sys.getenv('TASK_ID'))]
print(sample_id)

# spe = readRDS(here(Dr,"spe.rds"))
# speb = spe[,which(spe$sample_id==sample_id)]
# 
# counts = assays(speb)$counts
# rownames(counts) = rowData(speb)$gene_id
# coords = as.data.frame(spatialCoords(speb))
# colnames(coords) = c("ycoord","xcoord")
# nUMI <- colSums(counts) # In this case, total counts per pixel is nUMI
# 
# ### Create SpatialRNA object
# puck <- SpatialRNA(coords, counts, nUMI)
# 
# ### plotting
# barcodes <- colnames(puck@counts) # 
# dir.create(here(plots_out, cell_group, sample_id))
# pdf(here(plots_out, cell_group, sample_id, "UMIcount.pdf"), width = 10, height = 10)
# plot_puck_continuous(puck, barcodes, puck@nUMI, ylimit = c(0,round(quantile(puck@nUMI,0.9))), title ='plot of nUMI') 
# 
# myRCTD <- create.RCTD(puck, reference, max_cores = 1, MAX_MULTI_TYPES = 5)
# 
# dir.create(here(processed_out,cell_group,sample_id))
# saveRDS(myRCTD, here(processed_out,cell_group,sample_id,'myRCTD.rds'))