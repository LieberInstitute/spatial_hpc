setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")

suppressPackageStartupMessages(library("SingleCellExperiment"))
suppressPackageStartupMessages(library("SpatialExperiment"))
suppressPackageStartupMessages(library("sessioninfo"))
suppressPackageStartupMessages(library("here"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("spacexr"))

#  Paths
Dr <- here("processed-data","spot_deconvo","shared_utilities")
#cell_group = "broad"
cell_group = "layer"

# processed_out = here("processed-data","spot_deconvo","RCTD","2ndRun_newClass_RCTDmarkers")
# plots_out = here("plots","spot_deconvo","RCTD","2ndRun_newClass_RCTDmarkers")

processed_out = here("processed-data","spot_deconvo","RCTD","3rdRun_newClass_deconvoMarkers")
plots_out = here("plots","spot_deconvo","RCTD","3rdRun_newClass_deconvoMarkers")

spaceranger_dirs = read.csv(file.path(here::here("code","VistoSeg","code","samples.txt")), header = FALSE, sep = '\t', stringsAsFactors = FALSE, col.names = c('SPpath','sample_id','brain'))
# spaceranger_dirs = spaceranger_dirs[1:36,]#for HE
spaceranger_dirs = spaceranger_dirs[37:44,]#for IF
sample_ids = spaceranger_dirs$sample_id

sample_id = sample_ids[as.numeric(Sys.getenv("SGE_TASK_ID"))]
print(sample_id)

#spe = readRDS(here(Dr,"spe.rds"))
spg = readRDS(here(Dr,"spg.rds")) #for IF
# markers = readLines(here(Dr,"markers_broad_class.txt"))
# speb = spe[rowData(spe)$gene_id %in% markers,spe$sample_id==sample_id]

#speb = spe[,which(spe$sample_id==sample_id)]
speb = spg[,which(spg$sample_id==sample_id)] #for IF

counts = assays(speb)$counts
rownames(counts) = rowData(speb)$gene_id
coords = as.data.frame(spatialCoords(speb))
colnames(coords) = c("ycoord","xcoord")
nUMI <- colSums(counts) # In this case, total counts per pixel is nUMI

### Create SpatialRNA object
puck <- SpatialRNA(coords, counts, nUMI)

### plotting
barcodes <- colnames(puck@counts) #
dir.create(here(plots_out, cell_group, sample_id))
pdf(here(plots_out, cell_group, sample_id, "UMIcount.pdf"), width = 10, height = 10)
plot_puck_continuous(puck, barcodes, puck@nUMI, ylimit = c(0,round(quantile(puck@nUMI,0.9))), title ='plot of nUMI')
reference = readRDS(here(processed_out,cell_group,'SCRef.rds'))
myRCTD <- create.RCTD(puck, reference, max_cores = 1, MAX_MULTI_TYPES = 5, UMI_min = 2)
#myRCTD <- create.RCTD(puck, reference, max_cores = 1, MAX_MULTI_TYPES = 5)

dir.create(here(processed_out,cell_group,sample_id))
saveRDS(myRCTD, here(processed_out,cell_group,sample_id,'myRCTD.rds'))