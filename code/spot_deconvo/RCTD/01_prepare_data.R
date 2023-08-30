setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")

suppressPackageStartupMessages(library("SingleCellExperiment"))
suppressPackageStartupMessages(library("SpatialExperiment"))
suppressPackageStartupMessages(library("sessioninfo"))
suppressPackageStartupMessages(library("here"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("spacexr"))
  
#  Paths
Dr <- here("processed-data","spot_deconvo","shared_utilities")

#   Load objects
sce = readRDS(here(Dr,"sce.rds"))
spe = readRDS(here(Dr,"spe.rds"))

# reference data
counts = assays(sce)$counts
meta_data = data.frame(barcode = unique(rownames(colData(sce))), cluster = colData(sce)$broad.type, nUMI = colData(sce)$sum)
cell_types <- meta_data$cluster; names(cell_types) <- meta_data$barcode # create cell_types named list
cell_types <- as.factor(cell_types) # convert to factor data type
nUMI <- meta_data$nUMI; names(nUMI) <- meta_data$barcode # create nUMI named list

reference <- Reference(counts, cell_types, nUMI)

print(dim(reference@counts)) #observe Digital Gene Expression matrix
#[1] 36601 40664
table(reference@cell_types) #number of occurences for each cell type
# ExcN     InhN     Glia   Immune      CSF Vascular 
# 10000    10000    10000     3940     3587     3137 

## Save RDS object (optional)
saveRDS(reference, file.path(refdir,'SCRef.rds'))

## build puck 
speb = spe[,which(spe$sample_id=="V11L05-333_D1")]

counts = assays(speb)$counts
rownames(counts) = rowData(speb)$gene_id
coords = as.data.frame(spatialCoords(speb))
colnames(coords) = c("ycoord","xcoord")
nUMI <- colSums(counts) # In this case, total counts per pixel is nUMI

### Create SpatialRNA object
puck <- SpatialRNA(coords, counts, nUMI)

### plotting
barcodes <- colnames(puck@counts) # 
plot_puck_continuous(puck, barcodes, puck@nUMI, ylimit = c(0,round(quantile(puck@nUMI,0.9))), title ='plot of nUMI') 

myRCTD <- create.RCTD(puck, reference, max_cores = 1)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'multi')

results <- myRCTD@results
norm_weights = normalize_weights(results$weights) 
cell_type_names <- myRCTD@cell_type_info$info[[2]] #list of cell type names
spatialRNA <- myRCTD@spatialRNA
resultsdir <- 'RCTD_Plots' ## you may change this to a more accessible directory on your computer.
dir.create(resultsdir)

plot_all_cell_types(results$results_df, spatialRNA@coords, cell_type_names, resultsdir) 

