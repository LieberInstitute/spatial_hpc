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
sce = readRDS(here(Dr,"sce.rds"), verbose = TRUE)
spe = readRDS(here(Dr,"spe.rds"), verbose = TRUE)

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

