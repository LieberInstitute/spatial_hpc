setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")

suppressPackageStartupMessages(library("SingleCellExperiment"))
suppressPackageStartupMessages(library("SpatialExperiment"))
suppressPackageStartupMessages(library("sessioninfo"))
suppressPackageStartupMessages(library("here"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("spacexr"))

#  Paths
Dr <- here("processed-data","spot_deconvo","shared_utilities")
cell_group = "broad"
subtype = "_class"
cell_type_var = 'broad.class'

# cell_group = "layer"
# subtype = "_celltype_class1_noHATAGABAAmy"
# cell_type_var = 'cell.class'

processed_out = here("processed-data","spot_deconvo","RCTD","2ndRun_newClass")
#   Load objects
sce = readRDS(here(Dr,"sce_class.rds"))
# sce = readRDS(here(Dr,"sce_class_noHATAGABA.rds"))

cell_types = colData(sce)$broad.class
# cell_types = colData(sce)$cell.class

# reference data
counts = assays(sce)$counts
#meta_data = data.frame(barcode = unique(rownames(colData(sce))), cluster = colData(sce)$broad.class, nUMI = colData(sce)$sum)
#cell_types <- meta_data$cluster; 
#nUMI <- meta_data$nUMI;
#names(cell_types) <- meta_data$barcode # create cell_types named list
#names(nUMI) <- meta_data$barcode # create nUMI named list

names(cell_types) = unique(rownames(colData(sce)))
nUMI = colData(sce)$sum

cell_types <- as.factor(cell_types) # convert to factor data type
names(nUMI) <- unique(rownames(colData(sce)))

reference <- Reference(counts, cell_types, nUMI)

print(dim(reference@counts)) #observe Digital Gene Expression matrix
#[1] 36601 40664
#[1] 36601 43411
table(reference@cell_types) #number of occurences for each cell type
# ExcN     InhN     Glia   Immune      CSF Vascular 
# 10000    10000    10000     3940     3587     3137 

# Astro           CSF          ExcN          InhN Micro_Macro_T 
# 4298          3587         10000         10000          3940 
# Oligo           OPC      Vascular 
# 6787          1662          3137 

## Save RDS object (optional)
saveRDS(reference, here(processed_out,'SCRef.rds'))
