setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")

suppressPackageStartupMessages(library("SingleCellExperiment"))
suppressPackageStartupMessages(library("SpatialExperiment"))
suppressPackageStartupMessages(library("sessioninfo"))
suppressPackageStartupMessages(library("here"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("spacexr"))

#  Paths
 Dr <- here("processed-data","spot_deconvo","shared_utilities")
# cell_group = "broad"
# subtype = "_class"
# cell_type_var = 'broad.class'

cell_group = "layer"
subtype = "_celltype_class1_noHATAGABAAmy"
cell_type_var = 'cell.class'

# processed_out = here("processed-data","spot_deconvo","RCTD","2ndRun_newClass_RCTDmarkers",cell_group)
processed_out = here("processed-data","spot_deconvo","RCTD","3rdRun_newClass_deconvoMarkers",cell_group)

#   Load objects
 sce = readRDS(here(Dr,"sce_class.rds"))
# sce = readRDS(here(Dr,"sce_class1_noHATAGABAAmy.rds"))

markers = readLines(here(Dr,"markers_broad_class.txt"))
#markers = readLines(here(Dr,"markers_sce_class1_noHATAGABAAmy.txt"))
sceb = sce[rowData(sce)$gene_id %in% markers,]
sce = sceb

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

#layer
#[1] 36601 68528
table(reference@cell_types) #number of occurences for each cell type
# ExcN     InhN     Glia   Immune      CSF Vascular 
# 10000    10000    10000     3940     3587     3137 

# Astro           CSF          ExcN          InhN Micro_Macro_T 
# 4298          3587         10000         10000          3940 
# Oligo           OPC      Vascular 
# 6787          1662          3137 

# layer
# Astro      CA1_ProS         CA2-4         Cajal       Choroid 
# 4298          4493          5319            40          3166 
# Ependy      GABA.CGE    GABA.LAMP5      GABA.MGE            GC 
# 421          2573          3547          3776         10000 
# L2_3.Prs.Ent  L2_3.PrS.PaS            L5         L6_6b Micro_Macro_T 
# 4899          2729          2059          3756          3940 
# Oligo           OPC         Sub.1         Sub.2          Thal 
# 6787          1662           994           897            35 
# Vascular 
# 3137 

## Save RDS object (optional)
saveRDS(reference, here(processed_out,'SCRef.rds'))
