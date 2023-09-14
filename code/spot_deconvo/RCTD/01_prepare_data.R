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
plots_out = here("plots","spot_deconvo","RCTD","2ndRun_newClass")
#   Load objects
sce = readRDS(here(Dr,"sce_class.rds"))
# sce = readRDS(here(Dr,"sce_class_noHATAGABA.rds"))
spe = readRDS(here(Dr,"spe.rds"))

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

## build puck 
sample_ids = unique(spe$sample_id)
marker_genes = read.delim(here(Dr,paste0("markers_",cell_group,subtype,".txt")), header = FALSE, sep = "\t")
marker_genes = as.character(marker_genes[,1])

for (sample_id in sample_ids){ 
speb = spe[rownames(spe) %in% marker_genes ,which(spe$sample_id==sample_id)]

counts = assays(speb)$counts
rownames(counts) = rowData(speb)$gene_id 
coords = as.data.frame(spatialCoords(speb))
colnames(coords) = c("ycoord","xcoord")
nUMI <- colSums(counts) # In this case, total counts per pixel is nUMI

### Create SpatialRNA object
puck <- SpatialRNA(coords, counts, nUMI)

### plotting
barcodes <- colnames(puck@counts) # 
p1 = plot_puck_continuous(puck, barcodes, puck@nUMI, ylimit = c(0,round(quantile(puck@nUMI,0.9))), title ='plot of nUMI', size = 1) 

pdf(file.path(plots_out, paste0("UMI",cell_group,name,".pdf")),width = 35, height = 35)
print(plot_list)
dev.off()

myRCTD <- create.RCTD(puck, reference, max_cores = 1)

mkdir(here(processed_out,cell_group,sample_id))
## Save RDS object 
saveRDS(myRCTD, here(processed_out,cell_group, sample_id, paste0(sample_id,'.rds')))

}


