
#cd /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/
suppressPackageStartupMessages(library("here"))
#install.packages("SpatialExperiment")
suppressPackageStartupMessages(library("SpatialExperiment"))
suppressPackageStartupMessages(library("spatialLIBD"))
suppressPackageStartupMessages(library("rtracklayer"))
suppressPackageStartupMessages(library("lobstr"))
suppressPackageStartupMessages(library("sessioninfo"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("gridExtra"))
  
load(file=here::here("processed-data","pilot_data_checks","spe_Basic.Rdata"))

# Rotations by sample:
# V10B01-085_A1 = 90
# V10B01-085_B1 = 180
# V10B01-085_C1 = 180
# V10B01-085_D1 = 180
# V10B01-086_A1 = 0
# V10B01-086_B1 = 0
# V10B01-086_C1 = 0
# V10B01-086_D1 = 180
  
angle_list = c(180, 180, 180, 0, 0, 0, 180)
sampleID = c('V10B01-085_B1', 'V10B01-085_C1', 'V10B01-085_D1', 'V10B01-086_A1', 'V10B01-086_B1', 'V10B01-086_C1',  'V10B01-086_D1')
source(file = here::here("code","pilot_data_checks","transform_spe.R"))

pdf(file = here::here("plots", "pilot_data_checks", "ReferenceMapping.pdf"), h = 10, w = 20)

id = 'V10B01-085_A1'
spe = trans_geom(spe_basic, sample_id = id , degrees = 90)
spe = rotateImg(spe, sample_id = id, image_id = TRUE, degrees = 90)
orig_spe = spe_basic[, (colData(spe_basic)$in_tissue & colData(spe_basic)$sample_id == id)]
xnew = spe[, (colData(spe)$in_tissue & colData(spe)$sample_id == id)]

p1=plot_spatial_coords(orig_spe, paste0(id," original"))
p2=plot_spatial_coords(xnew, paste0(id," new"))
grid.arrange(p1,p2,ncol=2)

par(mfrow = c(1, 2))
plot(imgRaster(spe_basic[, colData(spe_basic)$sample_id == id]), main= 'Original')
plot(imgRaster(spe[, colData(spe)$sample_id == id]), main = 'New')

for (i in seq_along(angle_list)){
id = sampleID[i]
x = trans_geom(spe_basic, sample_id = id, degrees = angle_list[i])
spe = cbind(spe,x)
spe = rotateImg(spe, sample_id = id, image_id = TRUE, degrees = angle_list[i])

orig_spe = spe_basic[, (colData(spe_basic)$in_tissue & colData(spe_basic)$sample_id == id)]
xnew = spe[, (colData(spe)$in_tissue & colData(spe)$sample_id == id)]

p1 = plot_spatial_coords(orig_spe, 'Original Coords')
p2 = plot_spatial_coords(xnew, 'New Coords')
grid.arrange(p1,p2,ncol=2)

par(mfrow = c(1, 2))
plot(imgRaster(spe_basic[, colData(spe_basic)$sample_id == id]), main= 'Original')
plot(imgRaster(spe[, colData(spe)$sample_id == id]), main = 'New')
}

dev.off()

# ## Remove genes with no data
# no_expr <- which(rowSums(counts(spe)) == 0)
# length(no_expr)
# # [1] 8968
# length(no_expr) / nrow(spe) * 100
# # [1] 24.50206
# spe <- spe[-no_expr, ]
# 
# 
# ## Now drop the spots outside the tissue
# spe <- spe_raw[, spatialData(spe_raw)$in_tissue]
# 
# ## Size in Gb
# lobstr::obj_size(spe) / 1024 ^ 3
# # 1.305857
# dim(spe)
# # [1] 27633 28871
# 
# ## Remove spots without counts
# if (any(colSums(counts(spe)) == 0)) {
#   message("removing spots without counts for spe")
#   spe <- spe[, -which(colSums(counts(spe)) == 0)]
#   dim(spe)
# }
# 
# lobstr::obj_size(spe) / 1024 ^ 3
# # 1.305857
# dim(spe)
# # [1] 27633 28871
# 
# save(spe, file = here::here("processed-data", "pilot_data_checks", "spe.Rdata"))
# 
# ## Reproducibility information
# print('Reproducibility information:')
# Sys.time()
# proc.time()
# options(width = 120)
# session_info()
