#   Rotations by sample:
# V10B01-085_A1 = 90
# V10B01-085_B1 = 180
# V10B01-085_C1 = 180
# V10B01-085_D1 = 180
# 
# V10B01-086_A1 = 0
# V10B01-086_B1 = 0
# V10B01-086_C1 = 0
# V10B01-086_D1 = 180

#   The function here uses clockwise as the positive direction, so the code


load('/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/pilot_data_checks/spe_raw.Rdata')
spe <- spe_raw[, spatialData(spe_raw)$in_tissue]

angle_list = c(90, 180, 180, 180, 180)
names(angle_list) = c('V10B01-085_A1', 'V10B01-085_C1', 'V10B01-085_D1', 'V10B01-086_D1')

id = 'V10B01-085_A1'
orig_spe = spe[, colData(spe)$sample_id == id]
plot_spatial_coords(orig_spe, 'Original Coords')
x = trans_geom(spe, sample_id = id, degrees = 90)
plot_spatial_coords(x, 'New Coords')

## Remove genes with no data
no_expr <- which(rowSums(counts(spe)) == 0)
length(no_expr)
# [1] 8968
length(no_expr) / nrow(spe) * 100
# [1] 24.50206
spe <- spe[-no_expr, ]


## Now drop the spots outside the tissue
spe <- spe_raw[, spatialData(spe_raw)$in_tissue]

## Size in Gb
lobstr::obj_size(spe) / 1024 ^ 3
# 1.305857
dim(spe)
# [1] 27633 28871

## Remove spots without counts
if (any(colSums(counts(spe)) == 0)) {
  message("removing spots without counts for spe")
  spe <- spe[, -which(colSums(counts(spe)) == 0)]
  dim(spe)
}

lobstr::obj_size(spe) / 1024 ^ 3
# 1.305857
dim(spe)
# [1] 27633 28871

save(spe, file = here::here("processed-data", "pilot_data_checks", "spe.Rdata"))

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
