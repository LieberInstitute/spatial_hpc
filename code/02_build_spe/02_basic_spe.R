
#cd /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/
suppressPackageStartupMessages(library("here"))
suppressPackageStartupMessages(library("SpatialExperiment"))
suppressPackageStartupMessages(library("spatialLIBD"))
suppressPackageStartupMessages(library("rtracklayer"))
suppressPackageStartupMessages(library("lobstr"))
suppressPackageStartupMessages(library("sessioninfo"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("gridExtra"))
suppressPackageStartupMessages(library("ggspavis"))

load(file=here::here("processed-data","02_build_spe","spe_raw.Rdata"))

# Rotations by sample:

# V10B01-085_A1= 90, TR
# V10B01-085_B1= 180, TL
# V10B01-085_C1= 180, BR
# V10B01-085_D1= 180, BL
# V10B01-086_A1= 0, TL
# V10B01-086_B1= 0, BL
# V10B01-086_C1= 0, TR
# V10B01-086_D1= 180, BR
# V11L05-333_A1= 0, TL
# V11L05-333_B1= 0, TR
# V11L05-333_C1= 0, BL
# V11L05-333_D1= 0, BR
# V11L05-335_A1= 0, BL
# V11L05-335_B1= 0, TR
# V11L05-335_C1= 0, TL
# V11L05-335_D1= 0, BR
# V11U08-084_A1= 0, TL
# V11U08-084_B1= 0, TR
# V11U08-084_C1= 0, BL
# V11U08-084_D1= 0, BR
# V11A20-297_A1= 0, BL
# V11A20-297_B1= 0, BR
# V11A20-297_C1= 0, TL
# V11A20-297_D1= 0, TR
# V11U08-081_A1= 0, TL
# V11U08-081_B1= 0, TR
# V11U08-081_C1= 0, BL
# V11U08-081_D1= 0, BR
# V11L05-336_A1= 0, BR
# V11L05-336_B1= 0, BL
# V11L05-336_C1= 0, TR
# V11L05-336_D1= 0, TL

sampleID = unique(spe_raw$sample_id)
angle_list = c(90, 180, 180, 180, 0, 0, 0, 180, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
source(file = here::here("code","pilot_data_checks","transform_spe.R"))

pdf(file = here::here("plots", "ReferenceMapping.pdf"), h = 10, w = 20)
for (i in seq_along(angle_list)){
  id = sampleID[i]
  x = trans_geom(spe_raw, sample_id = id , degrees = angle_list[i])
  
  if (i==1) {
    spe = x
  } else {
    spe = cbind(spe,x)
  }
  
  spe = rotateImg(spe, sample_id = id, image_id = TRUE, angle_list[i])
  orig_spe = spe_raw[, (colData(spe_raw)$in_tissue & colData(spe_raw)$sample_id == id)]
  xnew = spe[, (colData(spe)$in_tissue & colData(spe)$sample_id == id)]
  p1 = plotVisium(orig_spe,spots = TRUE,y_reverse = TRUE)
  p2 = plotVisium(xnew,spots = TRUE,y_reverse = TRUE)
  grid.arrange(p1,p2,ncol=2)
}
dev.off()

slide_order = unique(spe$slide)
sample_order <- unlist(sapply(slide_order, function(i) {
  (unique(spe$sample_id)[grepl(i, unique(spe$sample_id))])
}))
sample_order

spe$position = factor(spe$position, levels = c("TL","TR","BR","BL"))
spe$sample_id = spe$samlple_id[order(spe$position)]

Newsample_order <- unlist(sapply(slide_order, function(i) {
  (unique(spe$sample_id)[grepl(i, unique(spe$sample_id))])
}))
Newsample_order