setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")
library("here")
library("jaffelab")
library("reshape")

load(here("processed-data","NMF","vspg_weights.rda"))
nmf = as.data.frame(proj)
nmf$key = rownames(nmf)
nmf$key <-sapply(strsplit(nmf$key,"_Br"), `[`, 1)

write.csv(nmf, file = here("processed-data","VSPG_image_stitching", "nmf.csv"))