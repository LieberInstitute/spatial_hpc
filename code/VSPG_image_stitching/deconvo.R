setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")
library("here")
library("jaffelab")
library("reshape")

midfile = read.csv(here("processed-data","VSPG_image_stitching", "broad2broad_deconvo.csv"), row.names =1)
split_midfile = split(midfile, midfile$tool) 
CART = split_midfile[[1]][c(1,3:7)]
colnames(CART)[2:6] <- paste('CART', colnames(CART)[2:6], sep = "_")
mid2broad_cell2location = split_midfile[[2]][c(1,3:7)]
mid2broad_cell2location <- mid2broad_cell2location[match(CART$key, mid2broad_cell2location$key), ]
colnames(mid2broad_cell2location) <- paste('mid2broad_cell2location', colnames(mid2broad_cell2location), sep = "_")
mid2broad_RCTD = split_midfile[[3]][c(1,3:7)]
mid2broad_RCTD <- mid2broad_RCTD[match(CART$key, mid2broad_RCTD$key), ]
colnames(mid2broad_RCTD) <- paste('mid2broad_RCTD', colnames(mid2broad_RCTD), sep = "_")
mid2broad_tangram = split_midfile[[4]][c(1,3:7)]
mid2broad_tangram <- mid2broad_tangram[match(CART$key, mid2broad_tangram$key), ]
colnames(mid2broad_tangram) <- paste('mid2broad_tangram', colnames(mid2broad_tangram), sep = "_")

finefile = read.csv(here("processed-data","VSPG_image_stitching", "layer2broad_deconvo.csv"), row.names =1)
split_finefile = split(finefile, finefile$tool) 

fine2broad_cell2location = split_finefile[[2]][c(1,3:7)]
fine2broad_cell2location <- fine2broad_cell2location[match(CART$key, fine2broad_cell2location$key), ]
colnames(fine2broad_cell2location) <- paste('fine2broad_cell2location', colnames(fine2broad_cell2location), sep = "_")
fine2broad_RCTD = split_finefile[[3]][c(1,3:7)]
fine2broad_RCTD <- fine2broad_RCTD[match(CART$key, fine2broad_RCTD$key), ]
colnames(fine2broad_RCTD) <- paste('fine2broad_RCTD', colnames(fine2broad_RCTD), sep = "_")
fine2broad_tangram = split_finefile[[4]][c(1,3:7)]
fine2broad_tangram <- fine2broad_tangram[match(CART$key, fine2broad_tangram$key), ]
colnames(fine2broad_tangram) <- paste('fine2broad_tangram', colnames(fine2broad_tangram), sep = "_")


df = cbind(CART,mid2broad_cell2location[,2:6], mid2broad_RCTD[,2:6], mid2broad_tangram[,2:6], fine2broad_cell2location[,2:6], fine2broad_RCTD[,2:6], fine2broad_tangram[,2:6])
write.csv(df, file = here("processed-data","VSPG_image_stitching", "deconvo.csv"))