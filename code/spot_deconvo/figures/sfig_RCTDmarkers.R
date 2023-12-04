setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")

suppressPackageStartupMessages(library("sessioninfo"))
suppressPackageStartupMessages(library("here"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("spacexr"))
suppressPackageStartupMessages(library("spatialLIBD"))
library('gridExtra')
library('viridis')
library(dplyr)

plot_dir <- here("plots", "spot_deconvo", "figures", "fig_RCTDmarkers")

Dr <- here("processed-data","spot_deconvo","RCTD","2ndRun_newClass_RCTDmarkers", "broad")
myRCTD2 = readRDS(here(Dr,"V12D07-332_B1/V12D07-332_B1_full.rds"))

## RCTD
barcodes <- colnames(myRCTD2@spatialRNA@counts)
weights <- myRCTD2@results$weights
norm_weights <- normalize_weights(weights)
celltypes = colnames(norm_weights)
coords = myRCTD2@spatialRNA@coords
dat = data.frame(celltype = c(), weights = c(), x = c(), y = c())
for (i in celltypes){
  temp = cbind(rep(i,4663),norm_weights[,i],coords)
  dat = rbind(dat,temp)
}
colnames(dat) = c("celltype", "weights", "x", "y")

png(here(plot_dir, "RCTD.png"), width = 2400, height = 1200, units = "px") 
ggplot(dat, aes(x = x, y = y , color = weights)) + 
  geom_point(size = 2.3) +
  scale_color_distiller(type = "seq",palette = rev('Greys'),direction=1, na.value = "#CCCCCC40")+
  theme_void() + labs(y = "RCTD markers" )+ 
  scale_y_reverse()+
  theme(text = element_text(size=24, color="black"),
        axis.title.y = element_text(angle = 90),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank())+
  facet_wrap(~celltype, nrow=2)
dev.off()


## deconvo
rm(temp)
rm(dat)
Dr <- here("processed-data","spot_deconvo","RCTD","3rdRun_newClass_deconvoMarkers", "broad")
myRCTD2 = readRDS(here(Dr,"V12D07-332_B1/V12D07-332_B1_full.rds"))
barcodes <- colnames(myRCTD2@spatialRNA@counts)
weights <- myRCTD2@results$weights
norm_weights <- normalize_weights(weights)
celltypes = colnames(norm_weights)
coords1 = myRCTD2@spatialRNA@coords
dat = data.frame(celltype = c(), weights = c(), x = c(), y = c())
differ = anti_join(coords,coords1)
for (i in celltypes){
  temp = cbind(celltype = rep(i,4373), weights = norm_weights[,i],coords1)
  temp1 = cbind(celltype = rep(i,290),weights = rep(NA,290),differ)
  dat = rbind(dat,temp,temp1)
}

png(here(plot_dir, "deconvo.png"), width = 2400, height = 1200, units = "px") 
ggplot(dat, aes(x = x, y = y , color = weights)) + 
  geom_point(size = 2.3) +
  scale_color_distiller(type = "seq",palette = rev('Greys'),direction=1, na.value = "red")+
  theme_void() + labs(y = "User defined markers" )+ 
  scale_y_reverse()+
  theme(text = element_text(size=24, color="black"),
        axis.title.y = element_text(angle = 90),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank())+
  facet_wrap(~celltype, nrow=2)
dev.off()

