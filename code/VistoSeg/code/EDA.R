
setwd('/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/')
  
library("here")
library("ggplot2")
library("scater")

final <- read.csv(file.path(here::here("processed-data", "Images", "VistoSeg", "Capture_areas", "combined_refineVNS_metric.csv")), header = TRUE, stringsAsFactors = FALSE)
final$name = sapply(strsplit(final$name, "_refine"), '[', 1)


final$AreaM = isOutlier(final$Area, nmads = 3, batch = final$name)

p <- ggplot(final, aes(x=name, y=Area, colour_by = AreaM)) + geom_violin() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylim(0,1000)
