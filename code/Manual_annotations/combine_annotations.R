
setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")
suppressPackageStartupMessages(library("here"))
library(plyr)

files = list.files(path = here::here("processed-data", "manual_annotation_csv"), full.names = TRUE)
files = files[1:15]
compiled_annotation <- ldply(files, read.csv, header=TRUE)

save(compiled_annotation, file = here::here("processed-data", "compiled_annotation.Rdata"))