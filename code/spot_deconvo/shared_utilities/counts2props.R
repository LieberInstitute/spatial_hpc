setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")
library("here")
library("ggplot2")
library("jaffelab")
library("tidyverse")
library("reshape2")
library("spatialLIBD")
library("cowplot")
library("sessioninfo")

tool =  "cell2location"
Dr = here("processed-data", "spot_deconvo",tool, "HE", "2ndRun_newClass", "layer")
csv_files = list.files(Dr, pattern = "clusters.csv" , recursive = TRUE, full.names = TRUE)

for (csv_file in csv_files) {
temp_df <- read.csv(csv_file)
celltypes = names(temp_df)[3:length(names(temp_df))]
df = data.frame(key = temp_df$key)
 for (celltype in celltypes) {
   df$count = rowSums(temp_df[,c(3:dim(temp_df)[2])])
   df[,celltype] = temp_df[,celltype]/df$count
 }
}

Dr1 <- here("processed-data","spot_deconvo","shared_utilities")
spe = readRDS(here(Dr1,"spe.rds"))
speb = spe[,spe$sample_id == "V11L05-333_D1"]
coords = data.frame(y = colData(speb)$array_row, x = colData(speb)$array_col)

plot_list <- lapply(celltypes, function(i) {
  ggplot(coords, aes(x = coords$y, y=coords$x , color = df[,i])) + labs(title = i, x="", y="") + 
    geom_point(size = 0.5)+scale_color_gradientn(colours = viridis(10, option = "magma"), limits = c(0,1)) +
    theme(legend.key.width = unit(0.1, "cm"), legend.title = element_blank())
})

gridplot = grid.arrange(grobs = plot_list, ncol = 4)
