setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")
library("here")
library("sessioninfo")
library("SpatialExperiment")
library('ggplot2')
library('dplyr')
library('gridExtra')

## CART
Dr = here("processed-data", "spot_deconvo", "groundTruth", "03_CART")
csv_files = data.frame(files = list.files(Dr, pattern = "cell_metrics.csv" , recursive = TRUE, full.names = TRUE))
csv_files$sample_id = sapply(strsplit(csv_files$files,"/"), `[`, 11)
counts_list <- lapply(seq_along(csv_files$files), function(i) 
{data <- read.csv(csv_files$files[i], row.names = NULL) 
data$sample_id <- csv_files$sample_id[i]
return(data)})
temp_df <- do.call(rbind, counts_list)
rownames(temp_df) <- NULL

sample_ids = unique(temp_df$sample_id)
plot_list <- lapply(sample_ids, function(i) {
  ggplot(data = temp_df[which(temp_df$sample_id == i), ], aes(x=y, y=x, color = cell_type))+
  geom_point(size = 1)+labs(title = i)+
  scale_color_manual(values = c("astrocyte" = "orange", "microglia" = "yellow", "neuron" = "green", "oligo" = "magenta", "other" = "blue"))+
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "black",colour = "black"))})
  #gridplot = grid.arrange(grobs = plot_list, nrow = length(celltypes))
  ggsave(here("plots","spot_deconvo","groundTruth", "03_CART","CART.pdf"), plot = marrangeGrob(plot_list, nrow=1, ncol=1),  width = 20, height = 20)
