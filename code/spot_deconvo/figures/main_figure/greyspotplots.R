library(spatialLIBD)

source(here("code","spot_deconvo","shared_utilities","plottingfunctions.R"))
source(here("code", "spot_deconvo", "shared_utilities", "shared_function.R"))

spg_in <- here("processed-data", "spot_deconvo", "shared_utilities" , "spg.rds")
spg = readRDS(spg_in)
na_color = "#CCCCCC40"

dat = as.data.frame(colData(spg)) %>% select("key", "array_row", "array_col", "sample_id")
gene <- rownames(spg)[match("PPFIA2", rowData(spg)$gene_name)]
spT_small <- spg[gene, spg$sample_id == "V12D07-332_B1"]
spT_small$gene = as.numeric(t(assays(spT_small)$counts))
datb = as.data.frame(colData(spT_small))   

png(here("plots","spot_deconvo","figures","main_figure", "PPFIA2.png"), width = 600, height = 600, units = "px") 
ggplot(data = datb, aes(x=array_row, y=array_col, color = gene))+
  geom_point(size = 2.3)+
  scale_color_distiller(type = "seq",palette = rev('Greys'),direction=1, na.value = na_color)+
  labs(title = "(i) PPFIA2\n GCL marker", color = "counts\nmin>0", x= "", y = "" )+ theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        axis.text = element_blank(),
        axis.ticks =element_blank(),
        legend.text = element_text(size = 12, colour = "black"),
        plot.title = element_text(size = 48, colour = "black", hjust = 0.5),
        legend.justification="right",
        legend.box.spacing = unit(-15, "pt"),
        plot.margin = unit(c(0,0,-10,-10), "pt"),
        axis.ticks.length = unit(0, "pt"))
dev.off()

Dr <- here("processed-data","spot_deconvo","shared_utilities")
marker_stats = readRDS(here(Dr, "marker_stats_layer_celltype_class1_noHATAGABAAmy.rds"))
markers <- marker_stats |>filter(cellType.target == "GC",rank_ratio <= 25,ratio > 1) |>pull(gene)

spT_small <- spg[markers, spg$sample_id == "V12D07-332_B1"]
spT_small$prop_nonzero_marker <- colMeans(assays(spT_small)$counts > 0)
datb = as.data.frame(colData(spT_small))   
png(here("plots","spot_deconvo","figures","main_figure", "top25genesGC.png"), width = 600, height = 600, units = "px") 
 ggplot(data = datb, aes(x=array_row, y=array_col, color = prop_nonzero_marker))+
  geom_point(size = 2.3)+ theme_bw() +
  scale_color_distiller(type = "seq",palette = rev('Greys'),direction=1, na.value = na_color)+
  labs(title = "(ii) Top 25 markers\n      GC celltype", color = "min>0", x = "", y = "" )+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        axis.text = element_blank(),
        axis.ticks =element_blank(),
        plot.title = element_text(size = 48, colour = "black", hjust = 0.5),
        legend.text = element_text(size = 12, colour = "black"),
        legend.justification="right",
        legend.box.spacing = unit(-15, "pt"),
        plot.margin = unit(c(0,0,-10,-10), "pt"),
        axis.ticks.length = unit(0, "pt"))
 dev.off()
 
