library(spatialLIBD)

source(here("code","spot_deconvo","shared_utilities","plottingfunctions.R"))
source(here("code", "spot_deconvo", "shared_utilities", "shared_function.R"))

spe = readRDS(here("processed-data", "spot_deconvo", "shared_utilities", "spe.rds"))
dat = as.data.frame(colData(spe)) %>% select("key", "cluster_collapsed", "array_row", "array_col", "sample_id")
levels(dat$cluster_collapsed)[levels(dat$cluster_collapsed)=="WM.1"] <- "WM"
levels(dat$cluster_collapsed)[levels(dat$cluster_collapsed)=="WM.2"] <- "WM"
levels(dat$cluster_collapsed)[levels(dat$cluster_collapsed)=="WM.3"] <- "WM"
na_color = "#CCCCCC40"

Dr <- here("processed-data","spot_deconvo","shared_utilities")
marker_stats = readRDS(here(Dr, "marker_stats_layer_celltype_class1_noHATAGABAAmy.rds"))
  #   Get markers for this cell type
markers <- marker_stats |>filter(cellType.target == "GC",rank_ratio <= 25,ratio > 1) |>pull(gene)

spT_small <- spe[markers, spe$sample_id == "V11L05-333_D1"]

spT_small$prop_nonzero_marker <- colMeans(assays(spT_small)$counts > 0)
datb = as.data.frame(colData(spT_small))   

png(here("plots","spot_deconvo","figures","main_figure", "top25genesGC.png"), width = 600, height = 600, units = "px") 
 ggplot(data = datb, aes(x=array_row, y=array_col, color = prop_nonzero_marker))+
  geom_point(size = 2.3)+
  scale_color_distiller(type = "seq",palette = rev('Greys'),direction=1, na.value = na_color)+
  labs(title = i, color = "min>0" )+ theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
 dev.off()
 
gene <- rownames(spe)[match("PPFIA2", rowData(spe)$gene_name)]
spT_small <- spe[gene, spe$sample_id == "V11L05-333_D1"]
spT_small$gene = as.numeric(t(assays(spT_small)$counts))
datb = as.data.frame(colData(spT_small))   

png(here("plots","spot_deconvo","figures","main_figure", "PPFIA2.png"), width = 600, height = 600, units = "px") 
ggplot(data = datb, aes(x=array_row, y=array_col, color = gene))+
  geom_point(size = 2.3)+
  scale_color_distiller(type = "seq",palette = rev('Greys'),direction=1, na.value = na_color)+
  labs(title = i, color = "counts\nmin>0" )+ theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
dev.off()

na_color = "#CCCCCC40"

sample_id = ""
for (sample_id in unique(dat1$sample_id)){
  datb = dat1[which(dat1$sample_id == sample_id), ]
  datb[, celltypes][datb[, celltypes] <= 0.05] <- NA
  plot_list <- lapply(celltypes, function(i){
    ggplot(data = datb, aes(x=array_row, y=array_col, color = datb[,i]))+
      geom_point(size = 2.3)+facet_wrap(~tool)+
      #scale_color_gradientn(colours = viridisLite::plasma(21), na.value = na_color)+
      scale_color_distiller(type = "seq",palette = rev('Greys'),direction=1, na.value = na_color)+
      labs(title = i, color = "min>0.05" )+ theme_bw() +
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  })
  #gridplot = grid.arrange(grobs = plot_list, nrow = length(celltypes))
  ggsave(here("plots","spot_deconvo","shared_utilities",group,paste0(sample_id,"_new.pdf")), plot = marrangeGrob(plot_list, nrow=1, ncol=1),  width = 24, height = 8)
  print(paste0("done ", sample_id))
}
