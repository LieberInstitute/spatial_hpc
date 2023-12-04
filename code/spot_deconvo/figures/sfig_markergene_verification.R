setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")
suppressPackageStartupMessages(library("SingleCellExperiment"))
suppressPackageStartupMessages(library("DeconvoBuddies"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("sessioninfo"))
suppressPackageStartupMessages(library("here"))
suppressPackageStartupMessages(library("spatialLIBD"))
suppressPackageStartupMessages(library("cowplot"))

source(here("code","spot_deconvo","shared_utilities","plottingfunctions.R"))
source(here("code", "spot_deconvo", "shared_utilities", "shared_function.R"))

Dr <- here("processed-data","spot_deconvo","shared_utilities")
plot_dir <- here("plots", "spot_deconvo", "figures", "fig_markergene_verification")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

#sce = readRDS(here(Dr,"sce_class1_noHATAGABAAmy.rds"))
spg = readRDS(here(Dr,"spg.rds"))
marker_stats = readRDS(here(Dr,paste0("marker_stats_layer_celltype_class1_noHATAGABAAmy.rds")))

rownames(spg) <- rowData(spg)$gene_id
## top 15
markers <- marker_stats |>filter(cellType.target == "GC",rank_ratio <= 15,ratio > 1) |>pull(gene)
spT <- spg[markers, spg$brnum == "Br3942_VSPG"]
spT$prop_nonzero_marker <- colMeans(assays(spT)$counts > 0)
datb = as.data.frame(colData(spT))   
png(here(plot_dir, "top15.png"), width = 2400, height = 600, units = "px") 
ggplot(data = datb, aes(x=array_row, y=array_col, color = prop_nonzero_marker))+
  geom_point(size = 2.3)+
  scale_color_distiller(type = "seq",palette = rev('Greys'),direction=1, na.value = "#CCCCCC40")+
  theme_void() + labs(y = "Top 15 markers", color = "min>0" )+ 
  theme(text = element_text(size=24, color="black"),
        axis.title.y = element_text(angle = 90),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank())+
  facet_wrap(~sample_id, nrow=1)
dev.off()

## top 25
markers <- marker_stats |>filter(cellType.target == "GC",rank_ratio <= 25,ratio > 1) |>pull(gene)
spT <- spg[markers, spg$brnum == "Br3942_VSPG"]
spT$prop_nonzero_marker <- colMeans(assays(spT)$counts > 0)
datb = as.data.frame(colData(spT))   
png(here(plot_dir, "top25.png"), width = 2400, height = 600, units = "px") 
ggplot(data = datb, aes(x=array_row, y=array_col, color = prop_nonzero_marker))+
  geom_point(size = 2.3)+
  scale_color_distiller(type = "seq",palette = rev('Greys'),direction=1, na.value = "#CCCCCC40")+
  theme_void() + labs(y = "Top 25 markers", color = "min>0" )+ 
  theme(text = element_text(size=24, color="black"),
        axis.title.y = element_text(angle = 90),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank())+
  facet_wrap(~sample_id, nrow=1)
dev.off()

## top 50
markers <- marker_stats |>filter(cellType.target == "GC",rank_ratio <= 50,ratio > 1) |>pull(gene)
spT <- spg[markers, spg$brnum == "Br3942_VSPG"]
spT$prop_nonzero_marker <- colMeans(assays(spT)$counts > 0)
datb = as.data.frame(colData(spT))   
png(here(plot_dir, "top50.png"), width = 2400, height = 600, units = "px") 
ggplot(data = datb, aes(x=array_row, y=array_col, color = prop_nonzero_marker))+
  geom_point(size = 2.3)+
  scale_color_distiller(type = "seq",palette = rev('Greys'),direction=1, na.value = "#CCCCCC40")+
  theme_void() + labs(y = "Top 50 markers", color = "min>0" )+ 
  theme(text = element_text(size=24, color="black"),
        axis.title.y = element_text(angle = 90),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank())+
  facet_wrap(~sample_id, nrow=1)
dev.off()

## celltype genes
celltypes = unique(marker_stats$cellType.target)
for (i in celltypes){
  markers <- marker_stats |>filter(cellType.target == i,rank_ratio <= 25,ratio > 1) |>pull(gene)
  spT <- spg[markers, spg$brnum == "Br3942_VSPG"]
  spT$prop_nonzero_marker <- colMeans(assays(spT)$counts > 0)
  datb = as.data.frame(colData(spT)) 
  png(here(plot_dir, paste0(i,".png")), width = 2400, height = 600, units = "px") 
  print(here(plot_dir, paste0(i,".png")))
  p = ggplot(data = datb, aes(x=array_row, y=array_col, color = prop_nonzero_marker))+
    geom_point(size = 2.3)+
    scale_color_distiller(type = "seq",palette = rev('Greys'),direction=1, na.value = "#CCCCCC40")+
    theme_void() + labs(y = i, color = "min>0" ) + 
    theme(text = element_text(size=24, color="black"),
          axis.title.y = element_text(angle = 90),
          panel.grid.minor = element_blank(), 
          panel.grid.major = element_blank())+
    facet_wrap(~sample_id, nrow=1)
  print(p)
  dev.off()
  print(i)
}

##### figure verification 2
classical_markers <- c("PPFIA2", "AMPH", "FNDC1", "GFRA1", "KRT17", "C5orf63", "GAD2", "MIF", "FABP7", "MAN1A2", "SFRP2", "MOBP", "MAG", "MTURN", "PHLDB1", "ACTA2",    "TTR",     "SLC17A6",      "NTRK2")
clusters <-c("GCL", "CA2_4", "CA1",   "SUB", "SUB.RHP",  "RHP",  "GABA", "SL_SR", "ML", "SR_SLM", "SLM_WM", "WM",   "WM1",  "WM2",    "WM3", "vascular", "choroid", "SUM.RHP_RHP", "SR_SLM_WM")

##domain genes
rownames(spg) <- rowData(spg)$gene_name
for (i in seq_along(classical_markers)){
  gene <- rownames(spg)[match(classical_markers[i], rowData(spg)$gene_name)]
  spT <- spg[gene, spg$brnum == "Br3942_VSPG"]
  spT$gene = as.numeric(t(assays(spT)$counts))
  datb = as.data.frame(colData(spT))   
  png(here(plot_dir, paste0(clusters[i],".png")), width = 2400, height = 600, units = "px") 
  p = ggplot(data = datb, aes(x=array_row, y=array_col, color = gene))+
    geom_point(size = 2.3)+
    scale_color_distiller(type = "seq",palette = rev('Greys'),direction=1, na.value = "#CCCCCC40")+
    theme_void() + labs(y = classical_markers[i], color = "counts\nmin>0" ) + 
    theme(text = element_text(size=24, color="black"),
          axis.title.y = element_text(angle = 90),
          panel.grid.minor = element_blank(), 
          panel.grid.major = element_blank())+
    facet_wrap(~sample_id, nrow=1)
  print(p)
  dev.off()
  print(i)
}


