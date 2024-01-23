setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")
suppressPackageStartupMessages(library("SingleCellExperiment"))
#suppressPackageStartupMessages(library("DeconvoBuddies"))
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

#spe = readRDS(here(Dr,"spe.rds"))
spe = load(here("processed-data","06_clustering", "PRECAST" , "spe_precast_final.rda"))
spg = spe[,which(spe$slide == "V12D07-332" | spe$slide == "V12D07-335")]
reducedDims(spg)$spatial <- spatialCoords(spg)
rownames(spg) <- rowData(spg)$gene_id

#load(here("processed-data", "06_clustering", "PRECAST", "spe_norm_with_domain.rda"))
#spg = spe_norm
marker_stats = readRDS(here(Dr,paste0("marker_stats_layer_celltype_class1_noHATAGABAAmy.rds")))

rownames(spg) <- rowData(spg)$gene_id
## top 15
markers <- marker_stats |>filter(cellType.target == "GC",rank_ratio <= 15,ratio > 1) |>pull(gene)
spT <- spg[markers, spg$brnum == "Br3942_VSPG"]
spT$prop_nonzero_marker <- colMeans(assays(spT)$counts > 0)
datb = as.data.frame(colData(spT))   
png(here(plot_dir, "top15.png"), width = 2400, height = 600, units = "px") 
ggplot(data = datb, aes(x=array_row, y=array_col, color = prop_nonzero_marker))+
  geom_point(size = 2.3)+ theme_void() + 
  scale_color_distiller(type = "seq",palette = rev('Greys'),direction=1, na.value = "#CCCCCC40")+
  labs(y = "Top 15 markers", color = "min>0" )+ 
  facet_wrap(~sample_id, nrow=1)+
  theme(text = element_text(size=48, color="black"),
        panel.border=element_blank(),
        axis.title.y = element_text(angle = 90),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        panel.spacing.x = unit(-15, "pt"),
        axis.ticks.length = unit(0, "pt"),
        legend.text = element_text(size=24))
  
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
  facet_wrap(~sample_id, nrow=1)+
  theme(text = element_text(size=48, color="black"),
        panel.border=element_blank(),
        axis.title.y = element_text(angle = 90),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        panel.spacing.x = unit(-15, "pt"),
        axis.ticks.length = unit(0, "pt"),
        legend.text = element_text(size=24))

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
  facet_wrap(~sample_id, nrow=1)+
  theme(text = element_text(size=48, color="black"),
        panel.border=element_blank(),
        axis.title.y = element_text(angle = 90),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        panel.spacing.x = unit(-15, "pt"),
        axis.ticks.length = unit(0, "pt"),
        legend.text = element_text(size=24))

dev.off()

## celltype genes
celltypes = unique(marker_stats$cellType.target)
for (i in celltypes){
  markers <- marker_stats |>filter(cellType.target == i,rank_ratio <= 25,ratio > 1) |>pull(gene)
  spT <- spg[markers, spg$brnum == "Br3942_VSPG"]
  spT$prop_nonzero_marker <- colMeans(assays(spT)$counts > 0)
  datb = as.data.frame(colData(spT)) 
  png(here(plot_dir, "top25forcelltypes", paste0(i,".png")), width = 2400, height = 600, units = "px") 
  print(here(plot_dir, paste0(i,".png")))
  p = ggplot(data = datb, aes(x=array_row, y=array_col, color = prop_nonzero_marker))+
    geom_point(size = 2.3)+
    scale_color_distiller(type = "seq",palette = rev('Greys'),direction=1, na.value = "#CCCCCC40")+
    theme_void() + labs(y = i, color = "min>0" ) + 
    facet_wrap(~sample_id, nrow=1) +
    theme(text = element_text(size=48, color="black"),
          panel.border=element_blank(),
          axis.title.y = element_text(angle = 90),
          panel.grid.minor = element_blank(), 
          panel.grid.major = element_blank(),
          panel.spacing.x = unit(-15, "pt"),
          axis.ticks.length = unit(0, "pt"),
          legend.text = element_text(size=24))
  
  print(p)
  dev.off()
  print(i)
}

##### figure verification 2
classical_markers <- c("PPFIA2", "AMPH", "FNDC1", "GFRA1", "KRT17", "C5orf63", "GAD2", "MIF", "FABP7", "MAN1A2", "SFRP2", "MOBP", "MAG", "MTURN", "PHLDB1", "ACTA2",    "TTR",     "SLC17A6",      "NTRK2")
          clusters <-c("GCL",   "CA2_4", "CA1",   "SUB",  "SUB.RHP",  "RHP",   "GABA", "SL_SR", "ML", "SR_SLM", "SLM_WM", "WM",   "WM1",  "WM2",    "WM3", "Vascular", "Choroid", "SUM.RHP_RHP", "SR_SLM_WM")

##domain genes
rownames(spg) <- rowData(spg)$gene_name
for (i in seq_along(classical_markers)){
  gene <- rownames(spg)[match(classical_markers[i], rowData(spg)$gene_name)]
  spT <- spg[gene, spg$brnum == "Br3942_VSPG"]
  spT$gene = as.numeric(t(assays(spT)$counts))
  datb = as.data.frame(colData(spT))   
  png(here(plot_dir, "domainMarkers", paste0(clusters[i],".png")), width = 2400, height = 600, units = "px") 
  p = ggplot(data = datb, aes(x=array_row, y=array_col, color = gene))+
    geom_point(size = 2.3)+
    scale_color_distiller(type = "seq",palette = rev('Greys'),direction=1, na.value = "#CCCCCC40")+
    theme_void() + labs(y = paste0(clusters[i], " : ", classical_markers[i]), color = "counts\nmin>0" ) + 
    facet_wrap(~sample_id, nrow=1) + 
    theme(text = element_text(size=48, color="black"),
          panel.border=element_blank(),
          axis.title.y = element_text(angle = 90),
          panel.grid.minor = element_blank(), 
          panel.grid.major = element_blank(),
          panel.spacing.x = unit(-15, "pt"),
          axis.ticks.length = unit(0, "pt"),
          legend.text = element_text(size=24))
  
  print(p)
  dev.off()
  print(i)
}



gene <- rownames(spg)[match(classical_markers[i], rowData(spg)$gene_name)]
spT <- spg[classical_markers, spg$brnum == "Br3942_VSPG"]
temp = as.data.frame(as.matrix(t(assays(spT)$counts)))
temp$key = rownames(temp)
datb = as.data.frame(colData(spT)) 
key = match(datb$key,temp$key)
temp = temp[key,]
datb = cbind(datb,temp)
datb = datb[,c("key", "array_row", "array_col", "sample_id", classical_markers)]
df = melt(datb, id = c("key", "array_row", "array_col", "sample_id"))
png(here(plot_dir, "classical_markers.png"), width = 2400, height = 11400, units = "px") 
p = ggplot(data = df, aes(x=array_row, y=array_col, color = value))+
  geom_point(size = 2.3)+
  scale_color_distiller(type = "seq",palette = rev('Greys'),direction=1, na.value = "#CCCCCC40")+
  theme_void() + labs(color = "counts\nmin>0", x = "", y = "" ) + 
  facet_grid(variable~sample_id, switch = "y") + 
  theme(text = element_text(size=48, color="black"),
        panel.border=element_blank(),
        axis.title.y = element_text(angle = 90),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        panel.spacing = unit(-15, "pt"),
        axis.ticks.length = unit(0, "pt"),
        legend.text = element_text(size=24),
        strip.text.y = element_text(angle = 90))


print(p)
dev.off()


## celltype genes violins

load(here("processed-data", "06_clustering", "PRECAST", "spe_norm_with_domain.rda"))
spg = spe_norm
rownames(spg) <- rowData(spg)$gene_id
# levels(spg$domain)[levels(spg$domain)=="WM.1"] <- "WM"
# levels(spg$domain)[levels(spg$domain)=="WM.2"] <- "WM"
# levels(spg$domain)[levels(spg$domain)=="WM.3"] <- "WM"

marker_stats = readRDS(here(Dr,paste0("marker_stats_layer_celltype_class1_noHATAGABAAmy.rds")))
celltypes = unique(marker_stats$cellType.target)

load('/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/plots/spatial_palettes.rda')
colors = srt.palette
#colors = unname(colors)

for (i in celltypes){
  markers <- marker_stats |>filter(cellType.target == i,rank_ratio <= 25,ratio > 1) |>pull(gene)
  spT <- spg[markers, ]
  spT$prop_nonzero_marker <- colMeans(assays(spT)$counts > 0)
 # spT$prop_nonzero_marker = as.numeric(t(assays(spT)$counts))
  datb = as.data.frame(colData(spT)) 
  png(here(plot_dir, "top25forcelltypes", paste0(i,"_box.png")), width = 1200, height = 600, units = "px") 
  p = ggplot(data = datb, aes(x=domain, y=prop_nonzero_marker, fill = domain))+
#   geom_violin() + scale_fill_manual(values = colors)+
    geom_boxplot(outlier.shape = NA) + scale_fill_manual(values = colors)+
    theme_bw() + labs(y = i, x = "", fill = "projected domain" ) + 
    ylim(0,0.75) +
    theme(text = element_text(size = 20, colour = "black"),
          axis.text = element_text(size = 24, colour = "black"),
          axis.text.x = element_text(angle = 90),
          axis.line = element_line(linewidth=2, colour = "black"),
          panel.grid.minor = element_blank(), 
          panel.grid.major = element_blank(),
          panel.border = element_blank(),
          axis.title.y = element_text(size = 40,  face = "bold"))
  
  print(p)
  dev.off()
  print(i)
}
