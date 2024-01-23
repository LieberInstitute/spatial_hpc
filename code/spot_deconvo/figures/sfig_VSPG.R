setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")
suppressPackageStartupMessages(library("SingleCellExperiment"))
suppressPackageStartupMessages(library("DeconvoBuddies"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("sessioninfo"))
suppressPackageStartupMessages(library("here"))
suppressPackageStartupMessages(library("spatialLIBD"))
suppressPackageStartupMessages(library("cowplot"))

Dr <- here("processed-data","spot_deconvo","shared_utilities")
#load('/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/06_clustering/PRECAST/spe_norm_with_domain.rda')
#spg = spe_norm
spg = readRDS(here(Dr,"spg.rds"))
plot_dir <- here("plots", "spot_deconvo", "figures", "fig_VSPG")

#brnum = "Br3942_VSPG"
brnum = "Br8325_VSPG"

#gen = "SNAP25"
gen = "MBP"

gene <- rownames(spg)[match(gen, rowData(spg)$gene_name)]
spT <- spg[gene, spg$brnum == brnum]
spT$gene = as.numeric(t(assays(spT)$logcounts))
datb = as.data.frame(colData(spT)) 

png(here(plot_dir, paste0(brnum,"_", gen,".png")), width = 1200, height = 1200, units = "px") 
ggplot(data = datb, aes(x=array_row, y=array_col, color = gene))+
  geom_point(size = 2.3)+ theme_void() + 
  scale_color_distiller(type = "seq",palette = rev('Greys'),direction=1, na.value = "#CCCCCC40")+
  labs(title = gen, color = gen)+ 
#facet_wrap(~factor(sample_id, levels = c("V12D07-332_A1", "V12D07-332_D1", "V12D07-332_C1", "V12D07-332_B1")), nrow=2)+
facet_wrap(~factor(sample_id, levels = c("V12D07-335_B1", "V12D07-335_D1", "V12D07-335_C1", "V12D07-335_A1")), nrow=2)+
  
  theme(text = element_text(size=18, color="black"),
        panel.border=element_blank(),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        panel.spacing.x = unit(-20, "pt"),
        panel.spacing.y = unit(-20, "pt"),
        axis.ticks.length = unit(0, "pt"),
        plot.title = element_text(size = 48, colour = "black", hjust = 0.5))

dev.off()
