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

load(here("processed-data", "06_clustering","PRECAST","spe_norm_with_domain.rda"))
#load(here("processed-data", "06_clustering", "PRECAST", "spe_norm_with_domain.rda"))
marker_stats = readRDS(here("processed-data/spot_deconvo/shared_utilities/marker_stats_broad_class.rds"))


load('/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/plots/spatial_palettes.rda')
colors = c("CA1" = "#00dc00", "WM" = "#7a007a", "Choroid" = "#00006a", "SL.SR" = "#000000", "SLM.SGZ" = "#dfa56e" ,
			"SR.SLM" = "#555555", "SUB.RHP" = "#61963d", "GABA" = "#5ffffb", "Vascular" = "#1e1eff", "RHP" =  "#99ff99",
			 "SUB" = "#add294", "CA2.4" = "#008000", "GCL" =  "#002800" , "ML" = "#a1a1a1") 
load(here("plots","snRNAseq_palettes.rda"))

 marker_stats = readRDS(here("processed-data/spot_deconvo/shared_utilities/marker_stats_broad_class.rds"))
  celltypes = unique(marker_stats$cellType.target)
  i = celltypes[1]
  Oligo = "WM"
  InhN = c("RHP", "GABA")
  ExcN = c("GCL", "CA2.4", "CA1", "SUB", "SUB.RHP", "RHP")
  Astro = c("SL.SR", "ML", "SR.SLM", "SLM.SGZ")
  Vascular = "Vascular"  
  CSF = "Choroid"

sn.mid.palette
       ExcN          InhN Micro/Macro/T         Astro         Oligo 
  "#009000"     "#33fff9"     "#FAFA33"     "#dfa56e"     "#a300a0" 
        OPC      Vascular           CSF 
  "#7070ff"     "#0000f4"     "#000030" 

  
  for (i in celltypes){

    markers <- marker_stats |>filter(cellType.target == i,rank_ratio <= 25,ratio > 1) |>pull(gene)
    spT <- spe_norm[match(markers, rowData(spe_norm)$gene_id), ]
    spT$prop_nonzero_marker <- colMeans(assays(spT)$counts > 0)
   # spT$prop_nonzero_marker = as.numeric(t(assays(spT)$counts))
    datb = as.data.frame(colData(spT)) 
	levels(datb$domain)[levels(datb$domain)=="WM.1"] <- "WM"
	levels(datb$domain)[levels(datb$domain)=="WM.2"] <- "WM"
	levels(datb$domain)[levels(datb$domain)=="WM.3"] <- "WM"
	levels(datb$domain)[levels(datb$domain)=="SLM.WM"] <- "SLM.SGZ"
	datb$domain1 <- ifelse(datb$domain %in% CSF, "pos", "other")
	
    png(here("plots", "revision_maddy", paste0(i,"_box.png")), width = 1200, height = 300, units = "px") 
    p = ggplot(data = datb, aes(x=domain, y=prop_nonzero_marker, fill = domain1))+
  #   geom_violin() + scale_fill_manual(values = colors)+
      geom_boxplot(outlier.shape = NA) + scale_fill_manual(values = c("pos" = "#000030" , "others" = "#CCCCCC40"))+
      theme_bw() + labs(y = i, x = "", fill = "domain" ) + 
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
  
  
  Dr <- here("processed-data","spot_deconvo","shared_utilities")
  spe = readRDS(here(Dr,"spe.rds"))
  spe_colData = as.data.frame(colData(spe))

