setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")
suppressPackageStartupMessages(library("SingleCellExperiment"))
suppressPackageStartupMessages(library("here"))
suppressPackageStartupMessages(library("spatialLIBD"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("escheR"))

Dr <- here("processed-data","spot_deconvo","shared_utilities")
spe = readRDS(here(Dr,"spe.rds"))

load('/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/plots/spatial_palettes.rda')
colors = srt.palette
names(colors) <- c("WM", "PCL-CA1", "SO", "SR", "SLM", "SUB", "ML", "CA4", "GCL", "SGZ", "PCL-CA3", "CP", "SL", "CTX", "THAL")

library(escheR)
#sample_id = "V10B01-086_D1"
#
#p = make_escheR(spe[, spe$sample_id == sample_id]) |> add_fill("nmf77") |> add_ground("ManualAnnotation") + scale_color_manual(values = colors) +
#      scale_fill_distiller(type = "seq",palette = rev('Greys'),direction=1, na.value = "#CCCCCC40")  +  labs(title = sample_id)
#
 pdf(here("plots", "revision", "nmfVsmanualannotions.pdf"), width = 10, height = 8)  # Adjust dimensions as needed

 # List of sample IDs
 sample_ids <- unique(spe$sample_id)

 # Loop over each sample_id and create the plots
 for (id in sample_ids) {

   # Create the plot
   p <- make_escheR(spe[, spe$sample_id == id]) |> 
         add_fill("nmf77") |> 
         add_ground("ManualAnnotation") + 
         scale_color_manual(values = colors) + 
         scale_fill_distiller(type = "seq", palette = rev('Greys'), direction = 1, na.value = "#CCCCCC40") +
         labs(title = id)

   # Print the plot to the PDF
   print(p)
 }

 # Close the PDF device
 dev.off()
	  

## compare to RCTD
Dr = here("processed-data", "spot_deconvo", "RCTD", "2ndRun_newClass_RCTDmarkers", "layer")
csv_files = data.frame(files = list.files(Dr, pattern = "clusters.csv" , recursive = TRUE, full.names = TRUE))
csv_files = data.frame(files = csv_files[1:32,])
counts_list <- lapply(csv_files$files, function(file_path) 
    {data <- read.csv(file_path, row.names = NULL)  
    return(data)})
temp_df <- do.call(rbind, counts_list)
rownames(temp_df) <- NULL
temp_df$key = temp_df$X
temp_df$key = sapply(strsplit(temp_df$key,"_Br"), `[`, 1)

merged_df <- merge(as.data.frame(colData(spe)),temp_df, by = "key")
RCTD_nmf = merged_df

##Oligo
cor(RCTD_nmf$nmf77, RCTD_nmf$Oligo)
[1] 0.9441788
cor(RCTD_nmf$nmf42, RCTD_nmf$Oligo)
[1] 0.2036639
cor(RCTD_nmf$nmf44, RCTD_nmf$Oligo)
[1] 0.5036317
cor(RCTD_nmf$nmf38, RCTD_nmf$Oligo)
[1] 0.4283612

##Astro
cor(RCTD_nmf$nmf79, RCTD_nmf$Astro)
[1] 0.1435998

cor(RCTD_nmf$nmf81, RCTD_nmf$Astro)
[1] 0.8849263

cor(merged_data$nmf81, merged_data$CART_astrocyte)
[1] 0.02924326

> cor(RCTD_nmf$nmf87, RCTD_nmf$Ependy)
[1] 0.8513248

> cor(RCTD_nmf$nmf36, RCTD_nmf$OPC)
[1] 0.7379624

> cor(RCTD_nmf$nmf40, RCTD_nmf$Sub.1)
[1] 0.8222314

> cor(RCTD_nmf$nmf54, RCTD_nmf$Sub.2)
[1] 0.924581

> cor(RCTD_nmf$nmf15, RCTD_nmf$CA1_ProS)
[1] 0.8925183

> cor(RCTD_nmf$nmf61, RCTD_nmf$CA2.4)
[1] 0.6903548

> cor(RCTD_nmf$nmf48, RCTD_nmf$Choroid)
[1] 0.5243216

> cor(RCTD_nmf$nmf35, RCTD_nmf$GABA.LAMP5)
[1] 0.6612972

> cor(RCTD_nmf$nmf47,RCTD_nmf$GABA.LAMP5)
[1] 0.6773249

> cor(RCTD_nmf$nmf26,RCTD_nmf$GC)
[1] 0.4238688

> cor(RCTD_nmf$nmf5,RCTD_nmf$GC)
[1] 0.9607683

> cor(RCTD_nmf$nmf84,RCTD_nmf$L2_3.Prs.Ent)
[1] 0.4079313

> cor(RCTD_nmf$nmf53,RCTD_nmf$L6_6b)
[1] 0.6210496

> cor(RCTD_nmf$nmf65,RCTD_nmf$L6_6b)
[1] 0.5938282

> cor(RCTD_nmf$nmf22,RCTD_nmf$L6_6b)
[1] 0.445467

> cor(RCTD_nmf$nmf51,RCTD_nmf$L5)
[1] 0.683681

> cor(RCTD_nmf$nmf68,RCTD_nmf$L5)
[1] 0.5910643

> cor(RCTD_nmf$nmf82,RCTD_nmf$Micro_Macro_T)
[1] 0.7334912

> cor(RCTD_nmf$nmf90,RCTD_nmf$Micro_Macro_T)
[1] 0.701673

colData(spg) <- as(colData(RCTD_nmf), "DataFrame")

library(patchwork)
sample_ids <- unique(spg$sample_id)
pdf(here("plots", "revision", "nmf77VsRCTDoligo.pdf"), width = 20, height = 8)  # Adjust dimensions as needed
for (id in sample_ids) {
subset_data <- RCTD_nmf[RCTD_nmf$sample_id == id, ]

p = make_escheR(subset_data) |> add_fill("nmf77") |> add_ground("domain") + scale_color_manual(values = colors) +
      scale_fill_distiller(type = "seq",palette = rev('Greys'),direction=1, na.value = "#CCCCCC40")  +  labs(title = id)

p1 = make_escheR(subset_data) |> add_fill("Oligo") |> add_ground("domain") + scale_color_manual(values = colors) +
      scale_fill_distiller(type = "seq",palette = rev('Greys'),direction=1, na.value = "#CCCCCC40")  +  labs(title = id)

combined_plot = p + p1+ plot_layout(ncol = 2)

  # Print combined plot to the PDF device
  print(combined_plot)
	  }

  # Close the PDF device
  dev.off()
  
  load('/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/plots/spatial_palettes.rda')
  colors = srt.palette  
  
  sample_ids <- unique(spg$sample_id)
  pdf(here("plots", "revision", "nmf54VsRCTDsub2.pdf"), width = 20, height = 8)  # Adjust dimensions as needed
  for (id in sample_ids) {
	  
 subset_data <- RCTD_nmf[RCTD_nmf$sample_id == id, ]
 p = ggplot(data = subset_data, aes(x = array_row, y = array_col)) +
    geom_point(aes(color = domain, fill = nmf54), size = 2.3, shape = 21, stroke = 0.5) +  # shape = 21 allows separate fill and color
    scale_fill_distiller(type = "seq", palette = rev('Greys'), direction = 1, na.value = "#CCCCCC40") +
    scale_color_manual(values = colors) + 
    labs(title = unique(subset_data$brnum)) + 
    #facet_wrap(~factor(sample_id), nrow = 2) +
    theme_void() 
	
    p1 = ggplot(data = subset_data, aes(x = array_row, y = array_col)) +
       geom_point(aes(color = domain, fill = Sub.2), size = 2.3, shape = 21, stroke = 0.5) +  # shape = 21 allows separate fill and color
       scale_fill_distiller(type = "seq", palette = rev('Greys'), direction = 1, na.value = "#CCCCCC40") +
       scale_color_manual(values = colors) + 
       labs(title = unique(subset_data$brnum)) + 
       #facet_wrap(~factor(sample_id), nrow = 2) +
       theme_void() 
	   
	   combined_plot = p + p1+ plot_layout(ncol = 2)

	     # Print combined plot to the PDF device
	     print(combined_plot)  
		 } 
  dev.off()
    
    theme(text = element_text(size=18, color="black"),
          panel.border=element_blank(),
          panel.grid.minor = element_blank(), 
          panel.grid.major = element_blank(),
          panel.spacing.x = unit(-20, "pt"),
          panel.spacing.y = unit(-20, "pt"),
          axis.ticks.length = unit(0, "pt"),
          plot.title = element_text(size = 48, colour = "black", hjust = 0.5))


 RCTD_nmf <- merge(merged_df,temp_df, by = "key", all.x = TRUE)
 RCTD_nmf <- RCTD_nmf[, !(colnames(RCTD_nmf) %in% c("X", "...1"))]
 colData(spg) <- cbind(colData(spg), RCTD_nmf)

#	  suppressPackageStartupMessages({
#	      library(SpatialExperiment)
#	      library(spatialLIBD)
#	      library(here)
#	      library(edgeR)
#	      library(scuttle)
#	      library(scater)
#	      library(scran)
#	      library(dplyr)
#	      library(PCAtools)
#	      library(sessioninfo)
#	      library(gridExtra)
#	      library(ggforce)
#	      library(pheatmap)
#	      library(scater)
#	      library(scran)
#	  })
#	  
#	  spe_pseudo <- aggregateAcrossCells(
#	      spe,
#	      DataFrame(
#	          domain = colData(spe)$cluster_collapsed,
#	          sample_id = colData(spe)$sample_id
#	      ))
#	  dim(spe_pseudo)
#	 
#	 
#	  unique_values <- colData(spe) %>%
#	    select(PRECAST_k16_nnSVG, cluster_collapsed) %>%
#	    distinct() %>%
#	    arrange(PRECAST_k16_nnSVG)
#		