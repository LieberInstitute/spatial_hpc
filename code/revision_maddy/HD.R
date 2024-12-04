library(raster)


for (i in genes) {
    vis_grid_gene(
        spe = spe,
        geneid = i,
        pdf = here::here("plots", "marker_genes", paste0(gsub("; ", "_", i), ".pdf")),
        assayname = "logcounts",
        minCount = 0,
        viridis = FALSE,
        point_size = 1.5
    )
}

vis_grid_gene(
    spe = sfe,
    geneid = 'SFRP2',
   # pdf = here::here("plots", "marker_genes", paste0(gsub("; ", "_", i), ".pdf")),
    assayname = "logcounts",
    minCount = 0,
    viridis = FALSE,
    point_size = 1.5)




	img_data <- imgData(sfe)
	print(img_data)  # Basic inspection
	str(img_data)    # Detailed structure


	image_object <- img_data$data[[1]]
	image_values <- image_object@values
	
	definition_str <- image_object@definition
	print(definition_str)

	# Extract values from the definition string if possible
	ncols <- as.numeric(sub(".*ncols=([0-9]+).*", "\\1", definition_str))
	nrows <- as.numeric(sub(".*nrows=([0-9]+).*", "\\1", definition_str))
	
	length(image_values) 
	image_array <- array(image_values, dim = c(nrows, ncols, 3))
	image_raster <- raster(image_array[, , 1])  #
	raster_matrix <- as.matrix(image_raster)
	
	raster_matrix_normalized <- raster_matrix / 255

	# Now, create the raster grob with the normalized matrix
	image_grob <- rasterGrob(raster_matrix_normalized, interpolate = TRUE)
	ggplot() +
	  annotation_custom(image_grob, 
	                    xmin = 0, xmax = max(plot_data$pxl_col_in_fullres),
	                    ymin = 0, ymax = max(plot_data$pxl_row_in_fullres)) +
	  geom_point(data = plot_data, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, color = expression), size = 0.5) +
	  scale_color_viridis_c(option = "magma", na.value = "transparent") +
	  coord_fixed() +
	  ggtitle("PROX1 Expression on H&E Image") +
	  theme_minimal()
	
	
	  plotSpots(sfe, annotate="SFRP2", assay_name="logcounts", point_size=0.01)+
	    #scale_color_viridis_c(option="F", direction=-1)+theme_void()+
	  theme(
	      plot.background = element_blank(), # Remove plot background
	      panel.background = element_blank(), # Remove panel background
	      legend.background = element_blank(), # Remove legend background
	      legend.key = element_blank() # Optional: remove background from legend keys
	    )
	
	
	
	ncols <- image_object@definition$ncols
	nrows <- image_object@definition$nrows

	# Reshape the image data into a 2D matrix (rows x columns)
	image_matrix <- matrix(image_values, nrow = nrows, ncol = ncols)
	image_raster <- raster(image_matrix)

	# Create a rasterGrob object for plotting
	image_grob <- rasterGrob(image_raster, interpolate = TRUE)
	
	coords <- spatialCoords(sfe)
	expression <- assay(sfe, "logcounts")["PROX1", ]

	plot_data <- data.frame(
	  pxl_col_in_fullres = coords[, 1],
	  pxl_row_in_fullres = coords[, 2],
	  expression = expression
	)

	ggplot(plot_data, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, color = expression)) +
	  geom_point(size = 0.5) +
	  scale_color_viridis_c(option = "magma", na.value = "transparent") +
	  coord_fixed() +
	  annotation_custom(image_grob, 
	                    xmin = 0, xmax = max(plot_data$pxl_col_in_fullres),
	                    ymin = 0, ymax = max(plot_data$pxl_row_in_fullres)) +
	  ggtitle("PROX1 Expression on H&E Image") +
	  theme_minimal()
	  
	 #######building SPE 
	 # run /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/code/revision_maddy/03_ficture_transcripts.sh to extract transcripts file
	  
	  library(here)
	  library(tidyverse)
	  library(scran)
	  library(spatialLIBD)
	  library(sessioninfo)
	  library(HDF5Array)

	  sample_id = 'HHPC_spe_016'
	  reference_gtf = '/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-GRCh38-2020-A/genes/genes.gtf'
	  
	  ## rename 'square_016um' to 'outs' in /dcs04/lieber/marmaypag/VisiumHD_HPC_pilot_LIBD4185/VisiumHD_HPC_pilot/processed-data/01_spaceranger/H1-6PGDCDB_A1/outs/binned_outputs/
	  sr_dir = '/dcs04/lieber/marmaypag/VisiumHD_HPC_pilot_LIBD4185/VisiumHD_HPC_pilot/processed-data/01_spaceranger/H1-6PGDCDB_A1/outs/binned_outputs/outs'
	  message(Sys.time(), ' | Building SpatialExperiment...')
	  spe <- read10xVisiumWrapper(
	      samples = sr_dir,
	      sample_id = sample_id,
	      type = "sparse",
	      data = "raw",
	      images = "lowres",
	      load = FALSE,
	      reference_gtf = reference_gtf
	  )
	  
	  
	  
 ##### plotting ####
 library(SpatialExperiment)
 library(scuttle)
 library(ggplot2)
 setwd('/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/')
 library(here)
 
 load(here('processed-data/revision/HPC_spe_016.Rdata'))
 spe <- logNormCounts(spe)
 plt_df <- data.frame(colData(spe), spatialCoords(spe))
 plt_df[['PROX1']] <- assay(spe, 'counts')['PROX1', ]
 plt_df[['PROX1']] <- scales::rescale(plt_df[['PROX1']], to = c(0, 1))
 
 # SFRP2, TNC, RELN, NPY
 var = 'RELN'
 plt_df[[var]] <- assay(spe, 'logcounts')[var, ]
 plt_df[[var]] <- scales::rescale(plt_df[[var]], to = c(0, 1))

 ## BLACK BACKGROUND
#plt_df$G <- plt_df[['PROX1']]
#plt_df$R <- plt_df[[var]]
#plt_df$B <- rep(0, nrow(plt_df))
##plt_df$RGB <- rgb(plt_df$R, plt_df$G, plt_df$B)
#plt_df$RGB <- rgb(plt_df$R, plt_df$G, plt_df$B, maxColorValue = 1)

## WHITE BACKGROUND
plt_df$G <- 1 - plt_df[['PROX1']]  # Scale PROX1 for green channel
plt_df$R <- 1 - plt_df[[var]]      # Scale NPY for red channel
plt_df$B <- rep(1, nrow(plt_df))   # Keep blue constant for white base

# Generate RGB values
plt_df$RGB <- rgb(plt_df$R, plt_df$G, plt_df$B, maxColorValue = 1)

	plot <- ggplot(plt_df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, color = RGB)) +
	  geom_point(size = 5) +  # Adjust point size as needed
	  scale_color_identity() 

	# Save the plot
	ggsave(here('plots', 'revision_maddy', paste0(var,'_HD.pdf')), plot = plot, width = 60, height = 60, dpi = 300, limitsize=FALSE)

	library(ggplot2)
	library(tidyr)
	library(dplyr)

	# Reshape the data to long format
	plt_long <- plt_df %>%
	  select(pxl_col_in_fullres, pxl_row_in_fullres, G, R) %>%
	  pivot_longer(cols = c(G, R), names_to = "Color", values_to = "Intensity")

	# Plot the data
	plot = ggplot(plt_long, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres)) +
	  geom_tile(aes(fill = Intensity)) +
	  facet_wrap(~Color) + # Separate panels for green and red
	  scale_fill_gradient(low = "white", high = "magenta", name = "G", na.value = "white") +
	  theme_minimal() +
	  labs(title = "Intensity Maps for C (Green) and P (Red)",
	       x = "Pixel Column",
	       y = "Pixel Row")

		 ggsave(here('plots', 'revision_maddy', 'p_legend.png'), plot = plot, width = 6, height = 10, dpi = 300, limitsize=FALSE)
 