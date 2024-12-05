setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")
suppressPackageStartupMessages(library("SingleCellExperiment"))
suppressPackageStartupMessages(library("here"))
suppressPackageStartupMessages(library("spatialLIBD"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("escheR"))

library(readr)

nmf_df = as.data.frame(read_csv(here("processed-data", "VSPG_image_stitching", "nmf.csv"),col_names = TRUE,show_col_types = FALSE))
deconvo_df = as.data.frame(read_csv(here("processed-data", "VSPG_image_stitching", "deconvoVSPG.csv"),col_names = TRUE,show_col_types = FALSE))
merged_data <- merge(nmf_df, deconvo_df, by = "key", all.x = TRUE)

######## CART vs RCTD #######
cor(merged_data$CART_oligo, merged_data$fine2broad_RCTD_oligo, use = "complete.obs")
[1] 0.4832567
cor(merged_data$CART_astrocyte, merged_data$fine2broad_RCTD_astrocyte, use = "complete.obs")
[1] 0.02601568

cor(merged_data$CART_oligo, merged_data$mid2broad_RCTD_oligo, use = "complete.obs")
[1] 0.4960319
cor(merged_data$CART_astrocyte, merged_data$mid2broad_RCTD_astrocyte, use = "complete.obs")
[1] 0.02973257


######## nmf vs CART #######
cor(merged_data$nmf77, merged_data$CART_oligo, use = "complete.obs")
[1] 0.4615806
cor(merged_data$nmf44, merged_data$CART_oligo, use = "complete.obs")
[1] 0.4390692

cor(merged_data$nmf79, merged_data$CART_astrocyte, use = "complete.obs")
[1] 0.03943824
cor(merged_data$nmf81, merged_data$CART_astrocyte, use = "complete.obs")
[1] 0.02840451

######## nmf vs RCTD #######
cor(merged_data$nmf44, merged_data$fine2broad_RCTD_oligo, use = "complete.obs")
[1] 0.8763121
cor(merged_data$nmf81, merged_data$fine2broad_RCTD_astrocyte, use = "complete.obs")
[1] 0.9038443

cor(merged_data$nmf44, merged_data$mid2broad_RCTD_oligo, use = "complete.obs")
[1] 0.8757857
cor(merged_data$nmf81, merged_data$mid2broad_RCTD_astrocyte, use = "complete.obs")
[1] 0.9060335

