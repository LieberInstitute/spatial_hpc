setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")
suppressPackageStartupMessages(library("SingleCellExperiment"))
suppressPackageStartupMessages(library("here"))
suppressPackageStartupMessages(library("spatialLIBD"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("escheR"))

## correlation	  
load(here("processed-data", "06_clustering","PRECAST","spe_norm_with_domain.rda"))
spg = spe_norm[,which(spe_norm$slide == "V12D07-332" | spe_norm$slide == "V12D07-335")]
#load(here('processed-data/NMF/spe_nmf_final.rda'))
#spg = spe[,which(spe$slide == "V12D07-332" | spe$slide == "V12D07-335")]
colData_df <- as.data.frame(colData(spg))
colData_df$key = sapply(strsplit(colData_df$key,"_Br"), `[`, 1)

library(readr)

nmf_df = read_csv(here("processed-data", "VSPG_image_stitching", "nmf.csv"),col_names = TRUE,show_col_types = FALSE)
merged_df <- merge(colData_df,nmf_df, by = "key", all.x = TRUE)

deconvo_df = read_csv(here("processed-data", "VSPG_image_stitching", "deconvo.csv"),col_names = TRUE,show_col_types = FALSE)
merged_data <- merge(merged_df, deconvo_df, by = "key", all.x = TRUE)


colData(spg) <- as(colData(merged_data), "DataFrame")


cor_nmf77_CART_oligo <- cor(merged_data$nmf77, merged_data$CART_oligo, use = "complete.obs")
[1] 0.5220068

cor_nmf79_CART_astrocyte <- cor(merged_data$nmf79, merged_data$CART_astrocyte, use = "complete.obs")
cat("Correlation between nmf77 and CART_oligo:", cor_nmf77_CART_oligo, "\n")
cat("Correlation between nmf79 and CART_astrocyte:", cor_nmf79_CART_astrocyte, "\n")
Correlation between nmf77 and CART_oligo: 0.5220068 
Correlation between nmf79 and CART_astrocyte: 0.03595106 

cor(merged_data$nmf81, merged_data$CART_astrocyte)
[1] 0.02924326

cor(merged_data$nmf7, merged_data$CART_neuron)
cor(merged_data$nmf13, merged_data$CART_neuron)

cor(merged_data$nmf7, merged_data$mid2broad_RCTD_neuron)
cor(merged_data$nmf13, merged_data$mid2broad_RCTD_neuron)

cor(merged_data$nmf44, merged_data$CART_oligo)
cor(merged_data$nmf89, merged_data$CART_astro)
