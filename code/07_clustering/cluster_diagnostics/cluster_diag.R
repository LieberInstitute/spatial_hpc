load("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/02_build_spe/spe_bayes_clus.Rdata")
print("spe loaded")
setwd('/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/')
suppressPackageStartupMessages({
    library(dplyr)
    library(purrr)
    library(Seurat)
    library(SpatialExperiment)
    library(PRECAST)
    library(spatialLIBD)
    library(ggplot2)
    library(gridExtra)
    library("here")
    library(bluster)
})
Sys.time()
# ###Run silhouette boxplot code
# print('setting up data for boxplots')
# 
# clustering_columns <- colData(spe)[,c(20,50:74)]
# column_names <- colnames(colData(spe))[c(20,50:74)]
# 
# sil.data.list <- list()
# 
# for (i in seq_along(clustering_columns)) {
#     current_clustering <- clustering_columns[, i]
#     current_colname <- column_names[i]
#     
#     # Remove NAs from the current_clustering and corresponding reduced dimensions
#     non_na_indices <- !is.na(current_clustering)
#     current_clustering_no_na <- current_clustering[non_na_indices]
#     reduced_dim_no_na <- reducedDim(spe, "HARMONY")[non_na_indices,]
#     
#     sil.approx <- approxSilhouette(reduced_dim_no_na, clusters = current_clustering_no_na)
#     sil.data <- as.data.frame(sil.approx)
#     sil.data$closest <- factor(ifelse(sil.data$purity > 0, current_clustering_no_na, sil.data$other))
#     sil.data$cluster <- current_clustering_no_na
#     
#     sil.data.list[[current_colname]] <- sil.data
# }
# 
# avg_purity <- data.frame(clustering = character(), cluster = integer(), avg_purity = numeric())
# for (colname in names(sil.data.list)) {
#     sil.data <- sil.data.list[[colname]]
#     avg_purity_col <- aggregate(sil.data$purity, by = list(sil.data$cluster), FUN = mean)
#     colnames(avg_purity_col) <- c("cluster", "avg_purity")
#     avg_purity_col$clustering <- colname
#     avg_purity <- rbind(avg_purity, avg_purity_col)
# }
# 
# avg_purity$clustering <- factor(avg_purity$clustering, levels = names(sil.data.list))
# # Add a column for the algorithm
# avg_purity$algorithm <- factor(ifelse(grepl("bayes", avg_purity$clustering), "BayesSpace", 
#                                       ifelse(grepl("PRECAST", avg_purity$clustering), "PRECAST",
#                                              "ManualAnnotation")))
# 
# save(avg_purity,file=here::here("plots","06_clustering","cluster_diagnostics","sil_boxplot.rda")
#      boxplot <- ggplot(avg_purity, aes(x = clustering, y = avg_purity,fill=algorithm)) +
#          geom_boxplot() +
#          theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
#          facet_wrap(~ algorithm, scales = "free_x", strip.position = "bottom")+
#          scale_fill_manual(values = c("BayesSpace" = "#f88379", "PRECAST" = "#afeeee",
#                                       "ManualAnnotation" = "#f9e09c"))+
#          ggtitle("Average cluster silhouette purity for different clusterings")
#      pdf(here::here('plots', '06_clustering', 'cluster_diagnostics','avg_silhouette_boxplot.pdf'))
#      print(boxplot)
#      dev.off()
#      
#      Sys.time()
#      
     ##purity
     
     print("setting up purity boxplots")
     clustering_columns <- colData(spe)[,c(20,50:74)]
     column_names <- colnames(colData(spe))[c(20,50:74)]
     
     purity.data.list <- list()
     
     for (i in seq_along(clustering_columns)) {
         current_clustering <- clustering_columns[, i]
         current_colname <- column_names[i]
         
         # Remove NAs from the current_clustering and corresponding reduced dimensions
         non_na_indices <- !is.na(current_clustering)
         current_clustering_no_na <- current_clustering[non_na_indices]
         reduced_dim_no_na <- reducedDim(spe, "HARMONY")[non_na_indices,]
         
         purity.approx <- neighborPurity(reduced_dim_no_na, clusters = current_clustering_no_na)
         purity.data <- as.data.frame(purity.approx)
         purity.data$maximum <- factor(purity.data$maximum)
         purity.data$cluster <- current_clustering_no_na
         
         purity.data.list[[current_colname]] <- purity.data
     }
     
     avg_purity <- data.frame(clustering = character(), cluster = integer(), avg_purity = numeric())
     for (colname in names(purity.data.list)) {
         purity.data <- purity.data.list[[colname]]
         avg_purity_col <- aggregate(purity.data$purity, by = list(purity.data$cluster), FUN = mean)
         colnames(avg_purity_col) <- c("cluster", "avg_purity")
         avg_purity_col$clustering <- colname
         avg_purity <- rbind(avg_purity, avg_purity_col)
     }
     
     avg_purity$clustering <- factor(avg_purity$clustering, levels = names(purity.data.list))
     # Add a column for the algorithm
     avg_purity$algorithm <- factor(ifelse(grepl("bayes", avg_purity$clustering), "BayesSpace", 
                                    ifelse(grepl("PRECAST", avg_purity$clustering), "PRECAST",
                                           "ManualAnnotation")))

save(avg_purity,file=here::here("plots","06_clustering","cluster_diagnostics","pure_boxplot.rda"))
     
     
     boxplot <- ggplot(avg_purity, aes(x = clustering, y = avg_purity,fill=algorithm)) +
         geom_boxplot() +
         theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
         #facet_wrap(~ algorithm, strip.position = "bottom")+
         scale_fill_manual(values = c("BayesSpace" = "#f88379", "PRECAST" = "#afeeee",
                                      "ManualAnnotation" = "#f9e09c"))+
         ggtitle("Average cluster purity for different clusterings")
     pdf(here::here('plots', '06_clustering', 'cluster_diagnostics','avg_purityhouette_boxplot.pdf'))
     print(boxplot)
     dev.off()