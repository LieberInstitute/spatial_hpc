
setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")
suppressPackageStartupMessages({
  library("here")
  library("spatialLIBD")
  library("ggplot2")
  library("patchwork")
  library("scater")
  library("harmony")
  library("BayesSpace")
  library("scran")
})


load(file=here::here("processed-data","pilot_data_checks","spe_harmony.Rdata"))
# OFFSET
# The spatial locations for each sample need to be offset so that spots of different samples are not neighbors. 
# summary(spatialData(speH)$array_row)
# summary(spatialData(speH)$array_col)
auto_offset_row <- as.numeric(factor(unique(speH$sample_id))) * 100
names(auto_offset_row) <-unique(speH$sample_id)
speH$row <- spatialData(speH)$array_row + auto_offset_row[speH$sample_id]
speH$col <- spatialData(speH)$array_col
# summary(colData(speH)$row)
# summary(colData(speH)$col)

pdf(file=here::here("processed-data", "pilot_data_checks", "plots", "hpc_BayesSpace_OffsetCheck.pdf"))
clusterPlot(speH, "subject", color = NA) + #make sure no overlap between samples
  labs(fill = "Subject", title = "Offset check")
dev.off()

Sys.time() 
speB = spatialCluster(speH, use.dimred = "HARMONY", q = 7, nrep = 10000) #use HARMONY
Sys.time()


pdf(file=here::here("processed-data", "pilot_data_checks", "plots", "hpc_BayesSpace_clusterPlot_10k.pdf"))
clusterPlot(speB, color = NA) + #plot clusters
  labs(title = "BayesSpace joint clustering")
dev.off()

save(speB, file=here::here("processed-data","pilot_data_checks", "spe_bayesSpace_10k.Rdata"))

