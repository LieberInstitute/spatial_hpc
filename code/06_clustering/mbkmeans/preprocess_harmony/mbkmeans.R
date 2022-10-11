setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")

suppressPackageStartupMessages({
  library("SpatialExperiment")
  library("mbkmeans")
  library("fasthplus")
  library("here")
  library("sessioninfo")
  library("spatialLIBD")
  library("Polychrome")
  library("gridExtra")
})

temp = "OSCApreprocess_harmony_brain" #"spatialPreprocess_harmony", "OSCApreprocess_harmony_captureArea", "OSCApreprocess_harmony_brain"
## load data
load(file = here::here("processed-data", "05_preprocess_batchCorrection", paste0(temp ,"_spe.Rdata")), verbose = TRUE)
dim(spe)

set.seed(610)
k_list <- seq(5, 30)

message("Apply mbkmeans from 5:30 - ", Sys.time())
km_res <- lapply(k_list, function(k) {
  message("k=", k)
  mbkmeans(spe,
           clusters = k,
           batch_size = 500,
           reduceMethod = "HARMONY",
           calc_wcss = TRUE
  )
})

names(km_res[[1]])
save(km_res, file = here("processed-data", "06_clustering", "mbkmeans", "preprocess_harmony", paste0(temp, "_mbkmeans.Rdata")))

wcss <- sapply(km_res, function(x) sum(x$WCSS_per_cluster))

pdf(here("plots", "06_clustering", "mbkmeans", "preprocess_harmony", paste0(temp, "_mbkmeans_wcss.pdf")))
plot(k_list, wcss, type = "b")
abline(v = 15, lty = 2, col = "red")
dev.off()

for (k in c(5:30)) {
  cols <- Polychrome::palette36.colors(k)
  spe$kmeans <- paste0("mbk", km_res[[which(k_list==k)]]$Clusters)
  names(cols) <- sort(unique(spe$kmeans))
  pdf(here("plots", "06_clustering", "mbkmeans", "preprocess_harmony", paste0(temp, "_mbkmeans", k, ".pdf")), width = 21, height = 20)
  brains = unique(spe$brnum)
  for (i in seq_along(brains)){
    speb <- spe[, which(spe$brnum == brains[i])]
    samples <- unique(speb$sample_id)
    samples
    
    if (length(samples) == 1){
      p1 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "kmeans", colors = cols, point_size = 2, ... = paste0("_", brains[i]))
    } else if (length(samples) == 2){
      p1 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "kmeans", colors = cols, point_size = 3, ... = paste0("_", brains[i]))
      p2 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "kmeans", colors = cols, point_size = 3, ... = paste0("_", brains[i]))
      grid.arrange(p1, p2, nrow = 2)
    } else if (length(samples) == 3){
      p1 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "kmeans", colors = cols, point_size = 2, ... = paste0("_", brains[i]))
      p2 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "kmeans", colors = cols, point_size = 2, ... = paste0("_", brains[i]))
      p3 <- vis_clus(spe = speb, sampleid = samples[3], clustervar = "kmeans", colors = cols, point_size = 2, ... = paste0("_", brains[i]))
      grid.arrange(p1, p2, p3, nrow = 2)
    } else if (length(samples) == 4){
      p1 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "kmeans", colors = cols, point_size = 2, ... = paste0("_", brains[i]))
      p2 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "kmeans", colors = cols, point_size = 2, ... = paste0("_", brains[i]))
      p3 <- vis_clus(spe = speb, sampleid = samples[3], clustervar = "kmeans", colors = cols, point_size = 2, ... = paste0("_", brains[i]))
      p4 <- vis_clus(spe = speb, sampleid = samples[4], clustervar = "kmeans", colors = cols, point_size = 2, ... = paste0("_", brains[i]))
      grid.arrange(p1, p2, p3, p4, nrow = 2)
    } else if (length(samples) == 5){
      p1 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "kmeans", colors = cols, point_size = 2, ... = paste0("_", brains[i]))
      p2 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "kmeans", colors = cols, point_size = 2, ... = paste0("_", brains[i]))
      p3 <- vis_clus(spe = speb, sampleid = samples[3], clustervar = "kmeans", colors = cols, point_size = 2, ... = paste0("_", brains[i]))
      p4 <- vis_clus(spe = speb, sampleid = samples[4], clustervar = "kmeans", colors = cols, point_size = 2, ... = paste0("_", brains[i]))
      p5 <- vis_clus(spe = speb, sampleid = samples[5], clustervar = "kmeans", colors = cols, point_size = 2, ... = paste0("_", brains[i]))
      grid.arrange(p1, p2, p3, p4, p5, nrow = 2)}
  }
  
  dev.off()}