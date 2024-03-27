## adapted from https://www.stephaniehicks.com/biocdemo/articles/Demo.html
setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")

library("SpatialExperiment")
library("scater")
library("mbkmeans")
library("fasthplus")
library("here")
library("sessioninfo")
library("jaffelab")
library("ggplot2")
library("spatialLIBD")
library("Polychrome")
suppressPackageStartupMessages(library("gridExtra"))


## load data
load(file = here::here("processed-data", "06_clustering", "BayesSpace", "1st_run", "spe_modify.Rdata"), verbose = TRUE)

## neron paper returned 19 clusters for DLPFC, try 5:50
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
save(km_res, file = here("processed-data", "06_clustering", "mbkmeans", "1st_run", "mbkmeans.Rdata"))

# load(here("processed-data", "03_build_sce","km_res.Rdata"),verbose = TRUE)
## get wcss
wcss <- sapply(km_res, function(x) sum(x$WCSS_per_cluster))

pdf(here("plots", "06_clustering", "mbkmeans", "1st_run", "mbkmeans_wcss.pdf"))
plot(k_list, wcss, type = "b")
abline(v = 15, lty = 2, col = "red")
dev.off()

brains = c("Br6423","Br6432","Br2743","Br8325","Br3942","Br6471","Br8667","Br8492","Br6522")

## Visualize mbkmeans results
for (k in c(5:30)) {
cols <- Polychrome::palette36.colors(k)
spe$kmeans <- paste0("mbk", km_res[[which(k_list==k)]]$Clusters)
names(cols) <- sort(unique(spe$kmeans))
#
# vis_grid_clus(
#   spe = spe,
#   clustervar = "kmeans",
#   pdf_file = here("plots", "06_Clustering", "mbkmeans", paste0("vis_grid_clus_mbkmeans_k",k,".pdf")),
#   sort_clust = FALSE,
#   colors = cols,
#   spatial = FALSE,
#   point_size = 1,
#   sample_order = unique(spe$sample_id)
# )
# }

pdf(here("plots", "06_Clustering", "mbkmeans", "1st_run", paste0("mbkmeans_", k, ".pdf")), width = 21, height = 20)
# QC plot of tissue spots discarded

ii <- 1
speb <- spe[, which(spe$brnum == brains[ii])]
samples <- unique(speb$sample_id)
samples
p1 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "kmeans", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p2 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "kmeans", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p3 <- vis_clus(spe = speb, sampleid = samples[3], clustervar = "kmeans", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p4 <- vis_clus(spe = speb, sampleid = samples[4], clustervar = "kmeans", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))

grid.arrange(p1, p2, p3, p4, nrow = 2)

##
ii <- 2
speb <- spe[, which(spe$brnum == brains[ii])]
samples <- unique(speb$sample_id)
samples

p1 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "kmeans", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p2 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "kmeans", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p3 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "kmeans", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p4 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "kmeans", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))

grid.arrange(p1, p3, p2, p4, nrow = 2)

##
ii <- 3
speb <- spe[, which(spe$brnum == brains[ii])]
samples <- unique(speb$sample_id)
samples

p1 <- vis_clus(spe = speb, sampleid = samples[3], clustervar = "kmeans", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p2 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "kmeans", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p3 <- vis_clus(spe = speb, sampleid = samples[3], clustervar = "kmeans", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p4 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "kmeans", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))

grid.arrange(p1, p2, p3, p4, nrow = 2)

p1 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "kmeans", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p2 <- vis_clus(spe = speb, sampleid = samples[4], clustervar = "kmeans", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p3 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "kmeans", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p4 <- vis_clus(spe = speb, sampleid = samples[4], clustervar = "kmeans", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))

grid.arrange(p1, p3, p2, p4, nrow = 2)

##
ii <- 4
speb <- spe[, which(spe$brnum == brains[ii])]
samples <- unique(speb$sample_id)
samples

p1 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "kmeans", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p2 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "kmeans", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p3 <- vis_clus(spe = speb, sampleid = samples[5], clustervar = "kmeans", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p4 <- vis_clus(spe = speb, sampleid = samples[3], clustervar = "kmeans", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p5 <- vis_clus(spe = speb, sampleid = samples[4], clustervar = "kmeans", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
grid.arrange(p1, p2, p3, p4, p5, nrow = 2)

##
ii <- 5
speb <- spe[, which(spe$brnum == brains[ii])]
samples <- unique(speb$sample_id)
samples

p1 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "kmeans", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p2 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "kmeans", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p3 <- vis_clus(spe = speb, sampleid = samples[3], clustervar = "kmeans", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p4 <- vis_clus(spe = speb, sampleid = samples[4], clustervar = "kmeans", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))

grid.arrange(p1, p2, p3, p4, nrow = 2)

##
ii <- 6
speb <- spe[, which(spe$brnum == brains[ii])]
samples <- unique(speb$sample_id)
samples

p1 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "kmeans", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p2 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "kmeans", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p3 <- vis_clus(spe = speb, sampleid = samples[3], clustervar = "kmeans", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))

grid.arrange(p1, p2, p3, nrow = 2)

##
ii <- 7
speb <- spe[, which(spe$brnum == brains[ii])]
samples <- unique(speb$sample_id)
samples

p1 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "kmeans", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p2 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "kmeans", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p3 <- vis_clus(spe = speb, sampleid = samples[4], clustervar = "kmeans", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p4 <- vis_clus(spe = speb, sampleid = samples[3], clustervar = "kmeans", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))

grid.arrange(p1, p2, p3, p4, nrow = 2)

##
ii <- 8
speb <- spe[, which(spe$brnum == brains[ii])]
samples <- unique(speb$sample_id)
samples

p1 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "kmeans", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p2 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "kmeans", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p3 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "kmeans", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p4 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "kmeans", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))

grid.arrange(p1, p2, p3, p4, nrow = 2)

##
ii <- 9
speb <- spe[, which(spe$brnum == brains[ii])]
samples <- unique(speb$sample_id)
samples

p1 <- vis_clus(spe = speb, sampleid = samples[1], clustervar = "kmeans", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p2 <- vis_clus(spe = speb, sampleid = samples[2], clustervar = "kmeans", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p3 <- vis_clus(spe = speb, sampleid = samples[3], clustervar = "kmeans", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))
p4 <- vis_clus(spe = speb, sampleid = samples[4], clustervar = "kmeans", colors = cols, point_size = 2, ... = paste0("_", brains[ii]))

grid.arrange(p1, p2, p3, p4, nrow = 2)
dev.off()

}


#### Use fasthplus + wcss to refine k ####
# hpb estimate. t = pre-bootstrap sample size, D = reduced dimensions matrix, L = cluster labels, r = number of bootstrap iterations
l <- lapply(km_res, `[[`, 5)
length(l)

find_t <- function(L, proportion = 0.05) {
    initial_t <- floor(length(L) * proportion)
    # message("init. t = ", initial_t)
    smallest_cluster_size <- min(table(L))
    n_labels <- length(unique(L))
    message("smallest cluster: ", smallest_cluster_size, ", n lables: ", n_labels)
    ifelse(smallest_cluster_size > (initial_t / n_labels), initial_t, smallest_cluster_size * n_labels)
}

message("Find fasthplus for clusters - ", Sys.time())
fasthplus <- lapply(l, function(li) {
    message(Sys.time())
    initial_t <- find_t(L = li, proportion = 0.01)
    h <- hpb(
        D = reducedDims(spe)$HARMONY,
        L = li,
        t = initial_t,
        r = 30
    )
})

fasthplus <- unlist(fasthplus)


km_metrics <- data.frame(k = k_list, wcss = wcss, fasthplus = fasthplus)
write.csv(km_metrics, file = here("processed-data", "06_Clustering", "mbkmeans", "1st_run", "mbkmeans_metrics.csv"))

#### Plot metrics to select best k ####

pdf(here("plots", "06_Clustering", "mbkmeans", "1st_run", "mbkmeans_fastH.pdf"))
plot(k_list, fasthplus, type = "b")
abline(v = 15, lty = 2, col = "red")
dev.off()

## Save data
save(km_res, km_metrics, file = here("processed-data", "06_Clustering", "mbkmeans", "1st_run", "km_res.Rdata"))

# sgejobs::job_single('cluster_mb_kmeans', create_shell = TRUE, queue= 'bluejay', memory = '25G', command = "Rscript cluster_mb_kmeans.R")
## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
