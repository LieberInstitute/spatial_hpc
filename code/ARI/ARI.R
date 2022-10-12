setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")
suppressPackageStartupMessages({
library("SpatialExperiment")
library("mclust")
library("spatialLIBD")
library("tidyr")
library("ggplot2")
library("ggpubr")
library("viridis")
library("here")
})

# load spe wih mbkmeans
load(file = here::here("processed-data", "05_preprocess_batchCorrection", "OSCApreprocess_harmony_captureArea_spe.Rdata"))
load(file = here("processed-data", "06_clustering", "mbkmeans", "preprocess_harmony", "spatialPreprocess_harmony_mbkmeans.Rdata"))
spe$mbkmeans_spatialPreprocess <- km_res[[13]]$Clusters
s1 = table(spe$sample_id,spe$mbkmeans_spatialPreprocess)
colSums(s1)

load(file = here("processed-data", "06_clustering", "mbkmeans", "preprocess_harmony", "OSCApreprocess_harmony_brain_mbkmeans.Rdata"))
spe$mbkmeans_OSCApreprocess_brain <- km_res[[13]]$Clusters
s2 = table(spe$sample_id,spe$mbkmeans_OSCApreprocess_brain)
colSums(s2)

load(file = here("processed-data", "06_clustering", "mbkmeans", "preprocess_harmony", "OSCApreprocess_harmony_captureArea_mbkmeans.Rdata"))
spe$mbkmeans_OSCApreprocess_captureArea <- km_res[[13]]$Clusters
s3 = table(spe$sample_id,spe$mbkmeans_OSCApreprocess_captureArea)
colSums(s3)

dim(spe)

load(file = here::here("processed-data", "manual_annotation_csv", "compiled_annotation_before_match.Rdata"))
csv2$sample_id = NULL
temp = merge(colData(spe), csv2, by.x = "key", by.y = "spot_name")
temp = as.data.frame(temp)

speA = temp[which(temp$mbkmeans_spatialPreprocess != "10"),]
speA = speA[which(speA$mbkmeans_spatialPreprocess != "6"),]
speB = temp[which(temp$mbkmeans_OSCApreprocess_brain != "5"),]
speB = speB[which(speB$mbkmeans_OSCApreprocess_brain != "17"),]
speC = temp[which(temp$mbkmeans_OSCApreprocess_captureArea != "3"),]
speC = speC[which(speC$mbkmeans_OSCApreprocess_captureArea != "9"),]

sample_ids <- unique(temp$sample_id)
ari.df <- data.frame(matrix(ncol = 4, nrow = 32))
row.names(ari.df) <- sample_ids
colnames(ari.df) <- c("sample_id", "spatialPreprocess_harmony", "OSCApreprocess_harmony_brain", "OSCApreprocess_harmony_captureArea")

for (i in seq_along(sample_ids)) {
  spe_subA <- temp[which(speA$sample_id == sample_ids[i]),]
  spe_subB <- temp[which(speB$sample_id == sample_ids[i]),]
  spe_subC <- temp[which(speC$sample_id == sample_ids[i]),]
  
  ari.df$sample_id <- sample_ids[i]
  ari.df[sample_ids[i], "spatialPreprocess_harmony"] <- adjustedRandIndex(spe_subA$ManualAnnotation.y, spe_subA$mbkmeans_spatialPreprocess)
  ari.df[sample_ids[i], "OSCApreprocess_harmony_brain"] <- adjustedRandIndex(spe_subB$ManualAnnotation.y, spe_subB$mbkmeans_OSCApreprocess_brain)
  ari.df[sample_ids[i], "OSCApreprocess_harmony_captureArea"] <- adjustedRandIndex(spe_subC$ManualAnnotation.y, spe_subC$mbkmeans_OSCApreprocess_captureArea)
}

library(reshape)
meltData <- melt(ari.df)
boxplot(data=meltData, value~variable)

spe$annotations = "none"
spe$annotations[which(spe$mbkmeans == "1")]="SLM"
spe$annotations[which(spe$mbkmeans == "2")]="SLM"
spe$annotations[which(spe$mbkmeans == "3")]="CA1"
spe$annotations[which(spe$mbkmeans == "4")]="crap"
spe$annotations[which(spe$mbkmeans == "5")]="SR"
spe$annotations[which(spe$mbkmeans == "6")]="choroid"
spe$annotations[which(spe$mbkmeans == "7")]="CA2_3"
spe$annotations[which(spe$mbkmeans == "8")]="WM"
spe$annotations[which(spe$mbkmeans == "9")]="choroid"
spe$annotations[which(spe$mbkmeans == "10")]="CA1"
spe$annotations[which(spe$mbkmeans == "11")]="ML"
spe$annotations[which(spe$mbkmeans == "12")]="crap"
spe$annotations[which(spe$mbkmeans == "13")]="SO"
spe$annotations[which(spe$mbkmeans == "14")]="SLM"
spe$annotations[which(spe$mbkmeans == "15")]="GCL"
spe$annotations[which(spe$mbkmeans == "16")]="SUB/SL"
spe$annotations[which(spe$mbkmeans == "17")]="CA4"

spe$annotations[which(spe$mbkmeans == "1")]="SLM"
spe$annotations[which(spe$mbkmeans == "2")]="ML"
spe$annotations[which(spe$mbkmeans == "3")]="PCL_CA1"
spe$annotations[which(spe$mbkmeans == "4")]="crap"
spe$annotations[which(spe$mbkmeans == "5")]="SR"
spe$annotations[which(spe$mbkmeans == "6")]="CP"
spe$annotations[which(spe$mbkmeans == "7")]="PCL_CA3"
spe$annotations[which(spe$mbkmeans == "8")]="WM"
spe$annotations[which(spe$mbkmeans == "9")]="CP"
spe$annotations[which(spe$mbkmeans == "10")]="PCLCA1"
spe$annotations[which(spe$mbkmeans == "11")]="ML"
spe$annotations[which(spe$mbkmeans == "12")]="crap"
spe$annotations[which(spe$mbkmeans == "13")]="WM"
spe$annotations[which(spe$mbkmeans == "14")]="SLM"
spe$annotations[which(spe$mbkmeans == "15")]="GCL"
spe$annotations[which(spe$mbkmeans == "16")]="SO"
spe$annotations[which(spe$mbkmeans == "17")]="CA4"


spe$mbkmeans = paste0(spe$mbkmeans,"_",spe$annotations)
sample_ids <- unique(spe$sample_id)
ari.df <- data.frame(matrix(ncol = 1, nrow = 32))
row.names(ari.df) <- sample_ids
colnames(ari.df) <- c("mbkmeans")

for (i in seq_along(sample_ids)) {
  spe_sub <- spe[, colData(spe)$sample_id == sample_ids[i]]
  ari.df$sample_id <- sample_ids[i]
  ari.df[sample_ids[i], "mbkmeans"] <- adjustedRandIndex(spe_sub$ManualAnnotation, spe_sub$mbkmeans)
 }

ari.df.long <- gather(ari.df, method, ari, SNN_k10_k7:bayesSpace, factor_key = TRUE)
save(ari.df.long, file = here::here("processed-data", "rdata", "pilot_dlpfc_data", "05_ARI", "pilot_ari_clustering_across.Rdata"))


#### make plots for fiure 2a, clustering across pilot data
# pdf(here::here("plots","05_ARI","pilot_data_ARI_clustering_across_2.pdf"))
# ggplot(ari.df.long, aes(x = method, y=ari)) +
#   geom_boxplot()+
#   theme_bw()+
#   geom_jitter(color="black", size=0.4, alpha=0.9)+
#   ylim(0,0.6)+
#   theme(text = element_text(size = 40))
# dev.off()



#### spaGCN ARI for pilot data

# load old ARI results
# load(file = here::here("processed-data", "rdata", "pilot_dlpfc_data","05_ARI", "pilot_ari_clustering_across.Rdata"))
# ari.df.long <-ari.df.long[1:48,]

# load pilot data
# load(file = here::here("processed-data","rdata","pilot_dlpfc_data","spe_pilot_bayesSpace_batch_corr_sampleID.Rdata"))

# sample_ids <- unique(spe$sample_id)
# df <- data.frame()
# for(i in seq_along(sample_ids)){
#   x <-read.csv(file = here::here("..","spython","spagcn","processed-data","03-our_data_analysis",sample_ids[i],"7_clusters","clusters.csv"))
#   df <- rbind(df,x)
# }
#
# write.csv(df,file = here::here("processed-data","rdata","pilot_dlpfc_data","clustering_results","spaGCN_k7","clusters.csv"))

# import spaGCN clusters
spe <- cluster_import(
  spe,
  cluster_dir = here::here("processed-data", "rdata", "pilot_dlpfc_data", "clustering_results"),
  prefix = "spaGCN_"
)

ari.df <- data.frame(matrix(ncol = 2, nrow = 12))
row.names(ari.df) <- sample_ids
colnames(ari.df) <- c("sample_id", "spaGCN_refined_cluster")

for (i in seq_along(sample_ids)) {
  spe_sub <- spe[, colData(spe)$sample_id == sample_ids[i]]
  ari.df$sample_id <- sample_ids[i]
  ari.df[sample_ids[i], "spaGCN_refined_cluster"] <- adjustedRandIndex(spe_sub$layer_guess_reordered, spe_sub$spaGCN_refined_cluster)
}
# append to ari.df.long
ari.df.long.2 <- gather(ari.df, method, ari, spaGCN_refined_cluster, factor_key = TRUE)

dim(ari.df.long)
# 48  3

dim(ari.df.long.2)
# 12  3

ari.df.long <- rbind(ari.df.long, ari.df.long.2)

dim(ari.df.long)
# 60  3

levels(ari.df.long$method) <- c(levels(ari.df.long$method), "Graph-Based", "Graph-based(BC)", "BayesSpace", "BayesSpace(BC)", "SpaGCN")
ari.df.long$method[ari.df.long$method == "SNN_k10_k7"] <- "Graph-Based"
ari.df.long$method[ari.df.long$method == "batch_corr_SNN_k10_k7"] <- "Graph-based(BC)"
ari.df.long$method[ari.df.long$method == "bayesSpace_pc"] <- "BayesSpace"
ari.df.long$method[ari.df.long$method == "bayesSpace"] <- "BayesSpace(BC)"
ari.df.long$method[ari.df.long$method == "spaGCN_refined_cluster"] <- "SpaGCN"

ari.df.long$general_method[c(1:12)] <- "GB"
ari.df.long$general_method[c(25:48)] <- "BS"
ari.df.long$general_method[c(49:60)] <- "SG"

save(ari.df.long, file = here::here("processed-data", "rdata", "pilot_dlpfc_data", "05_ARI", "pilot_ari_clustering_across.Rdata"))

# level_order <- c("Graph-Based","Graph-based(BC)","BayesSpace","BayesSpace(BC)","SpaGCN")
# pdf(here::here("plots","05_ARI","pilot_data_ARI_clustering_across.pdf"))
# ggplot(ari.df.long, aes(x = factor(method,level =  level_order), y=ari)) +
#   geom_boxplot(outlier.shape = NA)+
#   theme_bw()+
#   geom_jitter(color="black", size=1.0, alpha=0.9)+
#   ylim(0,0.6)+
#   theme(axis.text.x = element_text(size = 20,angle = 90, vjust = 0.5, hjust = 1, colour = c("blue","blue","red","red","green")),text = element_text(size = 30),axis.title = element_text(size = 30))+
#   ylab("Adjusted Rand Index")+
#   xlab("Clustering Method")
# dev.off()

pdf(here::here("plots", "05_ARI", "ggboxplot_pilot_data_ARI_clustering_across.pdf"))
ggboxplot(ari.df.long, x = "method", y = "ari", color = "general_method",
          palette = viridis(3), add = "jitter", repel = TRUE, 
          font.label = list(size = 10), legend = "none", ggtheme = theme_pubr(base_size = 20), 
          ylab = "Adjusted Rand Index", xlab = "Clustering Method", size = 1) +
  font("xy.text", size = 11) +
  font("xlab", size = 16) +
  font("ylab", size = 16)
dev.off()