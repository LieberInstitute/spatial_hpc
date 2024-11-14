#testing BayesSpace
library(SpatialExperiment)
library(scater)
set.seed(123)

#loading the file that was used as input for PRECAST
load(file = here::here("processed-data", "05_preprocess_batchCorrection", "spe_norm_final.rda"))
#there are V-SPG samples here but they are subset out in Eric's code by filtering to the first 10 brnum
table(colData(spe)[,c("brnum","VSPG")])

spe<-spe[,spe$brnum %in% levels(spe$brnum)[1:10]]
spe$brnum<-droplevels(spe$brnum)
table(colData(spe)[,c("brnum","VSPG")])

#loading nnSVG list
load(file = here::here("processed-data","nnSVG","nnSVG_gene_lists.rda"))
length(nnSVG) #2k genes

#redo PCA to make sure it is on nnSVGs
spe <- runPCA(spe, subset_row=nnSVG)
#check for need to correct for batches
spe <- runUMAP(spe, dimred="PCA")
plotReducedDim(spe, dimred="UMAP", color_by = "brnum", point_size=.3)
plotReducedDim(spe, dimred="UMAP", color_by = "brnum", point_size=.3)+facet_wrap(vars(colour_by))


#batch correction (brain donor)
mnn <- batchelor::reducedMNN(reducedDim(spe,'PCA'), batch=spe$brnum, k=20)
reducedDim(spe,'MNN') <- mnn$corrected
spe <- runUMAP(spe, dimred="MNN", name="UMAP_MNN")
plotReducedDim(spe, dimred="UMAP_MNN", color_by = "brnum", point_size=.3)
plotReducedDim(spe, dimred="UMAP_MNN", color_by = "brnum", point_size=.3)+facet_wrap(vars(colour_by))

spe <- runUMAP(spe, dimred="MNN", name="UMAP_MNN_20", n_dimred=20)
plotReducedDim(spe, dimred="UMAP_MNN_20", color_by = "brnum", point_size=.3)
plotReducedDim(spe, dimred="UMAP_MNN_20", color_by = "brnum", point_size=.3)+facet_wrap(vars(colour_by))

spe_save <- spe
dim(spe_save)

#PRECAST domains
load(here::here('processed-data','06_clustering','PRECAST','spe_precast_HE_domain.rda'))
dim(spe)
identical(colnames(spe), colnames(spe_save))
colData(spe_save) = cbind(colData(spe_save), colData(spe)[,c("PRECAST_k18","cluster","broad.domain","domain")])
spe <- spe_save
save(spe, file="processed-data/06_clustering/BayesSpace/revision/spe_bayes.Rdata")

#BayesSpace clustering
setwd('/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/')
suppressPackageStartupMessages({
  library(SpatialExperiment)
  library(BayesSpace)
})
# vignette: https://www.bioconductor.org/packages/release/bioc/vignettes/BayesSpace/inst/doc/BayesSpace.html
### needs log transformed data
### needs PCs
set.seed(123)
load("processed-data/06_clustering/BayesSpace/revision/spe_bayes.Rdata")
cat("\n",format(Sys.time()),"spatialPreprocess\n")
spe$row <- spe$array_row
spe$col <- spe$array_col
spe <- spatialPreprocess(spe, platform="Visium", skip.PCA=TRUE, log.normalize=FALSE)
cat(format(Sys.time()),"spatialCluster\n")
spe <- spatialCluster(spe, q=18, use.dimred="MNN", d=20, #start with using 20 MNN dims
               init.method="kmeans", #initial clustering method options = mclust or kmeans
               model="t", #Using t-distributed errors in likelihood/bayesian models is a different way to obtain robust methods, as the t-distribution has heavier tails than the normal.
               nrep=10000 #default nrep is 50k, vignette recommends min 10k so will start there
               )
table(colData(spe)[,c("cluster","spatial.cluster")])
cat("\n\n",format(Sys.time()),"save results\n")
save(spe, file="processed-data/06_clustering/BayesSpace/revision/spe_bayes_k18-kmeans-10k.Rdata")

## Reproducibility information
cat("\n\nReproducibility information:\n")
format(Sys.time())
proc.time()
options(width = 120)
sessionInfo()

#compare results
library(ggspavis)
library(dplyr)

load("processed-data/06_clustering/BayesSpace/revision/spe_bayes_k18-kmeans-10k.Rdata")
load(file=here::here('plots','spatial_palette_final.rda'))

fix_order = distinct(as.data.frame(colData(spe)), slide, array, brnum, sample_id, position, sex) %>% 
  arrange(slide, array)
sub4 = fix_order$sample_id[c(14,16,
                             20,21)]
spe_sub4 = spe[,spe$sample_id %in% sub4]

plotSpots(spe_sub4, sample_id = "sample_id", annotate = "domain")+scale_color_manual(values=spatial.palette)


prop.df = group_by(as.data.frame(colData(spe)), cluster, spatial.cluster) %>% tally() %>%
  group_by(spatial.cluster) %>% mutate(spatial.total=sum(n)) %>%
  mutate(prop.spatial=n/spatial.total) %>% ungroup() %>%
  filter(prop.spatial>.1) %>% group_by(cluster) %>%
  slice_max(n=3, prop.spatial)
print(n=54, prop.df)

spe_sub4$spatial.cluster.annot = factor(spe_sub4$spatial.cluster,
                                        levels=c(13,
                                                4,15,
                                                8,10,
                                                3,
                                                18,
                                                5,
                                                16,
                                                11,#SL.SR
                                                14,
                                                9,
                                                1,#11,#SLM.SGZ
                                                7,17,12,
                                                6,
                                                2),
                                        labels=c("GCL",
                                                 "CA2.4","CA2.4",
                                                 "CA1","CA1",
                                                 "SUB",
                                                 "SUB.RHP",
                                                 "RHP",
                                                 "GABA",
                                                 "SL.SR",
                                                 "ML",
                                                 "SR.SLM",
                                                 "SLM.SGZ",#"SLM.SGZ",
                                                 "WM.1","WM.2","WM.3",
                                                 "Vascular",
                                                 "Choroid"))

#spatial.clusters with non-obvious match to precast clusters: 6,1,11,10
### 1 = 40% CA2.4, ~30% neuropil (SL.SR/SLM.SGZ) --> SLM/SGZ
### 6 = 34% vascular, 22% SR.SLM, 14% SLM.SGZ --> vascular
### 10 = 44% SL.SR, 25% CA1, 17% SR.SLM --> SL.SR
### 11 = 43% SR.SLM, 17% SLM.SGZ, 10% WM.1 --> SLM.SGZ


group_by(as.data.frame(colData(spe)), cluster, spatial.cluster) %>% tally() %>%
  group_by(spatial.cluster) %>% mutate(spatial.total=sum(n)) %>%
  mutate(prop.spatial=n/spatial.total) %>% filter(spatial.cluster %in% c(6,1,11,10)) %>%
  slice_max(n=3, prop.spatial)


plotSpots(spe_sub4, sample_id = "sample_id", annotate = "spatial.cluster.annot")+
  scale_color_manual(values=spatial.palette2)






