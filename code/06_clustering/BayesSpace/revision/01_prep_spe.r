#performed in interactive mode
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

#PRECAST domain labels (for comparison down the line)
load(here::here('processed-data','06_clustering','PRECAST','spe_precast_HE_domain.rda'))
dim(spe)
identical(colnames(spe), colnames(spe_save))
colData(spe_save) = cbind(colData(spe_save), colData(spe)[,c("PRECAST_k18","cluster","broad.domain","domain")])
spe <- spe_save

save(spe, file="processed-data/06_clustering/BayesSpace/revision/spe_bayes.Rdata")
