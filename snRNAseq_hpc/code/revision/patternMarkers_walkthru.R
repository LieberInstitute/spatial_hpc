#NMF markers and CoGAPS exploration for lab meeting 2024/12/10
setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc")

library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(gridExtra)
set.seed(123)

#load sce for avg counts
load(file='snRNAseq_hpc/processed-data/sce/sce_final.rda')
avg.expr = rowMeans(logcounts(sce))

#load nmfs
load(file='snRNAseq_hpc/processed-data/NMF/nmf_final.rda')
loads<-x@w
no_expr <- which(rowSums(loads) == 0)
loads <- loads[-no_expr, ]
avg.expr = avg.expr[rownames(loads)]

df1 = cbind.data.frame(loads, avg.expr)

# scatter plot of gene weight ~ avg expr
plist1 <- lapply(c("nmf44","nmf37","nmf55"), function(x) {
  tmp = df1[,c(x,"avg.expr")]
  colnames(tmp) = c("nmf","avg.expr")
  ggplot(tmp, aes(x=avg.expr, y=nmf))+
    geom_point(size=.3)+theme_bw()+
    labs(title=x, y="gene weights", x="avg. expr")+
    theme(text=element_text(size=10), axis.text.y=element_blank(),
          aspect.ratio=1)
})
do.call(grid.arrange, c(plist1, ncol=3))

##########################################################################
######################### CoGAPS pattern markers ######################### 

#rewrite for matrix: patternMarkers_m()
source("snRNAseq_hpc/code/revision/patternMarkers_matrix.R")

##############
### validate that patternMarkers rewrite didn't change outcome
##############

# cogaps example for comparison
cogapsresult_url <- "https://zenodo.org/record/7709664/files/cogapsresult.Rds"
bfc <- BiocFileCache::BiocFileCache()
cogapsresult <- BiocFileCache::bfcrpath(bfc, cogapsresult_url)
cogapsresult <- readRDS(cogapsresult)

test1 = CoGAPS::patternMarkers(cogapsresult, axis=1, threshold="all")
# all genes marker for an NMF
setdiff(rownames(cogapsresult@featureLoadings), unlist(test1$PatternMarkers)) 
### no genes in W matrix excluded

# a gene can only be a marker for 1 NMF pattern
length(unlist(test1$PatternMarkers))==length(unique(unlist(test1$PatternMarkers)))
### no gene in the pattern markers list repeated

test2 = patternMarkers_m(cogapsresult@featureLoadings, cogapsresult@sampleFactors, rescale=TRUE, 
                         threshold="all", axis=1)
# for cogaps example, rescale parameter diesn't matter because the feature loadings are already rescaled in this way
test3 = patternMarkers_m(cogapsresult@featureLoadings, cogapsresult@sampleFactors, rescale=FALSE, 
                         threshold="all", axis=1)

identical(test1$PatternScores, test2$PatternScores)
identical(test1$PatternScores, test3$PatternScores)

identical(test1$PatternMarkers, test2$PatternMarkers)
identical(test1$PatternMarkers, test3$PatternMarkers)

plot(test1$PatternScores, test1$PatternRanks)

##############
### run patternMarkers for HPC matrix or load pre-run results
##############

all.markers1 <- patternMarkers_m(x@w, x@h, "all", rescale=FALSE, axis=1)
saveRDS(all.markers1, "snRNAseq_hpc/processed-data/revision/patternMarkers_CoGAPS_all-rescale-FALSE.rda")
all.markers1 <- readRDS("snRNAseq_hpc/processed-data/revision/patternMarkers_CoGAPS_all-rescale-FALSE.rda")

plot(all.markers1$PatternScores, all.markers1$PatternRanks)
#didn't remove genes that had non-zero weights in all NMF patterns before running
plot(all.markers1$PatternScores[-no_expr,], all.markers1$PatternRanks[-no_expr,])

# weights ~ expr scatter plot with patternMarkers results labeled
plist2 <- lapply(c("nmf44","nmf37","nmf55"), function(x) {
  tmp = df1[,c(x,"avg.expr")]
  colnames(tmp) = c("nmf","avg.expr")
  tmp$is_top = rownames(tmp) %in% all.markers1[[1]][[x]]
  ggplot(tmp, aes(x=avg.expr, y=nmf, color=is_top))+
    geom_point(size=.3)+theme_bw()+
    geom_point(data=filter(tmp, is_top==T), size=.5)+
    scale_color_manual(values=c("black","red"))+
    labs(title=x, y="gene weights", x="avg. expr")+
    theme(text=element_text(size=10), axis.text.y=element_blank(), 
          legend.position="none", aspect.ratio=1)
})

# density plot of weights for patternMarker results genes vs other genes
plist3 <- lapply(c("nmf44","nmf37","nmf55"), function(x) {
  tmp = df1[,c(x,"avg.expr")]
  colnames(tmp) = c("nmf","avg.expr")
  tmp$is_top = rownames(tmp) %in% all.markers1[[1]][[x]]
  ggplot(tmp, aes(x=nmf, color=is_top))+
    geom_density(linewidth=.7)+theme_bw()+
    scale_color_manual(values=c("black","red"))+
    scale_x_log10()+
    labs(title=x, x="gene weights")+
    theme(text=element_text(size=10), axis.text.y=element_blank(),
          legend.position="none", aspect.ratio=1)
})
do.call(grid.arrange, c(c(plist2, plist3), ncol=3))

# additional (poorer performing) examples
plist2.1 <- lapply(c("nmf10","nmf48"), function(x) {
  tmp = df1[,c(x,"avg.expr")]
  colnames(tmp) = c("nmf","avg.expr")
  tmp$is_top = rownames(tmp) %in% all.markers1[[1]][[x]]
  ggplot(tmp, aes(x=avg.expr, y=nmf, color=is_top))+
    geom_point(size=.3)+theme_bw()+
    geom_point(data=filter(tmp, is_top==T), size=.5)+
    scale_color_manual(values=c("black","red"))+
    labs(title=x, y="gene weights", x="avg. expr")+
    theme(text=element_text(size=10), axis.text.y=element_blank(), 
          legend.position="none", aspect.ratio=1)
})
plist3.1 <- lapply(c("nmf10","nmf48"), function(x) {
  tmp = df1[,c(x,"avg.expr")]
  colnames(tmp) = c("nmf","avg.expr")
  tmp$is_top = rownames(tmp) %in% all.markers1[[1]][[x]]
  ggplot(tmp, aes(x=nmf, color=is_top))+
    geom_density(linewidth=.7)+theme_bw()+
    scale_color_manual(values=c("black","red"))+
    scale_x_log10()+
    labs(title=x, x="gene weights")+
    theme(text=element_text(size=10), axis.text.y=element_blank(),
          legend.position="none", aspect.ratio=1)
})
do.call(grid.arrange, c(c(plist2.1, plist3.1), ncol=2))

##############
### using their rescale step to normalize W by H max doesn't help
##############

all.markers2 <- patternMarkers_m(x@w, x@h, "all", rescale=TRUE, axis=1)
saveRDS(all.markers2, "snRNAseq_hpc/processed-data/revision/patternMarkers_CoGAPS_all.rda")
all.markers2 <- readRDS("snRNAseq_hpc/processed-data/revision/patternMarkers_CoGAPS_all.rda")

# reduces number of marker genes in most cases but increases in other cases
plot(sapply(all.markers1[[1]], length), sapply(all.markers2[[1]], length), 
     xlab="# markers per NMF (no rescale)", ylab="# markers per NMF (rescale)")
### remove nmf1 which has >15k marker genes
plot(sapply(all.markers1[[1]][2:100], length), sapply(all.markers2[[1]][2:100], length), 
     xlab="# markers per NMF (no rescale)", ylab="# markers per NMF (rescale)")
abline(a=0, b=1)

# weights ~ expr scatter plot with patternMarkers results labeled
plist4 <- lapply(c("nmf44","nmf37","nmf55"), function(x) {
  tmp = df1[,c(x,"avg.expr")]
  colnames(tmp) = c("nmf","avg.expr")
  tmp$is_top = rownames(tmp) %in% all.markers2[[1]][[x]]
  ggplot(tmp, aes(x=avg.expr, y=nmf, color=is_top))+
    geom_point(size=.3)+theme_bw()+
    geom_point(data=filter(tmp, is_top==T), size=.5)+
    scale_color_manual(values=c("black","red"))+
    labs(title=x, y="gene weights", x="avg. expr")+
    theme(text=element_text(size=10), axis.text.y=element_blank(), 
          legend.position="none", aspect.ratio=1)
})

# density plot of weights for patternMarker results genes vs other genes
plist5 <- lapply(c("nmf44","nmf37","nmf55"), function(x) {
  tmp = df1[,c(x,"avg.expr")]
  colnames(tmp) = c("nmf","avg.expr")
  tmp$is_top = rownames(tmp) %in% all.markers2[[1]][[x]]
  ggplot(tmp, aes(x=nmf, color=is_top))+
    geom_density(linewidth=.7)+theme_bw()+
    scale_color_manual(values=c("black","red"))+
    scale_x_log10()+
    labs(title=x, x="gene weights")+
    theme(text=element_text(size=10), axis.text.y=element_blank(),
          legend.position="none", aspect.ratio=1)
})
do.call(grid.arrange, c(c(plist4, plist5), ncol=3))


##############
### using the "cut" threshold is also not helpful
##############

cut.markers <- patternMarkers_m(x@w, x@h, "cut", rescale=FALSE, axis=1)
saveRDS(cut.markers, "snRNAseq_hpc/processed-data/revision/patternMarkers_CoGAPS_cut-rescale-FALSE.rda")
cut.markers <- readRDS("snRNAseq_hpc/processed-data/revision/patternMarkers_CoGAPS_cut-rescale-FALSE.rda")

# produces waaay fewer markers
plot(sapply(all.markers1[[1]], length), sapply(cut.markers[[1]], length), 
     xlab="# markers per NMF ('all')", ylab="# markers per NMF ('cut')")
### remove nmf1 with >15k markers
plot(sapply(all.markers1[[1]][2:100], length), sapply(cut.markers[[1]][2:100], length), 
     xlab="# markers per NMF ('all')", ylab="# markers per NMF ('cut')")
abline(a=0, b=1)

# weights ~ expr scatter plot with patternMarkers results labeled
plist6 <- lapply(c("nmf44","nmf37","nmf55"), function(x) {
  tmp = df1[,c(x,"avg.expr")]
  colnames(tmp) = c("nmf","avg.expr")
  tmp$is_top = rownames(tmp) %in% cut.markers[[1]][[x]]
  ggplot(tmp, aes(x=avg.expr, y=nmf, color=is_top))+
    geom_point(size=.3)+theme_bw()+
    geom_point(data=filter(tmp, is_top==T), size=.5)+
    scale_color_manual(values=c("black","red"))+
    labs(title=x, y="gene weights", x="avg. expr")+
    theme(text=element_text(size=10), axis.text.y=element_blank(), 
          legend.position="none", aspect.ratio=1)
})

# density plot of weights for patternMarker results genes vs other genes
plist7 <- lapply(c("nmf44","nmf37","nmf55"), function(x) {
  tmp = df1[,c(x,"avg.expr")]
  colnames(tmp) = c("nmf","avg.expr")
  tmp$is_top = rownames(tmp) %in% cut.markers[[1]][[x]]
  ggplot(tmp, aes(x=nmf, color=is_top))+
    geom_density(linewidth=.7)+theme_bw()+
    scale_color_manual(values=c("black","red"))+
    scale_x_log10()+
    labs(title=x, x="gene weights")+
    theme(text=element_text(size=10), axis.text.y=element_blank(),
          legend.position="none", aspect.ratio=1)
})
do.call(grid.arrange, c(c(plist6, plist7), ncol=3))

###########################################################################
######################### rank by weight  markers ######################### 

n_top = 100
# weight ~ expr scatter
plist8 <- lapply(c("nmf44","nmf37","nmf55"), function(x) {
  tmp = df1[,c(x,"avg.expr")]
  colnames(tmp) = c("nmf","avg.expr")
  tmp$rank = (1+nrow(loads))-rank(tmp$nmf)
  tmp$is_top = tmp$rank<= n_top
  ggplot(tmp, aes(x=avg.expr, y=nmf, color=is_top))+
    geom_point(size=.3)+theme_bw()+
    geom_point(data=filter(tmp, n_top==T), size=.5)+
    scale_color_manual(values=c("black","red"))+
    labs(title=paste(x,"- top", n_top), y="gene weights", x="avg. expr")+
    theme(text=element_text(size=10), axis.text.y=element_blank(), 
          legend.position="none", aspect.ratio=1)
})

# density
plist9 <- lapply(c("nmf44","nmf37","nmf55"), function(x) {
  tmp = df1[,c(x,"avg.expr")]
  colnames(tmp) = c("nmf","avg.expr")
  tmp$rank = (1+nrow(loads))-rank(tmp$nmf)
  tmp$is_top = tmp$rank<= n_top
  ggplot(tmp, aes(x=nmf, color=is_top))+
    geom_density(linewidth=.7)+theme_bw()+
    scale_color_manual(values=c("black","red"))+
    scale_x_log10()+
    labs(title=paste(x,"- top", n_top), x="gene weights")+
    theme(text=element_text(size=10), axis.text.y=element_blank(),
          legend.position="none", aspect.ratio=1)
})

do.call(grid.arrange, c(c(plist8, plist9), ncol=3))

# additional examples
# weight ~ expr scatter
plist8.1 <- lapply(c("nmf10","nmf48"), function(x) {
  tmp = df1[,c(x,"avg.expr")]
  colnames(tmp) = c("nmf","avg.expr")
  tmp$rank = (1+nrow(loads))-rank(tmp$nmf)
  tmp$is_top = tmp$rank<= n_top
  ggplot(tmp, aes(x=avg.expr, y=nmf, color=is_top))+
    geom_point(size=.3)+theme_bw()+
    geom_point(data=filter(tmp, n_top==T), size=.5)+
    scale_color_manual(values=c("black","red"))+
    labs(title=paste(x,"- top", n_top), y="gene weights", x="avg. expr")+
    theme(text=element_text(size=10), axis.text.y=element_blank(), 
          legend.position="none", aspect.ratio=1)
})

# density
plist9.1 <- lapply(c("nmf10","nmf48"), function(x) {
  tmp = df1[,c(x,"avg.expr")]
  colnames(tmp) = c("nmf","avg.expr")
  tmp$rank = (1+nrow(loads))-rank(tmp$nmf)
  tmp$is_top = tmp$rank<= n_top
  ggplot(tmp, aes(x=nmf, color=is_top))+
    geom_density(linewidth=.7)+theme_bw()+
    scale_color_manual(values=c("black","red"))+
    scale_x_log10()+
    labs(title=paste(x,"- top", n_top), x="gene weights")+
    theme(text=element_text(size=10), axis.text.y=element_blank(),
          legend.position="none", aspect.ratio=1)
})

do.call(grid.arrange, c(c(plist8.1, plist9.1), ncol=2))
