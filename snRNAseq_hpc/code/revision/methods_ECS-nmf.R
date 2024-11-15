library(SingleCellExperiment)
library(scater)
library(pheatmap)
library(dplyr)
library(ggplot2)

set.seed(123)

sce.ms = readRDS("snRNAseq_hpc/processed-data/sce/sce_mouse_ECS.rds")
# scaling nmf nuclei weights like he did with our snRNAseq and SRT spots
tmp = as.matrix(colData(sce.ms)[,paste0("V",1:100)])
tmp<-apply(tmp,2,function(x){x/sum(x)})
colnames(tmp) <- paste0("nmf",1:100)
#colSums(tmp)
colData(sce.ms) <- cbind(colData(sce.ms), tmp)

#GCs auto split into act and non act so must combine
sce.ms$cell.type.mouse<-sce.ms$cellType
levels(sce.ms$cell.type.mouse)[1]<-'GC'
levels(sce.ms$cell.type.mouse)[2]<-'GC'
levels(sce.ms$cell.type.mouse)[2:5]<-'CA2-4'
levels(sce.ms$cell.type.mouse)[4:6]<-'PS/Sub'

#load nmfs
load(file=here::here('snRNAseq_hpc','processed-data','NMF','nmf_final.rda'))
loads<-x@w
no_expr <- which(rowSums(loads) == 0)
loads <- loads[-no_expr, ]

#de analysis
non0.nmf91 = rownames(loads)[loads[,"nmf91"]>0]
sce.ms.gcs = sce.ms[intersect(rownames(sce.ms), non0.nmf91),sce.ms$cell.type.mouse=="GC"]
exc.gc.markers = scran::findMarkers(sce.ms.gcs, groups=sce.ms.gcs$condition, test.type="binom")
nmf91.volcano = tibble::rownames_to_column(cbind.data.frame(exc.gc.markers[[2]], 
                                                            nmf91=loads[rownames(exc.gc.markers[[2]]),"nmf91"]), var="gene")
write.csv(nmf91.volcano, "snRNAseq_hpc/processed-data/revision/nmf91_sham-ecs_results.csv", row.names=F)

non0.nmf20 = rownames(loads)[loads[,"nmf20"]>0]
sce.ms.gcs = sce.ms[intersect(rownames(sce.ms), non0.nmf20),sce.ms$cell.type.mouse=="GC"]
exc.gc.markers = scran::findMarkers(sce.ms.gcs, groups=sce.ms.gcs$condition, test.type="binom")
nmf20.volcano = tibble::rownames_to_column(cbind.data.frame(exc.gc.markers[[2]], 
                                                            nmf20=loads[rownames(exc.gc.markers[[2]]),"nmf20"]), var="gene")
write.csv(nmf20.volcano, "snRNAseq_hpc/processed-data/revision/nmf20_sham-ecs_results.csv", row.names=F)

non0.nmf55 = rownames(loads)[loads[,"nmf55"]>0]
sce.ms.gcs = sce.ms[intersect(rownames(sce.ms), non0.nmf55),sce.ms$cell.type.mouse=="GC"]
exc.gc.markers = scran::findMarkers(sce.ms.gcs, groups=sce.ms.gcs$condition, test.type="binom")
nmf55.volcano = tibble::rownames_to_column(cbind.data.frame(exc.gc.markers[[2]], 
                                                            nmf55=loads[rownames(exc.gc.markers[[2]]),"nmf55"]), var="gene")
write.csv(nmf55.volcano, "snRNAseq_hpc/processed-data/revision/nmf55_sham-ecs_results.csv", row.names=F)

non0.nmf10 = rownames(loads)[loads[,"nmf10"]>0]
sce.ms.gcs = sce.ms[intersect(rownames(sce.ms), non0.nmf10),sce.ms$cell.type.mouse=="GC"]
exc.gc.markers = scran::findMarkers(sce.ms.gcs, groups=sce.ms.gcs$condition, test.type="binom")
nmf10.volcano = tibble::rownames_to_column(cbind.data.frame(exc.gc.markers[[2]], 
                                                            nmf10=loads[rownames(exc.gc.markers[[2]]),"nmf10"]), var="gene")
write.csv(nmf10.volcano, "snRNAseq_hpc/processed-data/revision/nmf10_sham-ecs_results.csv", row.names=F)

non0.nmf14 = rownames(loads)[loads[,"nmf14"]>0]
sce.ms.gcs = sce.ms[intersect(rownames(sce.ms), non0.nmf14),sce.ms$cell.type.mouse=="GC"]
exc.gc.markers = scran::findMarkers(sce.ms.gcs, groups=sce.ms.gcs$condition, test.type="binom")
nmf14.volcano = tibble::rownames_to_column(cbind.data.frame(exc.gc.markers[[2]], nmf14=loads[rownames(exc.gc.markers[[2]]),"nmf14"]), var="gene")
write.csv(nmf14.volcano, "snRNAseq_hpc/processed-data/revision/nmf14_sham-ecs_results.csv", row.names=F)

