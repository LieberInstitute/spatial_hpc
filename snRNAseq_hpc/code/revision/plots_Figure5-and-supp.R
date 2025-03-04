library(SingleCellExperiment)
library(scater)
library(pheatmap)
library(dplyr)
library(ggplot2)
library(ggrepel)

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

###some duplicate barcodes, just paste sample to uniquify
# this step was included here: https://github.com/LieberInstitute/spatial_hpc/blob/main/snRNAseq_hpc/code/nmf/cross_species_NMF.R
colnames(sce.ms)<-paste(sce.ms$Sample,colnames(sce.ms),sep='_')

## standard QC approach 
non0.nuc = colSums(as.matrix(colData(sce.ms)[,paste0("nmf",1:100)])>0)
plot(ecdf(log10(non0.nuc)), xlim=c(0,5), xlab="# non-zero weighted nuclei per NMF")

remove.nmf = names(non0.nuc[non0.nuc<1000])
length(remove.nmf)
plot(ecdf(log10(non0.nuc)), xlim=c(0,5), xlab="# non-zero weighted spots per NMF",
     col=c(rep("red",13), rep("black",100-13)), 
     main="NMF with <1000 spots removed")

remove.nmf[1] <- c("nmf2")
remove.nmf = c(remove.nmf, "nmf37", "nmf28")


nmf.ordered = c("nmf37", "nmf9", "nmf3", "nmf1", "nmf55", "nmf56", "nmf57", "nmf58", "nmf59", "nmf75", "nmf80", "nmf91", "nmf92", "nmf94", "nmf24",#broad/non-specific
                "nmf2", "nmf6", "nmf8", "nmf18", "nmf20", "nmf71", "nmf72", #nrn not glia
                "nmf12", "nmf28", 'nmf7',#inhb neurons
                "nmf31","nmf4", "nmf13", "nmf16", "nmf21", "nmf25", "nmf34", "nmf49", #exc nrn
                "nmf26", "nmf10", "nmf14", "nmf5", "nmf66", #GCs
                "nmf52", "nmf11", "nmf63", "nmf61", #MC, CA3, CA2
                "nmf15", "nmf32", #CA1, ProS,
                "nmf40", "nmf54", #Sub.1, Sub.2
                "nmf65", "nmf22", "nmf53", "nmf68", "nmf51", #L6, L5
                "nmf45", "nmf84", "nmf27", "nmf17", "nmf78", #L2/3
                "nmf62", "nmf69", "nmf29", "nmf43", "nmf64", #amy
                "nmf23", #thal, cajal
                "nmf67", "nmf73", "nmf83", "nmf50", "nmf35", "nmf47", "nmf88", "nmf46", "nmf60", "nmf74", "nmf86", "nmf93", #GABA 
                "nmf81", "nmf19", "nmf76", "nmf79", #astro
                "nmf42", "nmf44", "nmf38", "nmf77", "nmf33", "nmf36", #oligo
                "nmf90", "nmf39", "nmf98", "nmf82", "nmf96", "nmf100", #micro immune
                "nmf87", "nmf30", "nmf41", "nmf48", #ependy, CP
                "nmf97", "nmf89", "nmf70", "nmf85", "nmf99", "nmf95" #endo
)
nmf.ordered.keep = setdiff(nmf.ordered, remove.nmf)
nmf.ordered.remove = intersect(nmf.ordered, remove.nmf)
nmf.ordered.remove = nmf.ordered.remove[c(1,3,2,4:16)]

#################
## dotplot ecs nmf filter
#################
seed1 = as.matrix(colData(sce.ms)[,paste0("nmf",1:100)])
seed1 = seed1>0
d1 = cbind.data.frame(cell.type.mouse=as.data.frame(colData(sce.ms))[,"cell.type.mouse"],
                      seed1) %>% 
  group_by(cell.type.mouse) %>% add_tally(name="total") %>%
  group_by(cell.type.mouse, total) %>%
  summarise_at(nmf.ordered.keep, sum) %>%
  tidyr::pivot_longer(nmf.ordered.keep, values_to="n", names_to="nmf") %>%
  mutate(prop=n/total)

seed2 = as.matrix(colData(sce.ms)[,paste0("nmf",1:100)])
seed2 = apply(seed2, 2, scale)
d2 = cbind.data.frame(cell.type.mouse=as.data.frame(colData(sce.ms))[,"cell.type.mouse"],
                      seed2) %>% 
  group_by(cell.type.mouse) %>%
  summarise_at(nmf.ordered.keep, mean) %>% 
  tidyr::pivot_longer(nmf.ordered.keep, values_to="scaled.avg", names_to="nmf")


dot.df = left_join(d1[,c("cell.type.mouse","nmf","prop")], 
                   d2[,c("cell.type.mouse","nmf","scaled.avg")]) %>%
  mutate(cell.type.mouse=factor(cell.type.mouse, 
                                levels=levels(cell.type.mouse)[c(1:4,6:7,5,8:12)]))

dot.df$nmf_f = factor(dot.df$nmf, levels=nmf.ordered.keep)

ggplot(dot.df, aes(x=nmf_f, y=cell.type.mouse, size=prop, color=scaled.avg))+
  geom_count()+theme_bw()+
  scale_size(range=c(0,3))+scale_color_viridis_c(option="F", direction=-1)+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5))

#################
## dotplot look at only ones that are well mapped to spots as well
#################
load(file=here::here('processed-data','NMF','spe_nmf_final.rda'))
non0.spots = colSums(as.matrix(colData(spe)[,paste0("nmf",1:100)])>0)
remove.nmf.spe = names(non0.spots[non0.spots<1050])
remove.nmf.spe[c(2,3,6)] <- c("nmf2","nmf3","nmf16")

remove.nmf = union(remove.nmf, remove.nmf.spe)
nmf.ordered.keep = setdiff(nmf.ordered, remove.nmf)

d3 = cbind.data.frame(as.data.frame(colData(sce.ms))[,c("cell.type.mouse","condition")],
                      seed1) %>% 
  group_by(cell.type.mouse, condition) %>% add_tally(name="total") %>%
  group_by(cell.type.mouse, condition, total) %>%
  summarise_at(nmf.ordered.keep, sum) %>%
  tidyr::pivot_longer(nmf.ordered.keep, values_to="n", names_to="nmf") %>%
  mutate(prop=n/total)

d4 = cbind.data.frame(as.data.frame(colData(sce.ms))[,c("cell.type.mouse","condition")],
                      seed2) %>% 
  group_by(cell.type.mouse, condition) %>%
  summarise_at(nmf.ordered.keep, mean) %>% 
  tidyr::pivot_longer(nmf.ordered.keep, values_to="scaled.avg", names_to="nmf")


dot.df2 = left_join(d3[,c("cell.type.mouse","condition","nmf","prop")], 
                   d4[,c("cell.type.mouse","condition","nmf","scaled.avg")]) %>%
  mutate(ylab=paste(cell.type.mouse, condition), 
         ylab_f = factor(ylab, levels=paste(rep(levels(sce.ms$cell.type.mouse)[c(1:4,6:7,5,8:12)],each=2), 
                                            c("Sham","ECS"))))

dot.df2$nmf_f = factor(dot.df2$nmf, levels=nmf.ordered.keep)

ggplot(dot.df2, aes(x=nmf_f, y=ylab_f, size=prop, color=scaled.avg))+
  geom_count()+theme_bw()+
  scale_size(range=c(0,3))+scale_color_viridis_c(option="F", direction=-1)+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5))


#################
## small dotplot of only patterns we highlight 
#################
seed3 = as.matrix(colData(sce.ms)[,paste0("nmf",1:100)])
seed3 = seed3>0
d3 = cbind.data.frame(as.data.frame(colData(sce.ms))[,c("cell.type.mouse","condition")],
                      seed3) %>% 
  group_by(cell.type.mouse, condition) %>% add_tally(name="total") %>%
  group_by(cell.type.mouse, condition, total) %>%
  summarise_at(all_of(c("nmf14","nmf10","nmf55","nmf91","nmf20")), sum) %>%
  tidyr::pivot_longer(c("nmf14","nmf10","nmf55","nmf91","nmf20"), values_to="n", names_to="nmf") %>%
  mutate(prop=n/total)

seed4 = as.matrix(colData(sce.ms)[,paste0("nmf",1:100)])
seed4 = apply(seed4, 2, scale)
d4 = cbind.data.frame(as.data.frame(colData(sce.ms))[,c("cell.type.mouse","condition")],
                      seed4) %>% 
  group_by(cell.type.mouse, condition) %>%
  summarise_at(c("nmf14","nmf10","nmf55","nmf91","nmf20"), mean) %>% 
  tidyr::pivot_longer(c("nmf14","nmf10","nmf55","nmf91","nmf20"), values_to="scaled.avg", names_to="nmf")


dot.df2 = left_join(d3[,c("cell.type.mouse","condition","nmf","prop")], 
                    d4[,c("cell.type.mouse","condition","nmf","scaled.avg")]) %>%
  mutate(ylab=paste(cell.type.mouse, condition), 
         ylab_f = factor(ylab, levels=paste(rep(levels(sce.ms$cell.type.mouse)[c(1:4,6:7,5,8:12)],each=2), 
                                            c("Sham","ECS"))))

dot.df2$nmf_f = factor(dot.df2$nmf, levels=c("nmf91","nmf20","nmf55","nmf10","nmf14"))

dot.df2$ylab_c = as.numeric(dot.df2$ylab_f)
dot.df2$ylab_c = ifelse(dot.df2$condition=="ECS", dot.df2$ylab_c-.25, dot.df2$ylab_c)

p1 <- ggplot(filter(dot.df2, nmf_f!="nmf55"), aes(x=nmf_f, y=ylab_c, size=prop, color=scaled.avg))+
  geom_count()+theme_bw()+
  scale_size(range=c(0,3))+scale_color_viridis_c(option="F", direction=-1)+
  scale_y_continuous(breaks=unique(dot.df2$ylab_c), labels=levels(dot.df2$ylab_f), minor_breaks = NULL)+
  labs(x="",y="")+
  theme(text=element_text(size=14), axis.text.x=element_text(angle=90, hjust=1, vjust=.5))

pdf(file = "snRNAseq_hpc/plots/revision/Fig5_small-dotplot.pdf",
    width=4, height=7)
p1
dev.off()

#################
## volcanoes
#################

nmf91.volcano = read.csv("snRNAseq_hpc/processed-data/revision/nmf91_sham-ecs_results.csv")
p1 <- ggplot(nmf91.volcano, aes(x=summary.logFC, y=-log10(FDR), color=nmf91))+
  geom_point(size=2)+xlim(-6.6,6.6)+
  geom_point(data=slice_max(nmf91.volcano, n=50, nmf91), size=2)+
  geom_text_repel(data=filter(nmf91.volcano, nmf91>.0015, summary.logFC>1, -log10(FDR)>30),
                  aes(label=gene), color="black", min.segment.length = 0, force = 5,
                  size=4.5)+
  theme_bw()+scale_color_viridis_c(option="F", direction=-1, labels=function(x) format(x, scientific=T, digits=1))+
  labs(x="log2FC (sham GCs vs ECS GCs)")+theme(text=element_text(size=12), legend.key.size = unit(10,"pt"), legend.box.spacing=unit(6,"pt"))

nmf20.volcano = read.csv("snRNAseq_hpc/processed-data/revision/nmf20_sham-ecs_results.csv")
p2 <- ggplot(nmf20.volcano, aes(x=summary.logFC, y=-log10(FDR), color=nmf20))+
  geom_point(size=2)+xlim(-6,6)+
  geom_point(data=slice_max(nmf20.volcano, n=150, nmf20), size=2)+
  geom_text_repel(data=filter(nmf20.volcano, nmf20>.00065, summary.logFC>1, -log10(FDR)>30),
                  #data=filter(nmf20.volcano, gene %in% c("SORCS3","PDE10A","HOMER1")),
                  aes(label=gene), color="black", min.segment.length = 0, force = 5,
                  size=4.5)+
  theme_bw()+scale_color_viridis_c(option="F", direction=-1, labels=function(x) format(x, scientific=T, digits=1))+
  labs(x="log2FC (sham GCs vs ECS GCs)")+theme(text=element_text(size=12), legend.key.size = unit(10,"pt"), legend.box.spacing=unit(6,"pt"))

nmf55.volcano = read.csv("snRNAseq_hpc/processed-data/revision/nmf55_sham-ecs_results.csv")
p3 <- ggplot(nmf55.volcano, aes(x=summary.logFC, y=-log10(FDR), color=nmf55))+
  geom_point(size=2)+xlim(-4.2,4.2)+
  geom_point(data=slice_max(nmf55.volcano, n=150, nmf55), size=2)+
  geom_hline(aes(yintercept=-log10(.05)), lty=2, color="grey30")+
  geom_text_repel(data=filter(nmf55.volcano, nmf55>.0025, summary.logFC>0.5, FDR<.05),
                  aes(label=gene), color="black", min.segment.length = 0, force = 5, 
                  max.overlaps=Inf, size=4.5)+
  theme_bw()+scale_color_viridis_c(option="F", direction=-1, labels=function(x) format(x, scientific=T, digits=1))+
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10))+
  labs(x="log2FC (sham GCs vs ECS GCs)")+theme(text=element_text(size=12), legend.key.size = unit(10,"pt"), legend.box.spacing=unit(6,"pt"))

nmf10.volcano = read.csv("snRNAseq_hpc/processed-data/revision/nmf10_sham-ecs_results.csv")
p4 <- ggplot(nmf10.volcano, aes(x=summary.logFC, y=-log10(FDR), color=nmf10))+
  geom_point(size=2)+xlim(-6,6)+
  geom_point(data=slice_max(nmf10.volcano, n=100, nmf10), size=2)+
  geom_point(data=filter(nmf10.volcano, nmf10>.00065, summary.logFC< -1, -log10(FDR)>30), size=2)+
  geom_text_repel(data=filter(nmf10.volcano, nmf10>.00065, summary.logFC< -1, -log10(FDR)>30),
                  aes(label=gene), color="black", min.segment.length = 0, force = 5,
                  size=4.5)+
  theme_bw()+scale_color_viridis_c(option="F", direction=-1, labels=function(x) format(x, scientific=T, digits=1))+
  labs(x="log2FC (sham GCs vs ECS GCs)")+theme(text=element_text(size=12), legend.key.size = unit(10,"pt"), legend.box.spacing=unit(6,"pt"))

nmf14.volcano <- read.csv("snRNAseq_hpc/processed-data/revision/nmf14_sham-ecs_results.csv")
p5 <- ggplot(nmf14.volcano, aes(x=summary.logFC, y=-log10(FDR), color=nmf14))+
  geom_point(size=2)+xlim(-6,6)+
  geom_point(data=slice_max(nmf14.volcano, n=100, nmf14), size=2)+
  geom_point(data=filter(nmf14.volcano, gene %in% c("ACVR1C","BDNF","SGK1")), size=2)+
  geom_text_repel(data=filter(nmf14.volcano, nmf14>.0005, summary.logFC>1, -log10(FDR)>30),
                  aes(label=gene), color="black", min.segment.length = 0, force = 5,
                  size=4.5)+
  theme_bw()+scale_color_viridis_c(option="F", direction=-1, labels=function(x) format(x, scientific=T, digits=1))+
  labs(x="log2FC (sham GCs vs ECS GCs)")+theme(text=element_text(size=12), legend.key.size = unit(10,"pt"), legend.box.spacing=unit(6,"pt"))

pdf(file = "snRNAseq_hpc/plots/revision/Fig5_volcanoes.pdf",
    width=13, height=7)
grid.arrange(p1, p2, p3, p4, p5, ncol=3)
dev.off()

##################################
######## GC subgroups plots
##################################

load(file=here::here('snRNAseq_hpc','processed-data','sce','sce_final.rda'))

plotGroupedHeatmap(sce, features=c("FAT4","CDH10","CHST9","CNTN1","ADAM22","ST6GALNAC5"), group="gc_labels",
                   scale=T, center=T)

sce$gc_labels = ifelse(sce$fine.cell.class=="GC", as.character(sce$superfine.cell.class), "other")

plotReducedDim(sce, dimred="UMAP", color_by="nmf10", point_size=.1)+
  scale_color_viridis_c(option="F", direction=-1, labels=function(x) format(x, scientific=T, digits=1))+
  ggtitle("nmf10")+theme(plot.title=element_text(size=18, hjust=.5))
plotReducedDim(sce, dimred="UMAP", color_by="nmf14", point_size=.1)+
  scale_color_viridis_c(option="F", direction=-1, labels=function(x) format(x, scientific=T, digits=1))+
  ggtitle("nmf14")+theme(plot.title=element_text(size=18, hjust=.5))

plotReducedDim(sce, dimred="UMAP", color_by="CHST9", point_size=.1)+
  scale_color_viridis_c(option="F", direction=-1)+ggtitle("CHST9")+
  theme(plot.title=element_text(size=18, hjust=.5))
plotReducedDim(sce, dimred="UMAP", color_by="SORCS3", point_size=.1)+
  scale_color_viridis_c(option="F", direction=-1)+ggtitle("SORCS3")+
  theme(plot.title=element_text(size=18, hjust=.5))

#subset to gcs only
sce_gc = sce[,sce$gc_labels!="other"]
dim(sce_gc) #36601 10508
umap.emb = reducedDim(sce_gc, "UMAP")
remove.extra = union(rownames(umap.emb)[umap.emb[,"UMAP1"]<5],rownames(umap.emb)[umap.emb[,"UMAP2"]< -3])
sce_gc = sce_gc[,!colnames(sce_gc) %in% remove.extra]
dim(sce_gc) #36601 10506

p1 <- plotReducedDim(sce, dimred="UMAP", color_by="gc_labels", point_size=.1)+
  scale_color_manual(values=c(RColorBrewer::brewer.pal(n=5,"Set1"),"grey"))+
  theme_void()+theme(legend.position="none", axis.line=element_line(color="black", linewidth=.5), 
                     plot.margin=unit(c(5,5,5,5), units="pt"))
p2 <- plotReducedDim(sce_gc, dimred="UMAP", color_by="superfine.cell.class", point_size=.3)+
  scale_color_brewer(palette="Set1")+
  theme_void()+theme(legend.key.width=unit(5, units="pt"),
                     legend.text = element_text(size=8, margin=margin(l=3,unit="pt")), legend.title=element_blank(), 
                     axis.line=element_line(color="black", linewidth=.5), plot.margin=unit(c(5,5,5,5), units="pt"))

ggsave("snRNAseq_hpc/plots/revision/Figure5_gc-cluster-umap.png", gridExtra::grid.arrange(p1, p2, ncol=2), 
       bg="white", width=4, height=2, units="in")


p3 <- plotReducedDim(sce_gc, dimred="UMAP", color_by="nmf10", point_size=.3)+
  scale_color_viridis_c("nmf10", option="F", direction=-1, labels=function(x) format(x, scientific=T, digits=1),
			breaks=c(0, 1e-4, 2e-4))+
  theme_void()+theme(legend.key.width=unit(5, units="pt"), legend.title=element_text(size=8),
                     legend.text = element_text(size=7, margin=margin(l=2,unit="pt")),
                     axis.line=element_line(color="black", linewidth=.5), plot.margin=unit(c(5,0,5,5), units="pt"))

p4 <- plotReducedDim(sce_gc, dimred="UMAP", color_by="nmf14", point_size=.3)+
  scale_color_viridis_c("nmf14", option="F", direction=-1, labels=function(x) format(x, scientific=T, digits=1),
			breaks=c(0, 2e-4, 4e-4))+
  theme_void()+theme(legend.key.width=unit(5, units="pt"), legend.title=element_text(size=8),
                     legend.text = element_text(size=7, margin=margin(l=2,unit="pt")),
                     axis.line=element_line(color="black", linewidth=.5), plot.margin=unit(c(5,0,5,5), units="pt"))


p5 <- plotReducedDim(sce_gc, dimred="UMAP", color_by="CHST9", point_size=.1)+
  scale_color_gradient(low="grey90", high="black")+labs(color="CHST9")+
  theme_void()+theme(legend.key.width=unit(5, units="pt"), legend.title=element_text(size=8),
                     legend.text = element_text(size=8, margin=margin(l=3,unit="pt")),
                     axis.line=element_line(color="black", linewidth=.5), plot.margin=unit(c(5,5,5,5), units="pt"))

p6 <- plotReducedDim(sce_gc, dimred="UMAP", color_by="SORCS3", point_size=.1)+
  scale_color_gradient(low="grey90", high="black")+labs(color="SORCS3")+
  theme_void()+theme(legend.key.width=unit(5, units="pt"), legend.title=element_text(size=8),
                     legend.text = element_text(size=8, margin=margin(l=3,unit="pt")),
                     axis.line=element_line(color="black", linewidth=.5), plot.margin=unit(c(5,5,5,5), units="pt"))

ggsave("snRNAseq_hpc/plots/revision/Figure5_nmf10-chst9-umap.png", gridExtra::grid.arrange(p3, p5, ncol=2), 
       bg="white", width=4, height=2, units="in")

ggsave("snRNAseq_hpc/plots/revision/Figure5_nmf14-sorcs3-umap.png", gridExtra::grid.arrange(p4, p6, ncol=2), 
       bg="white", width=4, height=2, units="in")
