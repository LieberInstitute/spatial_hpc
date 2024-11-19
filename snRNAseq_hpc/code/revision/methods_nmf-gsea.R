library(dplyr)
library(ggplot2)
library(fgsea)
library(reactome.db)
library(org.Hs.eg.db)

set.seed(123)

#load nmfs
load(file="snRNAseq_hpc/processed-data/NMF/nmf_final.rda")
loads<-x@w
no_expr <- which(rowSums(loads) == 0)
loads <- loads[-no_expr, ]

# prep reactome
xx <- as.list(reactomePATHID2NAME)
reactome.h <- xx[grep("^Homo",xx)]
x <- as.list(reactomePATHID2EXTID)
reactome.h = reactome.h[intersect(names(reactome.h), names(x))]
x.h <- x[names(reactome.h)]
identical(names(x.h), names(reactome.h))
reactome.h = gsub("Homo sapiens: ","",reactome.h)
names(x.h) = reactome.h


################# GSEA

#### oligodendrocyte
non0.nmf44 = rownames(loads)[loads[,"nmf44"]>0]
non0.44.id = mapIds(org.Hs.eg.db, keys=non0.nmf44, keytype="SYMBOL", column="ENTREZID", multiVals = "first")
names(non0.44.id) = non0.nmf44
non0.44.id = non0.44.id[!is.na(non0.44.id)]
pathways.44 <- reactomePathways(non0.44.id)
pathways.44 <- x.h[names(pathways.44)]

nmf44.stats = loads[names(non0.44.id),"nmf44"]
names(nmf44.stats) = non0.44.id
#in order to get reproducible results, must call set.seed every time before you run the function
### the function has an internal set seed that makes this necessary
set.seed(123)
olig.results.44 = fgseaMultilevel(pathways.44, stats=nmf44.stats, scoreType="pos", minSize=15, maxSize=500)
olig.results.44$leadingEdge2 = sapply(olig.results.44$leadingEdge, paste, collapse="/")
write.csv(olig.results.44[,c(1:7,9)], "snRNAseq_hpc/processed-data/revision/nmf44_reactome_results_3-sig.csv")



#### astrocyte
non0.nmf81 = rownames(loads)[loads[,"nmf81"]>0]
non0.81.id = mapIds(org.Hs.eg.db, keys=non0.nmf81, keytype="SYMBOL", column="ENTREZID", multiVals = "first")
names(non0.81.id) = non0.nmf81
non0.81.id = non0.81.id[!is.na(non0.81.id)]
pathways.81 <- reactomePathways(non0.81.id)
pathways.81 <- x.h[names(pathways.81)]

nmf81.stats = loads[names(non0.81.id),"nmf81"]
names(nmf81.stats) = non0.81.id
#in order to get reproducible results, must call set.seed every time before you run the function
### the function has an internal set seed that makes this necessary
set.seed(123)
astro.results.81 = fgseaMultilevel(pathways.81, stats=nmf81.stats, scoreType="pos", minSize=15, maxSize=500)
astro.results.81$leadingEdge2 = sapply(astro.results.81$leadingEdge, paste, collapse="/")
write.csv(astro.results.81[,c(1:7,9)], "snRNAseq_hpc/processed-data/revision/nmf81_reactome_results.csv")



#### exc postsynaptic
non0.nmf13 = rownames(loads)[loads[,"nmf13"]>0]
non0.13.id = mapIds(org.Hs.eg.db, keys=non0.nmf13, keytype="SYMBOL", column="ENTREZID", multiVals = "first")
names(non0.13.id) = non0.nmf13
non0.13.id = non0.13.id[!is.na(non0.13.id)]
pathways.13 <- reactomePathways(non0.13.id)
pathways.13 <- x.h[names(pathways.13)]

nmf13.stats = loads[names(non0.13.id),"nmf13"]
names(nmf13.stats) = non0.13.id
#in order to get reproducible results, must call set.seed every time before you run the function
### the function has an internal set seed that makes this necessary
set.seed(123)
nrn.results.13 = fgseaMultilevel(pathways.13, stats=nmf13.stats, scoreType="pos", minSize=15, maxSize=500)
nrn.results.13$leadingEdge2 = sapply(nrn.results.13$leadingEdge, paste, collapse="/")
write.csv(nrn.results.13[,c(1:7,9)], "snRNAseq_hpc/processed-data/revision/nmf13_reactome_results.csv")



#### inhb postsynaptic
non0.nmf7 = rownames(loads)[loads[,"nmf7"]>0]
non0.7.id = mapIds(org.Hs.eg.db, keys=non0.nmf7, keytype="SYMBOL", column="ENTREZID", multiVals = "first")
names(non0.7.id) = non0.nmf7
non0.7.id = non0.7.id[!is.na(non0.7.id)]
pathways.7 <- reactomePathways(non0.7.id)
pathways.7 <- x.h[names(pathways.7)]

nmf7.stats = loads[names(non0.7.id),"nmf44"]
names(nmf7.stats) = non0.7.id
#in order to get reproducible results, must call set.seed every time before you run the function
### the function has an internal set seed that makes this necessary
set.seed(123)
nrn.results.7 = fgseaMultilevel(pathways.7, stats=nmf7.stats, scoreType="pos", minSize=15, maxSize=500)
nrn.results.7$leadingEdge2 = sapply(nrn.results.7$leadingEdge, paste, collapse="/")
write.csv(nrn.results.7[,c(1:7,9)], "snRNAseq_hpc/processed-data/revision/nmf7_reactome_results.csv")


################# GSEA TABLES

nmf44.terms = c("L1CAM interactions","Nervous system development","Axon guidance")
plotGseaTable(x.h[nmf44.terms],
              nmf44.stats, olig.results.44, 
              gseaParam = 0.5)

nmf81.terms = c("Signaling by Receptor Tyrosine Kinases","Extracellular matrix organization","RHO GTPase cycle")
plotGseaTable(x.h[nmf81.terms],
              nmf81.stats, astro.results.81, 
              gseaParam = 0.5)

nmf13.terms = c("Transmission across Chemical Synapses","Neurotransmitter receptors and postsynaptic signal transmission","Protein-protein interactions at synapses")
plotGseaTable(x.h[nmf13.terms],
              nmf13.stats, nrn.results.13, 
              gseaParam = 0.5)

nmf7.terms = c("Transmission across Chemical Synapses","Neurotransmitter receptors and postsynaptic signal transmission","RHO GTPase Effectors")
plotGseaTable(x.h[nmf7.terms],
              nmf7.stats, nrn.results.7, 
              gseaParam = 0.5)


################# HEATMAP OF SELECT MARKER GENES
l1cam.genes = strsplit(filter(olig.results.44, pathway=="L1CAM interactions")$leadingEdge, split="/")[[1]]
l1cam.genes2 = mapIds(org.Hs.eg.db, keys=l1cam.genes, keytype="ENTREZID", column="SYMBOL", multiVals = "first")
rank.olig = (1+nrow(loads))-rank(loads[,"nmf44"])
top.olig = names(rank.olig)[rank.olig<=50]
intersect(top.olig, l1cam.genes2) #ANK3, DNM3, DLG1, SHTN1

ecm.genes = strsplit(filter(astro.results.81, pathway=="Extracellular matrix organization")$leadingEdge2, split="/")[[1]]
ecm.genes2 = mapIds(org.Hs.eg.db, keys=ecm.genes, keytype="ENTREZID", column="SYMBOL", multiVals = "first")
rank.astro = (1+nrow(loads))-rank(loads[,"nmf81"])
top.astro = names(rank.astro)[rank.astro<=100]
intersect(top.astro, ecm.genes2) #NRXN1 (top50), NCAN (top50), VCAN, BCAN, PLEC

n13.genes = filter(nrn.results.13, pathway=="Protein-protein interactions at synapses")$leadingEdge[[1]]
n13.genes2 = mapIds(org.Hs.eg.db, keys=n13.genes, keytype="ENTREZID", column="SYMBOL", multiVals = "first")
rank.n13 = (1+nrow(loads))-rank(loads[,"nmf13"])
top.n13 = names(rank.n13)[rank.n13<=50]
intersect(top.n13, n13.genes2) 
#top50 in nt receptors: "NRG1"    "RPS6KA2" "PRKAG2"  "ADCY9"   "CAMKK1"
#Top 50 in chem synaps trans: "NRG1"    "RPS6KA2" "PRKAG2"  "ADCY9"   "CAMKK1" 
#Top50 in p-p interact at synaps: "NRXN3"  "NTRK3"  "APBA2"  "LRFN2"  "DLGAP3"


n7.genes = filter(nrn.results.7, pathway=="RHO GTPase Effectors")$leadingEdge[[1]]
n7.genes2 = mapIds(org.Hs.eg.db, keys=n7.genes, keytype="ENTREZID", column="SYMBOL", multiVals = "first")
rank.n7 = (1+nrow(loads))-rank(loads[,"nmf7"])
top.n7 = names(rank.n7)[rank.n7<=50]
intersect(top.n7, n7.genes2) 
#Top50 in nt receptors: "CALM1"   "CALM3"   "GABRA1"  "PRKAR1A"
#Top50 in chem synaps trans: "CALM1"   "VAMP2"   "CALM3"   "GABRA1"  "PRKAR1A"
#Top50 in rho gtpase: "CALM1"  "CALM3"  "YWHAG"  "KIF5A"  "DYNLL2"

################# HEATMAP OF SELECT MARKER GENES
select.nmfs = c("nmf44","nmf81","nmf13","nmf7")
nmf.genes = c("ANK3", "DNM3", "DLG1", "SHTN1",
              "NRXN1","NCAN","PLEC",
              "CAMKK1","DLGAP3","LRFN2",
              #"PRKAG2","ADCY9","LRFN2","NTRK3",
              #"CALM1","CALM3","PRKAR1A","KIF5A",
              "GABRA1","KIF5A","DYNLL2")


m1 = loads[nmf.genes,select.nmfs]
pheatmap::pheatmap(m1, scale="row", cluster_rows = F, cluster_cols = F,
                   angle_col=0)
