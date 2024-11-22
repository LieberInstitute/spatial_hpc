library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(scater)

set.seed(123)

load(file=here::here('snRNAseq_hpc','processed-data','sce','sce_final.rda'))

sce_subset = sce[,sce$mid.cell.class=="ExcN"]
sce_subset$superfine.cell.class = droplevels(sce_subset$superfine.cell.class)

sce_subset$plot.cell.class = as.character(sce_subset$superfine.cell.class)
sce_subset$plot.cell.class = ifelse(sce_subset$fine.cell.class %in% c("HATA","Amy","Thal","Cajal","GC","CA2-4"), 
                                    as.character(sce_subset$fine.cell.class), sce_subset$plot.cell.class)
sce_subset$new.cell.class = factor(sce_subset$plot.cell.class,
                                   levels=c("GC","CA2-4","HATA","Amy","Thal","Cajal","L2/3.5",
                                            "CA1","ProS",
                                            "Sub.1","Sub.2","L6.1","L2/3.1",
                                            "L6b","L6.2","L5.2","L5.1",
                                            "L2/3.2","L2/3.3","L2/3.6","L2/3.4"),
                                   labels=c("GC","CA2-4","HATA","Amy","Thal","Cajal","L2/3.5",
                                            "CA1","ProS",
                                            "Sub.1","Sub.2","Sub.3","PreS",
                                            "RHP.L6b","RHP.L6","RHP.CBLN2+","ENT.L5",
                                            "ENT.sup3","ENT.sup2b","ENT.sup2a","ENT.sup1")
)

#subiculum DE
#sce_sub = sce_subset[,sce_subset$superfine.cell.class %in% c("CA1","Sub.1","Sub.2","ProS","L6.1","L2/3.1")]
sce_sub = sce_subset[,sce_subset$new.cell.class %in% c("CA1","Sub.1","Sub.2","ProS","Sub.3","PreS")]
sce_sub$new.cell.class = droplevels(sce_sub$new.cell.class)
sce_sub = sce_sub[rowSums(counts(sce_sub)>0)>100,] #22572  8783
sub.results = scran::findMarkers(sce_sub, group=sce_sub$new.cell.class, test="binom")

saveRDS(sub.results, "snRNAseq_hpc/processed-data/revision/subiculum-DE_sub1-sub2-pros.rda")

sub.1.list <- list("CA1"=rownames(filter(as.data.frame(sub.results[['Sub.1']]), -log10(FDR)>30 & logFC.CA1>1)),
                   "ProS"=rownames(filter(as.data.frame(sub.results[['Sub.1']]), -log10(FDR)>30 & logFC.ProS>1)),
                   "Sub.2"=rownames(filter(as.data.frame(sub.results[['Sub.1']]), -log10(FDR)>30 & logFC.Sub.2>1)),
                   "L6.1"=rownames(filter(as.data.frame(sub.results[['Sub.1']]), -log10(FDR)>30 & logFC.Sub.3>1)),
                   "L2/3.1"=rownames(filter(as.data.frame(sub.results[['Sub.1']]), -log10(FDR)>30 & logFC.PreS>1)))
sapply(sub.1.list, length)
#897    666    670   1198   1730 

t1 = table(unlist(sub.1.list))
sub1.unique = names(t1)[t1==length(sub.1.list)]
length(sub1.unique) #146

avg.expr = mutate(as.data.frame(sub.results[['Sub.1']][sub1.unique,]), 
                  avg.sub1.expr = rowMeans(logcounts(sce_sub)[sub1.unique,sce_sub$new.cell.class=="Sub.1"]))
length(rownames(filter(avg.expr, avg.sub1.expr>1))) #10
#[1] "AC007368.1" "AL138694.1" "ATP6V1C2"   "COL21A1"    "EBF4"      
#[6] "FN1"        "NDST4"      "PARD3B"     "PRKCH"      "RAPGEF3"
plot.list = list("Sub.1"=rownames(filter(avg.expr, avg.sub1.expr>1)))


sub.2.list <- list("CA1"=rownames(filter(as.data.frame(sub.results[['Sub.2']]), -log10(FDR)>30 & logFC.CA1>1)),
                   "ProS"=rownames(filter(as.data.frame(sub.results[['Sub.2']]), -log10(FDR)>30 & logFC.ProS>1)),
                   "Sub.1"=rownames(filter(as.data.frame(sub.results[['Sub.2']]), -log10(FDR)>30 & logFC.Sub.1>1)),
                   "L6.1"=rownames(filter(as.data.frame(sub.results[['Sub.2']]), -log10(FDR)>30 & logFC.Sub.3>1)),
                   "L2/3.1"=rownames(filter(as.data.frame(sub.results[['Sub.2']]), -log10(FDR)>30 & logFC.PreS>1)))
sapply(sub.2.list, length)
#1276   1071    717    798   1512

t2 = table(unlist(sub.2.list))
sub2.unique = names(t2)[t2==length(sub.2.list)]
length(sub2.unique) #184

avg.expr = mutate(as.data.frame(sub.results[['Sub.2']][sub2.unique,]), 
                  avg.sub2.expr = rowMeans(logcounts(sce_sub)[sub2.unique,sce_sub$new.cell.class=="Sub.2"]))
length(rownames(filter(avg.expr, avg.sub2.expr>1))) #7
#"GDNF-AS1" "LHFPL3"   "MAMDC2"   "PCED1B"   "RORB"     "SULF1"    "TRPC3"
plot.list$Sub.2 = rownames(filter(avg.expr, avg.sub2.expr>1))


pros.list <- list("CA1"=rownames(filter(as.data.frame(sub.results[['ProS']]), -log10(FDR)>30 & logFC.CA1>1)),
                  "Sub.2"=rownames(filter(as.data.frame(sub.results[['ProS']]), -log10(FDR)>30 & logFC.Sub.2>1)),
                  "Sub.1"=rownames(filter(as.data.frame(sub.results[['ProS']]), -log10(FDR)>30 & logFC.Sub.1>1)),
                  "L6.1"=rownames(filter(as.data.frame(sub.results[['ProS']]), -log10(FDR)>30 & logFC.Sub.3>1)),
                  "L2/3.1"=rownames(filter(as.data.frame(sub.results[['ProS']]), -log10(FDR)>30 & logFC.PreS>1)))
sapply(pros.list, length)
#218    583    321    807   1048

t3 = table(unlist(pros.list))
pros.unique = names(t3)[t3==length(pros.list)]
length(pros.unique) #17

avg.expr = mutate(as.data.frame(sub.results[['ProS']][pros.unique,]), 
                  avg.pros.expr = rowMeans(logcounts(sce_sub)[pros.unique,sce_sub$new.cell.class=="ProS"]))
length(rownames(filter(avg.expr, avg.pros.expr>1))) #0
plot.list$ProS = rownames(filter(avg.expr, avg.pros.expr>1)) 


#deep DE
sce_subd = sce_subset[,sce_subset$new.cell.class %in% c("Sub.1","Sub.2","Sub.3","RHP.L6b","RHP.L6")]
sce_subd$new.cell.class = droplevels(sce_subd$new.cell.class)
sce_subd = sce_subd[rowSums(counts(sce_subd)>0)>100,] #21037  5647
subd.results = scran::findMarkers(sce_subd, group=sce_subd$new.cell.class, test="binom")

saveRDS(subd.results, "snRNAseq_hpc/processed-data/revision/subiculum-DE_sub3-deep.rda")

sub.3.list <- list("RHP.L6b"=rownames(filter(as.data.frame(subd.results[['Sub.3']]), -log10(FDR)>30 & logFC.RHP.L6b>1)),
                 "Sub.2"=rownames(filter(as.data.frame(subd.results[['Sub.3']]), -log10(FDR)>30 & logFC.Sub.2>1)),
                 "Sub.1"=rownames(filter(as.data.frame(subd.results[['Sub.3']]), -log10(FDR)>30 & logFC.Sub.1>1)),
                 "RHP.L6"=rownames(filter(as.data.frame(subd.results[['Sub.3']]), -log10(FDR)>30 & logFC.RHP.L6>1)))
sapply(sub.3.list, length)
#155   119   198   184 

t4 = table(unlist(sub.3.list))
sub3.unique = names(t4)[t4==length(sub.3.list)]
length(sub3.unique) #68

avg.expr = mutate(as.data.frame(subd.results[['Sub.3']][sub3.unique,]), 
                  avg.sub3.expr = rowMeans(logcounts(sce_subd)[sub3.unique,sce_subd$new.cell.class=="Sub.3"]))
length(rownames(filter(avg.expr, avg.sub3.expr>1))) #21
#[1] "AC007100.1" "AC010967.1" "AC023503.1" "AC046195.2" "AL356295.1"
#[6] "CD36"       "COL4A1"     "COL4A2"     "DISC1"      "FSTL5"     
#[11] "GUCA1C"     "LINC01194"  "LINC01239"  "LINC01821"  "NR2F2-AS1" 
#[16] "PLEKHG1"    "RASGEF1B"   "SCN7A"      "SCUBE1"     "SNCAIP"    
#[21] "VEGFC" 
plot.list$Sub.3 = rownames(filter(avg.expr, avg.sub3.expr>1)) 

#superficial DE
sce_subs = sce_subset[,sce_subset$new.cell.class %in% c("Sub.1","Sub.2","ProS","PreS","ENT.sup1","ENT.sup2a","ENT.sup2b","ENT.sup3")]
sce_subs$new.cell.class = droplevels(sce_subs$new.cell.class)
sce_subs = sce_subs[rowSums(counts(sce_subs)>0)>100,] #22708  9835
subs.results = scran::findMarkers(sce_subs, group=sce_subs$new.cell.class, test="binom")

filter(as.data.frame(subs.results[['PreS']]), -log10(FDR)>30 & logFC.Sub.1>1) %>% nrow() #446
filter(as.data.frame(subs.results[['PreS']]), -log10(FDR)>30 & logFC.Sub.2>1) %>% nrow() #419
filter(as.data.frame(subs.results[['PreS']]), -log10(FDR)>30 & logFC.ProS>1) %>% nrow() #463
filter(as.data.frame(subs.results[['PreS']]), -log10(FDR)>30 & logFC.ENT.sup1>1) %>% nrow() #351
filter(as.data.frame(subs.results[['PreS']]), -log10(FDR)>30 & logFC.ENT.sup2a>1) %>% nrow() #418
filter(as.data.frame(subs.results[['PreS']]), -log10(FDR)>30 & logFC.ENT.sup2b>1) %>% nrow() #321
filter(as.data.frame(subs.results[['PreS']]), -log10(FDR)>30 & logFC.ENT.sup3>1) %>% nrow() #362

saveRDS(subs.results, "snRNAseq_hpc/processed-data/revision/subiculum-DE_pres-superficial.rda")


pres.list <- list("ProS"=rownames(filter(as.data.frame(subs.results[['PreS']]), -log10(FDR)>30 & logFC.ProS>1)),
                 "Sub.2"=rownames(filter(as.data.frame(subs.results[['PreS']]), -log10(FDR)>30 & logFC.Sub.2>1)),
                 "Sub.1"=rownames(filter(as.data.frame(subs.results[['PreS']]), -log10(FDR)>30 & logFC.Sub.1>1)),
                 "ENT.sup1"=rownames(filter(as.data.frame(subs.results[['PreS']]), -log10(FDR)>30 & logFC.ENT.sup1>1)),
                 "ENT.sup2a"=rownames(filter(as.data.frame(subs.results[['PreS']]), -log10(FDR)>30 & logFC.ENT.sup2a>1)),
                 "ENT.sup2b"=rownames(filter(as.data.frame(subs.results[['PreS']]), -log10(FDR)>30 & logFC.ENT.sup2b>1)),
                 "ENT.sup3"=rownames(filter(as.data.frame(subs.results[['PreS']]), -log10(FDR)>30 & logFC.ENT.sup3>1)))

t5 = table(unlist(pres.list))
pres.unique = names(t5)[t5==length(pres.list)]
length(pres.unique) #38

avg.expr = mutate(as.data.frame(subs.results[['PreS']][pres.unique,]), 
                  avg.pres.expr = rowMeans(logcounts(sce_subs)[pres.unique,sce_subs$new.cell.class=="PreS"]))
length(rownames(filter(avg.expr, avg.pres.expr>1))) #5
#[1] "AC008662.1" "AL161629.1" "FSTL5"      "MDFIC"      "WSCD1" 

plot.list$PreS = rownames(filter(avg.expr, avg.pres.expr>1))

saveRDS(plot.list, "snRNAseq_hpc/processed-data/revision/subiculum-DE_sig-genes-list.rda")
