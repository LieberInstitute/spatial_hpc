library(SpatialExperiment)
library(dplyr)
library(ggplot2)
set.seed(123)

ldsc.new <- read.csv("snRNAseq_hpc/processed-data/revision/ldsc_results_toprank.csv", row.names=1)
nrow(filter(ldsc.new, FDR<.05)) #538

###
load(file=here::here('processed-data','NMF','spe_nmf_final.rda'))

non0.spots = colSums(as.matrix(colData(spe)[,paste0("nmf",1:100)])>0)
remove.nmf = names(non0.spots[non0.spots<1050])
remove.nmf[c(2,3,6)] <- c("nmf2","nmf3","nmf16")
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

#ldsc sig results in keep
sig.ldsc = filter(ldsc.new, cell %in% nmf.ordered.keep, FDR<.05) 
#433 significant results

unique(sig.ldsc$trait)
sig.ldsc$trait.ordered = factor(sig.ldsc$trait, levels=c("BMI","Height","Education Years","Drinks per week",
                                                         "Smoking initiation","Smoking cessation","Age of smoking","Cigarettes per day",
                                                         "Intelligence","Neuroticism",
                                                         "Epilepsy","Schizophrenia",
                                                         "Bipolar","Depression","Autism","ADHD","Anorexia",
                                                         "Alzheimer Disease","Parkinson Disease"))
sig.ldsc$nmf_f= factor(sig.ldsc$cell, levels=nmf.ordered.keep)

palette1= RColorBrewer::brewer.pal(n=7,"RdYlBu")[7:1]
ggplot(sig.ldsc, aes(x=nmf_f, y=trait.ordered, size=-log10(FDR), color=Coefficient_z.score))+
  geom_count()+theme_bw()+
  scale_size(range=c(1,4))+
  #scale_color_gradient2(low=palette1[1], mid="lightgrey", high=palette1[7])+
  scale_color_gradientn(colors=c(palette1[1:2],"lightgrey",palette1[6:7]),
                        values=c(0,.2,.4,.8,1))+
  theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5))

