library(dplyr)
library(ggplot2)

results <- read.csv("code/enrichment_analysis/project_all/ldsc_results.csv", row.names=1)
colnames(results)
sn.res <- results[grep("t_stat", results$cell),]
nmf.res <- results[grep("nmf", results$cell),]
srt.res <- anti_join(results, bind_rows(sn.res, nmf.res))
unique(srt.res$cell)

sn.res$cluster = substr(sn.res$cell, start=8, stop=30)
sort(unique(sn.res$cluster))

sn.res$cluster_f = factor(sn.res$cluster,
                                  levels=c("GC.3","GC.1","GC.2","GC.4","GC.5", #checked
                                           "MC","CA3.1","CA3.2","CA2","CA1","ProS","Sub.1","Sub.2",#checked
                                           "L6.1","L6.2","L6b","L5.2","L5.1",#checked
                                           "L2.3.3","L2.3.6","L2.3.4","L2.3.2",#checked
                                           "L2.3.1","L2.3.5",#checked
                                           "HATA",#"AHi.1","AHi.2","AHi.3","AHi.4","Thal", #checked
                                           #"Cajal",
                                           "CXCL14","HTR3A","VIP", #checked
                                           "LAMP5.CGE","LAMP5.MGE","CRABP1", #checked
                                           "C1QL1","PV.FS","SST","CORT",#"PENK", #checked
                                           "Astro.1","Astro.2","Astro.3", #checked
                                           "Oligo.1","Oligo.2","OPC","COP",#checked
                                           "Micro.1","Micro.2","Macro.Tcell",#checked
                                           "Ependy","CP.1","CP.2","CP.3","Endo.2","Endo.1","PC.SMC","VLMC"), #checked
                                  labels=c("GC.3","GC.1","GC.2","GC.4","GC.5", #checked
                                           "MC","CA3.1","CA3.2","CA2","CA1","ProS","Sub.1","Sub.2",#checked
                                           "L6.1","L6.2","L6b","L5.2","L5.1",#checked
                                           "L2/3.3","L2/3.6","L2/3.4","L2/3.2",#checked
                                           "L2/3.1","L2/3.5",#checked
                                           "HATA",#"AHi.1","AHi.2","AHi.3","AHi.4","Thal", #checked
                                           #"Cajal",
                                           "CXCL14","HTR3A","VIP", #checked
                                           "LAMP5.CGE","LAMP5.MGE","CRABP1", #checked
                                           "C1QL1","PV.FS","SST","CORT",#"PENK", #checked
                                           "Astro.1","Astro.2","Astro.3", #checked
                                           "Oligo.1","Oligo.2","OPC","COP",#checked
                                           "Micro.1","Micro.2","Macro/Tcell",#checked
                                           "Ependy","CP.1","CP.2","CP.3","Endo.2","Endo.1","PC/SMC","VLMC"))
sn.res$y_order = factor(as.character(sn.res$cluster_f), levels=rev(levels(sn.res$cluster_f)))
sn.res$x_order = factor(sn.res$trait, levels=c("Type 2 Diabetes","BMI","Height","Education Years","Drinks per week",
                                               "Smoking initiation","Smoking cessation","Age of smoking","Cigarettes per day",
                                               "Intelligence","Neuroticism",
                                               "Epilepsy","Schizophrenia",
                                               "Bipolar","Depression","Autism","ADHD","Anorexia",
                                               "Alzheimer Disease","Parkinson Disease"))


p1 <- ggplot(sn.res, aes(y=y_order, x=x_order, color=Coefficient))+
  geom_point(data=filter(sn.res, FDR>.05), size=1, color="grey")+
  geom_point(data=filter(sn.res, FDR<.05), size=3.2, color="black")+
  geom_point(data=filter(sn.res, FDR<.05), size=3)+
  scale_color_gradientn(colors=RColorBrewer::brewer.pal(n=7,"RdYlBu")[7:1], limits=c(-4e-8, 6e-8),
                        values=c(0,.15,.35,.4,.45,.7,1))+
  labs(title="snRNA-seq S-LDSC results", x="trait", y="snRNA-seq cluster")+
  theme_minimal()+theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5),
                        axis.text.y=element_text(margin=margin(0,0,0,0,"pt")),
			legend.key.size=unit(10,"pt"), legend.text=element_text(size=7),
                        plot.title.position = "plot",
                        legend.box.spacing=unit(0,"pt"), legend.box.margin=margin(0,0,0,0,"pt"))


load("plots/spatial_palette_final.rda")
srt.res$domain_f = factor(gsub("CA1_srt","CA1",srt.res$cell), levels=names(spatial.palette))
srt.res$y_order = factor(gsub("CA1_srt","CA1",srt.res$cell), levels=rev(names(spatial.palette)))
srt.res$x_order = factor(srt.res$trait, levels=c("Type 2 Diabetes","BMI","Height","Education Years","Drinks per week",
                                               "Smoking initiation","Smoking cessation","Age of smoking","Cigarettes per day",
                                               "Intelligence","Neuroticism",
                                               "Epilepsy","Schizophrenia",
                                               "Bipolar","Depression","Autism","ADHD","Anorexia",
                                               "Alzheimer Disease","Parkinson Disease"))

p2 <- ggplot(srt.res, aes(y=y_order, x=x_order, color=Coefficient))+
  geom_point(data=filter(srt.res, FDR>.05), size=1, color="grey")+
  geom_point(data=filter(srt.res, FDR<.05), size=3.2, color="black")+
  geom_point(data=filter(srt.res, FDR<.05), size=3)+
  scale_color_gradientn(colors=RColorBrewer::brewer.pal(n=7,"RdYlBu")[7:1], limits=c(-3e-8, 6e-8),
                        values=c(0,.15,.3,.35,.4,.7,1), breaks=c(-2e-8,0,2e-8,4e-8,6e-8))+
  labs(title="SRT S-LDSC results", x="trait", y="SRT domain")+
  theme_minimal()+theme(axis.text.x=element_text(angle=90, hjust=1, vjust=.5),
                        legend.key.size=unit(10,"pt"), legend.text=element_text(size=7),
                        plot.title.position = "plot",
                        legend.box.spacing=unit(0,"pt"), legend.box.margin=margin(0,0,0,0,"pt"))

pdf(file="plots/revision/supp_ldsc-results_dotplot.pdf", height=10, width=8)
gridExtra::grid.arrange(p1, p2, layout_matrix=cbind(rep(1,9), c(rep(2,4),rep(NA,5))))
dev.off()

