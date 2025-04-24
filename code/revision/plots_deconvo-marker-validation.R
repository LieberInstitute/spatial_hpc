library(SpatialExperiment)
library(dplyr)
library(ggplot2)


load("processed-data/06_clustering/PRECAST/spe_norm_with_domain.rda"))
spg_norm = spe_norm[,spe_norm$brnum %in% c("Br3942_VSPG","Br8325_VSPG")]
dim(spg_norm)
## confirm same dims with spg.rds
#spg = readRDS(here("processed-data/spot_deconvo/shared_utilities/spg.rds"))
#dim(spg) #yes same dims
table(spg_norm$domain)
rm(spe_norm)
colData(spg_norm)[spg_norm$domain=="SLM.WM","domain"] <- "SLM.SGZ"
spg_norm$domain = droplevels(spg_norm$domain)
table(spg_norm$domain)


tbl = read.csv("plots/spot_deconvo/markergenes_stats.csv", row.names=1)
fine.markers = filter(tbl, cellTypeResolution=="Fine")
table(fine.markers$cellType.target) #yes 21 
#table(fine.markers$cellType)
fine.markers$fine_f = factor(as.character(fine.markers$cellType.target), 
                             levels=c("GC","CA2-4","CA1_ProS","Sub.1","Sub.2",
                                      "L6_6b","L5","L2_3.PrS.PaS","L2_3.Prs.Ent",
                                      #"HATA","Amy",
                                      "Thal","Cajal",#"GABA.PENK",
                                      "GABA.MGE","GABA.LAMP5","GABA.CGE",
                                      "Micro_Macro_T","Astro","Oligo","OPC",
                                      "Ependy","Choroid","Vascular"),
                             labels=c("GC","CA2-4","CA1/ProS","Sub.1","Sub.2",
                                      "L6/6b","L5/6","L2/3.1","L2/3.2",
                                      #"HATA","Amy",
                                      "Thal","Cajal",#"GABA.PENK",
                                      "GABA.MGE","GABA.LAMP5","GABA.CGE",
                                      "Micro/Macro/T","Astro","Oligo","OPC",
                                      "Ependy","Choroid","Vascular"))

load(file=here::here('plots','spatial_palette_final.rda'))

df = mutate(as.data.frame(colData(spg_norm)[,c("sample_id","key","domain")]), 
            domain_f=factor(domain, levels=names(spatial.palette)))


unique.df = distinct(fine.markers, cellType.target, fine_f)

for(i in levels(fine.markers$fine_f)) {
  encode = filter(unique.df, fine_f==i)$cellType.target
  genes = filter(fine.markers, fine_f==i)$symbol
  df[[encode]] <- colSums(logcounts(spg_norm)[genes,]>0)
}

df.long = pivot_longer(df, all_of(unique.df$cellType.target), names_to="cellType.target", values_to="n_markers_nonzero") %>%
  mutate(fine_f = factor(cellType.target, levels=c("GC","CA2-4","CA1_ProS","Sub.1","Sub.2",
                                                   "L6_6b","L5","L2_3.PrS.PaS","L2_3.Prs.Ent",
                                                   #"HATA","Amy",
                                                   "Thal","Cajal",#"GABA.PENK",
                                                   "GABA.MGE","GABA.LAMP5","GABA.CGE",
                                                   "Micro_Macro_T","Astro","Oligo","OPC",
                                                   "Ependy","Choroid","Vascular"),
                         labels=c("GC","CA2-4","CA1/ProS","Sub.1","Sub.2",
                                  "L6/6b","L5/6","L2/3.1","L2/3.2",
                                  #"HATA","Amy",
                                  "Thal","Cajal",#"GABA.PENK",
                                  "GABA.MGE","GABA.LAMP5","GABA.CGE",
                                  "Micro/Macro/T","Astro","Oligo","OPC",
                                  "Ependy","Choroid","Vascular")))

dim(df.long)

p1 <- ggplot(df.long, aes(x=domain_f, y=n_markers_nonzero, fill=domain_f))+
  geom_boxplot(outlier.size=.3)+scale_fill_manual(values=spatial.palette)+
  facet_wrap(vars(fine_f), ncol=3)+
  theme_minimal()+theme(legend.position="none", axis.text.x=element_text(angle=90, hjust=1, vjust=.5),
                        panel.grid.minor=element_blank(), panel.grid.major.x=element_blank(),
                        strip.text=element_text(size=12), axis.title.x=element_blank())
ggsave(file="plots/revision/Supp_deconvo_fine-marker-validation.png", p1, bg="white", width=8, height=11)


mid.markers = filter(tbl, cellTypeResolution=="Mid")
table(mid.markers$cellType.target) #yes 21 
mid.markers$mid_f = factor(mid.markers$cellType.target, levels=c("ExcN","InhN","Micro_Macro_T","Astro","Oligo","OPC","CSF","Vascular"),
                           labels=c("ExcN","InhN","Micro/Macro/T","Astro","Oligo","OPC","CSF","Vascular"))

df2 = mutate(as.data.frame(colData(spg_norm)[,c("sample_id","key","domain")]), 
             domain_f=factor(domain, levels=names(spatial.palette)))
unique.df2 = distinct(mid.markers, cellType.target, mid_f)
for(i in levels(mid.markers$mid_f)) {
  encode = filter(unique.df2, mid_f==i)$cellType.target
  genes = filter(mid.markers, mid_f==i)$symbol
  df2[[encode]] <- colSums(logcounts(spg_norm)[genes,]>0)
}

df2.long = pivot_longer(df2, all_of(unique.df2$cellType.target), names_to="cellType.target", values_to="n_markers_nonzero") %>%
  mutate(mid_f = factor(cellType.target, levels=c("ExcN","InhN","Micro_Macro_T","Astro","Oligo","OPC","CSF","Vascular"),
                        labels=c("ExcN","InhN","Micro/Macro/T","Astro","Oligo","OPC","CSF","Vascular")))

p2 <- ggplot(df2.long, aes(x=domain_f, y=n_markers_nonzero, fill=domain_f))+
  geom_boxplot(outlier.size=.3)+scale_fill_manual(values=spatial.palette)+
  facet_wrap(vars(mid_f), ncol=1)+
  theme_minimal()+theme(legend.position="none", axis.text.x=element_text(angle=90, hjust=1, vjust=.5),
                        panel.grid.minor=element_blank(), panel.grid.major.x=element_blank(),
                        strip.text=element_text(size=12), axis.title.x=element_blank())
ggsave(file="plots/revision/Supp_deconvo_mid-marker-validation.png", p2, bg="white", width=3.5, height=11)
