library(SpatialExperiment)
library(scater)
library(dplyr)
library(ggplot2)

load("processed-data/revision/spe_broad-domain-JT.Rdata")
dim(spe_pseudo)
#[1] 15404   317

#reproduce boxplots from Fig2 and see how they'd look as violins
#from this code: https://github.com/LieberInstitute/spatial_hpc/blob/main/code/08_pseudobulk/PRECAST/volcanoes_boxplots_heatmap_fig2.R

#had to search for the broad palette and rename the colors
load("plots/spatial_palettes.rda")
names(broad.srt.palette) <- c("Neuron","Neuropil","WM","Vasc_CSF")
"#008000" "#a1a1a1" "#ff3ffc" "#00006a" 


rownames(spe_pseudo)[rowData(spe_pseudo)$gene_name=="CLSTN3"] #"ENSG00000139182"
#confirm that the boxplot i'm getting looks like his since i had to recreate the pseudobulk object
ggcells(spe_pseudo, mapping=aes(x=domain, y=ENSG00000139182,fill=broad.domain)) +
  geom_boxplot()+theme(axis.text.x = element_text(angle = 90),
                       text=element_text(size = 13,colour='black'))+
  theme(legend.position='none',
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))+
  scale_fill_manual(values=broad.srt.palette)


p1 <- ggcells(spe_pseudo, mapping=aes(x=domain, y=ENSG00000139182,fill=broad.domain)) +
  #ggbeeswarm::geom_quasirandom(width=.25)+
  geom_violin(trim=F, scale="width", width=.5)+ #super flat tops with normal violin
  labs(y="CLSTN3")+
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=.5),
        text=element_text(size = 13,colour='black'))+
  theme(legend.position='none',
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))+
  scale_fill_manual(values=broad.srt.palette)

rownames(spe_pseudo)[rowData(spe_pseudo)$gene_name=="SLC1A3"] #"ENSG00000079215"
p2 <- ggcells(spe_pseudo, mapping=aes(x=domain, y=ENSG00000079215,fill=broad.domain)) +
  #ggbeeswarm::geom_quasirandom(width=.25)+
  geom_violin(trim=F, scale="width", width=.5)+ #super flat tops with normal violin
  labs(y="SLC1A3")+
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=.5),
        text=element_text(size = 13,colour='black'))+
  theme(legend.position='none',
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))+
  scale_fill_manual(values=broad.srt.palette)

rownames(spe_pseudo)[rowData(spe_pseudo)$gene_name=="SHTN1"] #"ENSG00000187164"
p3 <- ggcells(spe_pseudo, mapping=aes(x=domain, y=ENSG00000187164,fill=broad.domain)) +
  #ggbeeswarm::geom_quasirandom(width=.25)+
  geom_violin(trim=F, scale="width", width=.5)+ #super flat tops with normal violin
  labs(y="SHTN1")+
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=.5),
        text=element_text(size = 13,colour='black'))+
  theme(legend.position='none',
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))+
  scale_fill_manual(values=broad.srt.palette)

rownames(spe_pseudo)[rowData(spe_pseudo)$gene_name=="TPM2"] #"ENSG00000198467"
p4 <- ggcells(spe_pseudo, mapping=aes(x=domain, y=ENSG00000198467,fill=broad.domain)) +
  #ggbeeswarm::geom_quasirandom(width=.25)+
  geom_violin(trim=F, scale="width", width=.5)+ #super flat tops with normal violin
  labs(y="TPM2")+scale_y_continuous(breaks=c(0,3,6,9))+
  theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=.5),
        text=element_text(size = 13,colour='black'))+
  theme(legend.position='none',
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"))+
  scale_fill_manual(values=broad.srt.palette)

ggsave("plots/revision/Fig2_broad-domain-violins.pdf", 
       gridExtra::grid.arrange(p1, p2, p3, p4, ncol=4),
       width=12, height=3.5)


