setwd('/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/')
suppressPackageStartupMessages(library("scater"))
suppressPackageStartupMessages(library("scran"))
suppressPackageStartupMessages(library("SingleCellExperiment"))
suppressPackageStartupMessages(library("gridExtra"))
suppressPackageStartupMessages(library("here"))
suppressPackageStartupMessages(library("ggrastr"))

load(here("snRNAseq_hpc",'processed-data','sce','sce_no_empty_droplets.Rdata'), verbose = TRUE)

# ## low sum/detected genes
pdf(file=here::here('plots','figures','supp_figures','QC_plots_snRNAseq.pdf'),w=11.7,h=15)
p1<-plotColData(sce, x = "Sample", y = "sum", colour_by = "low_library",point_alpha=1,point_size=1) +
  scale_y_log10() + ylab('library size')+xlab('sample')+facet_wrap(~sce$sort,scales = 'free_x')+
  ggtitle("Total UMIs") +theme(axis.text.x = element_text(size=8, angle = 90))+theme(legend.position = 'bottom')
p2<-plotColData(sce, x = "Sample", y = "sum", colour_by = "discard_semiauto",point_alpha=1,point_size=1) +
  scale_y_log10() + ylab('library size')+xlab('sample')+facet_wrap(~sce$sort,scales = 'free_x')+
  ggtitle("Total UMIs") +theme(axis.text.x = element_text(size=8, angle = 90))+theme(legend.position = 'bottom')
p3<-plotColData(sce, x = "Sample", y = "detected", colour_by = "low_genes",point_alpha=1,point_size=1) +
    scale_y_log10() + ylab('detected genes')+xlab('sample')+facet_wrap(~sce$sort,scales = 'free_x')+
    ggtitle("Detected genes") +theme(axis.text.x = element_text(size=8, angle = 90))+theme(legend.position = 'bottom')
p4<-plotColData(sce, x = "Sample", y = "detected", colour_by = "discard_semiauto",point_alpha=1,point_size=1) +
    scale_y_log10() + ylab('detected genes')+xlab('sample')+facet_wrap(~sce$sort,scales = 'free_x')+
    ggtitle("Detected genes") +theme(axis.text.x = element_text(size=8, angle = 90))+theme(legend.position = 'bottom')
p5<-plotColData(sce, x = "Sample", y = "subsets_Mito_percent", colour_by = "high_mito",point_alpha=1,point_size=1) +
    ylab('mito expression percent')+xlab('sample')+facet_wrap(~sce$sort,scales = 'free_x')+
    ggtitle("Mitochondrial gene expression percentage") +theme(axis.text.x = element_text(size=8, angle = 90))+theme(legend.position = 'bottom')
p6<-plotColData(sce, x = "Sample", y = "subsets_Mito_percent", colour_by = "discard_semiauto",point_alpha=1,point_size=1) +
   ylab('mito expression percent')+xlab('sample')+facet_wrap(~sce$sort,scales = 'free_x')+
    ggtitle("Mitochondrial gene expression percentage") +theme(axis.text.x = element_text(size=8, angle = 90))+theme(legend.position = 'bottom')
grid.arrange(p1, p2,p3,p4,p5,p6,ncol=2)
dev.off()

##barplot showing  nuclei by sample, colored by discard status
# Convert SingleCellExperiment to dataframe
df <- as.data.frame(colData(sce))

# Calculate proportions
df <- df %>%
    group_by(discard_semiauto,Sample, sort) %>%
    summarise(count = n(), .groups = 'drop') %>%
    mutate(proportion = count / sum(count)) %>%
    ungroup() %>%
    group_by(Sample) %>%
    mutate(total = sum(proportion)) %>%
    mutate(proportion = proportion / total)

# Create barplot
pdf(file=here::here('plots','figures','supp_figures','barplot_preDiscard.pdf'),h=4,w=8)
ggplot(df, aes(x = Sample, y = count, fill = discard_semiauto)) +
    geom_bar(stat = "identity", position = "stack") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          text=element_text(colour='black',size=16),
          axis.text.x = element_text(angle = 90),
          legend.position='bottom') +
    ylab("nuclei") +
    xlab("sample") +
    labs(fill = "discard")+facet_wrap(~sort,scales='free_x')+
    geom_hline(yintercept=6000,linetype=2,color='magenta')+
    geom_hline(yintercept=9000,linetype=2,color='black')
dev.off()
