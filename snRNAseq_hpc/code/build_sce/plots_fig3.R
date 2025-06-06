####figure 3 plots
library("scater")
library("scran")
library("SingleCellExperiment")
library("pheatmap")
library("here")
library("ggrastr")
library("dplyr")
library("magrittr")

##load data
load(file=here::here("snRNAseq_hpc",'processed-data','sce','sce_final.rda'))

##load palettes
load(here::here('plots','snRNAseq_palettes.rda'))

################UMAP###############
pdf(file=here::here('plots','figures','figure_3','cluster_UMAP_final.pdf'),h=4,w=4)
plot<-plotUMAP(sce,text_by='fine.cell.class',colour_by='fine.cell.class',
               point_size=0.001,add_legend=F,point_alpha=1,text_size=6)+
    theme(axis.text.x=element_blank(),
                             axis.ticks.x=element_blank(),
                             axis.text.y=element_blank(),
                             axis.ticks.y=element_blank(),
                             text = element_text(colour='black',size=11))+
          scale_colour_manual(values=as.vector(sn.fine.palette))
rasterize(plot,layers='Point',dpi=500)
dev.off()


################Barplots###############
#####fine cell class#####
# Obtain the column metadata from sce
df <- as.data.frame(colData(sce))
palette<-sn.fine.palette


# Create a data frame with counts of each cluster for each sample (brnum)
count_df <- df %>%
  group_by(brnum,fine.cell.class) %>%
  summarise(count = n(), .groups = 'drop') %>%
    #count(fine.cell.class, brnum) %>%
    #mutate(proportion = n / sum(n)) %>%
    ungroup()

# Set the factor levels of sample_id based on the order in neuron_order
neuron_order<-c("Br2720", "Br6522", "Br8492", "Br8325", "Br3942",
                "Br8667", "Br6471", "Br6432", "Br6423", "Br2743")

count_df$brnum <- factor(count_df$brnum, levels = neuron_order)
# Aggregate data
agg_df <- df %>%
  group_by(fine.cell.class) %>%
  summarise(count = n(), .groups = 'drop') %>%
  #count(fine.cell.class, brnum) %>%
  #mutate(proportion = n / sum(n)) %>%
  ungroup()

# Generate plot
pdf(file=here::here('plots','figures','figure_3','barplot_total_finecellclass_final.pdf'),
    h=2.5,w=4)
ggplot(agg_df, aes(x = "Total", y = count, fill = fine.cell.class)) +
    geom_bar(stat = 'identity') +
    labs(x = "", y = 'Overall Proportion', fill = 'fine.cell.class') +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
          text = element_text(colour='black')) +
    scale_fill_manual(values = palette)
dev.off()


pdf(file=here::here('plots','figures','figure_3','barplot_donor_finecellclass_final.pdf'),
    w=2.5,h=4)
ggplot(count_df, aes(x = brnum, y = count, fill = fine.cell.class)) +
    geom_bar(stat = 'identity', position = 'fill') +
    labs(x = 'donor', y = 'proportion', fill = 'cell type') +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
          text = element_text(colour='black'),
          legend.position='none') +
    scale_fill_manual(values=palette)
dev.off()

# Obtain the column metadata from sce
df <- as.data.frame(colData(sce))

# Create a data frame with counts of each cluster for each sample (brnum)
count_df <- df %>%
  group_by(sort,fine.cell.class) %>%
  summarise(count = n(), .groups = 'drop') %>%
  #count(fine.cell.class, brnum) %>%
  #mutate(proportion = n / sum(n)) %>%
  ungroup()


pdf(file=here::here('plots','figures','figure_3','barplot_finecellclass_sort_final.pdf'),
    w=2.5,h=4)
ggplot(count_df, aes(x = sort, y = count, fill = fine.cell.class)) +
    geom_bar(stat = 'identity', position = 'fill') +
    labs(x = 'sort', y = 'proportion', fill = 'cell type') +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
          text = element_text(colour='black'),
          legend.position='none') +
    scale_fill_manual(values=palette)
dev.off()


#####mid cell class#####
palette=sn.mid.palette
# Obtain the column metadata from sce
df <- as.data.frame(colData(sce))

# Create a data frame with counts of each cluster for each sample (brnum)
count_df <- df %>%
  group_by(brnum,mid.cell.class) %>%
  summarise(count = n(), .groups = 'drop') %>%
  #count(mid.cell.class, brnum) %>%
  #mutate(proportion = n / sum(n)) %>%
  ungroup()

# Set the factor levels of sample_id based on the order in neuron_order
neuron_order<-c("Br2720", "Br6522", "Br8492", "Br8325", "Br3942",
                "Br8667", "Br6471", "Br6432", "Br6423", "Br2743")

count_df$brnum <- factor(count_df$brnum, levels = neuron_order)
# Aggregate data
agg_df <- df %>%
  group_by(mid.cell.class) %>%
  summarise(count = n(), .groups = 'drop') %>%
  #count(mid.cell.class, brnum) %>%
  #mutate(proportion = n / sum(n)) %>%
  ungroup()

# Generate plot
pdf(file=here::here('plots','figures','figure_3','barplot_total_midcellclass_final.pdf'),
    h=2.5,w=4)
ggplot(agg_df, aes(x = "Total", y = count, fill = mid.cell.class)) +
  geom_bar(stat = 'identity') +
  labs(x = "", y = 'Overall Proportion', fill = 'mid.cell.class') +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        text = element_text(colour='black')) +
  scale_fill_manual(values = palette)
dev.off()


pdf(file=here::here('plots','figures','figure_3','barplot_donor_midcellclass_final.pdf'),
    w=2.5,h=4)
ggplot(count_df, aes(x = brnum, y = count, fill = mid.cell.class)) +
  geom_bar(stat = 'identity', position = 'fill') +
  labs(x = 'donor', y = 'proportion', fill = 'cell type') +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        text = element_text(colour='black'),
        legend.position='none') +
  scale_fill_manual(values=palette)
dev.off()

# Obtain the column metadata from sce
df <- as.data.frame(colData(sce))

# Create a data frame with counts of each cluster for each sample (brnum)
count_df <- df %>%
  group_by(sort,mid.cell.class) %>%
  summarise(count = n(), .groups = 'drop') %>%
  #count(mid.cell.class, brnum) %>%
  #mutate(proportion = n / sum(n)) %>%
  ungroup()




pdf(file=here::here('plots','figures','figure_3','barplot_midcellclass_sort_final.pdf'),
    w=2.5,h=4)
ggplot(count_df, aes(x = sort, y = count, fill = mid.cell.class)) +
  geom_bar(stat = 'identity', position = 'fill') +
  labs(x = 'sort', y = 'proportion', fill = 'cell type') +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        text = element_text(colour='black'),
        legend.position='none') +
  scale_fill_manual(values=palette)
dev.off()

#####broad cell class#####

palette=sn.broad.palette
# Obtain the column metadata from sce
df <- as.data.frame(colData(sce))

# Create a data frame with counts of each cluster for each sample (brnum)
count_df <- df %>%
  group_by(brnum,broad.cell.class) %>%
  summarise(count = n(), .groups = 'drop') %>%
  #count(broad.cell.class, brnum) %>%
  #mutate(proportion = n / sum(n)) %>%
  ungroup()

# Set the factor levels of sample_id based on the order in neuron_order
neuron_order<-c("Br2720", "Br6522", "Br8492", "Br8325", "Br3942",
                "Br8667", "Br6471", "Br6432", "Br6423", "Br2743")

count_df$brnum <- factor(count_df$brnum, levels = neuron_order)
# Aggregate data
agg_df <- df %>%
  group_by(broad.cell.class) %>%
  summarise(count = n(), .groups = 'drop') %>%
  #count(broad.cell.class, brnum) %>%
  #mutate(proportion = n / sum(n)) %>%
  ungroup()

# Generate plot
pdf(file=here::here('plots','figures','figure_3','barplot_total_broadcellclass_final.pdf'),
    h=2.5,w=4)
ggplot(agg_df, aes(x = "Total", y = count, fill = broad.cell.class)) +
  geom_bar(stat = 'identity') +
  labs(x = "", y = 'Overall Proportion', fill = 'broad.cell.class') +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        text = element_text(colour='black')) +
  scale_fill_manual(values = palette)
dev.off()


pdf(file=here::here('plots','figures','figure_3','barplot_donor_broadcellclass_final.pdf'),
    w=2.5,h=4)
ggplot(count_df, aes(x = brnum, y = count, fill = broad.cell.class)) +
  geom_bar(stat = 'identity', position = 'fill') +
  labs(x = 'donor', y = 'proportion', fill = 'cell type') +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        text = element_text(colour='black'),
        legend.position='none') +
  scale_fill_manual(values=palette)
dev.off()

# Obtain the column metadata from sce
df <- as.data.frame(colData(sce))

# Create a data frame with counts of each cluster for each sample (brnum)
count_df <- df %>%
  group_by(sort,broad.cell.class) %>%
  summarise(count = n(), .groups = 'drop') %>%
  #count(broad.cell.class, brnum) %>%
  #mutate(proportion = n / sum(n)) %>%
  ungroup()




pdf(file=here::here('plots','figures','figure_3','barplot_broadcellclass_sort_final.pdf'),
    w=2.5,h=4)
ggplot(count_df, aes(x = sort, y = count, fill = broad.cell.class)) +
  geom_bar(stat = 'identity', position = 'fill') +
  labs(x = 'sort', y = 'proportion', fill = 'cell type') +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        text = element_text(colour='black'),
        legend.position='none') +
  scale_fill_manual(values=palette)
dev.off()


################Violins###############
pdf(file=here::here('plots','figures','figure_3','marker_violins_final.pdf'),w=6,h=8)
plot<-plotExpression(sce,features=c('SV2B','GAD2','C3','ETNPPL',
                                    'MOBP','PDGFRA',
                                    'DNAH11','EBF1'),
                     x="fine.cell.class", colour_by="mid.cell.class", point_alpha=1, point_size=.2,add_legend=F,ncol=1)+
    stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar",
                 width = 0.25)+
    theme(axis.text.x = element_text(angle = 90))+
         # text=element_text(size = 12))+
    labs(x='fine.cell.class')+scale_color_manual(values=sn.mid.palette)
rasterize(plot,layers = c('Point'),dpi=350)
dev.off()

pdf(file=here::here('plots','figures','figure_3','EBF1_axistext.pdf'),w=6,h=1.5)
ggcells(sce,mapping=aes(x=cell.class, y=EBF1,fill=broad.class)) +
    geom_boxplot()+theme(axis.text.x = element_text(angle = 90,size=12),
                         text=element_text(colour='black'))+
    theme(legend.position='none',
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))+scale_fill_manual(values=palette2)
dev.off()

############Heatmap#############
sce<-sce[,!sce$fine.cell.class %in% c('Amy','Thal','GABA.PENK')]
sce$fine.cell.class<-droplevels(sce$fine.cell.class)
sce$superfine.cell.class<-droplevels(sce$superfine.cell.class)
palette<-sn.fine.palette
palette<-palette[match(levels(sce$fine.cell.class),names(palette))]
palette2<-sn.mid.palette
palette2<-palette2[match(levels(sce$mid.cell.class),names(palette2))]
palette3<-sn.broad.palette
palette3<-palette3[match(levels(sce$broad.cell.class),names(palette3))]


sce$superfine.cell.class<-factor(sce$superfine.cell.class,
                      levels=c("GC.1", "GC.2", "GC.3", "GC.4", "GC.5", "MC",
                               "CA3.1", "CA3.2", "CA2", "CA1", "ProS", "Sub.1",
                               "Sub.2",
                               "L6.2", "L6.1", "L6b",
                               "L5.1", "L5.2",
                               "L2/3.1","L2/3.5", "L2/3.2", "L2/3.3", "L2/3.4",
                               "L2/3.6",
                               "HATA", "Cajal", "PENK", "SST", "CORT", "PV.FS",
                               "CRABP1", "C1QL1", "LAMP5.MGE", "LAMP5.CGE",
                               "CXCL14", "HTR3A", "VIP",
                               "Micro.1", "Micro.2", "Macro/Tcell",
                               "Astro.1", "Astro.2", "Astro.3",
                               "Oligo.1", "Oligo.2",
                               "OPC", "COP",
                               "Ependy",
                               "CP.1", "CP.2", "CP.3",
                               "Endo.1", "Endo.2",
                               "PC/SMC", "VLMC"))


###set up annotations
anno<-data.frame('superfine'=sce$superfine.cell.class,
                 'fine'=sce$fine.cell.class,
                 'mid'=sce$mid.cell.class,
                 'broad'=sce$broad.cell.class)
anno<-anno[!duplicated(anno$superfine),]
rownames(anno)=anno$superfine
anno$superfine=NULL
anno<-anno[match(levels(sce$superfine.cell.class),rownames(anno)),]

palettes<-list(palette,palette2,palette3)
names(palettes)<-c('fine','mid','broad')

###features to plot:
features<-c('PROX1','TSPAN18','CARTPT','ZNF208','ST18','SLC9A2','COL5A2',
            'FNDC1','CHRM5','GFRA1','CYP26B1',
            'TLE4','SATB2','CUX2','ESR1','LHX1','PENK','LHX6','SST','BTBD11',
            'LAMP5','KIT','CXCL14','HTR3A','VIP',
            'C3','SKAP1',
            'ETNPPL','MOBP','PDGFRA','GPR17','CFAP73','PRLR',
            'VWF','ACTA2','CEMIP')

##make heatmap!!!
pdf(file=here::here('plots','figures','figure_3','heatmap_grayscale.pdf'),h=5.6,w=5.6)
plotGroupedHeatmap(sce,features=features,group='superfine.cell.class',
                   color=colorRampPalette(c('white','black'))(25),
                   cluster_rows=F,cluster_cols=F,zlim=c(0,5),annotation_col=anno,
                   annotation_colors=palettes,annotation_legend=F,legend=F,fontsize=8)
dev.off()
