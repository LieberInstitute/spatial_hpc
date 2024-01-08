####figure 3 plots
###Fig 3A--UMAP colored by cluster:
levels(sce$cell.class)
# [1] "Amy"           "Astro"         "CA1/ProS"      "CA2-4"         "Cajal"         "Choroid"       "Ependy"
# [8] "GABA.CGE"      "GABA.LAMP5"    "GABA.MGE"      "GABA.PENK"     "GC"            "HATA"          "L2/3.Prs.Ent"
# [15] "L2/3.PrS.PaS"  "L5"            "L6/6b"         "Micro/Macro/T" "Oligo"         "OPC"           "Sub.1"
# [22] "Sub.2"         "Thal"          "Vascular"
sce$cell.class<-factor(ifelse(sce$fine.type %in% c('L6.1','L6.2','L6b'),'L6/6b',
                              ifelse(sce$fine.type %in% c('L5.1','L5.2'),'L5/6',
                                     ifelse(sce$fine.type %in% 'Sub.1','Sub.1',
                                            ifelse(sce$fine.type %in% 'Sub.2','Sub.2',
                                                   ifelse(sce$fine.type %in% c('CA1','ProS'),'CA1/ProS',
                                                          ifelse(sce$cell.type %in% c('Choroid'),'Choroid',
                                                                 ifelse(sce$cell.type %in% c('Ependy'),'Ependy',
                                                                        ifelse(sce$fine.type %in% c('L2/3.1','L2/3.5'),'L2/3.PrS.PaS',
                                                                               ifelse(sce$fine.type %in% c('L2/3.2','L2/3.4','L2/3.6','L2/3.3'),'L2/3.PrS.Ent',
                                                                                      ifelse(sce$fine.type %in% c('HATA'),'HATA',
                                                                                             ifelse(sce$fine.type %in% c('AHi.1','AHi.2','AHi.3','AHi.4'),'Amy',
                                                                                                    ifelse(sce$fine.type %in% c('PENK'),'GABA.PENK',
                                                                                                           ifelse(sce$fine.type %in% c('LAMP5.CGE','LAMP5.MGE','CXCL14'),'GABA.LAMP5',
                                                                                                                  ifelse(sce$fine.type %in% c('HTR3A','VIP'),'GABA.CGE',
                                                                                                                         ifelse(sce$fine.type %in% c('PV.FS','SST','CORT','CRABP1','C1QL1'),'GABA.MGE',
                                                                                                                                as.character(sce$cell.type)))))))))))))))))
sce$cell.class<-factor(sce$cell.class,levels=c('GC','CA2-4','CA1/ProS','Sub.1','Sub.2',
                                               'L6/6b','L5/6','L2/3.PrS.PaS','L2/3.PrS.Ent',
                                               'HATA','Amy','Thal','Cajal','GABA.PENK','GABA.MGE','GABA.LAMP5','GABA.CGE',
                                               'Oligo','OPC','Micro/Macro/T','Astro','Ependy','Choroid','Vascular'))

sce$broad.class<-factor(
    ifelse(sce$cell.class %in% levels(sce$cell.class)[1:11],'ExcN',
           ifelse(sce$cell.class %in% levels(sce$cell.class)[12:15],'InhN',
                  ifelse(sce$cell.class %in% levels(sce$cell.class)[16],'Oligo',
                         ifelse(sce$cell.class %in% levels(sce$cell.class)[17],'OPC',
                                ifelse(sce$cell.class %in% levels(sce$cell.class)[18],'Micro/Macro/T',
                                       ifelse(sce$cell.class %in% levels(sce$cell.class)[19],'Astro',
                                              ifelse(sce$cell.class %in% levels(sce$cell.class)[20:21],'CSF',
                                                     'Vascular'))))))),levels=c('ExcN','InhN','Oligo','OPC','Micro/Macro/T',
                                                                                'Astro','Vascular','CSF')

)

palette<-c("Choroid"="#000030",
           "Ependy"="#00008a",
          "Vascular"="#0000f4",
           "Micro/Macro/T"='#FAFA33',
               'Oligo'="#a300a0",
           'OPC'='#7070ff',
           'Astro'="#dfa56e",
           'GABA.LAMP5'='#00d4cd',
           'GABA.CGE'='#5fdbff',
           'GABA.MGE'='#33fff9',
           'GABA.PENK'='#00a39e',
           'GC'='#003800',
           'CA2-4'='#008000',
           'CA1/ProS'='#00bb00',
           'Sub.1'='#add294',
           'Sub.2'='#add294',
           'L6/6b'='#61963d',
           'L5/6'='#61963d',
           'L2/3.PrS.PaS'='#2fec00',
           'L2/3.PrS.Ent'='#2fec00',
           'Amy'='#99ff99',
           'HATA'='#99ff99',
           'Thal'='#676666',
           'Cajal'='#333333'
)


palette<-palette[match(levels(sce$cell.class),names(palette))]

palette2<-c('ExcN'='#009000',
            'InhN'='#33fff9',
            'Oligo'='#a300a0',
            'OPC'='#7070ff',
            "Micro/Macro/T"='#FAFA33',
            'Astro'='#dfa56e',
            "Vascular"='#0000f4',
            "CSF"="#000030"
)

palette3<-c('Neuron'='#009000',
          #  'InhN'='#33fff9',
          "Micro"='#FAFA33',
          'Astro'='#dfa56e',
            'Oligo'='#a300a0',
             'Other'='#7070ff'
            #"Vascular"='#0000f4',
            #"CSF"="#000030"
)
palette2<-palette2[match(levels(sce$broad.class),names(palette2))]

spe$domain<-spe$cluster
spe$broad.domain<-spe$broad
levels(spe$domain)[2]<-'CA2.4'
levels(spe$domain)[3]<-'CA2.4'
levels(spe$domain)[3]<-'CA1'
levels(spe$domain)[4]<-'CA1'


spatial.palette <- c(
    "#002800",  # Level 1 Brighter Forest Green
    "#008000",  # Level 5 Thistle (light muted purple)
    '#00dc00',
    '#add294',
    '#61963d',
    '#99ff99',  # Level 4 Purple (darker muted purple)
    '#5ffffb',  # Level 9 Mint Green (mintier)
    "#000000",  # Level 10 Lighter Black
    "#a1a1a1",  # Level 11 Lighter Gray
    "#555555",  # Level 12 Medium Gray (lighter)
    "#dfa56e",  # Level 13 Lighter Deep Navy
    "#ff3ffc",  # Level 14 Sky Blue
    "#7a007a",  # Level 15 Steel Blue
    "#ff80fe",  # Level 16 Lighter Dark Blue
    "#1e1eff",  # Level 17 Sandy Brown
    "#00006a"   # Level 18 Saddle Brown
)
names(spatial.palette)<-levels(spe$domain)




spatial.palette2<-c(
        "#008000",  # Level 1 Brighter Forest Green
        "#a1a1a1",  # Level 2 Bright Hot Pink
        "#ff3ffc",  # Level 3 Medium Pink
        "#00006a"  # Level 4 Purple (darker muted purple)
)
names(spatial.palette2)<-levels(spe$broad.domain)






pdf(file=here::here('plots','figures','figure_3','cluster_UMAP.pdf'),h=3.5,w=3.75)
plot<-plotUMAP(sce,text_by='cell.class',colour_by='cell.class',
               point_size=0.001,add_legend=T,point_alpha=1)+
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          text = element_text(colour='black',size=11))+
    scale_colour_manual(values=as.vector(palette))
rasterize(plot,layers='Point',dpi=1000)
dev.off()

# Obtain the column metadata from sce
df <- as.data.frame(colData(sce))

# Create a data frame with counts of each cluster for each sample (brnum)
count_df <- df %>%
    count(cell.class, brnum) %>%
    group_by(brnum) %>%
    mutate(proportion = n / sum(n)) %>%
    ungroup()
# Set the factor levels of sample_id based on the order in neuron_order
neuron_order<-c("Br2720", "Br6522", "Br8492", "Br8325", "Br3942",
                "Br8667", "Br6471", "Br6432", "Br6423", "Br2743")

count_df$brnum <- factor(count_df$brnum, levels = neuron_order)
# Aggregate data by broad
agg_df <- count_df %>%
    group_by(cell.class) %>%
    summarise(n = sum(n), proportion = sum(n) / sum(count_df$n))

# Generate plot
pdf(file=here::here('plots','figures','figure_3','barplot_proportion_total_cellclass.pdf'),
    h=2.5,w=4)
ggplot(agg_df, aes(x = "Total", y = proportion, fill = cell.class)) +
    geom_bar(stat = 'identity') +
    labs(x = "", y = 'Overall Proportion', fill = 'broad') +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
          text = element_text(colour='black')) +
    scale_fill_manual(values = palette)
dev.off()

pdf(file=here::here('plots','figures','figure_3','barplot_proportion_donor_cellclass.pdf'),
    w=2.5,h=4)
ggplot(count_df, aes(x = brnum, y = proportion, fill = cell.class)) +
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
    count(cell.class, sort) %>%
    group_by(sort) %>%
    mutate(proportion = n / sum(n)) %>%
    ungroup()

# Aggregate data by broad
agg_df <- count_df %>%
    group_by(cell.class) %>%
    summarise(n = sum(n), proportion = sum(n) / sum(count_df$n))


pdf(file=here::here('plots','figures','figure_3','barplot_proportion_cellclass_sort.pdf'),
    w=2.5,h=4)
ggplot(count_df, aes(x = sort, y = proportion, fill = cell.class)) +
    geom_bar(stat = 'identity', position = 'fill') +
    labs(x = 'sort', y = 'proportion', fill = 'cell type') +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
          text = element_text(colour='black'),
          legend.position='none') +
    scale_fill_manual(values=palette)
dev.off()



# Create a data frame with counts of each cluster for each sample (sort)
broad_df <- df %>%
    count(broad.class, sort, brnum) %>%
    group_by(sort) %>%
    mutate(proportion = n / sum(n)) %>%
    ungroup()

# Set the factor levels of sample_id based on the order in neuron_order
neuron_order<-c("Br2720", "Br6522", "Br8492", "Br8325", "Br3942",
                "Br8667", "Br6471", "Br6432", "Br6423", "Br2743")

broad_df$brnum <- factor(broad_df$brnum, levels = neuron_order)
# Aggregate data by broad
agg_df <- broad_df %>%
    group_by(broad.class) %>%
    summarise(n = sum(n), proportion = sum(n) / sum(broad_df$n))

# Generate plot
pdf(file=here::here('plots','figures','figure_3','barplot_proportion_total_broadclass.pdf'),
    h=2.5,w=4)
ggplot(agg_df, aes(x = "Total", y = proportion, fill = broad.class)) +
    geom_bar(stat = 'identity') +
    labs(x = "", y = 'Overall Proportion', fill = 'broad') +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
          text = element_text(colour='black')) +
    scale_fill_manual(values = palette2)
dev.off()

pdf(file=here::here('plots','figures','figure_3','barplot_proportion_donor_broadclass.pdf'),
    w=2.5,h=4)
ggplot(broad_df, aes(x = brnum, y = proportion, fill = broad.class)) +
    geom_bar(stat = 'identity', position = 'fill') +
    labs(x = 'donor', y = 'proportion', fill = 'cell type') +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
          text = element_text(colour='black'),
          legend.position='none') +
    scale_fill_manual(values=palette2)
dev.off()

# Create a data frame with counts of each cluster for each sample (brnum)
sort_df <- df %>%
    count(broad.class, sort) %>%
    group_by(sort) %>%
    mutate(proportion = n / sum(n)) %>%
    ungroup()

# Aggregate data by broad
agg_df <- sort_df %>%
    group_by(broad.class) %>%
    summarise(n = sum(n), proportion = sum(n) / sum(sort_df$n))


pdf(file=here::here('plots','figures','figure_3','barplot_proportion_broadclass_sort.pdf'),
    w=2.5,h=4)
ggplot(sort_df, aes(x = sort, y = proportion, fill = broad.class)) +
    geom_bar(stat = 'identity', position = 'fill') +
    labs(x = 'sort', y = 'proportion', fill = 'cell type') +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
          text = element_text(colour='black'),
          legend.position='none') +
    scale_fill_manual(values=palette2)
dev.off()


pdf(file=here::here('plots','figures','figure_3','marker_violins.pdf'),w=6,h=8)
plot<-plotExpression(sce,features=c('CLSTN3','SV2B','GAD2',
                                    'MOBP','PDGFRA','CSF1R','ETNPPL',
                                    'DNAH11','EBF1'),
                     x="cell.class", colour_by="broad.class", point_alpha=1, point_size=.2,add_legend=F,ncol=1)+
    stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar",
                 width = 0.25)+
    theme(axis.text.x = element_text(angle = 90))+
    # text=element_text(size = 12))+
    labs(x='cell class')+scale_color_manual(values=palette2)
rasterize(plot,layers = c('Point'),dpi=100)
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

###main fig1 heatmap
sce<-sce[,!sce$cell.class %in% c('Amy','Thal')]
sce$cell.class<-droplevels(sce$cell.class)


palette<-c("Choroid"="#4A0404",
           "Vascular"='#A0522D',
           "Ependy"='#D2B48C',
           "Micro/Macro/T"='#71797E',
           'Oligo'='#0000FF',
           'OPC'='#28282B',
           'Astro'='#D3D3D3',
           'GABA.LAMP5'='#40B5AD',
           'GABA.CGE'='#7FFFD4',
           'GABA.MGE'='#AFE1AF',
           'GABA.PENK'='#DFFF00',
           'GC'='#4CBB17',
           'CA2-4'='#F8C8DC',
           'CA1/ProS'='#DA70D6',
           'Sub.1'='#FFFF8F',
           'Sub.2'='#FFD700',
           'L6/6b'='#FFAA33',
           'L5/6'='#FF7518',
           'L2/3.PrS.PaS'='#FA5F55',
           'L2/3.PrS.Ent'='#800080',
           #'Amy'='#953553',
           'HATA'='#F33A6A',
           #'Thal'='#355E3B',
           'Cajal'='#4682B4'
)
palette<-palette[match(levels(sce$cell.class),names(palette))]
sce$fine.type<-droplevels(sce$fine.type)

sce$fine.type<-factor(sce$fine.type,
                      levels=c("GC.1", "GC.2", "GC.3", "GC.4", "GC.5", "MC",
                               "CA3.1", "CA3.2", "CA2", "CA1", "ProS", "Sub.1",
                               "Sub.2",
                               "L6.2", "L6.1", "L6b",
                               "L5.1", "L5.2",
                               "L2/3.1","L2/3.5", "L2/3.2", "L2/3.3", "L2/3.4",
                               "L2/3.6",
                               "HATA", "Cajal", "PENK", "SST", "CORT", "PV.FS",
                               "CRABP1", "C1QL1", "LAMP5.MGE", "LAMP5.CGE",
                               "CXCL14", "HTR3A", "VIP", "Oligo.1", "Oligo.2",
                               "OPC", "COP",
                               "Micro.1", "Micro.2", "Macro/Tcell",
                               "Astro.1", "Astro.2", "Astro.3",
                               "Ependy",
                               "CP.1", "CP.2", "CP.3",
                               "Endo.1", "Endo.2",
                               "PC/SMC", "VLMC"))


###set up annotations
anno<-data.frame('fine'=sce$fine.type,'cell.class'=sce$cell.class,'broad.class'=sce$broad.class)
anno<-anno[!duplicated(anno$fine),]
rownames(anno)=anno$fine
anno$fine=NULL
anno<-anno[match(levels(sce$fine.type),rownames(anno)),]

###features to plot:
features<-c('PROX1','TSPAN18','CARTPT','ZNF208','ST18','SLC9A2','COL5A2',
            'FNDC1','CHRM5','GFRA1','CYP26B1',
            'TLE4','SATB2','CUX2','ESR1','LHX1','PENK','LHX6','SST','BTBD11',
            'LAMP5','KIT','CXCL14','HTR3A','VIP',
            'MOBP','VCAN','GPR17','C3','SKAP1',
            'ETNPPL','CFAP73','PRLR',
            'VWF','ACTA2','CEMIP')

##make heatmap!!!
pdf(file=here::here('plots','figures','figure_3','heatmap_grayscale.pdf'),h=10,w=13)
plotGroupedHeatmap(sce,features=features,group='fine.type',
                   color=colorRampPalette(c('white','black'))(25),
                   cluster_rows=F,cluster_cols=F,zlim=c(0,5),annotation_col=anno,annotation_colors=palettes)
dev.off()
