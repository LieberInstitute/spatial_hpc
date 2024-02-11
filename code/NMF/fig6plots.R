#######fig 6 plots#####
brains <- unique(spe$brnum)


speb <- spe[, (colData(spe)$brnum == brains[9])]
speb$sample_id <- droplevels(speb$sample_id)
speb$sample_id <- as.character(speb$sample_id)
samples <- unique(speb$sample_id)
speb$sample_id <- factor(speb$sample_id, levels = samples)
samples
speb$brnum <- droplevels(speb$brnum)

pdf(here::here('plots','figures','figure_6','spotplot1_hpc_pats.pdf'),h=8,w=8)
plotVisiumRGB(speb,vars=c('nmf11','nmf15','nmf52','nmf61'),image=F,highlight='neuron_cell_body',values=c('gray50','black'))
dev.off()

pdf(here::here('plots','figures','figure_6','spotplot2_sub_pats.pdf'),h=8,w=8)
plotVisiumRGB(speb,vars=c('dnmf54','dnmf17','dnmf32','dnmf40'),image=F,highlight='neuron_cell_body',values=values)
dev.off()

data<-as.data.frame(colData(sce)[,colData(sce)$])

sce_pyr<-sce[,sce$superfine.cell.class %in% c('MC','CA3.1','CA3.2','CA2',
                                              'CA1','ProS','Sub.1','Sub.2',
                                              'L6.2','L6.1','L6b','L5.1','L5.2',
                                              'L2/3.1','L2/3.5','L2/3.2','L2/3.3',
                                              'L2/3.4','L2/3.6','HATA')]

sce_pyr<-sce[,sce$superfine.cell.class %in% c('MC','CA3.1','CA3.2','CA2',
                                        'CA1','ProS','Sub.1','Sub.2',
                                        'L6.2','L6.1','L6b','L5.1','L5.2',
                                        'L2/3.1','L2/3.5','L2/3.2','L2/3.3',
                                        'L2/3.4','L2/3.6','HATA')]

loads<-colData(sce_pyr)[,c('nmf52','nmf11','nmf63','nmf61','nmf15',
                           'nmf32','nmf40','nmf54','nmf22','nmf65',
                           'nmf53','nmf51','nmf68','nmf17','nmf78','nmf27','nmf45',
                           'nmf84','nmf62','nmf69')]
loads<-loads[,colnames(loads) %in% c('nmf52','nmf11','nmf63','nmf61','nmf15',
               'nmf32','nmf40','nmf54','nmf22','nmf65',
               'nmf53','nmf51','nmf68','nmf17','nmf78','nmf27','nmf45',
               'nmf84','nmf62','nmf69')]

loads<-as.data.frame(loads)

data<-as.data.frame(colData(sce_pyr))
heat<-aggregate(data[,c('nmf52','nmf11','nmf63','nmf61','nmf15',
                        'nmf32','nmf40','nmf54','nmf22','nmf65',
                        'nmf53','nmf51','nmf68','nmf17','nmf78','nmf27','nmf45',
                        'nmf84','nmf62','nmf69')],
                            by=list(data$superfine.cell.class),FUN=sum)
rownames(heat)<-heat$Group.1
heat<-heat[,-1]
pdf(here::here('plots','figures','figure_6','heatmap_og_pyr.pdf'),h=4,w=7)
pheatmap(heat,color=colorRampPalette(brewer.pal(n = 7, name =
                                                        "Greys"))(100),cluster_cols=F,cluster_rows=F)
dev.off()

spe_pyr<-spe[,spe$domain %in% levels(spe$domain)[2:6]]
spe_pyr$domain<-droplevels(spe_pyr$domain)

data2<-as.data.frame(colData(spe_pyr))
heat2<-aggregate(data2[,c('dnmf52','dnmf11','dnmf63','dnmf61','dnmf15',
                        'dnmf32','dnmf40','dnmf54','dnmf22','dnmf65',
                        'dnmf53','dnmf51','dnmf68','dnmf17','dnmf78','dnmf27','dnmf45',
                        'dnmf84','dnmf62','dnmf69')],
                by=list(data2$domain),FUN=sum)
rownames(heat2)<-heat2$Group.1
heat2<-heat2[,-1]
colnames(heat2)<-gsub(colnames(heat2),pattern='dnmf',replacement='nmf')

pdf(here::here('plots','figures','figure_6','heatmap_proj_pyr.pdf'),h=3,w=7)
pheatmap(heat2,color=colorRampPalette(brewer.pal(n = 7, name =
                                                    "Greys"))(100),cluster_cols=F,cluster_rows=F)
dev.off()

library(reshape2)
data<-as.data.frame(sce_pyr$superfine.cell.class)
colnames(data)<-'superfine.cell.class'
onehot_superfine.cell.class <-  dcast(data = data, rownames(data) ~ superfine.cell.class, length)
rownames(onehot_superfine.cell.class)<-onehot_superfine.cell.class[,1]
onehot_superfine.cell.class[,1]<-as.numeric(onehot_superfine.cell.class[,1])
onehot_superfine.cell.class<-onehot_superfine.cell.class[order(onehot_superfine.cell.class[,1],decreasing=F),]
onehot_superfine.cell.class[,1]<-NULL

weights<-as.data.frame(colData(sce_pyr)[,c('nmf52','nmf11','nmf63','nmf61','nmf15',
                                           'nmf32','nmf40','nmf54','nmf22','nmf65',
                                           'nmf53','nmf51','nmf68','nmf17','nmf78','nmf27','nmf45',
                                           'nmf84','nmf62','nmf69')])

#####Get all the projections in line!

###1: traced projections from Cembrowski et al:
colData(mouse)[,c(8:17)]<-NULL
mouse$projection<-factor(mouse$note2)
levels(mouse$projection)<-c('Amygdala','Nucleus accumbens','Prefrontal cortex')
##
# Translate from one species to the other using the orthology
names <- orthology[orthology$Column3 %in% rownames(mouse),]

names <- names[match(rownames(mouse), names$Column3),]

setdiff(names$Column3, rownames(mouse))

rownames(mouse) <- names$Column1


set.seed(1029)
i<-intersect(rownames(mouse),rownames(x@w))
loadings<-x@w
loadings<-loadings[rownames(loadings) %in% i,]
mouse2<-mouse[rownames(mouse) %in% i,]
loadings<-loadings[match(rownames(mouse2),rownames(loadings)),]
proj<-RcppML::project(loadings,logcounts(mouse2),L1=0)

col_sums<-colSums(t(proj))

# Rescale each column
proj<-t(proj)
proj2 <- apply(proj,1, function(row) row / col_sums)

# Check if columns now sum to 1
print(rowSums(proj2))

colData(mouse)<-cbind(colData(mouse),t(proj2))

##get aggregate projected weights
data<-as.data.frame(colData(mouse))

heat<-aggregate(data[,c('nmf52','nmf11','nmf63','nmf61','nmf15',
                        'nmf32','nmf40','nmf54','nmf22','nmf65',
                        'nmf53','nmf51','nmf68','nmf17','nmf78','nmf27','nmf45',
                        'nmf84','nmf62','nmf69')],
                by=list(data$projection),FUN=sum)
rownames(heat)<-heat$Group.1
heat<-heat[,-1]
pheatmap(heat,color=colorRampPalette(brewer.pal(n = 7, name =
                                                    "Greys"))(100),cluster_cols=F,cluster_rows=F)


###2: ding et al
sce.ding$subclass<-ifelse(sce.ding$cluster_id %in% c(1:3),'Sub-pyramidal',
                          ifelse(sce.ding$cluster_id %in% c(6:8),'ProS-pyramidal_deep',
                                 ifelse(sce.ding$cluster_id %in% c(10:14),'ProS-pyramidal_sup',
                                        ifelse(sce.ding$cluster_id %in% c(17:19),'HATA-pyramidal',
                                               ifelse(sce.ding$cluster_id %in% c(20:23),'ProS-pyramidal_mostSup',
                                                      ifelse(sce.ding$cluster_id %in% c(24:26),'Sub/ProS-polymorphic',
                                                             ifelse(sce.ding$cluster_id %in% c(28:29),'Sub/ProS/HATA-L6',
                                                                    ifelse(sce.ding$cluster_id %in% c(4:5),'ProS-Sub_border',
                                                                           ifelse(sce.ding$cluster_id %in% c(9,15,16),'CA1-pyramidal',
                                                                                  'PreS')))))))))
sce.ding$subclass<-factor(sce.ding$subclass)
sce.ding$subclass<-factor(sce.ding$subclass,
                          levels=levels(sce.ding$subclass)[c(1,7,6,4,
                                                             8,9,3,10,2,5)])
col_sums<-colSums(
    as.data.frame(colData(sce.ding)[,c(16:115)]))

# Rescale each column
colData(sce.ding)[,c(16:115)] <- t(apply(as.data.frame(colData(sce.ding)[,c(16:115)]),
                                         1, function(row) row / col_sums))

# Check if columns now sum to 1
print(colSums(as.data.frame(colData(sce.ding)[,c(16:115)])))

##get aggregate projected weights
data<-as.data.frame(colData(sce.ding))

heat<-aggregate(data[,c('nmf52','nmf11','nmf63','nmf61','nmf15',
                        'nmf32','nmf40','nmf54','nmf22','nmf65',
                        'nmf53','nmf51','nmf68','nmf17','nmf78','nmf27','nmf45',
                        'nmf84','nmf62','nmf69')],
                by=list(data$subclass),FUN=sum)
rownames(heat)<-heat$Group.1
heat<-heat[,-1]
pheatmap(heat,color=colorRampPalette(brewer.pal(n = 7, name =
                                                    "Greys"))(100),cluster_cols=F,cluster_rows=F)

Generally, seven major subclasses can be identified in the Sub-PS-HA region
(Figure 1C; Table S1): Spy (RC1–RC3), deep PSpy (RC6–RC8), superficial PSpy
(RC10–RC14), HApy (RC17–RC19), most superficial PSpy (RC20–RC23), Spo and PSpo
(RC24–RC26), and L6 of the Sub and HA (RC28 and RC29). Cells from two small
clusters (RC4 and RC5) lie in the PSpy at the border between the ventral Sub
and PS. Therefore, taking away adjoining CA1 (RC9, RC15, and RC16) and PrSd
(i.e., PoS; RC27) and adding L6b of the Sub and PS from the OC (OC48 and OC49),
a total of 27 clusters or cell types were revealed in the Sub-PS-HA region.

###3: Allen Institute
##rhistory stuff

##set up the anndata

##once you've got the keys set up:

##once you've got the onehots set up:
p<-cor(onehot_class,mat)
pheatmap(p,cluster_cols=F,cluster_rows=F)

p2<-heat[c('017 CA3 Glut','016 CA1-ProS Glut','022 L5 ET CTX Glut','023 SUB-ProS Glut',
           '032 L5 NP CTX Glut','033 NP SUB Glut','031 CT SUB Glut','030 L6 CT CTX Glut',
           '028 L6b/CT ENT Glut','005 L5 IT CTX Glut','006 L4/5 IT CTX Glut','004 L6 IT CTX Glut',
           '019 L2/3 IT PPP Glut','018 L2 IT PPP-APr Glut','008 L2/3 IT ENT Glut','007 L2/3 IT CTX Glut',
           '011 L2 IT ENT-po Glut','014 LA-BLA-BMA-PA Glut','012 MEA Slc17a7 Glut'),]


##get aggregate projected weights
data<-as.data.frame(colData(mch2))

heat<-aggregate(data[,c('nmf52','nmf11','nmf63','nmf61','nmf15',
                        'nmf32','nmf40','nmf54','nmf22','nmf65',
                        'nmf53','nmf51','nmf68','nmf17','nmf78','nmf27','nmf45',
                        'nmf84','nmf62','nmf69')],
                by=list(data$Target),FUN=sum)
rownames(heat)<-heat$Group.1
heat<-heat[,-1]
pheatmap(heat,color=colorRampPalette(brewer.pal(n = 7, name =
                                                  "Greys"))(100),cluster_cols=F,cluster_rows=F)


##########umaps#############
sce<-makeVisiumRGB(sce,vars=c('nmf11','nmf15','nmf52','nmf61'))

pdf(file=here::here('plots','figures','figure_6','UMAP_rgb_hpc.pdf'),h=5,w=5)
p<-plotUMAP(sce,colour_by='RGB',text_by='fine.cell.class',point_size=0.01)+scale_color_identity()
ggrastr::rasterize(p,dpi=500,layers='Point')
dev.off()

sce<-makeVisiumRGB(sce,vars=c('nmf54','nmf17','nmf32','nmf40'))

pdf(file=here::here('plots','figures','figure_6','UMAP_rgb_sub.pdf'),h=5,w=5)
p<-plotUMAP(sce,colour_by='RGB',text_by='fine.cell.class',point_size=0.01)+scale_color_identity()
ggrastr::rasterize(p,dpi=500,layers='Point')
dev.off()



###marker gene heatmap
heat<-x@w[rownames(x@w) %in% pick,c('nmf52','nmf11','nmf63','nmf61','nmf15',
                                 'nmf32','nmf40','nmf54','nmf22','nmf65',
                                 'nmf53','nmf51','nmf68','nmf17','nmf78','nmf27','nmf45',
                                 'nmf84','nmf62','nmf69')]

pheatmap(heat,color=colorRampPalette(brewer.pal(n = 7, name =
                                                  "Greys"))(100),cluster_cols=F,cluster_rows=T)



glut<-c('GRIK1','GRIK3','GRIK4','GRIA4','GRID2','GRM8','GRM3')
gad<-c('GABRB2','GABRA2','GABRA1')
htr<-c('HTR2C','HTR4','HTR7','HTR1E','HTR2A')

nt<-c(glut,gad,htr)

genes<-genes[names(genes) %in% c('nmf52','nmf11','nmf63','nmf61','nmf15',
  'nmf32','nmf40','nmf54','nmf22','nmf65',
  'nmf53','nmf51','nmf68','nmf17','nmf78','nmf27','nmf45',
  'nmf84','nmf62','nmf69')]



#genes in diff cats:
ion<-c('')
ntr<-
np<-

####get pyr-spec markers
pyr<-c('nmf52','nmf11','nmf63','nmf61','nmf15',
       'nmf32','nmf40','nmf54','nmf22','nmf65',
       'nmf53','nmf51','nmf68','nmf17','nmf78','nmf27','nmf45',
       'nmf84','nmf62','nmf69')
###get markers for ORA after filtering mito genes and non-protein coding genes
loads<-x@w
loads<-loads[,colnames(loads) %in% pyr]
##get rid of genes with no weights in any factor
no_expr <- which(rowSums(loads) == 0)
# length(no_expr)
# length(no_expr) / nrow(loads) * 100
loads <- loads[-no_expr, ]
##now filter mito genes and non-protein coding genes
protein<-rownames(sce)[rowData(sce)$gene_type=='protein_coding']
loads<-loads[rownames(loads) %in% protein,]
mito<-rownames(sce)[which(seqnames(sce) == "chrM")]
loads<-loads[!rownames(loads) %in% mito,]
##now get markers
marks<-patternMarkers(loads,x@h,'all',1,100)



markers<-marks$PatternMarkers[c(52,11,63,61,15,32,40,54,22,65,53,
                                51,68,17,78,27,45,84,62,69)]

features<-list()
for(i in 1:length(markers)){
    features[i]<-markers[[i]][1]
}

features[2]<-'ZNF208'
features[4]<-'NRIP3'
features[7]<-'GFRA1'
features[8]<-'NPSR1-AS1'
features[9]<-'ZNF804B'
features[10]<-'AL391117.1'
features[11]<-'RXFP1'
features[12]<-'SORCS1'
features[13]<-'AGBL1'
features[14]<-'VIPR2'
features[15]<-'NIPAL2'
features[16]<-'RASGEF1B'
features[17]<-'GNAL'
features[18]<-'ADAMTS19'
features[19]<-'ETV1'
features[20]<-'HS6ST2'



heat<-x@w[features,c('nmf52','nmf11','nmf63','nmf61','nmf15',
                     'nmf32','nmf40','nmf54','nmf22','nmf65',
                     'nmf53','nmf51','nmf68','nmf17','nmf78','nmf27','nmf45',
                     'nmf84','nmf62','nmf69')]

pheatmap(heat,cluster_cols=F,cluster_rows=F,color=colorRampPalette(brewer.pal(n = 7, name =
                                                                                  "Greys"))(100))
