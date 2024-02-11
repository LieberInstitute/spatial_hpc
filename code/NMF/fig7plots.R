###FIG7 plots
library(SingleCellExperiment)
library(scater)
library(scran)
library(edgeR)
library(clusterProfiler)
library(org.Hs.eg.db)
library(RcppML)
###projection
# Translate from one species to the other using the orthology
names <- orthology[orthology$Column3 %in% rownames(sce.subset),]

names <- names[match(rownames(sce.subset), names$Column3),]

setdiff(names$Column3, rownames(sce.subset))

rownames(sce.subset) <- names$Column1

spe <- sce.subset

# Project Erik's human snRNAseq NMF patterns onto mouse snRNAseq

set.seed(1029)
i<-intersect(rownames(spe),rownames(x2@w))
loadings<-x2@w
loadings<-loadings[rownames(loadings) %in% i,]
spe2<-spe[rownames(spe) %in% i,]
loadings<-loadings[match(rownames(spe2),rownames(loadings)),]
proj<-project(loadings,logcounts(spe2),L1=0)


###make new sce object with the patterns as "genes"
discard<-c('nmf37','nmf28','nmf18','nmf16','nmf94','nmf92','nmf2')
proj<-proj[!rownames(proj) %in% discard,]
#proj_rescaled <- proj / rowSums(proj)
###remove very lowly expressed patterns
# Assuming your data frame is named 'col'

# Calculate the proportion of non-zero values for each column
col<-proj
col<-t(col)
proportion_nonzero <- colSums(col > 0) / nrow(col)

# Identify columns where this proportion is greater than 10%
columns_gt_10_percent <- names(proportion_nonzero[proportion_nonzero > 0.20])

# Print the column names
print(columns_gt_10_percent)
proj<-proj[columns_gt_10_percent,]

sce.pats<-SingleCellExperiment(assays=list(counts=proj))
colData(sce.pats)<-colData(sce.subset)



sce.pats$cellType_DE<-as.character(sce.pats$cellType)
sce.pats$cellType_DE<-factor(ifelse(sce.pats$cellType_DE %in% c('GC.1','GC.2'),'GC',sce.pats$cellType_DE))
summed <- aggregateAcrossCells(sce.pats,
                               ids=colData(sce.pats)
                               [c('Sample','cellType_DE')])
#do DE analysis
de.dge <- pseudoBulkDGE(summed,
                        label=summed$cellType_DE,
                        design=~condition,
                        coef=2,
                        condition=summed$condition,
)
is.de <- decideTestsPerLabel(de.dge, threshold=0.1)
summarizeTestsPerLabel(is.de)

##Set up data for double-sided barplot for fig 1
x<-as.data.frame(summarizeTestsPerLabel(de.dge))



#aggregateAcrossCells to pseudobulk
summed <- aggregateAcrossCells(sce.pats,
                               ids=colData(sce.pats)
                               [c('Sample')])

#Filter out lowly expressed genes


#make DGElist
y <- DGEList(counts(summed), samples=colData(summed),
             genes=rowData(summed))
y$samples$condition<-factor(y$samples$condition)
y<-y[filterByExpr(y, group=summed$Sample),]

#make design matrix
design <- model.matrix(~condition,y$samples)

#normalization
y <- calcNormFactors(y)

#estimate dispersion parameter
y <- estimateDisp(y, design = design)

#run DE analysis and look at results
efit <- glmQLFit(y, design = design)
elrt <- glmQLFTest(efit, coef=2)
res<-topTags(elrt,n=Inf)
x<-as.data.frame(res)

library(dplyr)
library(tidyr)
library(ggplot2)

# Combine all data frames into one long-form data frame
long_df <- lapply(seq_along(de.dge), function(i) {
    data <- as.data.frame(de.dge[[i]])
    data$group <- names(de.dge)[i]
    data$nmf_pattern <- rownames(data)
    return(data)
}) %>%
    bind_rows()

# Add a column for -log10(FDR)
long_df$neg_log10_FDR <- -log10(long_df$FDR)

# Filter out rows where -log10(FDR) < 1
filtered_df <- long_df %>%
    filter(neg_log10_FDR >= 1)

# Get all unique group names (even those filtered out)
all_groups <- factor(long_df$group,levels=c('GC','CA4','CA3.1','CA3.2','CA1',
                                            'PS.1','PS.2','Sub','L2/3','L5/Po','L6/6b',
                                            'GABA.1','GABA.2','GABA.3','GABA.4','GABA.5'))
all_groups_levels <- levels(all_groups)
#Make second DF with only some factors
long_df2<-filtered_df[filtered_df$nmf_pattern %in% paste0('nmf',c(13,14,20,91,98,10)),]



# Create a dot plot with remaining data
pdf(file=here::here('plots','figures','figure_7','dotplot_de_analysis.pdf'),h=4,w=16)
ggplot(long_df2, aes(x = group, y = nmf_pattern, size = neg_log10_FDR, color = logFC)) +
    geom_point() +
    theme_minimal() +
    labs(x = "Group", y = "NMF Pattern", size = "-log10(FDR)", color = "logFC") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          text=element_text(size=24)) +
    scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    scale_size_continuous(range = c(1.2, max(long_df$neg_log10_FDR, na.rm = TRUE)))+
    scale_x_discrete(limits = all_groups_levels)+
    guides(size = guide_legend(reverse = TRUE)) # Reverse the size legend
dev.off()


###multicolor spotplots!
speb <- spe[, (colData(spe)$sample_id == 'V11L05-336_A1')]
speb$sample_id <- droplevels(speb$sample_id)
speb$sample_id <- as.character(speb$sample_id)
samples <- unique(speb$sample_id)
speb$sample_id <- factor(speb$sample_id, levels = samples)
samples
speb$brnum <- droplevels(speb$brnum)

pdf(file=here::here('plots','figures','figure_7','rgb_gcs.pdf'),h=3.9,w=3.9)
plotVisiumRGB(speb,green ='nmf10',blue ='nmf14',image=F,highlight='neuron_cell_body')
dev.off()

pdf(file=here::here('plots','figures','figure_7','rgb_nongcs.pdf'),h=3.9,w=3.9)
plotVisiumRGB(speb,yellow ='nmf91',cyan ='nmf13',pink='nmf20',image=F,highlight='neuron_cell_body')
dev.off()



###boxplots!
##mouse
colData(sce.subset)<-cbind(colData(sce.subset),t(proj_rescaled))
pdf(file=here::here('plots','figures','figure_7','mouse_nmf_10_boxplots.pdf'),h=5,w=11)
ggcells(sce.subset, mapping=aes(x=cell.type.mouse, y=nmf10,fill=condition)) +
    geom_boxplot(outlier.size = 0.05)+theme(axis.text.x = element_text(angle = 90),
                                            text=element_text(size = 36,colour='black'))+
    theme(legend.position='none',
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))#+scale_fill_manual(values=sn.mid.palette)
dev.off()

pdf(file=here::here('plots','figures','figure_7','mouse_nmf_14_boxplots.pdf'),h=5,w=11)
ggcells(sce.subset, mapping=aes(x=cell.type.mouse, y=nmf14,fill=condition)) +
    geom_boxplot(outlier.size = 0.05)+theme(axis.text.x = element_text(angle = 90),
                                            text=element_text(size = 36,colour='black'))+
    theme(legend.position='none',
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))#+scale_fill_manual(values=sn.mid.palette)
dev.off()


pdf(file=here::here('plots','figures','figure_7','mouse_nmf_13_boxplots.pdf'),h=5,w=11)
ggcells(sce.subset, mapping=aes(x=cell.type.mouse, y=nmf13,fill=condition)) +
    geom_boxplot(outlier.size = 0.05)+theme(axis.text.x = element_text(angle = 90),
                                            text=element_text(size = 36,colour='black'))+
    theme(legend.position='none',
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))#+scale_fill_manual(values=sn.mid.palette)
dev.off()

pdf(file=here::here('plots','figures','figure_7','mouse_nmf_91_boxplots.pdf'),h=5,w=11)
ggcells(sce.subset, mapping=aes(x=cell.type.mouse, y=nmf91,fill=condition)) +
    geom_boxplot(outlier.size = 0.05)+theme(axis.text.x = element_text(angle = 90),
                                            text=element_text(size = 36,colour='black'))+
    theme(legend.position='none',
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))#+scale_fill_manual(values=sn.mid.palette)
dev.off()

pdf(file=here::here('plots','figures','figure_7','mouse_nmf_98_boxplots.pdf'),h=5,w=11)
ggcells(sce.subset, mapping=aes(x=cell.type.mouse, y=nmf98,fill=condition)) +
    geom_boxplot(outlier.size = 0.05)+theme(axis.text.x = element_text(angle = 90),
                                            text=element_text(size = 36,colour='black'))+
    theme(legend.position='none',
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))#+scale_fill_manual(values=sn.mid.palette)
dev.off()

pdf(file=here::here('plots','figures','figure_7','mouse_nmf_20_boxplots.pdf'),h=5,w=11)
ggcells(sce.subset, mapping=aes(x=cell.type.mouse, y=PDGFB,fill=condition)) +
    geom_boxplot(outlier.size = 0.05)+theme(axis.text.x = element_text(angle = 90),
                                            text=element_text(size = 36,colour='black'))+
    theme(legend.position='none',
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))#+scale_fill_manual(values=sn.mid.palette)
dev.off()

##human
pdf(file=here::here('plots','figures','figure_7','human_nmf_10_boxplots.pdf'),h=6,w=11)
ggcells(sce.no, mapping=aes(x=fine.cell.class, y=nmf10,fill=mid.cell.class)) +
    geom_boxplot(outlier.size = 0.05)+theme(axis.text.x = element_text(angle = 90),
                                            text=element_text(size = 28,colour='black'))+
    theme(legend.position='none',
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))+scale_fill_manual(values=sn.mid.palette)+scale_y_continuous(labels = function(x) format(x, scientific = TRUE))
dev.off()

pdf(file=here::here('plots','figures','figure_7','human_nmf_14_boxplots.pdf'),h=6,w=11)
ggcells(sce.no, mapping=aes(x=fine.cell.class, y=nmf14,fill=mid.cell.class)) +
    geom_boxplot(outlier.size = 0.05)+theme(axis.text.x = element_text(angle = 90),
                                            text=element_text(size = 28,colour='black'))+
    theme(legend.position='none',
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))+scale_fill_manual(values=sn.mid.palette)
dev.off()

pdf(file=here::here('plots','figures','figure_7','human_nmf_13_boxplots.pdf'),h=6,w=11)
ggcells(sce.no, mapping=aes(x=fine.cell.class, y=nmf13,fill=mid.cell.class)) +
    geom_boxplot(outlier.size = 0.05)+theme(axis.text.x = element_text(angle = 90),
                                            text=element_text(size = 28,colour='black'))+
    theme(legend.position='none',
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))+scale_fill_manual(values=sn.mid.palette)+scale_y_continuous(labels = function(x) format(x, scientific = TRUE))
dev.off()

pdf(file=here::here('plots','figures','figure_7','human_nmf_20_boxplots.pdf'),h=6,w=11)
ggcells(spe, mapping=aes(x=brnum, y=nmf13,fill=brnum)) +
    geom_boxplot(outlier.size = 0.05)+theme(axis.text.x = element_text(angle = 90),
                                            text=element_text(size = 28,colour='black'))+
    theme(legend.position='none',
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))#+scale_fill_manual(values=sn.mid.palette)
dev.off()

pdf(file=here::here('plots','figures','figure_7','human_nmf_91_boxplots.pdf'),h=6,w=11)
ggcells(sce.no, mapping=aes(x=fine.cell.class, y=nmf91,fill=mid.cell.class)) +
    geom_boxplot(outlier.size = 0.05)+theme(axis.text.x = element_text(angle = 90),
                                            text=element_text(size = 28,colour='black'))+
    theme(legend.position='none',
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))+scale_fill_manual(values=sn.mid.palette)
dev.off()

pdf(file=here::here('plots','figures','figure_7','human_nmf_98_boxplots.pdf'),h=6,w=11)
ggcells(sce.no, mapping=aes(x=fine.cell.class, y=nmf98,fill=mid.cell.class)) +
    geom_boxplot(outlier.size = 0.05)+theme(axis.text.x = element_text(angle = 90),
                                            text=element_text(size = 28,colour='black'))+
    theme(legend.position='none',
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))+scale_fill_manual(values=sn.mid.palette)
dev.off()

##human GC only plots!
pdf(file=here::here('plots','figures','figure_7','human_nmf_10_gcs.pdf'),h=5,w=9.5)
ggcells(sce.gc, mapping=aes(x=superfine.cell.class, y=nmf10,fill=fine.cell.class)) +
    geom_boxplot(outlier.size = 0.05)+theme(axis.text.x = element_text(angle = 90),
                                            text=element_text(size = 36,colour='black'))+
    theme(legend.position='none',
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))+scale_fill_manual(values=sn.fine.palette)+scale_y_continuous(labels = function(x) format(x, scientific = TRUE))
dev.off()

pdf(file=here::here('plots','figures','figure_7','human_nmf_14_gcs.pdf'),h=5,w=9.5)
ggcells(sce.gc, mapping=aes(x=superfine.cell.class, y=nmf14,fill=fine.cell.class)) +
    geom_boxplot(outlier.size = 0.05)+theme(axis.text.x = element_text(angle = 90),
                                            text=element_text(size = 36,colour='black'))+
    theme(legend.position='none',
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))+scale_fill_manual(values=sn.fine.palette)+scale_y_continuous(labels = function(x) format(x, scientific = TRUE))
dev.off()

pdf(file=here::here('plots','figures','figure_7','dotplot_goterms.pdf'),h=12,w=8)
ggplot(merged2, aes(x = Cluster, y = Description, color = negFDR, size = geneRatio)) +
    geom_point() +
    theme_minimal() +
    labs(x = "Group", y = "NMF Pattern", color = "-log10(FDR)", size = "gene ratio") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1), text=element_text(size=30))+
              scale_color_distiller(type = "seq",
        palette = rev('Greys'),
        direction=1)
dev.off()

    #scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    #scale_size_continuous(range = c(1.2, max(long_df$neg_log10_FDR, na.rm = TRUE)))+
    #scale_y_discrete(limits = rev+
    #guides(size = guide_legend(reverse = TRUE)) # Reverse the size legend

pdf(file=here::here('plots','figures','figure_7','gene_heatmap.pdf'),h=2,w=6)
pheatmap(t(x@w[genesplot,c(10,14,91,20,13)]),
         cluster_rows=F,cluster_cols=F,breaks=seq(0,0.002,0.002/100),
         color=colorRampPalette(brewer.pal(n = 7, name ="Greys"))(100))
dev.off()

