library( "biomaRt" )
library(zellkonverter)
library(SingleCellExperiment)
library(RcppML)

mch<-readH5AD('rs2_mch_matrix.h5ad')

mart = useMart('ensembl')
# list all the ensembl database of organisms
listDatasets(mart)
#choose database of your interest
ensembl = useMart( "ensembl", dataset = "mmusculus_gene_ensembl" )
# choose attributes of your interest
listAttributes(ensembl)
gene <- getBM( attributes = c("ensembl_gene_id","external_gene_name"),values = rownames(mch),mart = ensembl)
#Macth your transcript id with ensembl_transcript_id
id <- match(test$Name , gene$ensembl_transcript_id)
#Add Gene symbol column in your data frame
test$Symbol <- gene$external_gene_name[id]
head(test)

mart <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
symb <- getBM(attributes = c("ensembl_gene_id","mgi_symbol"),
              filters = "ensembl_gene_id", values = rownames(mch),
              mart = mart)
symbs <- symb$mgi_symbol[match(rownames(mch), symb$ensembl_gene_id, nomatch = NA)]

rowData(mch)$gene_name<-symbs
rowData(mch)$start<-NULL
rowData(mch)$end<-NULL


mch<-mch[!is.na(rowData(mch)$gene_name),]
rownames(mch)<-rowData(mch)$gene_name


# Translate from one mchcies to the other using the orthology
names <- orthology[orthology$Column3 %in% rownames(mch),]

names <- names[match(rownames(mch), names$Column3),]

setdiff(names$Column3, rownames(mch))

rownames(mch) <- names$Column1

mch <- mch

# Project human snRNAseq NMF patterns onto mouse snRNAseq

###discard patterns associated with specific samples/rounds or technical/QC variables
discard<-c('nmf37','nmf28','nmf18','nmf16','nmf94','nmf92','nmf2')

###get markers for ORA after filtering mito genes and non-protein coding genes
loads<-x@w
#oads<-loads[,!colnames(loads) %in% discard]
##get rid of genes with no weights in any factor
no_expr <- which(rowSums(loads) == 0)
# length(no_expr)
# length(no_expr) / nrow(loads) * 100
loads <- loads[-no_expr, ]
dim(loads)
#[1] 19363    93

set.seed(1029)
i<-intersect(rownames(mch),rownames(loads))
loadings<-loads

loadings<-loadings[rownames(loadings) %in% i,]
mch2<-mch[rownames(mch) %in% i,]
loadings<-loadings[match(rownames(mch2),rownames(loadings)),]
proj<-project(loadings,assay(mch2,'zscores'),L1=0)


# Calculate the mean and standard deviation for each column
column_means <- apply(logcounts(mch), 2, mean)
column_sds <- apply(logcounts(mch), 2, sd)

# Convert the logcounts matrix to z-scores
z_scores <- sweep(logcounts(mch), 2, column_means, "-")  # Subtract the mean
z_scores <- sweep(z_scores, 2, column_sds, "/")  # Divide by the standard deviation

# Your matrix of z-scores is now in z_scores
dim(z_scores)


# subset mch matrix
index<-levels(mch$Cocluster)[c(58:81,132:175)]
mch.subset<-mch[,mch$Cocluster %in% index]

mch.subset$Source=droplevels(mch.subset$Source)
mch.subset$Target=droplevels(mch.subset$Target)
mch.subset$Cocluster=droplevels(mch.subset$Cocluster)
mch.subset$Subclass=droplevels(mch.subset$Subclass)

col_sums<-colSums(
    as.data.frame(colData(mch)[,c(19:118)]))

# Rescale each column
colData(mch)[,c(19:118)] <- t(apply(as.data.frame(colData(mch)[,c(19:118)]),
                                         1, function(row) row / col_sums))

# Check if columns now sum to 1
print(colSums(as.data.frame(colData(mch)[,c(19:118)])))



data<-as.data.frame(colData(mch))
heat<-aggregate(data[,c('nmf52','nmf11','nmf63','nmf61','nmf15',
                        'nmf32','nmf40','nmf54','nmf22','nmf65',
                        'nmf53','nmf51','nmf68','nmf17','nmf78','nmf27','nmf45',
                        'nmf84','nmf62','nmf69')],
                by=list(data$Subclass),FUN=sum)
rownames(heat)<-heat$Group.1
heat<-heat[,-1]

colnames(heat)<-c("nmf52 (CA4)", "nmf11 (CA3.1)", "nmf63 (CA3.2)", "nmf61 (CA2)",
                  "nmf15 (CA1)", "nmf32 (CA1/ProS", "nmf40 (Sub.1)", "nmf54 (Sub.2)",
                  "nmf22 (L6.2)", "nmf65 (L6.1/L6b)", "nmf53 (L6b)", "nmf51 (L5.1)",
                  "nmf68 (L5.1/L5.2)", "nmf17 (L2/3.1)", "nmf78 (L2/3.5)", "nmf27 (L2/3.2)",
                  "nmf45 (L2/3.2, L2/3.3, L2/3.6)","nmf84 (L2/3.4)", "nmf62 (HATA)", "nmf69 (HATA)")

pdf(here::here('plots','figures','figure_6','heatmap_allen_subclass.pdf'),h=4,w=7)
pheatmap(heat,color=colorRampPalette(brewer.pal(n = 7, name =
                                                    "Greys"))(100),cluster_cols=F,cluster_rows=T)
dev.off()

data<-as.data.frame(colData(mch))
heat<-aggregate(data[,c('nmf52','nmf11','nmf63','nmf61','nmf15',
                        'nmf32','nmf40','nmf54','nmf22','nmf65',
                        'nmf53','nmf51','nmf68','nmf17','nmf78','nmf27','nmf45',
                        'nmf84','nmf62','nmf69')],
                by=list(data$Cocluster),FUN=sum)
rownames(heat)<-heat$Group.1
heat<-heat[,-1]

colnames(heat)<-c("nmf52 (CA4)", "nmf11 (CA3.1)", "nmf63 (CA3.2)", "nmf61 (CA2)",
                  "nmf15 (CA1)", "nmf32 (CA1/ProS", "nmf40 (Sub.1)", "nmf54 (Sub.2)",
                  "nmf22 (L6.2)", "nmf65 (L6.1/L6b)", "nmf53 (L6b)", "nmf51 (L5.1)",
                  "nmf68 (L5.1/L5.2)", "nmf17 (L2/3.1)", "nmf78 (L2/3.5)", "nmf27 (L2/3.2)",
                  "nmf45 (L2/3.2, L2/3.3, L2/3.6)","nmf84 (L2/3.4)", "nmf62 (HATA)", "nmf69 (HATA)")

pdf(here::here('plots','figures','figure_6','heatmap_cocluster.pdf'),h=4,w=7)
pheatmap(heat,color=colorRampPalette(brewer.pal(n = 7, name =
                                                    "Greys"))(100),cluster_cols=F,cluster_rows=T)
dev.off()


data<-as.data.frame(colData(mch))
heat<-aggregate(data[,c('nmf52','nmf11','nmf63','nmf61','nmf15',
                        'nmf32','nmf40','nmf54','nmf22','nmf65',
                        'nmf53','nmf51','nmf68','nmf17','nmf78','nmf27','nmf45',
                        'nmf84','nmf62','nmf69')],
                by=list(data$Target),FUN=sum)
rownames(heat)<-heat$Group.1
heat<-heat[,-1]

colnames(heat)<-c("nmf52 (CA4)", "nmf11 (CA3.1)", "nmf63 (CA3.2)", "nmf61 (CA2)",
                  "nmf15 (CA1)", "nmf32 (CA1/ProS", "nmf40 (Sub.1)", "nmf54 (Sub.2)",
                  "nmf22 (L6.2)", "nmf65 (L6.1/L6b)", "nmf53 (L6b)", "nmf51 (L5.1)",
                  "nmf68 (L5.1/L5.2)", "nmf17 (L2/3.1)", "nmf78 (L2/3.5)", "nmf27 (L2/3.2)",
                  "nmf45 (L2/3.2, L2/3.3, L2/3.6)","nmf84 (L2/3.4)", "nmf62 (HATA)", "nmf69 (HATA)")

pdf(here::here('plots','figures','figure_6','heatmap_target.pdf'),h=4,w=7)
pheatmap(heat,color=colorRampPalette(brewer.pal(n = 7, name =
                                                    "Greys"))(100),cluster_cols=F,cluster_rows=T)
dev.off()

data<-as.data.frame(colData(mch))
heat<-aggregate(data[,c('nmf52','nmf11','nmf63','nmf61','nmf15',
                        'nmf32','nmf40','nmf54','nmf22','nmf65',
                        'nmf53','nmf51','nmf68','nmf17','nmf78','nmf27','nmf45',
                        'nmf84','nmf62','nmf69')],
                by=list(data$Source),FUN=sum)
rownames(heat)<-heat$Group.1
heat<-heat[,-1]

colnames(heat)<-c("nmf52 (CA4)", "nmf11 (CA3.1)", "nmf63 (CA3.2)", "nmf61 (CA2)",
                  "nmf15 (CA1)", "nmf32 (CA1/ProS", "nmf40 (Sub.1)", "nmf54 (Sub.2)",
                  "nmf22 (L6.2)", "nmf65 (L6.1/L6b)", "nmf53 (L6b)", "nmf51 (L5.1)",
                  "nmf68 (L5.1/L5.2)", "nmf17 (L2/3.1)", "nmf78 (L2/3.5)", "nmf27 (L2/3.2)",
                  "nmf45 (L2/3.2, L2/3.3, L2/3.6)","nmf84 (L2/3.4)", "nmf62 (HATA)", "nmf69 (HATA)")

pdf(here::here('plots','figures','figure_6','heatmap_target.pdf'),h=4,w=7)
pheatmap(heat,color=colorRampPalette(brewer.pal(n = 7, name =
                                                    "Greys"))(100),cluster_cols=F,cluster_rows=T)
dev.off()

