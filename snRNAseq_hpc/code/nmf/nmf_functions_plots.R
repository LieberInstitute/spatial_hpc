marker_genes <- function(W, n_top_genes) {
    max_other_patterns <- apply(W, 1, function(gene) max(gene[-which.max(gene)]))
    score <- W - max_other_patterns

    results <- apply(score, 2, function(pattern) {
        top_genes_indices <- order(pattern, decreasing = TRUE)[1:n_top_genes]
        list(
            gene_names = rownames(W)[top_genes_indices],
            scores = pattern[top_genes_indices]
        )
    })

    results
}

genes <- marker_genes(W = x$w, n_top_genes = 25)

pats<-t(x$h)
# pats<-colData(sce)[,c(54:153)]
data<-as.data.frame(sce$fine.type)
colnames(data)<-'fine.type'
onehot_fine.type <-  dcast(data = data, rownames(data) ~ fine.type, length)
rownames(onehot_fine.type)<-onehot_fine.type[,1]
onehot_fine.type[,1]<-as.numeric(onehot_fine.type[,1])
onehot_fine.type<-onehot_fine.type[order(onehot_fine.type[,1],decreasing=F),]
onehot_fine.type[,1]<-NULL
#pats<-colData(sce)[,c(109:183)]
heat<-cor(onehot_fine.type,as.data.frame(pats))
dist_matrix <- dist(heat)

# install the package if you haven't already
# install.packages("seriation")
library(seriation)

# Apply seriation on rows and columns separately

heat_ordered <- heat[get_order(seriate(heat, method = "BEA_TSP")),
                     get_order(seriate(t(heat), method = "BEA_TSP"))]

# Create the heatmap with the ordered data
heatmap(as.matrix(heat_ordered), Rowv = NA, Colv = NA, col = colorRampPalette(c("navy", "white", "firebrick3"))(256))


onehot_fine.type2<-onehot_fine.type[,c(4,5,29,30,31,11,32,27,48,16,15,51,44,13,14,
                                    17:20,25,8,9,7,6,47,46,41:43,40,39,34:38,
                                    21,1:3,9,49,
                                    45,12,33,28,24,23,50,22,26)]

reg<-registration_wrapper(
    sce,
    var_registration='fine.type',
    var_sample_id='brnum',
    covars = 'sex',
    gene_ensembl = 'gene_id',
    gene_name = 'gene_name',
    min_ncells = 10,
    pseudobulk_rds_file = NULL
)

reg<-registration_wrapper(
    spe,
    var_registration='cluster',
    var_sample_id='sample_id',
    covars = 'sex',
    gene_ensembl = 'gene_id',
    gene_name = 'gene_name',
    min_ncells = 50,
    pseudobulk_rds_file = NULL
)



# Define function
get_annotation <- function(df) {
    # Make sure df is a data frame
    if (!is.data.frame(df)) {
        stop("Input must be a dataframe")
    }

    # Initialize a data frame to store the factors and annotations
    key <- data.frame(factor = colnames(df), annotation = character(ncol(df)))

    # Loop through each column and find max
    for (i in 1:ncol(df)) {
        # If max is 0, assign 'other' to annotation
        if (max(df[[i]]) == 0) {
            key$annotation[i] <- 'Other'
        } else {
            # Assign rowname with max value to annotation
            key$annotation[i] <- rownames(df)[which.max(df[[i]])]
        }
    }

    # Return key dataframe
    return(key)
}

# Assuming heat_cor is your dataframe
key <- get_annotation(heat_cor)

# Print key
print(key)

key$annotation<-ifelse(key$annotation %in% c('PV','SST','LAMP5_LHX6','LAMP5_KIT',
+                              'VIP','HTR3A','MEIS2'),'GABA',key$annotation)
key$annotation<-factor(key$annotation,
                       levels=c('GC','MC','CA3','CA2','CA1',
                        'Sub.Sup','Sub.Deep',paste0('RHP.L2/3_',c(1:5)),
                         'RHP.L5','RHP.L6_1','RHP.L6_2','RHP.L6b',
                            'HATA','Amy','GABA','Oligo','Astro','Micro','OPC',
                        'Ependy','CP','Vascular','Tcell/Macro','Other'))

# Min-Max normalization
pats <- apply(pats, 2, function(x) (x - min(x)) / (max(x) - min(x)))

# Initialize new vector
nmfType <- rep(NA, nrow(pats))

# Subset 'pats' to only include columns in 'cols'
pats_subset <- pats[, colnames(pats) %in% cols]

# For each row in 'pats_subset', find the column with the highest value
for (i in 1:nrow(pats_subset)) {
    # Check if all values are 0
    if(all(pats_subset[i, ] == 0)) {
        nmfType[i] <- "Other"
    } else {
        # Find the column name for the max value
        max_col <- colnames(pats_subset)[which.max(pats_subset[i,])]

        # Find the matching annotation in 'key'
        nmfType[i] <- key$annotation[which(key$factor == max_col)]
    }
}


spe$nmfType<-factor(spe$nmfType,
       levels=c('GC','MC','CA3','CA2','CA1',
                'Sub.Sup','Sub.Deep','RHP.L2/3',
                'RHP.L5','RHP.L6/6b','HATA','Amy',
                'GABA','Oligo','Astro','Micro','OPC',
                'Ependy','CP','Vascular','Tcell/Macro'))


31 18 80 28 21
6 16 20 59

94
12
55

pats<-t(x$h)
# pats<-colData(sce)[,c(54:153)]
data2<-as.data.frame(sce$brnum)
colnames(data2)<-'brnum'
onehot_fine.type <-  dcast(data = data2, rownames(data2) ~ brnum, length)
rownames(onehot_fine.type)<-onehot_fine.type[,1]
onehot_fine.type[,1]<-as.numeric(onehot_fine.type[,1])
onehot_fine.type<-onehot_fine.type[order(onehot_fine.type[,1],decreasing=F),]
onehot_fine.type[,1]<-NULL
#pats<-colData(sce)[,c(109:183)]
heat2<-cor(onehot_fine.type,as.data.frame(pats))


set.seed(1029)
i<-intersect(rownames(spe),rownames(x$w))
loadings<-x$w
loadings<-loadings[rownames(loadings) %in% i,]
spe2<-spe[rownames(spe) %in% i,]
loadings<-loadings[match(rownames(spe2),rownames(loadings)),]

proj<-project(loadings,logcounts(spe2),L1=0)

index<-which(colMaxs(heat_pyr)>0.25)
pyrpats<-colnames(heat_pyr)[index]
heat_pyr<-heat_pyr[,index]

dim(t(heat_pyr))
heat_index<-t(heat_pyr)
rownames(heat_index)<-1:24
index<-c(13,1,5,18,16,2,8,9,15,4,20,14,21,12,11,3,6,24,23,17,22,7,10,19)
heat_index<-heat_index[index,]



pyr_cor<-layer_stat_cor(stats,spe_stats)
pyr_cor<-pyr_cor[colnames(heat_index),c(1:8)]
pyr_cor<-pyr_cor[colnames(heat_index),c(2:8)]
pdf(file='nmf_heatmap_pyr.pdf',w=10,h=6)
pheatmap(heat_index,cluster_cols=F,cluster_rows=F,color=viridis(15))
dev.off()

marks<-c('ABI3BP')
