keep_feature <- rowSums(counts(sce) > 0) > 0
sce <- sce[keep_feature, ]

ego_bp<-list()
for(i in 1:length(scoreList)){
    print(paste0('running for pattern ',i))
    print(Sys.time())
    ego_bp[[i]]<-gseGO(
        scoreList[[i]],
        ont = "BP",
        org.Hs.eg.db,
        keyType = "SYMBOL",
       # exponent = 1,
        minGSSize = 10,
        maxGSSize = 500,
        eps = 1e-100,
        pvalueCutoff = 0.001,
        pAdjustMethod = "BH",
        verbose = TRUE,
        seed = FALSE,
        by = "fgsea"    )
}

##make sure all genes have a weight for at least one pattern

loads<-x@w
loads<-loads[,!colnames(loads) %in% discard]
no_expr <- which(rowSums(loads) == 0)
length(no_expr)
# [1] 5118
length(no_expr) / nrow(loads) * 100
# [1] 13.98322
loads <- loads[-no_expr, ]
dim(loads)
# loads<-loads[,!colnames(loads) %in% nonSpec]
# dim(loads)
# no_expr <- which(rowSums(loads) == 0)
# length(no_expr)
# # [1] 5118
# length(no_expr) / nrow(loads) * 100
# # [1] 13.98322
# loads <- loads[-no_expr, ]
# dim(loads)
##now filter mito genes and non-protein coding genes plz
protein<-rownames(sce)[rowData(sce)$gene_type=='protein_coding']
loads<-loads[rownames(loads) %in% protein,]
mito<-rownames(sce)[which(seqnames(sce) == "chrM")]
loads<-loads[!rownames(loads) %in% mito,]
##now get markers
marks<-patternMarkers(loads,x@h,'all',1,100)
}

ego_bp<-list()
for(i in 1:length(genes)){
    print(paste0('running for pattern ',i))
    print(Sys.time())


    test<-loads[rownames(loads) %in% genes[[70]],]
    test<-test[order(test[,70],decreasing=T),]
    plot(test[,70])
    test<-as.data.frame(test)
    test$index=c(1:nrow(test))
    f2 <- lm(nmf75 ~ index, data = test)

    seg2 <- segmented(f2,
                      seg.Z = ~index,
                      npsi = 2
    )
    chosen<-rownames(test)[1:round(seg2$psi[2,2])]

    ego <- enrichGO(gene     = genes[[70]],
                            universe      = rownames(loads),
                            OrgDb         = org.Hs.eg.db,
                            #organism='hsa',
                            ont           = "BP",
                            pAdjustMethod = "BH",
                            pvalueCutoff  = 0.05,
                            qvalueCutoff  = 0.05,
                            readable      = TRUE,
                            keyType = 'SYMBOL')
}
keep_feature <- rowSums(counts(sce) > 0) > 0
sce <- sce[keep_feature, ]

colnames(colData(spe))[88:187]<-c('high_det1','high_det2',,'ExcN',,)

data<-as.data.frame(sce$brnum)
colnames(data)<-'brnum'
onehot_brnum <-  dcast(data = data, rownames(data) ~ brnum, length)
rownames(onehot_brnum)<-onehot_brnum[,1]
onehot_brnum[,1]<-as.numeric(onehot_brnum[,1])
onehot_brnum<-onehot_brnum[order(onehot_brnum[,1],decreasing=F),]
onehot_brnum[,1]<-NULL


#unknown patterns:6,16,

seu@reductions[['nbp']] <- new("DimReduc",
                                           cell.embeddings = t(x@h),
                                           feature.loadings = x@w,
                                           assay.used = 'logcounts',
                                           stdev = x@d,
                                           key = 'Nbp',
                                           global = FALSE)
                                           #misc = list("cv_data" = cv_data))


seu<-RunGSEA(
    seu,
    reduction = "nbp",
    species = "Homo sapiens",
    category = "C2",
    subcategory='CP:REACTOME',
    min.size = 10,
    max.size = 500,
    dims = NULL,
    verbose = TRUE,
    padj.sig = 0.01
)

seu<-RunGSEA(
    seu,
    reduction = "nbp",
    species = "Homo sapiens",
    category = "C5",
    min.size = 10,
    max.size = 500,
    dims = NULL,
    verbose = TRUE,
    padj.sig = 0.01
)

seu<-RunGSEA(
    seu,
    reduction = "nbp",
    species = "Homo sapiens",
    category = "C5",
    subcategory='GO:BP',
    min.size = 10,
    max.size = 500,
    dims = NULL,
    verbose = TRUE,
    padj.sig = 0.01
)

go_bp<-RunGSEA(
    seu,
    reduction = "nbp",
    species = "Homo sapiens",
    category = "C5",
    subcategory='GO:bp',
    min.size = 10,
    max.size = 500,
    dims = NULL,
    verbose = TRUE,
    padj.sig = 0.01
)

df <- seu@reductions[['nbp']]@misc$gsea$padj
df<-df[,colnames(df) %in% c('nbp76','nbp79','nbp19','nbp81','nbp3','nbp44','nbp43','nbp36')]
terms <- c()
for (i in 1:ncol(df)) {
    terms_i <- df[, i]
    idx <- terms_i > -log10(0.05)
    terms_i <- terms_i[idx]
    terms_j <- df[idx, i]
    v <- sort(terms_j, decreasing = TRUE)
    if (length(v) > max.terms.per.factor) {
        terms <- c(terms, names(v)[1:max.terms.per.factor])
    } else {
        terms <- c(terms, names(v))
    }
}
terms <- unique(terms)
df <- df[terms, ]
