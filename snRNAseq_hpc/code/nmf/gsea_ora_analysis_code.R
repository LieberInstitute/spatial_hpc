keep_feature <- rowSums(counts(sce) > 0) > 0
sce <- sce[keep_feature, ]

ego_bp<-list()
for(i in 1:length(pats2)){
    print(paste0('running for pattern ',i))
    print(Sys.time())
    ego_bp[[i]]<-gseGO(
        pats2[[i]],
        ont = "BP",
        org.Hs.eg.db,
        keyType = "SYMBOL",
       # exponent = 1,
        minGSSize = 15,
        maxGSSize = 500,
        eps = 1e-100,
        pvalueCutoff = 0.001,
        pAdjustMethod = "BH",
        verbose = TRUE,
        seed = FALSE,
        by = "fgsea"    )
}

stats2<-list()
for(i in 1:length(pats)){
    print(paste0('running for pattern ',i))
    print(Sys.time())
    stats2[[i]]<-fgsea(
        stats=pats[[i]],
        pathways=projections,
        minSize = 1,
        maxSize = 1000,
        gseaParam = 1)
}

ego_wp<-list()
for(i in 1:length(genes)){
    print(paste0('running for pattern ',i))
    print(Sys.time())
    ego_wp[[i]] <- enrichPathway(gene     = genes[[i]],
                            universe      = l,
                            OrgDb         = org.Hs.eg.db,
                            #organism='hsa',
                            #ont           = "ALL",
                            pAdjustMethod = "BH",
                            pvalueCutoff  = 0.05,
                            qvalueCutoff  = 0.1,
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
