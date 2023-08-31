keep_feature <- rowSums(counts(sce) > 0) > 0
sce <- sce[keep_feature, ]

ego_bp<-list()
for(i in 1:length(genes)){
    print(paste0('running for pattern ',i))
    print(Sys.time())
ego_bp[[i]] <- enrichGO(gene     = genes[[i]],
                universe      = scerows,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.1,
                readable      = TRUE,
                keyType = 'SYMBOL')
}

ego_cc<-list()
for(i in 1:length(genes)){
    print(paste0('running for pattern ',i))
    print(Sys.time())
    ego_cc[[i]] <- enrichGO(gene     = genes[[i]],
                         universe      = rownames(sce),
                         OrgDb         = org.Hs.eg.db,
                         ont           = "CC",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.05,
                         qvalueCutoff  = 0.1,
                         readable      = TRUE,
                         keyType = 'SYMBOL')
}

ego_mf<-list()
for(i in 1:length(genes)){
    print(paste0('running for pattern ',i))
    print(Sys.time())
    ego_mf[[i]] <- enrichGO(gene     = genes[[i]],
                            universe      = rownames(sce),
                            OrgDb         = org.Hs.eg.db,
                            ont           = "MF",
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

seu@reductions[['nmf']] <- new("DimReduc",
                                           cell.embeddings = t(x@h),
                                           feature.loadings = x@w,
                                           assay.used = 'logcounts',
                                           stdev = x@d,
                                           key = 'NMF',
                                           global = FALSE)
                                           #misc = list("cv_data" = cv_data))


seu<-RunGSEA(
    seu,
    reduction = "nmf",
    species = "Homo sapiens",
    category = "C2",
    subcategory='CP:KEGG',
    min.size = 10,
    max.size = 500,
    dims = NULL,
    verbose = TRUE,
    padj.sig = 0.01
)
