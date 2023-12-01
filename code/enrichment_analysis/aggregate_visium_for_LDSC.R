## Pseudo-bulk for PRECAST_k16
spe<-SingleCellExperiment(assays=list(counts=counts(spe)),colData=colData(spe),rowData=rowData(spe))
spe_pseudo <- aggregateAcrossCells(
    spe,
    DataFrame(
        cluster = colData(spe)$domain
        #sample_id = spe$sample_id
    ))

spe_pseudo$cluster <- factor(spe_pseudo$cluster)
spe_pseudo <- spe_pseudo[, spe_pseudo$ncells >= 50]
colData(spe_pseudo)<-colData(spe_pseudo)[,c(21,22,24,31:33,226:230)]
colData(spe_pseudo)$tissue.type<-factor(
    ifelse(spe_pseudo$cluster %in% levels(spe_pseudo$cluster)[1:9],'Neuron',
           ifelse(spe_pseudo$cluster %in% levels(spe_pseudo$cluster)[10:13],'Neuropil',
                  ifelse(spe_pseudo$cluster %in% levels(spe_pseudo$cluster)[14:16],'WM',
                         'Vasc/CSF'))))

dim(spe_pseudo)
# 31483   409

##
pdf(file = here::here("plots","08_pseudobulk", "PRECAST", "histogram_boxplot_precast16.pdf"), width = 14, height = 14)
hist(spe_pseudo$ncells, breaks = 200)
boxplot(ncells ~ spe_pseudo$cluster, data = colData(spe_pseudo))
dev.off()

#find a good expression cutoff using edgeR::filterByExpr
rowData(spe_pseudo)$high_expr_group_sample_id <- filterByExpr(spe_pseudo, group = spe_pseudo$sample_id)
rowData(spe_pseudo)$high_expr_group_cluster <- filterByExpr(spe_pseudo, group = spe_pseudo$cluster)

summary(rowData(spe_pseudo)$high_expr_group_sample_id)
# Mode   FALSE    TRUE
# logical   15907   15576


summary(rowData(spe_pseudo)$high_expr_group_cluster)
# Mode   FALSE    TRUE
# logical   17386   14097

with(rowData(spe_pseudo), table(high_expr_group_sample_id, high_expr_group_cluster))
#high_expr_group_sample_id FALSE  TRUE
#                    FALSE 14752     0
#                    TRUE   1479 14097

## Now filter
spe_pseudo <- spe_pseudo[rowData(spe_pseudo)$high_expr_group_cluster, ]
dim(spe_pseudo)
#15576   409

spe_pseudo<-spe_pseudo[rowData(spe_pseudo)$gene_type=='protein_coding',]
spe_pseudo<-spe_pseudo[!duplicated(rowData(spe_pseudo)$gene_name),]


lognormd <- edgeR::cpm(edgeR::calcNormFactors(spe_pseudo), log = TRUE, prior.count = 1)
normd <- edgeR::cpm(edgeR::calcNormFactors(spe_pseudo))

colnames(normd)<-spe_pseudo$cluster
colnames(lognormd)<-spe_pseudo$cluster

write.table(normd, 'visium_aggregated_cpm.tsv', na = "NA", col.names = TRUE,
            row.names = TRUE, sep = "\t")
#run voom
# x<-voom(counts(spe_pseudo), design = 'mod', lib.size = NULL,
#      block = spe_pseudo$batch, correlation = NULL, weights = NULL,
#      span = 0.5, plot = FALSE, save.plot = FALSE)
## Store the log normalized counts on the spe object
x <- edgeR::cpm(edgeR::calcNormFactors(spe_pseudo), log = TRUE, prior.count = 1)
#
## Verify that the gene order hasn't changed
stopifnot(identical(rownames(x), rownames(spe_pseudo)))
#
## Fix the column names. DGEList will have samples name as Sample1 Sample2 etc
dimnames(x) <- dimnames(spe_pseudo)
#
## Store the log normalized counts on the SingleCellExperiment object
logcounts(spe_pseudo) <- x
