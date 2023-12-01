setwd('/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/')
suppressPackageStartupMessages(library("SpatialExperiment"))
suppressPackageStartupMessages(library("scuttle"))
suppressPackageStartupMessages(library("scran"))
suppressPackageStartupMessages(library("scater"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("here"))
suppressPackageStartupMessages(library("sessioninfo"))
load("C:/Users/enelso40/spatial_hpc/sce_with_nmf110.rda")

## Pseudo-bulk for PRECAST_k16
sce_pseudo <- aggregateAcrossCells(
    sce,
    DataFrame(
        cluster = colData(sce)$fine.type
    ))

no_expr <- which(rowSums(counts(sce)) == 0)
length(no_expr)

length(no_expr) / nrow(counts(sce)) * 100

sce <- sce[-no_expr, ]

#find a good expression cutoff using edgeR::filterByExpr
rowData(sce_pseudo)$high_expr_group_sample_id <- filterByExpr(sce_pseudo, group = sce_pseudo$sample_id)
rowData(sce_pseudo)$high_expr_group_cluster <- filterByExpr(sce_pseudo, group = sce_pseudo$cluster)

summary(rowData(sce_pseudo)$high_expr_group_sample_id)
# Mode   FALSE    TRUE
# logical   15907   15576


summary(rowData(sce_pseudo)$high_expr_group_cluster)
# Mode   FALSE    TRUE
# logical   17386   14097

with(rowData(sce_pseudo), table(high_expr_group_sample_id, high_expr_group_cluster))
#high_expr_group_sample_id FALSE  TRUE
#                    FALSE 14752     0
#                    TRUE   1479 14097

## Now filter
sce_pseudo <- sce_pseudo[rowData(sce_pseudo)$high_expr_group_cluster, ]
dim(sce_pseudo)

sce_pseudo<-sce_pseudo[rowData(sce_pseudo)$gene_type=='protein_coding',]
sce_pseudo<-sce_pseudo[!duplicated(rowData(sce_pseudo)$gene_name),]

lognormd <- edgeR::cpm(edgeR::calcNormFactors(sce_pseudo), log = TRUE, prior.count = 1)
normd <- edgeR::cpm(edgeR::calcNormFactors(sce_pseudo))

colnames(normd)<-sce_pseudo$fine.type
colnames(lognormd)<-sce_pseudo$fine.type

write.table(normd, 'snRNAseq_aggregated_cpm.tsv', na = "NA", col.names = TRUE,
            row.names = TRUE, sep = "\t")

c

