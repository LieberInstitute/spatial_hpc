#supplemental tables
library(dplyr)
set.seed(123)


reg = read.csv("snRNAseq_hpc/processed-data/revision/sn_enrichment_stats_superfine.csv", row.names=1)

tstats.df = reg[,c(grep("t_stat",colnames(reg)),237:238)] %>% 
  tidyr::pivot_longer(cols=grep("t_stat",colnames(reg), value=T), names_to="superfine.cell.class", values_to="t_stat", names_prefix="t_stat_")
fdr.df = reg[,c(grep("fdr",colnames(reg)),237:238)] %>%
  tidyr::pivot_longer(cols=grep("fdr",colnames(reg), value=T), names_to="superfine.cell.class", values_to="fdr", names_prefix="fdr_")
logfc.df = reg[,c(grep("logFC",colnames(reg)),237:238)] %>%
  tidyr::pivot_longer(cols=grep("logFC",colnames(reg), value=T), names_to="superfine.cell.class", values_to="logfc", names_prefix="logFC_")

enrich.df = left_join(tstats.df, fdr.df, by=c("superfine.cell.class","ensembl","gene")) %>%
  left_join(logfc.df, by=c("superfine.cell.class","ensembl","gene"))

up.df = filter(enrich.df, fdr<.0001 & logfc>2)
colnames(up.df)
#[1] "ensembl"              "gene"                 "superfine.cell.class"
#[4] "t_stat"               "fdr"                  "logfc" 
dim(up.df) #33112     6
write.csv(up.df, "snRNAseq_hpc/processed-data/revision/snRNAseq_pseudobulk_DEG_fdr-0001-logfc-2.csv", row.names=F)
