#supplemental tables
library(dplyr)
set.seed(123)

#broad.domain
load("processed-data/08_pseudobulk/PRECAST/visiumHE_DE_stats_broad_final.rda")

broad.reg = stats$enrichment
colnames(broad.reg)
tstats.df2 = broad.reg[,c(grep("t_stat",colnames(broad.reg)),17:18)] %>% 
  tidyr::pivot_longer(cols=grep("t_stat",colnames(broad.reg), value=T), names_to="broad.domain", values_to="t_stat", names_prefix="t_stat_")
fdr.df2 = broad.reg[,c(grep("fdr",colnames(broad.reg)),17:18)] %>%
  tidyr::pivot_longer(cols=grep("fdr",colnames(broad.reg), value=T), names_to="broad.domain", values_to="fdr", names_prefix="fdr_")
logfc.df2 = broad.reg[,c(grep("logFC",colnames(broad.reg)),17:18)] %>%
  tidyr::pivot_longer(cols=grep("logFC",colnames(broad.reg), value=T), names_to="broad.domain", values_to="logfc", names_prefix="logFC_")

enrich.df2 = left_join(tstats.df2, fdr.df2, by=c("broad.domain","ensembl","gene")) %>%
  left_join(logfc.df2, by=c("broad.domain","ensembl","gene"))

up.df2 = filter(enrich.df2, fdr<.01 & abs(logfc)>1)
colnames(up.df2)
#[1] "ensembl"      "gene"         "broad.domain" "t_stat"       "fdr"         
#[6] "logfc"
dim(up.df2) # 7554    6
write.csv(up.df2, "processed-data/revision/broad-domain_pseudobulk_DEG_fdr-01-logfc-abs-1.csv", row.names=F)

#spatial domain
load("processed-data/08_pseudobulk/PRECAST/visiumHE_DE_stats_domain.rda")

domain.reg = stats$enrichment
colnames(domain.reg)
tstats.df3 = domain.reg[,c(grep("t_stat",colnames(domain.reg)),65:66)] %>% 
  tidyr::pivot_longer(cols=grep("t_stat",colnames(domain.reg), value=T), names_to="domain", values_to="t_stat", names_prefix="t_stat_")
fdr.df3 = domain.reg[,c(grep("fdr",colnames(domain.reg)),65:66)] %>%
  tidyr::pivot_longer(cols=grep("fdr",colnames(domain.reg), value=T), names_to="domain", values_to="fdr", names_prefix="fdr_")
logfc.df3 = domain.reg[,c(grep("logFC",colnames(domain.reg)),65:66)] %>%
  tidyr::pivot_longer(cols=grep("logFC",colnames(domain.reg), value=T), names_to="domain", values_to="logfc", names_prefix="logFC_")

enrich.df3 = left_join(tstats.df3, fdr.df3, by=c("domain","ensembl","gene")) %>%
  left_join(logfc.df3, by=c("domain","ensembl","gene"))

up.df3 = filter(enrich.df3, fdr<.01 & abs(logfc)>1)
colnames(up.df3)
#"ensembl" "gene"    "domain"  "t_stat"  "fdr"     "logfc" 
dim(up.df3) #21416     6
write.csv(up.df3, "processed-data/revision/domain_pseudobulk_DEG_fdr-01-logfc-abs-1.csv", row.names=F)

