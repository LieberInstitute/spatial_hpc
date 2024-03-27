is_top_10_percent <- function(column) {
  top_10_percent_threshold <- quantile(column, 0.9)
  as.numeric(column >= top_10_percent_threshold)
}

dat <- read.table("snRNAseq_aggregated_cpm.tsv",header=T)

dat.norm <- t(apply(dat, 1, function(x){x/sum(x)}))
res <- apply(dat.norm, 2, is_top_10_percent)
rownames(res) <- rownames(dat.norm)
write.csv(res,"snRNA_score.csv")

