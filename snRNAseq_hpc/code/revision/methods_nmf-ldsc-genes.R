#load nmfs
load(file=here::here('snRNAseq_hpc','processed-data','NMF','nmf_final.rda'))
loads<-x@w
no_expr <- which(rowSums(loads) == 0)
loads <- loads[-no_expr, ]

#load old score matrix
ldsc.score <- as.matrix(read.csv("code/enrichment_analysis/project_all/score/nmf_score.csv", 
                                 row.names=1))

#### new nmf_score_top-rank.csv
identical(rownames(loads), rownames(ldsc.score))
identical(colnames(loads), colnames(ldsc.score))
topN.mat = matrix(0, nrow=nrow(ldsc.score), ncol=ncol(ldsc.score))

rownames(topN.mat) <- rownames(ldsc.score)
colnames(topN.mat) <- colnames(ldsc.score)

for(i in colnames(ldsc.score)) { 
  tmp = loads[,i]
  rank1 = (1+length(tmp))-rank(tmp)
  tmp.bin = ifelse(rank1<=1500, 1, 0)
  stopifnot(identical(names(tmp.bin), rownames(topN.mat)))
  topN.mat[,i] = tmp.bin
}

write.csv(topN.mat, "code/enrichment_analysis/nmf_score_toprank.csv", row.names=T)


