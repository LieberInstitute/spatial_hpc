#load nmfs
load(file=here::here('snRNAseq_hpc','processed-data','NMF','nmf_final.rda'))
loads<-x@w
no_expr <- which(rowSums(loads) == 0)
loads <- loads[-no_expr, ]

top10 = lapply(colnames(loads), function(x) {
  tmp = sort(loads[,x], decreasing=T)[1:10]
  cbind.data.frame(gene=names(tmp), weights=tmp, nmf=rep(x,10))
})

top10.df = do.call(rbind, top10)
write.csv(top10.df, "snRNAseq_hpc/processed-data/revision/nmf_top-10-genes.csv", row.names=F)
