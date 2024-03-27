gene.anno.file <- "/dcs04/lieber/shared/statsgen/LDSC/base/gene_meta_hg19.txt"
gene.meta <- read.table(gene.anno.file,sep="\t",header=T)

cell <- read.csv("../score/srt_score.csv")
colnames(cell)[1] <- "geneName"
modules <- colnames(cell)[-1]
for(i in 1:length(modules)){
	idx <- which(colnames(cell) == modules[i])
	genes <- cell$geneName[cell[,idx]==1]
	genes <- intersect(genes,gene.meta$Gene.name)
	bed <- gene.meta[is.element(gene.meta$Gene.name,genes),3:5]
	bed[,1] <- paste0("chr",bed[,1])
	bed[,2] <- ifelse(bed[,2] - 100000 > 0, bed[,2] - 100000, 0)
	bed[,3] <- bed[,3] + 100000
	idx <- bed[,1] !="chrX" & bed[,1] !="chrY" & bed[,1] !="chrMT"
	bed <- bed[idx,]
	outfile <- paste0(modules[i],".bed")
	write.table(bed,outfile,row.names=F,col.names=F,sep="\t",quote=F)
}
