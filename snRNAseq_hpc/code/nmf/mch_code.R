library(biomaRt)
library(zellkonverter)
library(SingleCellExperiment)
library(RcppML)
library(here)

##load the data
mch<-readH5AD(file=here::here('snRNAseq_hpc','processed-data',
                              'NMF','rs2_mch_matrix.h5ad'))

##get gene names
mart <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
symb <- getBM(attributes = c("ensembl_gene_id","mgi_symbol"),
              filters = "ensembl_gene_id", values = rownames(mch),
              mart = mart)
symbs <- symb$mgi_symbol[match(rownames(mch), symb$ensembl_gene_id, nomatch = NA)]

rowData(mch)$gene_name<-symbs
rowData(mch)$start<-NULL
rowData(mch)$end<-NULL

##drop genes with no names (need these for ortholog matching)
mch<-mch[!is.na(rowData(mch)$gene_name),]
rownames(mch)<-rowData(mch)$gene_name


# Translate from one mchcies to the other using the orthology
orthology<-read.csv(file=here::here('snRNAseq_hpc','processed-data',
                                    'NMF','human_mouse_orthologs.csv'))
names <- orthology[orthology$Column3 %in% rownames(mch),]

names <- names[match(rownames(mch), names$Column3),]

setdiff(names$Column3, rownames(mch))

rownames(mch) <- names$Column1


############## Project human snRNAseq NMF patterns onto transformed mouse methylation data###########
##load nmf patterns
load(file=here::here('snRNAseq_hpc','processed-data','NMF','nmf_final.rda'))

###process W matrix for projection
loads<-x@w

no_expr <- which(rowSums(loads) == 0)
loads <- loads[-no_expr, ]
dim(loads)
#[1] 19363    100

set.seed(1029)
i<-intersect(rownames(mch),rownames(loads))
loadings<-loads
loadings<-loadings[rownames(loadings) %in% i,]
mch<-mch[rownames(mch) %in% i,]
#####We made this object with a lot of cortical regions that we actually don't need--let's subset
mch<-mch[,mch$Source %in% c('CAa','CAp','DGa','DGp','ENT')]
mch<-mch[,!is.na(mch$Subclass)]
mch$Source<-droplevels(mch$Source)
mch$Target<-droplevels(mch$Target)
mch$Subclass<-droplevels(mch$Subclass)
##now get rid of irrelevant cell types (just want pyramidal ones right now)
mch<-mch[,mch$Subclass %in% c("CA3 Glut", "CA2-FC-IG Glut", "CA1-ProS Glut",
                                 "SUB-ProS Glut", "NP SUB Glut", "CT SUB Glut",
                                 "L6 IT CTX Glut", "IT EP-CLA Glut", "L6b/CT ENT Glut",
                                 "L5 ET CTX Glut", "L5/6 IT TPE-ENT Glut", "L2/3 IT PPP Glut",
                                 "L2 IT PPP-APr Glut", "L2/3 IT ENT Glut", "L2 IT ENT-po Glut",
                                 "L2/3 IT PIR-ENTl Glut", "ENTmv-PA-COAp Glut", "COAp Grxcr2 Glut")]

######projection happens here######
loadings<-loadings[match(rownames(mch),rownames(loadings)),]
proj<-project(loadings,assay(mch,'X'),L1=0)
proj<-t(proj)

#######rescale patterns#########
proj<-apply(proj,2,function(x){x/sum(x)})

########add projected patterns to mch#####
colData(mch)<-cbind(colData(mch),proj)

######save mch#######
save(mch,file=here::here('snRNAseq_hpc','processed-data','NMF','mch.rda'))
