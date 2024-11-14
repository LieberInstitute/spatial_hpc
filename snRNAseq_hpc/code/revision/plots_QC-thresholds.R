library(scater)
library(dplyr)
library(ggplot2)
library(here)

## redo QC plots
## https://github.com/LieberInstitute/spatial_hpc/blob/main/snRNAseq_hpc/code/build_sce/droplet_qc.R
load(here("snRNAseq_hpc","processed-data","sce","sce_drops_removed.rda"))

#add mdata
sampleNum = c(1:2,10:27,32,33,36:39)
sampleNames = paste0(sampleNum, "c-scp")
tmp <- read.csv(here("snRNAseq_hpc","raw-data","sample_info",
                     "snRNAseq_U01_HPC_AllRounds_Master_Spreadsheet_04072023.csv"))
tmp = tmp[,1:5]
colnames(tmp) = c("sample","tissue","brnum","round","sorted")
mdata = filter(tmp, sample %in% paste0(sampleNum,"c_scp"))

donorinfo <- read.csv(here("raw-data","sample_info","demographicInfo_Geo.csv"))
mdata = left_join(mdata, donorinfo, by=c("brnum"))
mdata$sample = sampleNames

mdata = left_join(as.data.frame(colData(sce)), mdata, by=c("Sample"="sample"))
identical(sce$key, mdata$key)
colData(sce) <- cbind(colData(sce), mdata[,4:13])

#add qc metrics
is.mt = grep("^mt-|^MT-", rowData(sce)$Symbol)
sce <- scuttle::addPerCellQC(
  sce,
  subsets = list(Mito = rownames(sce)[is.mt])#which(seqnames(sce) == "chrM")),
)

sce$high_mito <- isOutlier(sce$subsets_Mito_percent, nmads = 3, type = "higher", batch = sce$Sample)
sce$low_lib <- isOutlier(sce$sum, log = TRUE, type = "lower", batch = sce$Sample)
sce$low_genes <- isOutlier(sce$detected, log = TRUE, type = "lower", batch = sce$Sample)

df = cbind.data.frame("low_lib"=attributes(sce$low_lib)$thresholds[1,],
                      "low_genes"=attributes(sce$low_genes)$thresholds[1,],
                      "high_mito"=attributes(sce$high_mito)$thresholds[2,]) %>%
  tibble::rownames_to_column(var="Sample") %>%
  left_join(distinct(as.data.frame(colData(sce)[,c("Sample","round","brnum","sorted")])))

### plot thresholds
ggplot(df, aes(x=sorted, y=low_lib))+
  geom_boxplot()+
  geom_point(aes(color=as.factor(round)), size=3)+
  scale_color_brewer(palette="Dark2")+
  labs(title="Min. library size", y="# UMI", x="", color="seq. rnd")+
  theme_bw()+theme(text=element_text(size=16))

ggplot(df, aes(x=sorted, y=low_genes))+
  geom_boxplot()+
  geom_point(aes(color=as.factor(round)), size=3)+
  scale_color_brewer(palette="Dark2")+
  labs(title="Min. detected genes", y="# genes", x="", color="seq. rnd")+
  theme_bw()+theme(text=element_text(size=16))

ggplot(df, aes(x=sorted, y=high_mito))+
  geom_boxplot()+
  geom_point(aes(color=as.factor(round)), size=3)+
  scale_color_brewer(palette="Dark2")+
  labs(title="Max. percent chrM genes", y="% reads", x="", color="seq. rnd")+
  theme_bw()+theme(text=element_text(size=16))

### number of outliers removed
ggplot(group_by(as.data.frame(colData(sce)), Sample, round, brnum, sorted, .drop=FALSE) %>% 
         summarise(rm_mito=sum(high_mito)), 
       aes(x=sorted, y=rm_mito))+
  geom_boxplot()+
  geom_point(aes(color=as.factor(round)), size=3)+
  scale_color_brewer(palette="Dark2")+
  labs(title="chrM outliers", y="# nuclei", x="", color="seq. rnd")+
  theme_bw()+theme(text=element_text(size=16))

ggplot(group_by(as.data.frame(colData(sce)), Sample, round, brnum, sorted, .drop=FALSE) %>% 
         summarise(rm_lib=sum(low_lib)), 
       aes(x=sorted, y=rm_lib))+
  geom_boxplot()+
  geom_point(aes(color=as.factor(round)), size=3)+
  scale_color_brewer(palette="Dark2")+
  labs(title="Library size outliers", y="# nuclei", x="", color="seq. rnd")+
  theme_bw()+theme(text=element_text(size=16))

ggplot(group_by(as.data.frame(colData(sce)), Sample, round, brnum, sorted, .drop=FALSE) %>% 
         summarise(rm_genes=sum(low_genes)), 
       aes(x=sorted, y=rm_genes))+
  geom_boxplot()+
  geom_point(aes(color=as.factor(round)), size=3)+
  scale_color_brewer(palette="Dark2")+
  labs(title="Detected genes outliers", y="# nuclei", x="", color="seq. rnd")+
  theme_bw()+theme(text=element_text(size=16))

# revised thresholds
good_samples<-colnames(attributes(sce$high_mito)$thresholds[,attributes(sce$high_mito)$thresholds[2,] < 5])
sce$strict_mito <- isOutlier(sce$subsets_Mito_percent,
                          type="higher", batch=sce$Sample,
                          subset=sce$Sample %in% good_samples)

x<-c('11c-scp','13c-scp','15c-scp','18c-scp','19c-scp','21c-scp','23c-scp','2c-scp','33c-scp','39c-scp')
y<-c('11c-scp','15c-scp','19c-scp','23c-scp','33c-scp','39c-scp')

sce$strict_genes <- ifelse(sce$Sample %in% x, sce$low_genes,
                           ifelse(sce$detected <=1000,T,F))

sce$strict_lib <- ifelse(sce$Sample %in% y, sce$low_lib,
                         ifelse(sce$sum <=1000,T,F))

df2 = cbind.data.frame("low_lib"=attributes(sce$low_lib)$thresholds[1,],
                      "low_genes"=attributes(sce$low_genes)$thresholds[1,],
                      "strict_mito"=attributes(sce$strict_mito)$thresholds[2,]) %>%
  tibble::rownames_to_column(var="Sample") %>%
  left_join(distinct(as.data.frame(colData(sce)[,c("Sample","round","brnum","sorted")]))) %>%
  left_join(filter(as.data.frame(colData(sce)), strict_lib==FALSE) %>% group_by(Sample) %>% summarise(strict_lib=min(sum))) %>%
  left_join(filter(as.data.frame(colData(sce)), strict_genes==FALSE) %>% group_by(Sample) %>% summarise(strict_genes=min(detected)))

ggplot(df2, aes(x=sorted, y=strict_mito))+
  geom_boxplot()+
  geom_point(aes(color=as.factor(round)), size=3)+ylim(0,4)+
  scale_color_brewer(palette="Dark2")+
  labs(title="Max. percent chrM genes", y="% reads", x="", color="seq. rnd")+
  theme_bw()+theme(text=element_text(size=16))


ggplot(df2, aes(x=sorted, y=strict_lib))+
  geom_boxplot()+
  geom_point(aes(color=as.factor(round)), size=3)+ylim(0,2000)+
  scale_color_brewer(palette="Dark2")+
  labs(title="Min. library size", y="# UMI", x="", color="seq. rnd")+
  theme_bw()+theme(text=element_text(size=16))

ggplot(df2, aes(x=sorted, y=strict_genes))+
  geom_boxplot()+
  geom_point(aes(color=as.factor(round)), size=3)+ylim(0,2250)+
  scale_color_brewer(palette="Dark2")+
  labs(title="Min. detected genes", y="# genes", x="", color="seq. rnd")+
  theme_bw()+theme(text=element_text(size=16))

### number of outliers removed
ggplot(group_by(as.data.frame(colData(sce)), Sample, round, brnum, sorted, .drop=FALSE) %>% 
         summarise(rm_mito=sum(strict_mito)), 
       aes(x=sorted, y=rm_mito))+
  geom_boxplot()+
  geom_point(aes(color=as.factor(round)), size=3)+ylim(0,5000)+
  scale_color_brewer(palette="Dark2")+
  labs(title="chrM outliers", y="# nuclei", x="", color="seq. rnd")+
  theme_bw()+theme(text=element_text(size=16))

ggplot(group_by(as.data.frame(colData(sce)), Sample, round, brnum, sorted, .drop=FALSE) %>% 
         summarise(rm_lib=sum(strict_lib)), 
       aes(x=sorted, y=rm_lib))+
  geom_boxplot()+
  geom_point(aes(color=as.factor(round)), size=3)+ylim(0,2250)+
  scale_color_brewer(palette="Dark2")+
  labs(title="Library size outliers", y="# nuclei", x="", color="seq. rnd")+
  theme_bw()+theme(text=element_text(size=16))

ggplot(group_by(as.data.frame(colData(sce)), Sample, round, brnum, sorted, .drop=FALSE) %>% 
         summarise(rm_genes=sum(strict_genes)), 
       aes(x=sorted, y=rm_genes))+
  geom_boxplot()+
  geom_point(aes(color=as.factor(round)), size=3)+ylim(0,3300)+
  scale_color_brewer(palette="Dark2")+
  labs(title="Detected genes outliers", y="# nuclei", x="", color="seq. rnd")+
  theme_bw()+theme(text=element_text(size=16))

# final QC filters
sce$check = sce$strict_mito | sce$strict_genes | sce$strict_lib
table(sce$check)
#FALSE  TRUE 
#86905 54875
## this is the number reported in the methods

ggplot(as.data.frame(colData(sce)) %>% mutate(check=factor(check, levels=c(TRUE, FALSE), labels=c("removed","kept")),
                                              brnum2= paste(round, brnum)),
       aes(x=brnum2, y=subsets_Mito_percent, fill=as.factor(round)))+
  geom_violin()+
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1),
                     breaks=c(.1,1,10,50), labels=format(c(.1,1,10,50), scientific=F))+
  scale_fill_brewer(palette="Dark2")+
  facet_grid(rows=vars(sorted), cols=vars(check))+
  labs(y="% reads", fill="seq. rnd", title="QC summary: Mitochondrial fraction")+
  theme_bw()+theme(text=element_text(size=16), axis.text.x=element_text(angle=90, hjust=1))

ggplot(as.data.frame(colData(sce)) %>% mutate(check=factor(check, levels=c(TRUE, FALSE), labels=c("removed","kept")),
                                             brnum2= paste(round, brnum)),
       aes(x=brnum2, y=sum, fill=as.factor(round)))+
  geom_violin()+scale_y_log10()+
  scale_fill_brewer(palette="Dark2")+
  facet_grid(rows=vars(sorted), cols=vars(check))+
  labs(y="UMI counts", fill="seq. rnd", title="QC summary: Library size")+
  theme_bw()+theme(text=element_text(size=16), axis.text.x=element_text(angle=90, hjust=1))

ggplot(as.data.frame(colData(sce)) %>% mutate(check=factor(check, levels=c(TRUE, FALSE), labels=c("removed","kept")),
                                              brnum2= paste(round, brnum)),
       aes(x=brnum2, y=detected, fill=as.factor(round)))+
  geom_violin()+scale_y_log10()+
  scale_fill_brewer(palette="Dark2")+
  facet_grid(rows=vars(sorted), cols=vars(check))+
  labs(y="# genes", fill="seq. rnd", title="QC summary: Genes detected")+
  theme_bw()+theme(text=element_text(size=16), axis.text.x=element_text(angle=90, hjust=1))


#final removal of 17c-scp with <5k detected genes
sce$discard<-ifelse(sce$Sample %in% "17c-scp" & sce$detected < 5000,T,F)
sce$check2 = sce$check | sce$discard
table(sce$check2)
#FALSE  TRUE 
#80594 61186
## thank god this is the number present here: https://github.com/LieberInstitute/spatial_hpc/blob/main/snRNAseq_hpc/code/build_sce/feature_selection_dimred_clustering.R

ggplot(as.data.frame(colData(sce)) %>% 
         mutate(check2=factor(check2, levels=c(TRUE, FALSE), labels=c("removed","kept")),
                brnum2= paste(round, brnum)), 
       aes(x=brnum2, fill=as.factor(round)))+
  geom_bar(stat="count", position="dodge")+
  scale_fill_brewer(palette="Dark2")+
  facet_grid(cols=vars(check2), rows=vars(sorted))+
  labs(fill="seq. rnd", y="# nuclei", title="QC summary: nuclei kept")+
  theme_bw()+theme(text=element_text(size=16), axis.text.x=element_text(angle=90, hjust=1))


ggplot(as.data.frame(colData(sce)) %>% 
         mutate(check2=factor(check2, levels=c(TRUE, FALSE), labels=c("removed","kept")),
                brnum2= paste(round, brnum)) %>%
         group_by(round, brnum2, sorted) %>% add_tally(name="total") %>% 
         count(round, brnum2, sorted, total, check2) %>% 
         mutate(prop=n/total), 
       aes(x=brnum2, y=prop, fill=as.factor(round)))+
  geom_bar(stat="identity")+ylim(0,1)+
  scale_fill_brewer(palette="Dark2")+
  facet_grid(rows=vars(sorted), cols=vars(check2))+
  labs(fill="seq. rnd", y="prop. of sample", title="QC summary: nuclei kept")+
  theme_bw()+theme(text=element_text(size=16), axis.text.x=element_text(angle=90, hjust=1))

#########
### compare to adata

sce_qc = sce[,sce$check2==FALSE]
sce_qc$key = paste(sce_qc$Sample, sce_qc$Barcode, sep="-")
dim(sce_qc)
# 36601 80594
adata_qc = read.csv("snRNAseq_hpc/python_analysis/processed-data/adata_qc-strict_scvi-2k-hdg_obs.csv", row.names=1)
dim(adata_qc)
# 81838    34
length(intersect(sce_qc$key, rownames(adata_qc)))
# 68298
length(union(sce_qc$key, rownames(adata_qc)))
# 94134
