library(scater)
library(scran)
library(dplyr)
library(ggplot2)
library(here)

load(here::here("snRNAseq_hpc","processed-data", "sce", "sce_post_qc.rda"))

set.seed(100)
clust <- quickCluster(sce)

#test clust to see if can justify the increased filter of sample 17-scp in another way
df = mutate(as.data.frame(colData(sce)), quickCluster=clust,
            syt1=counts(sce)["ENSG00000067715",])

ggplot(df %>% mutate(is_nrn= syt1>1) %>% group_by(quickCluster) %>% add_tally(name="total") %>%
         group_by(quickCluster, is_nrn, total) %>%
         tally() %>% filter(is_nrn==TRUE) %>% mutate(prop=n/total), 
       aes(x=quickCluster, y=prop))+
  geom_bar(stat="identity", fill="black")+
  labs(title="is_nrn = raw SYT1 > 1", y="prop. of nuclei")+
  theme_bw()

nrn_clus = c(1,2,4,10:12,14:16,18,22)

ggplot(df %>% filter(quickCluster %in% nrn_clus), 
       aes(x=quickCluster, y=sum))+
  geom_violin()+scale_y_log10()+
  labs(title="nrn clusters only")+
  theme_bw()

ggplot(df %>% filter(quickCluster %in% nrn_clus), 
       aes(x=quickCluster, y=detected))+
  geom_violin()+scale_y_log10()+
  labs(title="nrn clusters only")+
  theme_bw()

ggplot(df %>% filter(quickCluster %in% nrn_clus) %>%
         mutate(is_17= Sample=="17c-scp"), 
       aes(x=quickCluster, y=sum, fill=is_17))+
  geom_violin(position=position_dodge(width=.5))+scale_y_log10()+
  scale_fill_manual(values=c("grey","tomato3"))+
  labs(title="nrn clusters only")+
  theme_bw()

ggplot(df %>% filter(quickCluster %in% nrn_clus) %>%
         mutate(is_17= Sample=="17c-scp"), 
       aes(x=quickCluster, y=detected, fill=is_17))+
  geom_violin(position=position_dodge(width=.5))+
  scale_fill_manual(values=c("grey","tomato3"))+
  geom_hline(aes(yintercept=5000), lty=2)+
  labs(title="nrn clusters only")+
  theme_bw()

ggplot(filter(df, Sample=="17c-scp") %>% mutate(exclude= detected<5000),
       aes(x=exclude, y=subsets_Mito_percent))+
  geom_boxplot()+theme_bw()

df2 = mutate(df, discard = ifelse(Sample %in% "17c-scp" & detected < 5000,T,F)) %>%
  filter(discard==F)

ggplot(df2 %>% filter(quickCluster %in% nrn_clus) %>%
         mutate(is_17= Sample=="17c-scp"), 
       aes(x=quickCluster, y=sum, fill=is_17))+
  geom_violin(position=position_dodge(width=.5))+scale_y_log10()+
  scale_fill_manual(values=c("grey","tomato3"))+
  labs(title="nrn clusters only - low quality 17c-scp removed")+
  theme_bw()

ggplot(df2 %>% filter(quickCluster %in% nrn_clus) %>%
         mutate(is_17= Sample=="17c-scp"), 
       aes(x=quickCluster, y=detected, fill=is_17))+
  geom_violin(position=position_dodge(width=.5))+
  scale_fill_manual(values=c("grey","tomato3"))+
  geom_hline(aes(yintercept=5000), lty=2)+
  labs(title="nrn clusters only")+
  theme_bw()
