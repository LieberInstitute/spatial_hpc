library(SpatialExperiment)
library(dplyr)
library(ggplot2)
library(here)

load(here::here("processed-data", "04_QC", "spe_QC_allSamples.Rdata"))

keep_brnum = table(spe$brnum)[1:10]
spe_he = spe[,spe$brnum %in% names(keep_brnum)]

df = cbind.data.frame("low_lib"=attributes(spe_he$low_sum_id)$thresholds[1,],
                      "low_genes"=attributes(spe_he$low_detected_id)$thresholds[1,]) %>%
  tibble::rownames_to_column(var="sample_id") %>%
  left_join(distinct(as.data.frame(colData(spe_he)[,c("sample_id","brnum")]))) %>%
  filter(!is.na(brnum)) #because attributes are factor stored the VSPG thresholds are still here

### plot thresholds
ggplot(df, aes(x=brnum, y=low_lib))+
  geom_boxplot()+
  #geom_point(aes(color=as.factor(round)), size=3)+
  #scale_color_brewer(palette="Dark2")+
  labs(title="Min. library size", y="# UMI", x="")+
  theme_bw()+theme(text=element_text(size=16), axis.text.x=element_text(angle=90, hjust=1))

ggplot(df, aes(x=sorted, y=low_genes))+
  geom_boxplot()+
  geom_point(aes(color=as.factor(round)), size=3)+
  scale_color_brewer(palette="Dark2")+
  labs(title="Min. detected genes", y="# genes", x="", color="seq. rnd")+
  theme_bw()+theme(text=element_text(size=16))


### plot spot values

p1 <- ggplot(as.data.frame(colData(spe_he)), aes(x=brnum, y=sum, color=low_sum_id))+
  geom_boxplot(outlier.size=.5)+scale_color_manual(values=c("black","red3"))+
  scale_y_log10()+
  labs(x="donor", y="# UMI", color="spot\nremoved",
       title="QC summary: library size")+
  theme_bw()+theme(text=element_text(size=12), legend.title=element_text(size=10), legend.text=element_text(size=8),
                   legend.box.spacing= unit(2, "pt"), legend.margin=margin(0,0,0,2, "pt"))

p2 <- ggplot(as.data.frame(colData(spe_he)), aes(x=brnum, y=detected, color=low_detected_id))+
  geom_boxplot(outlier.size=.5)+scale_color_manual(values=c("black","red3"))+
  scale_y_log10()+
  labs(x="donor", y="# genes", color="outlier",
       title="QC summary: detected genes", color="spot\nremoved")+
  theme_bw()+theme(text=element_text(size=12), legend.title=element_text(size=10), legend.text=element_text(size=8),
                   legend.box.spacing= unit(2, "pt"), legend.margin=margin(0,0,0,2, "pt"))

### spots kept

p3 <- ggplot(as.data.frame(colData(spe_he)) %>% 
         mutate(check2=factor(discard_auto_id, levels=c(TRUE, FALSE), labels=c("removed","kept"))), 
       aes(x=brnum, group=as.factor(sample_id)))+
  geom_bar(stat="count", position=position_dodge2(width=.9, preserve="single"), 
           fill="grey", color="grey20", )+
  scale_y_log10()+
  facet_grid(rows=vars(check2))+
  labs(x="donor", y="# spots", title="QC summary: spots kept")+
  theme_bw()+theme(text=element_text(size=12))

pdf(file="plots/revision/supp_Vis-HE-qc.pdf", width=7, height=9)
gridExtra::grid.arrange(p1, p2, p3, ncol=1)
dev.off()

