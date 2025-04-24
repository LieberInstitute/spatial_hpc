library(SpatialExperiment)
library(dplyr)
library(ggplot2)
library(ggspavis)

load(here::here("processed-data", "04_QC", "spe_QC_allSamples.Rdata"))

keep_brnum = table(spe$brnum)[11:12]
spe_spg = spe[,spe$brnum %in% names(keep_brnum)]
dim(spe_spg)
#31483 38053

df = cbind.data.frame("low_lib"=attributes(spe$low_sum_id)$thresholds[1,],
                      "low_genes"=attributes(spe$low_detected_id)$thresholds[1,]) %>%
  tibble::rownames_to_column(var="sample_id") %>%
  left_join(distinct(as.data.frame(colData(spe_spg)[,c("sample_id","brnum")]))) %>%
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
p1 <- ggplot(as.data.frame(colData(spe_spg)), aes(x=brnum, y=sum, color=low_sum_id))+
  geom_boxplot(outlier.size=.3, linewidth=.3)+scale_color_manual(values=c("black","red3"))+
  scale_y_log10()+scale_x_discrete(labels=c("Br3942\n(VSPG)","Br8325\n(VSPG)"))+
  labs(x="donor", y="# UMI", color="spot\nremoved",
       title="Library size")+
  theme_bw()+theme(text=element_text(size=8), legend.title=element_text(size=7), legend.text=element_text(size=6),
                   legend.box.spacing= unit(2, "pt"), legend.margin=margin(0,0,0,2, "pt"),
                   legend.key.size=unit(8,"pt"), legend.position="bottom")


p2 <- ggplot(as.data.frame(colData(spe_spg)), aes(x=brnum, y=detected, color=low_detected_id))+
  geom_boxplot(outlier.size=.3, linewidth=.3)+scale_color_manual(values=c("black","red3"))+
  scale_y_log10()+scale_x_discrete(labels=c("Br3942\n(VSPG)","Br8325\n(VSPG)"))+
  labs(x="donor", y="# genes", 
       title="Detected genes", color="spot\nremoved")+
  theme_bw()+theme(text=element_text(size=8), legend.title=element_text(size=7), legend.text=element_text(size=6),
                   legend.box.spacing= unit(2, "pt"), legend.margin=margin(0,0,0,2, "pt"),
                   legend.key.size=unit(8,"pt"), legend.position="bottom")

### spots kept

p3 <- ggplot(as.data.frame(colData(spe_spg)) %>%
               mutate(check2=factor(discard_auto_id, levels=c(TRUE, FALSE), labels=c("removed","kept"))),
             aes(x=brnum, group=as.factor(sample_id)))+
  geom_bar(stat="count", position="dodge", fill="grey", color="grey20", linewidth=.3)+
  scale_y_log10()+scale_x_discrete(labels=c("Br3942\n(VSPG)","Br8325\n(VSPG)"))+
  facet_grid(rows=vars(check2))+
  labs(x="donor", y="# spots", title="QC summary")+
  theme_bw()+theme(text=element_text(size=8))

pdf(file="plots/revision/supp_Vis-SPG-qc.pdf", width=7, height=3)
gridExtra::grid.arrange(p1, p2, p3, ncol=3)
dev.off()


### spot plots
spe_spg$sample_id = factor(spe_spg$sample_id, levels=c("V12D07-332_A1","V12D07-332_D1","V12D07-332_C1","V12D07-332_B1",#br3942
                                                 "V12D07-335_B1","V12D07-335_A1","V12D07-335_C1","V12D07-335_D1"#br8325
                                                 ))

sp1 = plotSpots(spe_spg[,spe_spg$brnum=="Br3942_VSPG"], annotate="discard_auto_id", sample_id="sample_id",
                point_size=.2)+
  scale_color_manual(values=c("grey","red3"))+theme_void()+
  labs(title="Br3942 (VSPG)", color="spot\nremoved")+
  theme(strip.background=element_blank(), strip.text=element_blank(),
        plot.title.position = "panel", plot.title= element_text(hjust=.5))
sp2 = plotSpots(spe_spg[,spe_spg$brnum=="Br8325_VSPG"], annotate="discard_auto_id", sample_id="sample_id",
                point_size=.2)+
  scale_color_manual(values=c("grey","red3"))+theme_void()+
  labs(title="Br8325 (VSPG)", color="spot\nremoved")+
  theme(strip.background=element_blank(), strip.text=element_blank(),
        plot.title.position = "panel", plot.title= element_text(hjust=.5))

pdf(file="plots/revision/supp_Vis-SPG-qc_spot-plots.pdf", width=6, height=8)
gridExtra::grid.arrange(sp1, sp2, ncol=1)
dev.off()
