library(SpatialExperiment)
library(dplyr)
library(ggplot2)
library(ggspavis)
library(scater)
library(gridExtra)

set.seed(123)

load(file=here::here('processed-data','NMF','spe_nmf_final.rda'))
load("plots/spatial_palette_final.rda")

#Br3942, all domains, CA1, ProS, Sub1, Sub2
sub1 = spe[,spe$sample_id %in% c("V11L05-333_A1","V11L05-333_C1")] #
sub1$facet_row = ifelse(sub1$sample_id=="V11L05-333_A1","r1","r2")

p1 <- plotVisium(sub1, spots=T, annotate="nmf15", highlight="domain",
           facets="sample_id", image = FALSE, point_size=.8)+coord_cartesian()+
  facet_grid(rows=vars(facet_row), scales="free")+
  scale_y_continuous(expand=c(0,0))+scale_x_continuous(expand=c(0,0))+
  scale_color_manual(values=spatial.palette, guide="none")+guides(color=NULL)+
  scale_fill_gradient(low="white",high="black", limits=c(0,max(spe$nmf15)),
                      labels=function(x) format(x, scientific=T, digits=1))+
  theme(legend.text = element_text(size = 8), 
        strip.text=element_blank(), aspect.ratio=1,
        panel.spacing=unit(5,'points'),plot.margin=unit(c(3,5,3,0), "pt"))

p2 <- plotVisium(sub1, spots=T, annotate="nmf32", highlight="domain",
           facets="sample_id", image = FALSE, point_size=.8)+coord_cartesian()+
  facet_grid(rows=vars(facet_row), scales="free")+
  scale_y_continuous(expand=c(0,0))+scale_x_continuous(expand=c(0,0))+
  scale_color_manual(values=spatial.palette, guide="none")+guides(color=NULL)+
  scale_fill_gradient(low="white",high="black", limits=c(0,max(spe$nmf32)),
                      labels=function(x) format(x, scientific=T, digits=1))+
  theme(legend.text = element_text(size = 8), 
        strip.text=element_blank(), aspect.ratio=1,
        panel.spacing=unit(5,'points'),plot.margin=unit(c(3,5,3,0), "pt"))

p3 <- plotVisium(sub1, spots=T, annotate="nmf40", highlight="domain",
           facets="sample_id", image = FALSE, point_size=.8)+coord_cartesian()+
  facet_grid(rows=vars(facet_row), scales="free")+
  scale_y_continuous(expand=c(0,0))+scale_x_continuous(expand=c(0,0))+
  scale_color_manual(values=spatial.palette, guide="none")+guides(color=NULL)+
  scale_fill_gradient(low="white",high="black", limits=c(0,max(spe$nmf40)),
                      labels=function(x) format(x, scientific=T, digits=1))+
  theme(legend.text = element_text(size = 8), 
        strip.text=element_blank(), aspect.ratio=1,
        panel.spacing=unit(5,'points'),plot.margin=unit(c(3,5,3,0), "pt"))

p4 <- plotVisium(sub1, spots=T, annotate="nmf54", highlight="domain",
           facets="sample_id", image = FALSE, point_size=.8)+coord_cartesian()+
  facet_grid(rows=vars(facet_row), scales="free")+
  scale_y_continuous(expand=c(0,0))+scale_x_continuous(expand=c(0,0))+
  scale_color_manual(values=spatial.palette, guide="none")+guides(color=NULL)+
  scale_fill_gradient(low="white",high="black", limits=c(0,max(spe$nmf54)),
                      labels=function(x) format(x, scientific=T, digits=1))+
  theme(legend.text = element_text(size = 8), 
        strip.text=element_blank(), aspect.ratio=1,
        panel.spacing=unit(5,'points'),plot.margin=unit(c(3,5,3,0), "pt"))

p5<-plotVisium(sub1, spots=T, annotate="nmf17", highlight="domain",
           facets="sample_id", image = FALSE, point_size=.8)+coord_cartesian()+
  facet_grid(rows=vars(facet_row), scales="free")+
  scale_y_continuous(expand=c(0,0))+scale_x_continuous(expand=c(0,0))+
  scale_color_manual(values=spatial.palette, guide="none")+guides(color=NULL)+
  scale_fill_gradient(low="white",high="black", limits=c(0,max(spe$nmf17)),
                      labels=function(x) format(x, scientific=T, digits=1))+
  theme(legend.text = element_text(size = 8),
        strip.text=element_blank(), aspect.ratio=1,
        panel.spacing=unit(5,'points'),plot.margin=unit(c(3,5,3,0), "pt"))

pdf(file = "plots/revision/Fig6_Br3942_2-capture-area_domain-nmf.pdf",
    width=9, height=8)
grid.arrange(p1,p2,p3,p4,p5, ncol=3)
dev.off()

#Br3942 smaller
p1<-plotVisium(sub1, spots=T, annotate="nmf15", highlight="domain",
               facets="sample_id", image = FALSE, point_size=.5)+coord_cartesian()+
  facet_grid(rows=vars(facet_row), scales="free")+
  scale_y_continuous(expand=c(0,0))+scale_x_continuous(expand=c(0,0))+
  scale_color_manual(values=spatial.palette, guide="none")+guides(color=NULL)+
  scale_fill_gradient(low="white",high="black", limits=c(0,max(spe$nmf15)),
                      labels=function(x) format(x, scientific=T, digits=1))+
  theme(legend.text = element_text(size = 6), legend.key.size = unit(10,"pt"),
        legend.title = element_text(size=8),
        strip.text=element_blank(), aspect.ratio=1,
        panel.spacing=unit(5,'points'),plot.margin=unit(c(3,5,3,0), "pt"))

p2<-plotVisium(sub1, spots=T, annotate="nmf32", highlight="domain",
               facets="sample_id", image = FALSE, point_size=.5)+coord_cartesian()+
  facet_grid(rows=vars(facet_row), scales="free")+
  scale_y_continuous(expand=c(0,0))+scale_x_continuous(expand=c(0,0))+
  scale_color_manual(values=spatial.palette, guide="none")+guides(color=NULL)+
  scale_fill_gradient(low="white",high="black", limits=c(0,max(spe$nmf32)),
                      labels=function(x) format(x, scientific=T, digits=1))+
  theme(legend.text = element_text(size = 6), legend.key.size = unit(10,"pt"),
        legend.title = element_text(size=8),
        strip.text=element_blank(), aspect.ratio=1,
        panel.spacing=unit(5,'points'),plot.margin=unit(c(3,5,3,0), "pt"))

p3<-plotVisium(sub1, spots=T, annotate="nmf40", highlight="domain",
               facets="sample_id", image = FALSE, point_size=.5)+coord_cartesian()+
  facet_grid(rows=vars(facet_row), scales="free")+
  scale_y_continuous(expand=c(0,0))+scale_x_continuous(expand=c(0,0))+
  scale_color_manual(values=spatial.palette, guide="none")+guides(color=NULL)+
  scale_fill_gradient(low="white",high="black", limits=c(0,max(spe$nmf40)),
                      labels=function(x) format(x, scientific=T, digits=1))+
  theme(legend.text = element_text(size = 6), legend.key.size = unit(10,"pt"),
        legend.title = element_text(size=8),
        strip.text=element_blank(), aspect.ratio=1,
        panel.spacing=unit(5,'points'),plot.margin=unit(c(3,5,3,0), "pt"))

p4<-plotVisium(sub1, spots=T, annotate="nmf54", highlight="domain",
               facets="sample_id", image = FALSE, point_size=.5)+coord_cartesian()+
  facet_grid(rows=vars(facet_row), scales="free")+
  scale_y_continuous(expand=c(0,0))+scale_x_continuous(expand=c(0,0))+
  scale_color_manual(values=spatial.palette, guide="none")+guides(color=NULL)+
  scale_fill_gradient(low="white",high="black", limits=c(0,max(spe$nmf54)),
                      labels=function(x) format(x, scientific=T, digits=1))+
  theme(legend.text = element_text(size = 6), legend.key.size = unit(10,"pt"),
        legend.title = element_text(size=8),
        strip.text=element_blank(), aspect.ratio=1,
        panel.spacing=unit(5,'points'),plot.margin=unit(c(3,5,3,0), "pt"))

p5<-plotVisium(sub1, spots=T, annotate="nmf17", highlight="domain",
               facets="sample_id", image = FALSE, point_size=.5)+coord_cartesian()+
  facet_grid(rows=vars(facet_row), scales="free")+
  scale_y_continuous(expand=c(0,0))+scale_x_continuous(expand=c(0,0))+
  scale_color_manual(values=spatial.palette, guide="none")+guides(color=NULL)+
  scale_fill_gradient(low="white",high="black", limits=c(0,max(spe$nmf17)),
                      labels=function(x) format(x, scientific=T, digits=1))+
  theme(legend.text = element_text(size = 6), legend.key.size = unit(10,"pt"),
        legend.title = element_text(size=8),
        strip.text=element_blank(), aspect.ratio=1,
        panel.spacing=unit(5,'points'),plot.margin=unit(c(3,5,3,0), "pt"))

pdf(file = "plots/revision/Fig6_Br3942_2-capture-area_domain-nmf_small.pdf",
    width=6, height=5.35)
grid.arrange(p1,p2,p3,p4,p5, ncol=3)
dev.off()

#3-donor single capture area plots for center
spe$sub1 = spe$nmf40>((max(spe$nmf40)*.95)/5)
spe$sub2 = spe$nmf54>((max(spe$nmf54)*.95)/5)
spe$subiculum = spe$sub1|spe$sub2

spe$broad.domain2 = factor(ifelse(spe$subiculum==TRUE, "Subiculum", as.character(spe$broad.domain)),
                           levels=c("Subiculum","Neuron","Neuropil","WM","Vasc_CSF"))
broad.palette2= c("Subiculum"= "#00a000", "Neuron"="#add294", "Neuropil"="#dfa56e", "WM"="#ff80fe", "Vasc_CSF"="#00006a")

sub2 = mirrorObject(spe[,spe$sample_id=="V11U08-081_D1"], axis="v")
#sub.special = cbind(spe[,spe$sample_id %in% c("V10B01-085_B1","V10B01-086_A1")], sub2)
#sub.special$facet_row = as.character(factor(sub.special$brnum, levels=c("Br6423","Br2743","Br6432"), labels=c("r1","r2","r3")))

sup.list = c("nmf84","nmf45","nmf27","nmf51")
plist = lapply(sup.list, function(x) {
  plotVisium(sub2, spots=T, annotate=x, highlight="broad.domain2", image = FALSE, point_size=.5)+
    scale_color_manual(values=broad.palette2, guide="none")+guides(color=NULL)+
    scale_fill_gradient(low="white",high="black", 
                        labels=function(y) format(y, scientific=T, digits=1))+
    theme(strip.text=element_blank(), legend.text = element_text(size = 6, margin=margin(0,2,0,2,"pt")),
          legend.box.spacing = unit(0,"pt"), legend.key.size=unit(8,"pt"), legend.title=element_text(size=9))
})
ggsave("plots/revision/Figure6_ENT-spot-plots.pdf", do.call(gridExtra::grid.arrange, c(plist, ncol=2)),
       bg="white", height=4, width=4, unit="in")

dp.list = c("nmf68","nmf22","nmf53","nmf65")
plist2 = lapply(dp.list, function(x) {
  plotVisium(sub2, spots=T, annotate=x, highlight="broad.domain2", image = FALSE, point_size=.5)+
    scale_color_manual(values=broad.palette2, guide="none")+guides(color=NULL)+
    scale_fill_gradient(low="white",high="black", 
                        labels=function(y) format(y, scientific=T, digits=1))+
    theme(strip.text=element_blank(), legend.text = element_text(size = 6, margin=margin(0,2,0,2,"pt")),
          legend.box.spacing = unit(0,"pt"), legend.key.size=unit(8,"pt"), legend.title=element_text(size=9))
})
ggsave("plots/revision/Figure6_RHP-spot-plots.pdf", do.call(gridExtra::grid.arrange, c(plist2, ncol=2)),
       bg="white", height=4, width=4, unit="in")


####### supp plots
#binarizing domains
spe$ca1 = spe$nmf15>((max(spe$nmf15)*.95)/5)
spe$ca3.2 = spe$nmf63>((max(spe$nmf63)*.95)/5)
spe$ca3.1 = spe$nmf11>((max(spe$nmf11)*.95)/5)
spe$ca3 = spe$ca3.1|spe$ca3.2
spe$ENT_sup3 = spe$nmf27>((max(spe$nmf27)*.95)/5)
spe$ENT_sup2 = spe$nmf45>((max(spe$nmf45)*.95)/5)
spe$ENT_sup1 = spe$nmf84>((max(spe$nmf84)*.95)/5)
spe$ENT_sup = spe$ENT_sup1|spe$ENT_sup2|spe$ENT_sup3
spe$ENT_L5 =  spe$nmf51>((max(spe$nmf51)*.95)/5)

spe$binary_label = "other"
colData(spe)[spe$subiculum==T & spe$ca1==F & spe$ca3==F & spe$ENT_sup==F & spe$ENT_L5==F,"binary_label"] = "subiculum"
colData(spe)[spe$subiculum==F & spe$ca1==T & spe$ca3==F & spe$ENT_sup==F & spe$ENT_L5==F,"binary_label"] = "CA1"
colData(spe)[spe$subiculum==F & spe$ca1==F & spe$ca3==1 & spe$ENT_sup==F & spe$ENT_L5==F,"binary_label"] = "CA3"
colData(spe)[spe$subiculum==F & spe$ca1==F & spe$ca3==F & spe$ENT_sup==T & spe$ENT_L5==F,"binary_label"] = "ENT_sup"
colData(spe)[spe$subiculum==F & spe$ca1==F & spe$ca3==F & spe$ENT_sup==F & spe$ENT_L5==T,"binary_label"] = "ENT_L5"

#binarized threshold supp
#table
n.tbl = as.data.frame(colData(spe)[,c("subiculum","ca1","ca3","ENT_sup","ENT_L5")]) %>%
  mutate_all(as.numeric) %>%
  group_by(subiculum, ca1, ca3, ENT_sup, ENT_L5) %>% tally()
totals.row = c(sum(spe$subiculum),sum(spe$ca1),sum(spe$ca3),sum(spe$ENT_sup),sum(spe$ENT_L5),sum(n.tbl$n))
names(totals.row) = colnames(n.tbl)
n.tbl = rbind(n.tbl, totals.row)
write.csv(n.tbl, "processed-data/revision/binarized-nmf_table.csv", row.names=F)

#bar
ggplot(as.data.frame(colData(spe)), aes(y=domain, fill=factor(binary_label, levels=names(binary.palette))))+
  geom_bar(stat="count", position="stack")+theme_bw()+
  scale_fill_manual(values=binary.palette)+labs(fill="", x="# spots")+
  theme(text=element_text(size=14))


#ignore below: creates violin plots that was more confusing than helpful
cdata = as.data.frame(colData(spe))

plot.df = bind_rows(mutate(filter(cdata, sub1==T), xax="Sub.1"),
                   mutate(filter(cdata, sub2==T), xax="Sub.2"),
                   mutate(filter(cdata, ca1==T), xax="CA1"),
                   mutate(filter(cdata, ca3.1==T), xax="CA3.1"),
                   mutate(filter(cdata, ca3.2==T), xax="CA3.2"),
                   mutate(filter(cdata, ENT_sup1==T), xax="ENT.sup1"),
                   mutate(filter(cdata, ENT_sup2==T), xax="ENT.sup2"),
                   mutate(filter(cdata, ENT_sup3==T), xax="ENT.sup3"),
                   mutate(filter(cdata, ENT_L5==T), xax="ENT.L5"),
                   mutate(filter(cdata, binary_label=="other"), xax="other")) %>%
  mutate(xax=factor(xax, levels=c("CA3.1","CA3.2","CA1","Sub.1","Sub.2","ENT.sup1","ENT.sup2","ENT.sup3","ENT.L5","other")))

binary.palette = c("CA1"="#984ea3", "CA3"="#fec44f",
                       "subiculum"="#e41a1c",
                       "ENT_L5"="#377eb8",
                       "ENT_sup"="#61963d",#"L2/3.PrS.Ent"="#add294",
                       "other"="white")
library(ggbeeswarm)
nmf.list = c("nmf11","nmf63","nmf15","nmf40","nmf54","nmf84","nmf45","nmf27","nmf51")
vlist <- lapply(nmf.list, function(i) {
  tmp = plot.df
  colnames(tmp)[grep(i, colnames(tmp))] <- "plotME"
  ggplot(tmp, aes(x=xax, y=plotME, color=binary_label))+
    geom_violin(scale="width", width=.8, color="black")+
    geom_quasirandom(width=.4, bandwidth = 1, size=.3)+
    scale_color_manual(values=binary.palette)+
    scale_y_continuous(labels=function(x) format(x, scientific=T, digits=1))+
    labs(title=i, y=i)+theme_bw()+theme(axis.title.x=element_blank(), legend.position="none",
                                        axis.text.x=element_text(angle=45, hjust=1)) 
}
)
do.call(gridExtra::grid.arrange, c(vlist, ncol=3))

#binary spot plots for main figure
sub1 = spe[,spe$sample_id %in% c("V11L05-333_A1","V11L05-333_C1")] #
sub1$facet_row = ifelse(sub1$sample_id=="V11L05-333_A1","r1","r2")
sub3 = spe[,spe$sample_id %in% c("V10B01-085_B1","V10B01-085_D1")]
sub3$facet_row = ifelse(sub3$sample_id=="V10B01-085_B1","r1","r2")
sub4 = spe[,spe$sample_id %in% c("V11L05-336_A1","V11L05-336_B1")]

plotVisium(sub1, spots=T, annotate="nmf17", highlight="binary_label",
           facets="sample_id", image = FALSE)+coord_cartesian()+
  facet_grid(rows=vars(facet_row), scales="free")+
  scale_y_continuous(expand=c(0,0))+scale_x_continuous(expand=c(0,0))+
  scale_color_manual(values=binary.palette, guide="none")+guides(color=NULL)+
  scale_fill_gradient(low="grey90",high="black", limits=c(0,max(spe$nmf17)),
                      labels=function(x) format(x, scientific=T, digits=1))+
  theme(legend.text = element_text(size = 8), 
        strip.text=element_blank(), aspect.ratio=1,
        panel.spacing=unit(5,'points'),plot.margin=unit(c(0,0,0,0), "pt"))

