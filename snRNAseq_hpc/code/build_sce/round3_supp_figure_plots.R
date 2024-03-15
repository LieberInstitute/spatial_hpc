pdf(file=here::here('plots','figures','supp_figures','supp_figure_round_3','QC_plots_round3.pdf'),w=12,h=9)
p1<-plotColData(sce, x = "sample_id", y = "sum", colour_by = "discard_semiauto",point_alpha=1,point_size=0.1) +
    scale_y_log10() + ylab('library size')+xlab('sample')+facet_wrap(~sce$round,scales = 'free')+
    ggtitle("Total UMIs") +theme(axis.text.x = element_text(angle = 90),
                                 text=element_text(size=20),
                                 axis.text=element_text(size=20))+theme(legend.position = 'none')
print(rasterize(p1,dpi=300,layers='Point'))


p2<-plotColData(sce, x = "sample_id", y = "subsets_Mito_percent", colour_by = "discard_semiauto",point_alpha=1,point_size=0.1) +
    ylab('mitochondrial gene expression percentage')+xlab('sample')+facet_wrap(~sce$round,scales = 'free')+
    ggtitle("mitochondrial gene expression percentage") +theme(axis.text.x = element_text(angle = 90),
                                                               text=element_text(size=20),
                                                               axis.text=element_text(size=20))+theme(legend.position = 'none')
print(rasterize(p2,dpi=300,layers='Point'))

p3<-plotColData(sce, x = "sample_id", y = "detected", colour_by = "discard_semiauto",point_alpha=1,point_size=0.1) +
    ylab('detected genes')+xlab('sample')+facet_wrap(~sce$round,scales = 'free')+
    scale_y_log10()+ggtitle("Detected genes") +theme(axis.text.x = element_text(angle = 90),
                                                     text=element_text(size=20),
                                                         axis.text=element_text(size=20))+theme(legend.position = 'none')
print(rasterize(p3,dpi=300,layers='Point'))

dev.off()


##barplot showing  nuclei by sample, colored by discard status
# Convert SingleCellExperiment to dataframe
df <- as.data.frame(colData(sce))

# Calculate proportions
df <- df %>%
    group_by(discard_semiauto,sample_id, round) %>%
    summarise(count = n(), .groups = 'drop') %>%
    mutate(proportion = count / sum(count)) %>%
    ungroup() %>%
    group_by(sample_id) %>%
    mutate(total = sum(proportion)) %>%
    mutate(proportion = proportion / total)

# Create barplot
pdf(file=here::here('plots','figures','supp_figures','supp_figure_round_3','barplot_round3.pdf'),h=9,w=12)
ggplot(df, aes(x = sample_id, y = count, fill = discard_semiauto)) +
    geom_bar(stat = "identity", position = "stack") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"),
          text=element_text(colour='black',size=16),
          axis.text.x = element_text(angle = 90),
          legend.position='bottom') +
    ylab("nuclei") +
    xlab("sample") +
    labs(fill = "discard")+facet_wrap(~round,scales='free_x')+
    scale_fill_manual(values=c('orange','cornflowerblue'))
dev.off()
