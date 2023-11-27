set.seed(1029)
i<-intersect(rownames(spe),rownames(x@w))
loadings<-x@w
loadings<-loadings[rownames(loadings) %in% i,]
spe2<-spe[rownames(spe) %in% i,]
loadings<-loadings[match(rownames(spe2),rownames(loadings)),]

proj<-project(loadings,logcounts(spe2),L2=0.00005)


sce_hex <- make_hexbin(sce, nbins = 100,
                           dimension_reduction = "UMAP", use_dims=c(1,2))
plot_hexbin_meta(sce_hex, col="nmf20", action="median")+
    #theme(legend.position = 'bottom')+
    #legend.text = element_text(angle=90))+
    scale_fill_distiller(
        type = "seq",
        palette = rev('Greys'),
        direction=1,
        limits=c(0,2e-04))+labs(x='UMAP1',y='UMAP2',title='nmf20 (growth factor/MAPK/ERK signaling)')

plot_hexbin_meta(sce_hex, col="nmf77", action="median")+
   # theme(legend.position = 'bottom')+
          #legend.text = element_text(angle=90))+
    scale_fill_distiller(
        type = "seq",
        palette = rev('Greys'),
        direction=1,
        limits=c(0,2e-04))+labs(x='UMAP1',y='UMAP2',title='nmf77 (oligodendrocytes)')

pdf('nmf77_umap.pdf',h=5,w=4)
plotUMAP(sce,colour_by='nmf77',point_size=0.001)+
    theme(legend.position = 'bottom',
          legend.text = element_text(angle=90))+
    scale_color_distiller(
        type = "seq",
        palette = rev('Greys'),
        direction=1+
        limits=c(0,3e-04))
dev.off()

