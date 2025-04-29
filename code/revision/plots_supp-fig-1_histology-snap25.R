library(SpatialExperiment)
library(ggspavis)

#load(here("processed-data", "04_QC", "spe_QC_allSamples.Rdata"), verbose = TRUE)
load(file="processed-data/NMF/spe_nmf_final.rda")
#update sample_id order
spe$sample_id = factor(spe$sample_id, levels=c("V10B01-086_D1","V10B01-086_C1","V11U08-081_C1","V11U08-081_D1",
                                               "V11L05-333_A1","V11L05-333_B1","V11L05-333_C1","V11L05-333_D1",
                                               "V10B01-085_B1","V10B01-085_A1","V10B01-085_D1","V10B01-085_C1",
                                               "V10B01-086_A1","V10B01-086_B1",
                                               #"V11L05-335_C1","V11L05-335_B1","V11L05-335_A1",
                                               "V11L05-335_A1","V11L05-335_C1","V11L05-335_B1",
                                               "V11U08-084_A1","V11U08-084_B1","V11U08-084_C1","V11U08-084_D1",
                                               "V11A20-297_C1","V11A20-297_D1","V11A20-297_A1","V11L05-335_D1","V11A20-297_B1",
                                               "V11U08-081_A1","V11U08-081_B1",
                                               "V11L05-336_A1","V11L05-336_B1","V11L05-336_C1","V11L05-336_D1",
                                               "V12F14-051_C1","V12F14-051_D1","V12F14-051_A1","V12F14-051_B1"))
spe$facet_column = as.character(factor(spe$sample_id, levels=levels(spe$sample_id),
                                       labels=c("c1","c2","c1","c2",
                                                "c1","c2","c1","c2",
                                                "c1","c2","c1","c2",
                                                "c1","c1",
                                                #"c1","c2","c1",
                                                "c1","c2","c3",
                                                "c1","c2","c1","c2",
                                                "c1","c2","c1","c2","c2",
                                                "c1","c1",
                                                "c1","c2","c1","c2",
                                                "c1","c2","c1","c2")))
spe$facet_row = as.character(factor(spe$sample_id, levels=levels(spe$sample_id),
                                    labels=c("r1","r1","r2","r2",
                                             "r1","r1","r2","r2",
                                             "r1","r1","r2","r2",
                                             "r1","r2",
                                             #"r1","r1","r2",
                                             "r1","r1","r1",
                                             "r1","r1","r2","r2",
                                             "r1","r1","r2","r2","r3",
                                             "r1","r2",
                                             "r1","r1","r2","r2",
                                             "r1","r1","r2","r2")))

brains <- c("Br3942","Br6522","Br8667","Br2743","Br6423","Br2720","Br6432","Br8492","Br6471","Br8325")

mlist <- list(matrix(c(1,2,3,4), nrow=2, byrow=T),
              matrix(c(1,2,3,4), nrow=2, byrow=T),
              matrix(c(1,2,3,4), nrow=2, byrow=T),
              matrix(c(1,2), ncol=1),
              matrix(c(1,2,3), ncol=3),
              matrix(c(1,2,3,4), nrow=2, byrow=T),
              matrix(c(1,2,3,4,NA,5), nrow=3, byrow=T),
              matrix(c(1,2), ncol=1),
              matrix(c(1,2,3,4), nrow=2, byrow=T),
              matrix(c(1,2,3,4), nrow=2, byrow=T))
unique(spe$brnum[order(spe$sample_id)])
#fourth one should be 6432 or 8492
#fifth one should be 6471
names(mlist) <- unique(spe$brnum[order(spe$sample_id)])

for (j in brains) {
  speb <- spe[, (colData(spe)$brnum == j)]
  speb$sample_id <- droplevels(speb$sample_id)
  #save histology image
  for(k in unique(speb$sample_id)) {
    speb_k = speb[,speb$sample_id==k]
    speb_k_title = paste(unique(speb_k$brnum)[[1]], unique(speb_k$facet_row)[[1]], unique(speb_k)$facet_column[[1]], sep="-")
    ggsave(file=paste0("plots/revision/SF1/",speb_k_title,"_histo.png"), 
           plotVisium(speb_k, spots = F, image=T, image_ids="lowres")+
             theme(plot.title=element_blank(),strip.text=element_blank()), height=1, width=1)
  }
  #plot snap25
  p2 <- plotVisium(speb, spots = T, annotate="SNAP25", assay = "logcounts", facets = "sample_id", image=F, point_size=.3)+
    facet_grid(rows=vars(facet_row), cols=vars(facet_column))+
    scale_y_continuous(expand=expansion(add=4))+scale_x_continuous(expand=expansion(add=4))+
    scale_fill_gradient(low="grey90",high="black", guide="none")+
   guides(color=NULL, fill=NULL)+
    theme(#legend.text = element_text(size = 8),
      plot.title = element_blank(),
      strip.text=element_blank(),
      panel.spacing=unit(5,'points'),plot.margin=unit(c(0,5,0,0), "pt"),
    )
  img_dim = dim(mlist[[j]])
  print(img_dim)
  ggsave(file=paste0("plots/revision/SF1/",j,"_SNAP25.png"),
         p2, bg="white", height=img_dim[[1]], width=img_dim[[2]])
}
