#######fig 6 plots#####
library(SingleCellExperiment)
library(SpatialExperiment)
library(scater)
library(scran)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(ggspavis)
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)

###load patterns
##load nmf patterns
load(file=here::here('snRNAseq_hpc','processed-data','NMF','nmf_final.rda'))

##load spe object
load(file=here::here('processed-data','NMF','spe_nmf_final.rda'))

##load sce object
load(file=here::here('snRNAseq_hpc','processed-data','sce','sce_nmf_final.rda'))

##load palettes
load(here::here('plots','snRNAseq_palettes.rda'))

###get rewritten plotVisium()) script
source(file=here::here('code','NMF','plotVisium_rewrite.R'))

###get multicolor plotVisium()) script
source(file=here::here('code','NMF','plotvisiumRGB.R'))
#######spotplots#########

brains <- unique(spe$brnum)
speb <- spe[, (colData(spe)$brnum == brains[9])]
speb$sample_id <- droplevels(speb$sample_id)
speb$sample_id <- as.character(speb$sample_id)
samples <- unique(speb$sample_id)
speb$sample_id <- factor(speb$sample_id, levels = samples)
samples
speb$brnum <- droplevels(speb$brnum)

pdf(here::here('plots','figures','figure_6','spotplot1_hpc_pats.pdf'),h=8,w=8)
plotVisiumRGB(speb,green='nmf52',pink='nmf11',blue='nmf61',yellow='nmf15',image=F,highlight='neuron_cell_body',values=c('grey40','black'))
dev.off()

pdf(here::here('plots','figures','figure_6','spotplot2_sub_pats.pdf'),h=8,w=8)
plotVisiumRGB(speb,green='nmf32',pink='nmf54',blue='nmf40',yellow='nmf17',image=F,highlight='neuron_cell_body',values=c('grey40','black'))
dev.off()


###########dotplots################
sce_pyr<-sce[,sce$superfine.cell.class %in% c('MC','CA3.1','CA3.2','CA2',
                                              'CA1','ProS','Sub.1','Sub.2',
                                              'L6.2','L6.1','L6b','L5.1','L5.2',
                                              'L2/3.1','L2/3.5','L2/3.2','L2/3.3',
                                              'L2/3.4','L2/3.6','HATA')]



create_custom_dot_plot <- function(data, category_col, features_cols,
                                   plot_title, x_axis_title, y_axis_title,
                                   legend_size_title, legend_color_title) {

    # Ensuring that features_cols is a character vector
    features_cols <- as.character(features_cols)

    # Melting the data into long format
    long_data <- data %>%
        dplyr::select(!!sym(category_col), all_of(features_cols)) %>%
        pivot_longer(cols = -!!sym(category_col), names_to = "Feature", values_to = "Value")

    # Calculating sum and proportion of non-zero cells
    stats <- long_data %>%
        group_by(!!sym(category_col), Feature) %>%
        summarize(
            Sum = sum(Value),
            NonZeroProportion = sum(Value != 0) / n()  # Explicit proportion calculation
        ) %>%
        ungroup()

    # Adjusting Feature as a factor with the specified order and reversing the category_col order
    stats$Feature <- factor(stats$Feature, levels = features_cols)
    stats[[category_col]] <- fct_relevel(stats[[category_col]], rev)
    print(stats$NonZeroProportion)

    # Creating the plot
    ggplot(stats, aes(x = Feature, y = !!sym(category_col), size = NonZeroProportion, color = Sum)) +
        geom_point() +
        scale_size_continuous(range = c(0,10)) + # Set minimum size to zero
        scale_color_gradient(low = "white", high = "black") + # Greyscale color scale
        theme_minimal() +
        labs(
            title = plot_title,
            x = x_axis_title,
            y = y_axis_title,
            size = legend_size_title,
            color = legend_color_title
        ) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Adjust x-axis text angle for readability
    }

####make sce spotplot
data<-as.data.frame(colData(sce_pyr))
data$superfine.cell.class<-factor(data$superfine.cell.class,levels=c("MC",
                                                                     "CA3.1", "CA3.2", "CA2", "CA1", "ProS", "Sub.1",
                                                                     "Sub.2",
                                                                     "L6.2", "L6.1", "L6b",
                                                                     "L5.1", "L5.2",
                                                                     "L2/3.1","L2/3.5", "L2/3.2", "L2/3.3", "L2/3.4",
                                                                     "L2/3.6",
                                                                     "HATA"))
colnames(data)[match(c('nmf52','nmf11','nmf63','nmf61','nmf15',
                                  'nmf32','nmf40','nmf54','nmf22','nmf65',
                                  'nmf53','nmf51','nmf68','nmf17','nmf78','nmf27','nmf45',
                                  'nmf84','nmf62','nmf69'),colnames(data))]<-c(
                                      'nmf52 MC','nmf11 CA3.1','nmf63 CA3.2','nmf61 CA2','nmf15 CA1',
                                      'nmf32 ProS','nmf40 Sub.1','nmf54 Sub.2','nmf22 L6.2','nmf65 L6/6b',
                                      'nmf53 L6b','nmf51 L5/6.1','nmf68 L5/6','nmf17 L2/3.1','nmf78 L2/3.5','nmf27 L2/3.2','nmf45 L2/3',
                                      'nmf84 L2/3','nmf62 HATA','nmf69 HATA'
                                  )

indices<-c(
    'nmf52 MC','nmf11 CA3.1','nmf63 CA3.2','nmf61 CA2','nmf15 CA1',
    'nmf32 ProS','nmf40 Sub.1','nmf54 Sub.2','nmf22 L6.2','nmf65 L6/6b',
    'nmf53 L6b','nmf51 L5/6.1','nmf68 L5/6','nmf17 L2/3.1','nmf78 L2/3.5','nmf27 L2/3.2','nmf45 L2/3',
    'nmf84 L2/3','nmf62 HATA','nmf69 HATA'
)

pdf(file=here::here('plots','figures','figure_6','snRNAseq_dotplot_final.pdf'),h=10.5,w=16)
create_custom_dot_plot(data, "superfine.cell.class", indices, "snRNA-seq", "NMF pattern",
                                                        "superfine cell class", "proportion nuclei\nwith nonzero\nweight",
                                                        "aggregate\nnuclei-level\nweights")+
                                                        theme(axis.text=element_text(size=32,color='black'),text=element_text(size=32,color='black'))
dev.off()


####now make spatial spotplot
spe_pyr<-spe[,spe$domain %in% c('CA2.4','CA1','SUB','SUB.RHP','RHP')]
data<-as.data.frame(colData(spe_pyr))
colnames(data)[match(c('nmf52','nmf11','nmf63','nmf61','nmf15',
                       'nmf32','nmf40','nmf54','nmf22','nmf65',
                       'nmf53','nmf51','nmf68','nmf17','nmf78','nmf27','nmf45',
                       'nmf84','nmf62','nmf69'),colnames(data))]<-c(
                           'nmf52 MC','nmf11 CA3.1','nmf63 CA3.2','nmf61 CA2','nmf15 CA1',
                           'nmf32 ProS','nmf40 Sub.1','nmf54 Sub.2','nmf22 L6.2','nmf65 L6/6b',
                           'nmf53 L6b','nmf51 L5/6.1','nmf68 L5/6','nmf17 L2/3.1','nmf78 L2/3.5','nmf27 L2/3.2','nmf45 L2/3',
                           'nmf84 L2/3','nmf62 HATA','nmf69 HATA'
                       )

indices<-c(
    'nmf52 MC','nmf11 CA3.1','nmf63 CA3.2','nmf61 CA2','nmf15 CA1',
    'nmf32 ProS','nmf40 Sub.1','nmf54 Sub.2','nmf22 L6.2','nmf65 L6/6b',
    'nmf53 L6b','nmf51 L5/6.1','nmf68 L5/6','nmf17 L2/3.1','nmf78 L2/3.5','nmf27 L2/3.2','nmf45 L2/3',
    'nmf84 L2/3','nmf62 HATA','nmf69 HATA'
)
pdf(file=here::here('plots','figures','figure_6','spatial_dotplot_final.pdf'),h=8,w=17)
create_custom_dot_plot(data, "domain", indices, "SRT", "NMF pattern",
                       "domain", "proportion spots\nwith nonzero\nweight",
                       "aggregate\nspot-level\nweights")+
    theme(axis.text=element_text(size=32,color='black'),text=element_text(size=32,color='black'))
dev.off()

###now mch
load(file=here::here('snRNAseq_hpc','processed-data','NMF','mch.rda'))

data<-as.data.frame(colData(mch))
colnames(data)[match(c('nmf52','nmf11','nmf63','nmf61','nmf15',
                       'nmf32','nmf40','nmf54','nmf22','nmf65',
                       'nmf53','nmf51','nmf68','nmf17','nmf78','nmf27','nmf45',
                       'nmf84','nmf62','nmf69'),colnames(data))]<-c(
                           'nmf52 MC','nmf11 CA3.1','nmf63 CA3.2','nmf61 CA2','nmf15 CA1',
                           'nmf32 ProS','nmf40 Sub.1','nmf54 Sub.2','nmf22 L6.2','nmf65 L6/6b',
                           'nmf53 L6b','nmf51 L5/6.1','nmf68 L5/6','nmf17 L2/3.1','nmf78 L2/3.5','nmf27 L2/3.2','nmf45 L2/3',
                           'nmf84 L2/3','nmf62 HATA','nmf69 HATA'
                       )

indices<-c(
    'nmf52 MC','nmf11 CA3.1','nmf63 CA3.2','nmf61 CA2','nmf15 CA1',
    'nmf32 ProS','nmf40 Sub.1','nmf54 Sub.2','nmf22 L6.2','nmf65 L6/6b',
    'nmf53 L6b','nmf51 L5/6.1','nmf68 L5/6','nmf17 L2/3.1','nmf78 L2/3.5','nmf27 L2/3.2','nmf45 L2/3',
    'nmf84 L2/3','nmf62 HATA','nmf69 HATA'
)
pdf(file=here::here('plots','figures','figure_6','spatial_dotplot_final.pdf'),h=8,w=17)
create_custom_dot_plot(data, "Subclass", indices, "", "NMF pattern",
                       "Allen subclass", "proportion nuclei\nwith nonzero\nweight",
                       "aggregate\nnuclei-level\nweights")+
    theme(axis.text=element_text(size=32,color='black'),text=element_text(size=32,color='black'))
dev.off()

pdf(file=here::here('plots','figures','figure_6','spatial_dotplot_final.pdf'),h=8,w=17)
create_custom_dot_plot(data, "Target", indices, "", "NMF pattern",
                       "Target", "proportion nuclei\nwith nonzero\nweight",
                       "aggregate\nnuclei-level\nweights")+
    theme(axis.text=element_text(size=32,color='black'),text=element_text(size=32,color='black'))
dev.off()



###########heatmap###########
features<-c('EYS','ZNF208','HGF','TGFB2','FIBCD1','TG','FN1',
            'COL12A1','SULF1','SDK2','GSG1L','RORB','AGBL1','VIPR2','NIPAL2',
            'CRHR2','CBLN4','RSPO2','KLHL1','TACR3')



heat<-x@w[features,c('nmf52','nmf11','nmf63','nmf61','nmf15',
                     'nmf32','nmf40','nmf54','nmf22','nmf65',
                     'nmf53','nmf51','nmf68','nmf17','nmf78','nmf27','nmf45',
                     'nmf84','nmf62','nmf69')]

pdf(file=here::here('plots','figures','figure_6','marker_gene_heatmap.pdf'),h=4.5,w=4)
pheatmap(heat,cluster_cols=F,cluster_rows=F,
         color=colorRampPalette(brewer.pal(n = 7, name ="Greys"))(100),
         breaks=seq(0,0.001,0.00001
         ))
dev.off()

