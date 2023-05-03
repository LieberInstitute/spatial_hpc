library(SingleCellExperiment)
library(scater)
library(scran)
library(liana)
library(here)
library(dplyr)
library(magrittr)
load(file=here::here('snRNAseq_hpc','processed-data','sce','sce_clustered.rda'))
reducedDim(sce,'PCA')<-NULL
reducedDim(sce,'UMAP')<-NULL
reducedDim(sce,'HARMONY')<-NULL

##run liana
liana_test <- liana_wrap(sce,idents_col = 'k_50_label')
save(liana_test,file=here::here('snRNAseq_hpc','processed-data',
                                'cell_cell_communication','liana_test.rda'))
liana_ag <- liana_aggregate(liana_test)
save(liana_ag,file=here::here('snRNAseq_hpc','processed-data',
                                'cell_cell_communication','liana_test.rda'))

##make some plots
pdf(here::here("snRNAseq_hpc","plots","cell_cell_communication","liana_dots_GCs.pdf"),h=15,w=15)
liana_ag %>%
  liana_dotplot(source_groups = c("8"),
                target_groups = c("28", "6", "33", "25","26","32","19"),
                ntop = 20)
dev.off()

pdf(here::here("snRNAseq_hpc","plots","cell_cell_communication","liana_dots_CA3.pdf"),h=15,w=15)
liana_ag %>%
  liana_dotplot(source_groups = c("6"),
                target_groups = c("28", "8", "33", "25","26","32","19"),
                ntop = 20)
dev.off()

pdf(here::here("snRNAseq_hpc","plots","cell_cell_communication","liana_dots_CA1.pdf"),h=15,w=15)
liana_ag %>%
  liana_dotplot(source_groups = c("28"),
                target_groups = c("6", "8", "33", "25","26","32","19"),
                ntop = 25)
dev.off()

pdf(here::here("snRNAseq_hpc","plots","cell_cell_communication","liana_dots_CA1.pdf"),h=15,w=15)
liana_ag %>%
  liana_dotplot(source_groups = c("26"),
                target_groups = c("6", "8", "33", "25","28","32","19"),
                ntop = 25)
dev.off()

pdf(here::here("snRNAseq_hpc","plots","cell_cell_communication","liana_heat.pdf"),h=15,w=15)
liana_trunc <- liana_ag %>%
   # only keep interactions concordant between methods
  filter(aggregate_rank <= 0.01) # note that these pvals are already corrected

heat_freq(liana_trunc)

dev.off()

