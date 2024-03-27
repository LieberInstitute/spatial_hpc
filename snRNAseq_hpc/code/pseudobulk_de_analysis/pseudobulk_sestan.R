library(SingleCellExperiment)
library(scater)
library(dplyr)
library(spatialLIBD)
library(here)

load(file=here::here('snRNAseq_hpc','processed-data','sce','sce_sestan_all.rda')

sce_sestan_pseudo<-registration_pseudobulk(
    sce_sestan,
    var_registration='cluster',
    var_sample_id='samplename',
    covars = NULL,
    min_ncells = 10,
    pseudobulk_rds_file = NULL
)

save(sce_sestan_pseudo,file=here::here('snRNAseq_hpc','processed-data','sce','sce_sestan_pseudo.rda')

sce_sestan_pseudo$cluster<-droplevels(sce_sestan_pseudo$cluster)

mod<-registration_model(
    sce_sestan_pseudo,
    covars = NULL,
    var_registration = "cluster"
)

cors<-registration_block_cor(
    sce_sestan_pseudo,
    mod,
    var_sample_id = "samplename"
)

reg<-registration_stats_enrichment(
    sce_sestan_pseudo,
    block_cor=cors,
    covars = NULL,
    var_registration = "cluster",
    var_sample_id = "samplename",
    gene_ensembl = 'gene_id',
    gene_name = 'gene_name'
)

regp<-registration_stats_pairwise(
    sce_sestan_pseudo,
    block_cor=cors,
    registration_model=mod,
    var_registration = "cluster",
    var_sample_id = "samplename",
    gene_ensembl = 'gene_id',
    gene_name = 'gene_name'
)

rega<-registration_stats_anova(
    sce_sestan_pseudo,
    block_cor=cors,
    covars = NULL,
    var_registration = "cluster",
    var_sample_id = "samplename",
    gene_ensembl = 'gene_id',
    gene_name = 'gene_name'
)

save(rega,file=here::here('snRNAseq_hpc','processed-data','de_analysis','sce_anova_stats_sestan.rda'))
save(regp,file=here::here('snRNAseq_hpc','processed-data','de_analysis','sce_pairwise_stats_sestan.rda'))
save(reg,file=here::here('snRNAseq_hpc','processed-data','de_analysis','sce_enrichment_stats_sestan.rda'))