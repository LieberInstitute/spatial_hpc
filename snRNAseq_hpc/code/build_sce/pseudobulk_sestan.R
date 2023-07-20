library(SingleCellExperiment)
library(scater)
library(dplyr)
library(spatialLIBD)
library(here)

load(file="/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/snRNAseq_hpc/processed-data/sce/sce_sestan_all.rda")

sce_sestan_pseudo<-registration_pseudobulk(
    sce_sestan,
    var_registration='cluster',
    var_sample_id='samplename',
    covars = NULL,
    min_ncells = 10,
    pseudobulk_rds_file = NULL
)

save(sce_sestan_pseudo,file="/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/snRNAseq_hpc/processed-data/sce/sce_sestan_pseudo.rda")

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

rega<-registration_stats_anova(
    sce_sestan_pseudo,
    block_cor=cors,
    covars = NULL,
    var_registration = "cluster",
    var_sample_id = "samplename",
    gene_ensembl = 'gene_id',
    gene_name = 'gene_name'
)

save(reg,rega,file='/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/snRNAseq_hpc/processed-data/pseudobulk/sestan_stats.rda')
