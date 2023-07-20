library(SingleCellExperiment)
library(scater)
library(dplyr)
library(spatialLIBD)
library(here)

load(file="/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/snRNAseq_hpc/processed-data/sce/sce_final_clustered.rda")

sce_pseudo<-registration_pseudobulk(
    sce,
    var_registration='cell.type',
    var_sample_id='Sample',
    covars = NULL,
    min_ncells = 10,
    pseudobulk_rds_file = NULL
)

save(sce_pseudo,file="/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/snRNAseq_hpc/processed-data/sce/sce_pseudo.rda")

sce_pseudo$cell.type<-droplevels(sce_pseudo$cell.type)

mod<-registration_model(
    sce_pseudo,
    covars = NULL,
    var_registration = "cell.type"
)

cors<-registration_block_cor(
    sce_pseudo,
    mod,
    var_sample_id = "Sample"
)

reg<-registration_stats_enrichment(
    sce_pseudo,
    block_cor=cors,
    covars = NULL,
    var_registration = "cell.type",
    var_sample_id = "Sample",
    gene_ensembl = 'gene_id',
    gene_name = 'gene_name'
)

rega<-registration_stats_anova(
    sce_pseudo,
    block_cor=cors,
    covars = NULL,
    var_registration = "cell.type",
    var_sample_id = "Sample",
    gene_ensembl = 'gene_id',
    gene_name = 'gene_name'
)

save(reg,rega,file='/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/snRNAseq_hpc/processed-data/pseudobulk/sce_stats.rda')
