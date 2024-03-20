###############################
# spatial_HPC project
# Practice SCE obj creation
# Anthony Ramnauth, Dec 07 2023
###############################

suppressPackageStartupMessages({
    library(SingleCellExperiment)
    library(scuttle)
    library(DropletUtils)
    library(here)
    library(scater)
    library(uwot)
    library(rtracklayer)
    library(Seurat)
})

#########################################################
# Build SCE object from Zhao-Zhe Hao et al., 2022 dataset
# https://doi.org/10.1038/s41593-022-01073-x
#########################################################

# Load the sparse matrix

hao_sparse_mat <- Matrix::readMM(here::here("E_MTAB_10225", "fm_hippo_207785.mtx"))

dim(hao_sparse_mat)

# Wrong format, transpose?
hao_sparse_mat1 <- t(hao_sparse_mat)

# Load gene list and metadata

hao_genes <- read.csv(here::here("E_MTAB_10225", "fm_hippo_207785_gene.csv"), header = TRUE)

dim(hao_genes)

hao_barcode_info <- read.csv(here::here("E_MTAB_10225", "fm_hippo_207785_barcode.csv"), header = TRUE)

dim(hao_barcode_info)

# Add the rownames & colnames to the sparse matrix

rownames(hao_sparse_mat1) <- hao_genes$X
colnames(hao_sparse_mat1) <- hao_barcode_info$X

# Make temporary sce object from sestan lab data

sce_hao <- SingleCellExperiment(assays = list(counts = hao_sparse_mat1))

# Add coldata with metadata from sestan lab data

dim(sce_hao)

hao_coldata <- DataFrame(
    cell_name = hao_barcode_info[, "X"],
    sample_ID = hao_barcode_info[, "sample"],
    age = hao_barcode_info[, "age"],
    cellType = hao_barcode_info[, "label_20"]
)

stopifnot(identical(hao_coldata$cell_name, colnames(hao_sparse_mat1)))

colData(sce_hao) <- hao_coldata
rownames(colData(sce_hao)) <- sce_hao$cell_name
rowData(sce_hao) <- hao_genes
colnames(rowData(sce_hao)) <- c("gene_name", "gene_id")

# Save the Hao HPC  sce object

saveRDS(sce_hao, file = here::here("sce_objects", "sce_hao.rds"))

