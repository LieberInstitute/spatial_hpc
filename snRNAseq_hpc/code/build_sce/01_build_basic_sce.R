# qrsh -l bluejay,mem_free=80G,h_vmem=80G
# cd /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/

library("SingleCellExperiment")
library("DropletUtils")
library("here")
library("rtracklayer")
library("lobstr")
library("sessioninfo")
library("dplyr")

## Define some info for the samples
tmp <- read.delim(here("snRNAseq_hpc","raw-data","sample_info",
                  "snRNAseq_U01_HPC_AllRounds_Master_Spreadsheet_12012022.csv"),
                  header = T,sep=',')

##Remove excluded samples
tmp<-tmp[c(1:2,10:27),]

##set up sample data table
sample_info <- data.frame(
  sample_ID = tmp[,1],
  sort = tmp[,2],
  BrNum = tmp[,3],
  round = tmp[,4]
)

stopifnot(all(!duplicated(sample_info$Sample)))

sample_info$sample_path<-rep(NA,20) 
sample_info$sample_path[1:2]<- file.path(
  here::here("snRNAseq_hpc","processed-data", "01_align"),
  sample_info$sample_ID[1:2],
  "outs",
  "raw_feature_bc_matrix"
)
sample_info$sample_path[3:20]<- file.path(
  here::here("snRNAseq_hpc","code", "01_align"),
  sample_info$sample_ID[3:20],
  "outs",
  "raw_feature_bc_matrix"
)
#stopifnot(all(file.exists(sample_info$sample_path)))

## Get the rest of the donor info from the spatialDLPFC SPE object
## Load SPE data
load(here("processed-data","02_build_spe","spe_basic.Rdata"))
m <- match(sample_info$BrNum, spe$brnum)
sample_info$age <- spe$age[m]
sample_info$sex <- spe$sex[m]
rm(m, spe)

## Build basic SCE
message("Read 10x data and create sce - ", Sys.time())
sce <- read10xCounts(
  sample_info$sample_path,
  sample_info$sample_ID,
  type = "sparse",
  col.names = TRUE
)
message("RDone - ", Sys.time())

save(sce,file="/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/snRNAseq_hpc/processed-data/sce/sce_raw.rda")

## Use key similar to spe objects
sce$key <- paste0(sce$Barcode, "_", sce$sample_ID)

## Add the study design info
new_col <- merge(colData(sce), sample_info[, -which(colnames(sample_info) == "sample_path")])
## Fix order
new_col <- new_col[match(sce$key, new_col$key), ]
stopifnot(identical(sce$key, new_col$key))
rownames(new_col) <- colnames(sce)
colData(sce) <- new_col

## Use code from https://github.com/LieberInstitute/Visium_IF_AD/commit/08df3f7e4a3178563d6b4b1861b664b21466b395#diff-10cb35de98e2a3e5f4235cd88f6dabce5469eead2b2db1fd7121126849fcf585L100
## Read in the gene information from the annotation GTF file
gtf <-
  rtracklayer::import(
    "/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-GRCh38-2020-A/genes/genes.gtf"
  )
gtf <- gtf[gtf$type == "gene"]
names(gtf) <- gtf$gene_id

## Match the genes
match_genes <- match(rownames(sce), gtf$gene_id)
stopifnot(all(!is.na(match_genes)))

## Keep only some columns from the gtf
mcols(gtf) <- mcols(gtf)[, c("source", "type", "gene_id", "gene_version", "gene_name", "gene_type")]

## Add the gene info to our sce object
rowRanges(sce) <- gtf[match_genes]

## Inspect object
sce

## Note that no empty droplets have been filtered out yet!

## Save for later
if (!dir.exists(here("processed-data", "sce"))) dir.create(here("processed-data", "sce"))
save(sce, file = here("processed-data", "sce", "sce_raw.Rdata"))

## Size in Gb
lobstr::obj_size(sce)


# sgejobs::job_single('01_build_basic_sce', create_shell = TRUE, queue= 'bluejay', memory = '50G', command = "Rscript 01_build_basic_sce.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()


