library("spatialLIBD")
library("markdown") ## Hm... to avoid this error
# 2021-11-11T05:30:49.941401+00:00 shinyapps[5096402]: Listening on http://127.0.0.1:32863
# 2021-11-11T05:30:50.218127+00:00 shinyapps[5096402]: Warning: Error in loadNamespace: there is no package called ‘markdown’
# 2021-11-11T05:30:50.222437+00:00 shinyapps[5096402]:   111: <Anonymous>

## spatialLIBD uses golem
options("golem.app.prod" = TRUE)

## You need this to enable shinyapps to install Bioconductor packages
options(repos = BiocManager::repositories())

## Data setup

## Download to my laptop
# scp e:/dcs04/lieber/lcolladotor/with10x_LIBD001/Visium_IF_AD/processed-data/spe/spe.Rdata processed-data/spe/

## Create a soft link to the data, otherwise rsconnect::deployApp doesn't work
## Note that soft link has to be relative to work
# cd /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/code/03_shiny_app_basic/
# ln -s ../../../processed-data/02_build_spe/spe_basic.Rdata spe.Rdata

## Load the data
load(file="spe.Rdata", verbose = TRUE)

#speB$BayesSpace <- speB$spatial.cluster
#speB$BayesSpace_initial <- speB$cluster.init 
vars <- colnames(colData(spe))
lobstr::obj_size(spe) / 1024^3

imgData(spe) <- imgData(spe)[!imgData(spe)$image_id %in% c("hires", "detected", "aligned"),]
lobstr::obj_size(spe) / 1024^3

## Deploy the website
spatialLIBD::run_app(
  spe,
  sce_layer = NULL,
  modeling_results = NULL,
  sig_genes = NULL,
  title = "Visium HPC 2022",
  spe_discrete_vars = c(
    vars[grep("^10x_", vars)],
    "ManualAnnotation",
    "brnum",
    "sex",
    "race",
    "sample_id"
  ),
  spe_continuous_vars = c(
    "sum_umi",
    "sum_gene",
    "expr_chrM",
    "expr_chrM_ratio",
    "age"
  ),
  default_cluster = "10x_graphclust"
)
