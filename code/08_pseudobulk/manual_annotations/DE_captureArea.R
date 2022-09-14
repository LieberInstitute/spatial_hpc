setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")

library('SingleCellExperiment')
library('here')
library('jaffelab')
library('scater')
library('scran')
library('pheatmap')
library('readxl')
library('Polychrome')
library('cluster')
library('limma')
library('sessioninfo')
library('ggplot2') 
library('ggrepel')

## load spe data
load(file = here::here("processed-data", "08_pseudobulk", "manual_annotations", "spe_pseudo_captureArea_wo_CP-THAL-CTX_Fncells50_Fol.Rdata"))

## Extract the data
mat <- assays(spe_pseudo)$logcounts

# make mat_formula
var_oi <- "manual_annotations"
covars <- c("age", "sex")
mat_formula <- eval(str2expression(paste("~", "0", "+", var_oi, "+", paste(covars, collapse = " + "))))

# make sure everything is  a factor
colData(spe_pseudo)[[var_oi]] <- as.factor(colData(spe_pseudo)[[var_oi]])
colData(spe_pseudo)$age <- as.numeric(colData(spe_pseudo)$age)
colData(spe_pseudo)$sex <- as.factor(colData(spe_pseudo)$sex)
colData(spe_pseudo)$brnum <- as.factor(colData(spe_pseudo)$brnum)
colData(spe_pseudo)$sample_id <- as.factor(colData(spe_pseudo)$sample_id)

## Compute correlation
## Adapted from https://github.com/LieberInstitute/Visium_IF_AD/blob/7973fcebb7c4b17cc3e23be2c31ac324d1cc099b/code/10_spatial_registration/01_spatial_registration.R#L134-L150
mod <- model.matrix(mat_formula, data = colData(spe_pseudo))
message(Sys.time(), " running duplicateCorrelation()")
corfit <- duplicateCorrelation(mat, mod, block = spe_pseudo$sample_id)
message("Detected correlation: ", corfit$consensus.correlation)

######### ENRICHMENT t-stats ######################
## Adapted from https://github.com/LieberInstitute/HumanPilot/blob/7049cd42925e00b187c0866f93409196dbcdd526/Analysis/Layer_Guesses/layer_specificity.R#L1423-L1443

cluster_idx <- splitit(colData(spe_pseudo)[, var_oi])

message(Sys.time(), " running the enrichment model")
eb0_list <- lapply(cluster_idx, function(x) {
  res <- rep(0, ncol(spe_pseudo))
  res[x] <- 1
  res_formula <- paste("~", "res", "+", paste(covars, collapse = " + "))
  m <- with(
    colData(spe_pseudo),
    model.matrix(eval(str2expression(res_formula)))
  )
  eBayes(
    lmFit(
      mat,
      design = m,
      block = spe_pseudo$sample_id,
      correlation = corfit$consensus.correlation
    )
  )
})

save(eb0_list, file = here::here("processed-data", "08_pseudobulk", "manual_annotations", "DE_eb0_list_captureArea.Rdata"))