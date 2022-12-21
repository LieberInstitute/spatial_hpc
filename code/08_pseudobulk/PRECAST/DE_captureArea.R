###############################
# spatial_HPC project
# Pseudobulk DE by capture area
# Anthony Ramnauth, Dec 21 2022
###############################

setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")

suppressPackageStartupMessages({
library('SingleCellExperiment')
library('here')
library('jaffelab')
library('scater')
library('scran')
library('readxl')
library('Polychrome')
library('cluster')
library('limma')
library('sessioninfo')
library('ggplot2')
library('ggrepel')
})

## load spe data
load(file = here::here("processed-data", "08_pseudobulk", "PRECAST",
    "spe_pseudo_captureArea_wo_9-15-NA_Fncells50.Rdata"))

## Extract the data
mat <- assays(spe_pseudo)$logcounts

# make mat_formula
var_oi <- "PRECAST_k15"
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
# message("Detected correlation: ", corfit$consensus.correlation)
# Detected correlation: 0.0717245612782684

######### ENRICHMENT t-stats ######################
## Adapted from https://github.com/LieberInstitute/HumanPilot/blob/7049cd42925e00b187c0866f93409196dbcdd526/Analysis/Layer_Guesses/layer_specificity.R#L1423-L1443
# cluster_idx <- splitit(spe_pseudo$mbkmeans_harmony_9) #split by layers not path_grups
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

######### PAIRWISE t-stats ######################

fit <-
    lmFit(
        mat,
        design = mod,
        block = spe_pseudo$sample_id,
        correlation = corfit$consensus.correlation
    )
eb <- eBayes(fit)

## Define the contrasts for each group vs another one
cluster_combs <- combn(colnames(mod), 2)
cluster_constrats <- apply(cluster_combs, 2, function(x) {
    z <- paste(x, collapse = "-")
    makeContrasts(contrasts = z, levels = mod)
})
rownames(cluster_constrats) <- colnames(mod)
colnames(cluster_constrats) <-
    apply(cluster_combs, 2, paste, collapse = "-")
eb_contrasts <- eBayes(contrasts.fit(fit, cluster_constrats))

######### ANOVA t-stats ######################

## From layer_specificity.R
fit_f_model <- function(spe) {

    ## Extract the data
    mat <- assays(spe)$logcounts

    ## For dropping un-used levels
    colData(spe)[[var_oi]] <- as.factor(colData(spe)[[var_oi]])

    ## Build a group model
    # already made in beginning of script #remember to adjust for age or sex
    corfit <-
        duplicateCorrelation(mat, mod, block = spe$sample_id)
    message(paste(Sys.time(), "correlation:", corfit$consensus.correlation))
    fit <-
        lmFit(
            mat,
            design = mod,
            block = spe$sample_id,
            correlation = corfit$consensus.correlation
        )
    eb <- eBayes(fit)
    return(eb)
}

ebF_list <-
    lapply(list("full" = spe_pseudo), fit_f_model)

## Extract F-statistics
f_stats <- do.call(cbind, lapply(names(ebF_list), function(i) {
    x <- ebF_list[[i]]
    top <-
        topTable(
            x,
            coef = 2:ncol(x$coefficients),
            sort.by = "none",
            number = length(x$F)
        )
    # identical(p.adjust(top$P.Value, 'fdr'), top$adj.P.Val)
    res <- data.frame(
        "f" = top$F,
        "p_value" = top$P.Value,
        "fdr" = top$adj.P.Val,
        "AveExpr" = top$AveExpr,
        stringsAsFactors = FALSE
    )
    colnames(res) <- paste0(i, "_", colnames(res))
    return(res)
}))
f_stats$ensembl <- rownames(spe_pseudo)
f_stats$gene <- rowData(spe_pseudo)$gene_name
rownames(f_stats) <- NULL

head(f_stats)

save(f_stats, eb0_list, eb_contrasts, file = here::here("processed-data", "08_pseudobulk", "PRECAST",
    "DE_list_captureArea_adjBrnum.Rdata"))
