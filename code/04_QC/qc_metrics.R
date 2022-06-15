library("SpatialExperiment")
library("scuttle")
library("scran")
library("scater")
library("jaffelab")
library("tidyverse")
library("here")
library("sessioninfo")


#### Compute QC metrics ####
load(here("processed-data", "02_build_spe", "spe_basic.Rdata"), verbose = TRUE)
length(table(spe$sample_id))
# 32
length(table(spe$brnum))
# 9

spe <- scuttle::addPerCellQC(
  spe,
  subsets = list(Mito = which(seqnames(spe) == "chrM")),
  BPPARAM = BiocParallel::MulticoreParam(4)
)

#### Check for low quality spots ####

## High mito
spe$high_mito_id <- isOutlier(spe$subsets_Mito_percent, nmads = 3, type = "higher", batch = spe$sample_id)
spe$high_mito_br <- isOutlier(spe$subsets_Mito_percent, nmads = 3, type = "higher", batch = spe$brnum)
table(spe$high_mito_id)
#  FALSE   TRUE 
# 136202   1240 
table(spe$high_mito_br)
# FALSE   TRUE 
# 135392   2050 
table(spe$high_mito_id, spe$sample_id)
# V10B01-085_A1 V10B01-085_B1 V10B01-085_C1 V10B01-085_D1 V10B01-086_A1
# FALSE          3648          3287          3762          3332          4666
# TRUE             27             7            30            25            31
# 
# V10B01-086_B1 V10B01-086_C1 V10B01-086_D1 V11A20-297_A1 V11A20-297_B1
# FALSE          3946          2472          3415          4457          4622
# TRUE             31           146            46            12            44
# 
# V11A20-297_C1 V11A20-297_D1 V11L05-333_A1 V11L05-333_B1 V11L05-333_C1
# FALSE          3551          3510          4643          4955          4762
# TRUE            114            18           349            37            53
# 
# V11L05-333_D1 V11L05-335_A1 V11L05-335_B1 V11L05-335_C1 V11L05-335_D1
# FALSE          4969          4669          4765          4964          4493
# TRUE              0             0             1             1            25
# 
# V11L05-336_A1 V11L05-336_B1 V11L05-336_C1 V11L05-336_D1 V11U08-081_A1
# FALSE          4677          4458          4460          4681          4558
# TRUE              4            19             3             5            67
# 
# V11U08-081_B1 V11U08-081_C1 V11U08-081_D1 V11U08-084_A1 V11U08-084_B1
# FALSE          3728          3651          4366          4849          4600
# TRUE             37             3            81            24             0
# 
# V11U08-084_C1 V11U08-084_D1
# FALSE          4990          4296
# TRUE              0             0
table(spe$high_mito_br, spe$brnum)
# Br2743 Br3942 Br6423 Br6432 Br6471 Br6522 Br8325 Br8492 Br8667
# FALSE  13920  19057  13901   8525  14396  18759  20363   8214  18257
# TRUE     260    711    217    149      4      0    483    176     50

## low library size
spe$low_sum_id <- isOutlier(spe$sum, log = TRUE, type = "lower", batch = spe$sample_id)
table(spe$low_sum_id)
# FALSE   TRUE 
# 135695   1747 
spe$low_sum_br <- isOutlier(spe$sum, log = TRUE, type = "lower", batch = spe$brnum)
table(spe$low_sum_br)
# FALSE   TRUE 
# 135259   2183 


## low detected features
spe$low_detected_id <- isOutlier(spe$detected, log = TRUE, type = "lower", batch = spe$sample_id)
table(spe$low_detected_id)
# FALSE   TRUE 
# 135788   1654
spe$low_detected_br <- isOutlier(spe$detected, log = TRUE, type = "lower", batch = spe$brnum)
table(spe$low_detected_br)
# FALSE   TRUE 
# 135514   1928 

## are all low sum are also low detected?
table(spe$low_sum_br, spe$low_detected_br)
#        FALSE   TRUE
# FALSE 135215     44
# TRUE     299   1884
table(spe$low_sum_id, spe$low_detected_id)
#       FALSE   TRUE
# FALSE 135640     55
# TRUE     148   1599

## are all low sum are also high mito?
table(spe$low_sum_br, spe$high_mito_br)
#       FALSE   TRUE
# FALSE 133240   2019
# TRUE    2152     31
table(spe$low_sum_id, spe$high_mito_id)
#       FALSE   TRUE
# FALSE 134517   1178
# TRUE    1685     62

## Annotate spots to drop
spe$discard_auto_br <- spe$high_mito_br | spe$low_sum_br | spe$low_detected_br
spe$discard_auto_id <- spe$high_mito_id | spe$low_sum_id | spe$low_detected_id

table(spe$discard_auto_br)
# FALSE   TRUE 
# 133209   4233 
table(spe$discard_auto_id)
# FALSE   TRUE 
# 134466   2976 

## discard 3% of spots by brain
100 * sum(spe$discard_auto_br) / ncol(spe)
# [1] 3.079845
## discard 2% of spots by capture area
100 * sum(spe$discard_auto_id) / ncol(spe)
# [1] 2.165277

save(spe, file = here::here("processed-data", "04_QC", "spe_QC.Rdata"))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

