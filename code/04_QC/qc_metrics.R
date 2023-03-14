setwd('/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/')
suppressPackageStartupMessages(library("SpatialExperiment"))
suppressPackageStartupMessages(library("scuttle"))
suppressPackageStartupMessages(library("scran"))
suppressPackageStartupMessages(library("scater"))
suppressPackageStartupMessages(library("jaffelab"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("here"))
suppressPackageStartupMessages(library("sessioninfo"))


#### Compute QC metrics ####
load(here("processed-data", "02_build_spe", "spe_basic.Rdata"), verbose = TRUE)
dim(spe)
# [1]  30359 137442
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
table(spe$sample_id, spe$high_mito_id)
#                FALSE TRUE
# V10B01-085_A1  3648   27
# V10B01-085_B1  3287    7
# V10B01-085_C1  3762   30
# V10B01-085_D1  3332   25
# V10B01-086_A1  4666   31
# V10B01-086_B1  3946   31
# V10B01-086_C1  2472  146
# V10B01-086_D1  3415   46
# V11A20-297_A1  4457   12
# V11A20-297_B1  4622   44
# V11A20-297_C1  3551  114
# V11A20-297_D1  3510   18
# V11L05-333_A1  4643  349
# V11L05-333_B1  4955   37
# V11L05-333_C1  4762   53
# V11L05-333_D1  4969    0
# V11L05-335_A1  4669    0
# V11L05-335_B1  4765    1
# V11L05-335_C1  4964    1
# V11L05-335_D1  4493   25
# V11L05-336_A1  4677    4
# V11L05-336_B1  4458   19
# V11L05-336_C1  4460    3
# V11L05-336_D1  4681    5
# V11U08-081_A1  4558   67
# V11U08-081_B1  3728   37
# V11U08-081_C1  3651    3
# V11U08-081_D1  4366   81
# V11U08-084_A1  4849   24
# V11U08-084_B1  4600    0
# V11U08-084_C1  4990    0
# V11U08-084_D1  4296    0
table(spe$high_mito_br, spe$brnum)
#        Br2743 Br3942 Br6423 Br6432 Br6471 Br6522 Br8325 Br8492 Br8667
# FALSE  13920  19057  13901   8525  14396  18759  20363   8214  18257
# TRUE     260    711    217    149      4      0    483    176     50

#are high_mito_br spots and high_mito_id spots same? No
table(spe$high_mito_br,spe$high_mito_id)
#              id    id
#            FALSE   TRUE
# br FALSE  134496   896
# br TRUE    1706    344

## low library size
spe$low_sum_id <- isOutlier(spe$sum, log = TRUE, type = "lower", batch = spe$sample_id)
table(spe$low_sum_id)
# FALSE   TRUE
# 135695   1747
spe$low_sum_br <- isOutlier(spe$sum, log = TRUE, type = "lower", batch = spe$brnum)
table(spe$low_sum_br)
# FALSE   TRUE
# 135259   2183

#are low_sum_br spots and low_sum_id spots same? Mostly
table(spe$low_sum_br,spe$low_sum_id)
#              id    id
#           FALSE   TRUE
# br FALSE 135105    154
# br TRUE     590   1593

## low detected features
spe$low_detected_id <- isOutlier(spe$detected, log = TRUE, type = "lower", batch = spe$sample_id)
table(spe$low_detected_id)
# FALSE   TRUE
# 135788   1654
spe$low_detected_br <- isOutlier(spe$detected, log = TRUE, type = "lower", batch = spe$brnum)
table(spe$low_detected_br)
# FALSE   TRUE
# 135514   1928

#are low_detected_br spots and low_detected_br spots same? Yes
table(spe$low_detected_br,spe$low_detected_id)
#              id    id
#           FALSE   TRUE
# br FALSE 135514      0
# br TRUE       0   1928


## are all low sum are also low detected? Mostly
table(spe$low_sum_br, spe$low_detected_br)
#        FALSE   TRUE
# FALSE 135215     44
# TRUE     299   1884
table(spe$low_sum_id, spe$low_detected_id)
#       FALSE   TRUE
# FALSE 135640     55
# TRUE     148   1599

## are all low sum are also high mito? No
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
Sys.time()
proc.time()
options(width = 120)
session_info()
