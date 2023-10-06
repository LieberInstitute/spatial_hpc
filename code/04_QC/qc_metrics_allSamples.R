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
# [1]  31483 191136
length(table(spe$sample_id))
# 44
length(table(spe$brnum))
# 12

spe <- scuttle::addPerCellQC(
    spe,
    subsets = list(Mito = which(seqnames(spe) == "chrM")),
    BPPARAM = BiocParallel::MulticoreParam(4)
)

##### Check for low quality spots ####
#
### High mito
###NOTE:based on exploratory data analysis, discarding spots based on mitochondrial rate appears to discard large
###numbers of neuropil spots. This code has been commented. Mito rate is still stored in spe for analysis purposes
#spe$high_mito_id <- isOutlier(spe$subsets_Mito_percent, nmads = 3, type = "higher", batch = spe$sample_id)
#spe$high_mito_br <- isOutlier(spe$subsets_Mito_percent, nmads = 3, type = "higher", batch = spe$brnum)

## low library size
spe$low_sum_id <- isOutlier(spe$sum, log = TRUE, type = "lower", batch = spe$sample_id)
table(spe$low_sum_id)
# FALSE   TRUE
# 188834   2302
spe$low_sum_br <- isOutlier(spe$sum, log = TRUE, type = "lower", batch = spe$brnum)
table(spe$low_sum_br)
# FALSE   TRUE
# 1188440   2696

#are low_sum_br spots and low_sum_id spots same? Mostly
table(spe$low_sum_br,spe$low_sum_id)
#              id    id
#           FALSE   TRUE
# br FALSE 188212    228
# br TRUE     622   2074

## low detected features
spe$low_detected_id <- isOutlier(spe$detected, log = TRUE, type = "lower", batch = spe$sample_id)
table(spe$low_detected_id)
# FALSE   TRUE
# 188971   2165
spe$low_detected_br <- isOutlier(spe$detected, log = TRUE, type = "lower", batch = spe$brnum)
table(spe$low_detected_br)
# FALSE   TRUE
# 188735   2401



### are all low sum are also low detected? Mostly
#table(spe$low_sum_br, spe$low_detected_br)
## FALSE   TRUE
## FALSE 150522     44
## TRUE     328   2189
#table(spe$low_sum_id, spe$low_detected_id)
## FALSE   TRUE
## FALSE 150917     59
## TRUE     188   1919
#
### are all low sum are also high mito? No
#table(spe$low_sum_br, spe$high_mito_br)
## FALSE   TRUE
## FALSE 148499   2067
## TRUE    2485     32
#table(spe$low_sum_id, spe$high_mito_id)
## FALSE   TRUE
## FALSE 149779   1197
## TRUE    2044     63

## Annotate spots to drop
spe$discard_auto_br <- spe$low_sum_br | spe$low_detected_br
spe$discard_auto_id <- spe$low_sum_id | spe$low_detected_id

table(spe$discard_auto_br)
# FALSE   TRUE
# 188388   2748
table(spe$discard_auto_id)
# FALSE   TRUE
# 188762   2374

## discard 3% of spots by brain
100 * sum(spe$discard_auto_br) / ncol(spe)
# [1] 1.43772
## discard 2% of spots by capture area
100 * sum(spe$discard_auto_id) / ncol(spe)
# [1] 1.242048

save(spe, file = here::here("processed-data", "04_QC", "spe_QC_allSamples.Rdata"))

## Reproducibility information
Sys.time()
proc.time()
options(width = 120)
session_info()
