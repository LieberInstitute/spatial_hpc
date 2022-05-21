
# cd /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/
library("here")
suppressPackageStartupMessages(library("SpatialExperiment"))
suppressPackageStartupMessages(library("spatialLIBD"))
suppressPackageStartupMessages(library("rtracklayer"))
suppressPackageStartupMessages(library("lobstr"))
suppressPackageStartupMessages(library("sessioninfo"))

## Define some info for the samples
load(here::here("code", "REDCap", "REDCap_HPC.rda"))

sample_info <- data.frame(dateImg = as.Date(REDCap_HPC$date))
sample_info$experimenterImg <- as.factor(REDCap_HPC$experimenter_img)
sample_info$slide <- as.factor(REDCap_HPC$slide)
sample_info$array <- as.factor(REDCap_HPC$array)
sample_info$brnum <- as.factor(sapply(strsplit(REDCap_HPC$sample, "-"), `[`, 1))
sample_info$position <- as.factor(REDCap_HPC$adjacent)
sample_info$seqNum <- as.factor(REDCap_HPC$sample_number)
sample_info$experimenterSeq <- as.factor(REDCap_HPC$experimenter_seq)
sample_info$sample_id <- paste(sample_info$slide, sample_info$array, sep = "_")

sample_info$sample_path[sample_info$dateImg <= "2021-10-11"] <- file.path(here::here("processed-data", "spaceranger_novaseq"), sample_info$sample_id, "outs")
sample_info$sample_path[sample_info$dateImg > "2021-10-11"] <- file.path(here::here("processed-data", "spaceranger_2022-04-12_SPag033122"), sample_info$sample_id[sample_info$dateImg > "2021-10-11"], "outs")
stopifnot(all(file.exists(sample_info$sample_path)))

## Define the donor info using information from
donor_info <- read.csv(file.path(here::here("raw-data", "sample_info", "demographicInfo_Geo.csv")), header = TRUE, stringsAsFactors = FALSE)
donor_info <- donor_info[-1]

## Combine sample info with the donor info
sample_info <- merge(sample_info, donor_info)

## Build basic SPE
Sys.time()
spe <- read10xVisiumWrapper(
    sample_info$sample_path,
    sample_info$sample_id,
    type = "sparse",
    data = "raw",
    images = c("lowres", "hires", "detected", "aligned"),
    load = TRUE
)
Sys.time()

## Add the study design info
add_design <- function(spe) {
    new_col <- merge(colData(spe), sample_info)
    ## Fix order
    new_col <- new_col[match(spe$key, new_col$key), ]
    stopifnot(identical(new_col$key, spe$key))
    rownames(new_col) <- rownames(colData(spe))
    colData(spe) <-
        new_col[, -which(colnames(new_col) == "sample_path")]
    return(spe)
}
spe <- add_design(spe)

## basic SPE
spe_raw <- spe
# dir.create(here::here("processed-data", "pilot_data_checks"), showWarnings = FALSE)
save(spe_raw, file = here::here("processed-data", "02_build_spe", "spe_raw.Rdata"))

## Size in Gb
lobstr::obj_size(spe_raw) / 1024^3
# 5.141452 B
dim(spe_raw)
# 36601 159744

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
