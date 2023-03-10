
setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")
suppressPackageStartupMessages({
library("here")
library("SpatialExperiment")
library("spatialLIBD")
library("rtracklayer")
library("lobstr")
library("sessioninfo")
})

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

sample_info$sample_path[sample_info$dateImg <= "2021-10-11"] <- file.path(here::here("processed-data", "01_spaceranger", "spaceranger_novaseq"), sample_info$sample_id, "outs")
sample_info$sample_path[sample_info$dateImg > "2021-10-11"] <- file.path(here::here("processed-data", "01_spaceranger", "spaceranger_2022-04-12_SPag033122"), sample_info$sample_id[sample_info$dateImg > "2021-10-11"], "outs")
stopifnot(all(file.exists(sample_info$sample_path)))

## Define the donor info using information from
donor_info <- read.csv(file.path(here::here("raw-data", "sample_info", "demographicInfo_Geo.csv")), header = TRUE, stringsAsFactors = FALSE)
donor_info <- donor_info[-1]

## Combine sample info with the donor info
sample_info <- merge(sample_info, donor_info)

###We gotta make the round 9 info object too!
sample_info9<-data.frame(
  dateImg=rep(as.Date('2023-01-31'),4),
  experimenterImg=factor(rep("Stephanie Page",4)),
  slide=factor(rep('V12F14-051',4)),
  array=factor(c('A1','B1','C1','D1')),
  brnum=factor(rep('Br2720',4)),
  position=factor(c('BL','BR','TL','TR')),
  seqNum=factor(c("33-10v", "34v_scp", "35v_scp","36v_scp")),
  experimenterSeq=factor(rep("Stephanie Page",4)),
  sample_id = c(
    "Br2720_A1",
    "Br2720_B1",
    "Br2720_C1",
    "Br2720_D1"
  )
)

sample_info9$sample_path <-
  file.path(
    here::here("processed-data", "01_spaceranger", "spaceranger_2023-01-31_round9"),
    sample_info9$sample_id,
    "outs"
  )
stopifnot(all(file.exists(sample_info9$sample_path)))

donor_info9 <- data.frame(
  age = c(48.2, 48.2, 48.2, 48.2),
  sex = c("F", "F", "F", "F"),
  race = c("CAUC", "CAUC", "CAUC", "CAUC"),
  dx = c("Control", "Control", "Control", "Control"),
  pmi = c(25.5, 25.5, 25.5, 25.5)
)

sample_info9 <- cbind(sample_info9, donor_info9)

##match up the colnames
sample_info9<-sample_info9[,match(colnames(sample_info),colnames(sample_info9))]

##rbind the sample_infos
sample_info<-rbind(sample_info,sample_info9)

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

# dir.create(here::here("processed-data", "pilot_data_checks"), showWarnings = FALSE)
save(spe, file = here::here("processed-data", "02_build_spe", "spe_raw_allSamples.Rdata"))

## Size in Gb
lobstr::obj_size(spe)
# 6.16 GB
dim(spe)
# 36601 179712

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
