
#cd /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/
suppressPackageStartupMessages(library("here"))
suppressPackageStartupMessages(library("SpatialExperiment"))
suppressPackageStartupMessages(library("spatialLIBD"))
suppressPackageStartupMessages(library("rtracklayer"))
suppressPackageStartupMessages(library("lobstr"))
suppressPackageStartupMessages(library("sessioninfo"))

## Define some info for the samples
load(here::here("code","REDCap","REDCap_HPC.rda"))

sample_info = data.frame(dateImg = as.Date(REDCap_HPC$date))
sample_info$experimenterImg = as.factor(REDCap_HPC$experimenter_img)
sample_info$slide = as.factor(REDCap_HPC$slide)
sample_info$array = as.factor(REDCap_HPC$array)
sample_info$brnum = as.factor(sapply(strsplit(REDCap_HPC$sample,"-"), `[`, 1))
sample_info$position = as.factor(REDCap_HPC$adjacent)
sample_info$seqNum = as.factor(REDCap_HPC$sample_number)
sample_info$experimenterSeq = as.factor(REDCap_HPC$experimenter_seq)
sample_info$sample_id = paste(sample_info$slide,sample_info$array,sep ="_")

## Define the donor info using information from
donor_info = read.csv(file.path(here::here("raw-data","sample_info","demographicInfo_Geo.csv")),header=TRUE, stringsAsFactors=FALSE)
donor_info = donor_info[-1]
