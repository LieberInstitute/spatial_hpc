#https://bioconductor.org/packages/release/data/experiment/vignettes/humanHippocampus2024/inst/doc/humanHippocampus2024.html
BiocManager::install("humanHippocampus2024")

library(SummarizedExperiment)
library(SpatialExperiment)
library(humanHippocampus2024)

library(ExperimentHub)
ehub <- ExperimentHub()

myfiles <- query(ehub, "humanHippocampus2024")

spatial_hpc_spe <- myfiles[["EH9605"]]
table(spatial_hpc_spe$brnum)

srt.coldata = as.data.frame(colData(spatial_hpc_spe))
srt.coldata[1:4,1:4]
write.csv(srt.coldata, "processed-data/revision/humanHippocampus2024_SRT_colData.csv", row.names=T)


spatial_hpc_snrna_seq <- myfiles[["EH9606"]]

sn.coldata = as.data.frame(colData(spatial_hpc_snrna_seq))
sn.coldata[1:4,1:4]
write.csv(sn.coldata, "processed-data/revision/humanHippocampus2024_snRNAseq_colData.csv", row.names=T)




# SRT: H&E or SPG, orientation/ grid arrangement
# https://github.com/LieberInstitute/spatial_hpc/blob/main/code/02_build_spe/02_transform_spe_allSamples.R
# for rotation angle I need to follow some code because the data was reordered twice and then the
# rotation angle is an unnamed vector based on sample_id order, so I need the samples to be in the same order
library(dplyr)
load(file = "processed-data/02_build_spe/spe_raw_allSamples.Rdata")
spe_raw = spe
dim(spe_raw)
# [1]  36601 219648

spe_raw$brnum<-as.character(spe_raw$brnum)
spe_raw$brnum<-ifelse(spe_raw$slide %in% "V12D07-335",'Br8325_VSPG',spe_raw$brnum)
spe_raw$brnum<-ifelse(spe_raw$slide %in% "V12D07-332",'Br3942_VSPG',spe_raw$brnum)
spe_raw$brnum<-factor(spe_raw$brnum,levels=unique(spe_raw$brnum))

spe_raw$sample_id <- paste0(spe_raw$slide,"_",spe_raw$array)

## arrange by brain and array (A1, B1, C1, D1)
raw.coldata = as.data.frame(colData(spe_raw))
#spe_raw <- arrange(spe_raw, brnum, array)
raw.coldata <- arrange(raw.coldata, brnum , array)

samples <- unique(raw.coldata$sample_id)
samples #check against output in original code file

position_list <- c("TR", "BL", "TL", "BR",
                   "TL", "TR", "BL", "BR",
                   "TR", "TL", "BR", "BL",
                   "TL", "BL",
                   "BL", "TR", "TL",
                   "TL", "TR", "BL", "BR",
                   "BL", "BR", "TL", "TR", "TR",
                   "TL", "BL",
                   "TL", "TR", "BL", "BR",
                   "BL","BR","TL","TR",
                   "TL",'BR',"BL",'TR',
                   "BR",'TL',"BL",'TR')

for (i in seq_along(samples)){
  raw.coldata$position[raw.coldata$sample_id == samples[i]] = position_list[i]
}
table(raw.coldata$position, useNA="ifany")
raw.coldata$position <- factor(raw.coldata$position, levels = c("TL", "TR", "BL", "BR"))
table(raw.coldata$position, useNA="ifany")

## now arrange() to get the correct order
raw.coldata <- arrange(raw.coldata, brnum, position)


samples = unique(raw.coldata$sample_id)
samples #check against output in original code file


# Rotations by sample:
angle_list <- c(270, 90, 0, 0,
                0, 0, 0, 0,
                180, 90, 180, 180, 
                0, 0,
                0, 0, 0, 
                0, 0, 0, 0,
                0, 0, 0, 0, 0,
                90, 90, 
                0, 0, 0, 0,
                0,0,0,0,
                270,270,270,270,
                270,270,270,270)

srt.angles = data.frame("sample_id"=samples, "rotateImg_angle"=angle_list)


srt.tmp = group_by(srt.coldata, brnum, sample_id) %>% tally(name="n_spots")
srt.tmp$sample_id = factor(srt.tmp$sample_id, levels=c("V10B01-086_D1","V10B01-086_C1","V11U08-081_C1","V11U08-081_D1",
                                               "V11L05-333_A1","V11L05-333_B1","V11L05-333_C1","V11L05-333_D1",
                                               "V10B01-085_B1","V10B01-085_A1","V10B01-085_D1","V10B01-085_C1",
                                               "V10B01-086_A1","V10B01-086_B1",
                                               #"V11L05-335_C1","V11L05-335_B1","V11L05-335_A1",
                                               "V11L05-335_A1","V11L05-335_C1","V11L05-335_B1",
                                               "V11U08-084_A1","V11U08-084_B1","V11U08-084_C1","V11U08-084_D1",
                                               "V11A20-297_C1","V11A20-297_D1","V11A20-297_A1","V11L05-335_D1","V11A20-297_B1",
                                               "V11U08-081_A1","V11U08-081_B1",
                                               "V11L05-336_A1","V11L05-336_B1","V11L05-336_C1","V11L05-336_D1",
                                               "V12F14-051_C1","V12F14-051_D1","V12F14-051_A1","V12F14-051_B1"))
srt.tmp$grid_column = as.character(factor(srt.tmp$sample_id, levels=levels(srt.tmp$sample_id),
                                       labels=c("c1","c2","c1","c2",
                                                "c1","c2","c1","c2",
                                                "c1","c2","c1","c2",
                                                "c1","c1",
                                                #"c1","c2","c1",
                                                "c1","c2","c3",
                                                "c1","c2","c1","c2",
                                                "c1","c2","c1","c2","c2",
                                                "c1","c1",
                                                "c1","c2","c1","c2",
                                                "c1","c2","c1","c2")))
srt.tmp$grid_row = as.character(factor(srt.tmp$sample_id, levels=levels(srt.tmp$sample_id),
                                    labels=c("r1","r1","r2","r2",
                                             "r1","r1","r2","r2",
                                             "r1","r1","r2","r2",
                                             "r1","r2",
                                             #"r1","r1","r2",
                                             "r1","r1","r1",
                                             "r1","r1","r2","r2",
                                             "r1","r1","r2","r2","r3",
                                             "r1","r2",
                                             "r1","r1","r2","r2",
                                             "r1","r1","r2","r2")))
srt.tmp
table(srt.tmp$sample_id)

srt.tmp = left_join(srt.tmp, srt.angles)
table(srt.tmp$rotateImg_angle, useNA="ifany")
write.csv(srt.tmp, "processed-data/revision/humanHippocampus2024_SRT_sample-info.csv", row.names=F)


# snRNAseq: sort
colnames(sn.coldata)
head(sn.coldata[,1:25])

sn.tmp = group_by(sn.coldata, Sample, brnum, round, sort) %>% tally(name="n_nuclei") %>%
  mutate(round=as.numeric(as.character(round))) %>%
  arrange(round)
write.csv(sn.tmp, "processed-data/revision/humanHippocampus2024_snRNAseq_sample-info.csv", row.names=F)
