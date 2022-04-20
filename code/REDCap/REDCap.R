#cd /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/

suppressPackageStartupMessages(library("here"))

REDCap = read.csv(file.path(here::here("code","Visium_DATA_2022-04-20_1715.csv")),header=TRUE, stringsAsFactors=FALSE)
A1 = subset(REDCap, select = c("date","slide","experimenter","species_a1","sample_a1","serial_a1","adjacent_a1",
                               "region_a1","project_a1","sample_number1_a1", "master_sheet1_a1","experimenter1_a1"))
B1 = subset(REDCap, select = c("date","slide","experimenter","species_b1","sample_b1","serial_b1", "adjacent_b1",
                               "region_b1","project_b1","sample_number1_b1","master_sheet1_b1","experimenter1_b1"))
C1 = subset(REDCap, select = c("date","slide","experimenter", "species_c1","sample_c1","serial_c1", "adjacent_c1",                                              
                               "region_c1","project_c1","sample_number1_c1","master_sheet1_c1","experimenter1_c1"))
D1 = subset(REDCap, select = c("date","slide","experimenter","species_d1","sample_d1","serial_d1", "adjacent_d1",
                               "region_d1","project_d1","sample_number1_d1","master_sheet1_d1","experimenter1_d1"))
colnames(A1) = colnames(B1) = colnames(C1) = colnames(D1) = c("date","slide","experimenter_img","species","sample","serial","adjacent",
                                                              "region","project","sample_number", "master_sheet","experimenter_seq")
A1$array = "A1"
B1$array = "B1"
C1$array = "C1"
D1$array = "D1"

REDCap_table = rbind(A1,B1,C1,D1)
REDCap_table = REDCap_table[order(REDCap_table$slide),]

REDCap_HPC = REDCap_table[which(REDCap_table$project == "spatialHPC_LIBD4035"),]

Brain_nums = unique(REDCap_HPC$sample)
write.table(Brain_nums, file = (here::here("code","BrainList_for_Geo.txt")),row.names = FALSE, col.names = FALSE)

Samples = unique(paste0(REDCap_HPC$slide,"_",REDCap_HPC$array))
write.table(Samples,file = (here::here("code","samples.txt")),row.names = FALSE, col.names = FALSE)

save(REDCap_HPC, file = (here::here("code","REDCap_HPC.rda")))

##https://jhu-genomics.slack.com/archives/CR9NYA0BF/p1650383126365919
dir.create(here::here("raw-data","FASTQ","2022-04-12"))

samples = subset(REDCap_HPC, select = c("date","slide","array"))
samples$date = as.Date(samples$date)
samples[samples$date > "2021-10-06",]

sapply(here::here("raw-data","FASTQ","2022-04-12", paste(samples$slide,samples$array,sep = "_")), dir.create) 
