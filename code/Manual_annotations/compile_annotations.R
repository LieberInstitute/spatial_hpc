setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")
suppressPackageStartupMessages(library("here"))
library(plyr)

files = list.files(path = here::here("processed-data", "manual_annotation_csv"),
                   full.names = TRUE)
files = files[1:15]

##read in files as a list
csvList<-list()
for(i in 1:length(files)){
    csvList[[i]]<-read.csv(files[[i]])
}
names(csvList)<-files

##subset 297 SCP and EDN versions (Erik's B1 used, Stephanie's A1, C1, D1 used)
csvList[[10]]<-csvList[[10]][csvList[[10]]$sample_id=='V11A20-297_B1',]
csvList[[11]]<-[!csvList[[11]]$sample_id=='V11A20-297_B1',]

#compile
csv<-do.call(rbind,csvList)

#save OG compiled list
save(csv, file = here::here("processed-data", "manual_annotation_csv",
                             "compiled_annotation_before_match.Rdata"))

#now reorganize csv based on intersection between csv$spot_name and barcode in spe
## first have to remove repeated sample_id in spe$key
spe$real_key<-gsub('.{14}$', '', spe$key)

#now match csv$spot_name and spe$real_key
csv<-csv[match(csv$spot_name,spe$real_key),]

#now add manual annotation
spe$ManualAnnotation<-csv$ManualAnnotation

#save compiled list following reorg
save(csv, file = here::here("processed-data", "manual_annotation_csv",
                            "compiled_annotation_after_match.Rdata"))

