setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")
library("here")
library("sessioninfo")
library("SpatialExperiment")
library('ggplot2')
library('dplyr')
library('gridExtra')
library('scatterpie')

spg_in <- here("processed-data", "spot_deconvo", "shared_utilities" , "spg.rds")
spg = readRDS(spg_in)
dat = as.data.frame(colData(spg)) %>% select("key", "array_row", "array_col", "sample_id")
dat$key = sapply(strsplit(dat$key,"_Br"), `[`, 1)
celltypes = c("oligo","other","neuron", "microglia", "astrocyte")


## RCTD
tool =  "RCTD"
Dr = here("processed-data", "spot_deconvo",tool,"2ndRun_newClass_RCTDmarkers", "broad")
csv_files = data.frame(files = list.files(Dr, pattern = "clusters.csv" , recursive = TRUE, full.names = TRUE))
csv_files = data.frame(files = csv_files[c(33:40),])
csv_files$sample_id = sapply(strsplit(csv_files$files,"/"), `[`, 12)
counts_list <- lapply(seq_along(csv_files$files), function(i) 
{data <- read.csv(csv_files$files[i], row.names = NULL) 
data$sample_id <- csv_files$sample_id[i]
return(data)})
temp_df <- do.call(rbind, counts_list)
rownames(temp_df) <- NULL

temp_df$key = sapply(strsplit(temp_df$X,"_Br"), `[`, 1)
which(is.na(temp_df), arr.ind=TRUE)

temp_df$neuron = temp_df$ExcN + temp_df$InhN
temp_df$other = temp_df$CSF + temp_df$OPC + temp_df$Vascular
temp_df$microglia = temp_df$Micro_Macro_T
temp_df$astrocyte = temp_df$Astro
temp_df$oligo = temp_df$Oligo
temp_df$none = rowSums(temp_df[,celltypes])
temp_df$tool = tool
counts_match <- match(dat$key, temp_df$key)
temp_df <-temp_df[counts_match,]
which(is.na(temp_df), arr.ind=TRUE)

RCTD <- cbind(dat, temp_df)
RCTD <- RCTD[,c(celltypes,"sample_id", "tool", "array_row", "array_col")]
which(is.na(RCTD), arr.ind=TRUE)

sample_ids = as.character(unique(RCTD$sample_id))
plot_list <- lapply(sample_ids, function(i) {
  dat = RCTD[which(RCTD$sample_id == i),]
  ggplot(dat, aes(array_col, array_row)) +
    geom_scatterpie(aes(x=array_row, y=array_col, r=0.5),
                    data=dat, cols=c("oligo", "other", "neuron", "microglia", "astrocyte"),color=NA) +
    scale_fill_manual(values = c("magenta","blue", "green", "yellow", "orange")) +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "black",colour = "black"))+labs(title = i)
})

ggsave(here("plots","spot_deconvo","shared_utilities","RCTD_scatterpie.pdf"), plot = marrangeGrob(plot_list, nrow=1, ncol=1),  width = 20, height = 20)

sample_id = "V12D07−332_B1"
group = "broad"
## cell2location
tool =  "cell2location"
Dr = here("processed-data", "spot_deconvo",tool, "IF", "2ndRun_newClass", group)
csv_files = data.frame(files = list.files(Dr, pattern = "clusters.csv" , recursive = TRUE, full.names = TRUE))
csv_files$sample_id = sapply(strsplit(csv_files$files,"/"), `[`, 13)
counts_list <- lapply(seq_along(csv_files$files), function(i) 
{data <- read.csv(csv_files$files[i], row.names = NULL) 
data$sample_id <- csv_files$sample_id[i]
return(data)})
temp_df <- do.call(rbind, counts_list)
rownames(temp_df) <- NULL

celltypes1 = c("Oligo","Micro_Macro_T","InhN","ExcN","Astro","Vascular","OPC","CSF")
temp_df$count = rowSums(temp_df[,celltypes1])
temp_df$neuron = temp_df$ExcN + temp_df$InhN
temp_df$other = temp_df$CSF + temp_df$OPC + temp_df$Vascular
temp_df$microglia = temp_df$Micro_Macro_T
temp_df$astrocyte = temp_df$Astro
temp_df$oligo = temp_df$Oligo

temp_df[, celltypes] <- lapply(temp_df[, celltypes], function(col) col/temp_df$count)
temp_df$tool = tool
temp_df$key = sapply(strsplit(temp_df$key,"_Br"), `[`, 1)
counts_match <- match(dat$key, temp_df$key)
temp_df <-temp_df[counts_match,]
which(is.na(temp_df), arr.ind=TRUE)

cell2location <- cbind(dat, temp_df)
cell2location <- cell2location[,c(celltypes,"sample_id", "tool", "array_row", "array_col")]
which(is.na(cell2location), arr.ind=TRUE)


dat1 = cell2location[which(cell2location$sample_id == "V12D07-332_B1"),]
p = ggplot(dat1, aes(array_col, array_row)) +
    geom_scatterpie(aes(x=array_row, y=array_col, r=0.5),
                    data=dat1, cols=c("oligo", "other", "neuron", "microglia", "astrocyte"),color=NA) +
    scale_fill_manual(values = c("magenta","blue", "green", "yellow", "orange")) +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "black",colour = "black"))+labs(title = sample_id)
		 
ggsave(here("plots","revision_maddy","cell2location_scatterpie.png"), plot = p,  width = 20, height = 20)

## tangram
tool =  "tangram"
Dr = here("processed-data", "spot_deconvo",tool, "IF", "2ndRun_newClass", group)
csv_files = data.frame(files = list.files(Dr, pattern = "clusters.csv" , recursive = TRUE, full.names = TRUE))

csv_files$sample_id = sapply(strsplit(csv_files$files,"/"), `[`, 13)
counts_list <- lapply(seq_along(csv_files$files), function(i) 
{data <- read.csv(csv_files$files[i], row.names = NULL) 
data$sample_id <- csv_files$sample_id[i]
return(data)})
temp_df <- do.call(rbind, counts_list)
rownames(temp_df) <- NULL
temp_df$key = substr(temp_df$key, 1, nchar(temp_df$key) - 21)

temp_df$count = rowSums(temp_df[,celltypes1])
temp_df$neuron = temp_df$ExcN + temp_df$InhN
temp_df$other = temp_df$CSF + temp_df$OPC + temp_df$Vascular
temp_df$microglia = temp_df$Micro_Macro_T
temp_df$astrocyte = temp_df$Astro
temp_df$oligo = temp_df$Oligo

temp_df[, celltypes] <- lapply(temp_df[, celltypes], function(col) col/temp_df$count)
which(is.na(temp_df), arr.ind=TRUE)

temp_df[which(temp_df$count == 0),celltypes] = 0
which(is.na(temp_df), arr.ind=TRUE)

temp_df$tool = tool
counts_match <- match(dat$key, temp_df$key)
temp_df <-temp_df[counts_match,]
which(is.na(temp_df), arr.ind=TRUE)

tangram <- cbind(dat, temp_df)
tangram <- tangram[,c(celltypes,"sample_id", "tool", "array_row", "array_col")]
which(is.na(tangram), arr.ind=TRUE)


dat1 = tangram[which(tangram$sample_id == "V12D07-332_B1"),]
p = ggplot(dat1, aes(array_col, array_row)) +
    geom_scatterpie(aes(x=array_row, y=array_col, r=0.5),
                    data=dat1, cols=c("oligo", "other", "neuron", "microglia", "astrocyte"),color=NA) +
    scale_fill_manual(values = c("magenta","blue", "green", "yellow", "orange")) +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_rect(fill = "black",colour = "black"))+labs(title = sample_id)
		 
ggsave(here("plots","revision_maddy","tangram_scatterpie.png"), plot = p,  width = 20, height = 20)
