setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")
library("here")
library("sessioninfo")
library("SpatialExperiment")
library('ggplot2')
library('dplyr')
library('gridExtra')

spe = readRDS(here("processed-data", "spot_deconvo", "shared_utilities", "spe.rds"))
#dat = as.data.frame(colData(spe)) %>% select("key", "cluster_collapsed", "array_row", "array_col", "sample_id")
dat = as.data.frame(colData(spe)) %>% select("key", "broad2", "array_row", "array_col", "sample_id")

#group = "layer" 
#celltypes = c("Astro","CA1_ProS","CA2.4","Cajal","Choroid","Ependy","GABA.CGE","GABA.LAMP5","GABA.MGE"      "GC"            "L2_3.PrS.PaS"  "L2_3.Prs.Ent" 
#"L5","L6_6b","Micro_Macro_T","OPC","Oligo","Sub.1","Sub.2","Thal","Vascular") 

group = "broad"
celltypes = c("Oligo","Micro_Macro_T","InhN","ExcN","Astro","Vascular","OPC","CSF")      

## cell2location
tool =  "cell2location"
Dr = here("processed-data", "spot_deconvo",tool, "HE", "2ndRun_newClass", group)
csv_files = data.frame(files = list.files(Dr, pattern = "clusters.csv" , recursive = TRUE, full.names = TRUE))
csv_files$sample_id = sapply(strsplit(csv_files$files,"/"), `[`, 13)
counts_list <- list()
for (file_path in csv_files$files) {
  data <- read.csv(file_path, row.names = NULL)  # Read the CSV file
  counts_list[[file_path]] <- data  # Store the data frame in the list with the file path as the key
}

temp_df <- do.call(rbind, counts_list)
rownames(temp_df) <- NULL
#celltypes = names(temp_df)[3:length(names(temp_df))]
temp_df$count = rowSums(temp_df[,celltypes])
setdiff(temp_df$key, dat$key)
setdiff(dat$key,temp_df$key)

counts = temp_df %>% select("key")
for (celltype in celltypes) {
  counts[,celltype] = temp_df[,celltype]/temp_df$count
}

counts$tool = tool
counts_match <- match(dat$key, counts$key)
counts_info <-counts[counts_match,]
counts_info <- subset(counts_info, select = -key)
cell2location <- cbind(dat, counts_info)

which(is.na(cell2location), arr.ind=TRUE)

## tangram
tool =  "tangram"
Dr = here("processed-data", "spot_deconvo",tool, "HE", "2ndRun_newClass", group)
csv_files = data.frame(files = list.files(Dr, pattern = "clusters.csv" , recursive = TRUE, full.names = TRUE))

csv_files$sample_id = sapply(strsplit(csv_files$files,"/"), `[`, 13)
counts_list <- list()
for (file_path in csv_files$files) {
  data <- read.csv(file_path, row.names = NULL)  # Read the CSV file
  counts_list[[file_path]] <- data  # Store the data frame in the list with the file path as the key
}

temp_df <- do.call(rbind, counts_list)
rownames(temp_df) <- NULL
temp_df$key = substr(temp_df$key, 1, nchar(temp_df$key) - 14)
#celltypes = names(temp_df)[3:length(names(temp_df))]
temp_df$count = rowSums(temp_df[,celltypes])
counts = temp_df %>% select("key")

for (celltype in celltypes) {
  counts[,celltype] = temp_df[,celltype]/temp_df$count
}

#which(is.na(counts), arr.ind=TRUE)
counts <- replace(counts, is.na(counts), 0)
counts$key <- gsub("Br2720", "V12F14-051", counts$key)
setdiff(counts$key, dat$key)
setdiff(dat$key,counts$key)
which(is.na(counts), arr.ind=TRUE)

counts_match <- match(dat$key, counts$key)
counts_info <-counts[counts_match,]
which(is.na(counts_info), arr.ind=TRUE)

counts_info <- subset(counts_info, select = -key)
tangram <- cbind(dat, counts_info)
tangram$tool = tool

column_order <- colnames(cell2location)
tangram <- tangram[, column_order]
which(is.na(tangram), arr.ind=TRUE)

## RCTD
tool =  "RCTD"
Dr = here("processed-data", "spot_deconvo",tool,"2ndRun_newClass_RCTDmarkers", group)
csv_files = data.frame(files = list.files(Dr, pattern = "clusters.csv" , recursive = TRUE, full.names = TRUE))
csv_files = data.frame(files = csv_files[-c(33:40),])
csv_files$sample_id = sapply(strsplit(csv_files$files,"/"), `[`, 12)
counts_list <- list()
for (file_path in csv_files$files) {
  data <- read.csv(file_path, row.names = NULL)  # Read the CSV file
  counts_list[[file_path]] <- data  # Store the data frame in the list with the file path as the key
}

temp_df <- do.call(rbind, counts_list)
rownames(temp_df) <- NULL
#colnames(temp_df)[colnames(temp_df) == "CA2.4"] = "CA2-4"
colnames(temp_df)[colnames(temp_df) == "X"] = "key"
#celltypes = names(temp_df)[2:length(names(temp_df))]
which(is.na(temp_df), arr.ind=TRUE)
setdiff(temp_df$key, dat$key)

temp_df$key <- gsub("Br2720", "V12F14-051", temp_df$key)
setdiff(temp_df$key, dat$key)
setdiff(dat$key,temp_df$key)

rmv = which(dat$key %in% setdiff(dat$key,temp_df$key))
dat = dat[-rmv,]

#temp_df$slide = sapply(strsplit(temp_df$key,"_"), `[`, 2)
counts_match <- match(dat$key, temp_df$key)
counts_info <-temp_df[counts_match,]
which(is.na(counts_info), arr.ind=TRUE)

counts_info$tool = tool
counts_info <- subset(counts_info, select = -key)
RCTD <- cbind(dat, counts_info)

column_order <- colnames(cell2location)
RCTD <- RCTD[, column_order]
which(is.na(RCTD), arr.ind=TRUE)

cell2location = cell2location[-rmv,]
tangram = tangram[-rmv,]
dat1 = rbind(cell2location, tangram, RCTD)
# rmv = c("V11L05-335_B1", "V11U08-081_A1", "V11U08-081_C1", "V11U08-081_D1", "V11U08-084_D1")
# dat1 = dat1[!dat1$sample_id %in% rmv, ]  
which(is.na(dat1), arr.ind=TRUE)

## for box plots
plot_list <- lapply(celltypes, function(i){
  #ggplot(dat1, aes(x = dat1$cluster_collapsed, y = dat1[,i])) + geom_boxplot(aes(fill=tool), outlier.shape = NA)+labs(title = i)
  ggplot(dat1, aes(x = dat1$broad2, y = dat1[,i])) + geom_boxplot(aes(fill=tool), outlier.shape = NA)+labs(title = i)
})

ggsave(here("plots","spot_deconvo","shared_utilities",group,"boxplot.pdf"), plot = marrangeGrob(plot_list, nrow=1, ncol=1),  width = 24, height = 8)

## proportion plots
#celltypes = colnames(dat1)[6:26]
library(viridis)

for (sample_id in unique(dat1$sample_id)){
  datb = dat1[which(dat1$sample_id == sample_id), ]
  plot_list <- lapply(celltypes, function(i){
    ggplot(data = datb, aes(x=array_row, y=array_col, color = datb[,i]))+
      geom_point(size = 3)+facet_wrap(~tool)+scale_color_gradientn(colours = viridis(10, option = "magma"))+
      theme(legend.position = "none")+labs(title = i)
  })
  #gridplot = grid.arrange(grobs = plot_list, nrow = length(celltypes))
  ggsave(here("plots","spot_deconvo","shared_utilities",group,paste0(sample_id,".pdf")), plot = marrangeGrob(plot_list, nrow=1, ncol=1),  width = 24, height = 8)
  print(paste0("done ", sample_id))
}

