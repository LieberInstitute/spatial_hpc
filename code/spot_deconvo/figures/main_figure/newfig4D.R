setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")
library("here")
library("sessioninfo")
library("SpatialExperiment")
library('ggplot2')
library('dplyr')
library('gridExtra')

load(here("processed-data","06_clustering","PRECAST","spe_norm_with_domain.rda"))
dat = as.data.frame(colData(spe_norm)) %>% select("key", "domain", "array_row", "array_col", "sample_id", "brnum")
levels(dat$domain)[levels(dat$domain)=="WM.1"] <- "WM"
levels(dat$domain)[levels(dat$domain)=="WM.2"] <- "WM"
levels(dat$domain)[levels(dat$domain)=="WM.3"] <- "WM"
levels(dat$domain)[levels(dat$domain)=="SLM.WM"] <- "SLM.SGZ"
dat <- dat[dat$brnum %in% c("Br3942_VSPG", "Br8325_VSPG"), ]
#colors = load(here("plots","palettes.rda"))
#dat = as.data.frame(colData(spe)) %>% select("key", "broad2", "array_row", "array_col", "sample_id")

group = "layer"
celltypes = c("Astro","CA1_ProS","CA2.4","Cajal","Choroid","Ependy","GABA.CGE","GABA.LAMP5","GABA.MGE","GC","L2_3.PrS.PaS","L2_3.Prs.Ent",
              "L5","L6_6b","Micro_Macro_T","OPC","Oligo","Sub.1","Sub.2","Thal","Vascular")

# group = "broad"
# celltypes = c("Oligo","Micro_Macro_T","InhN","ExcN","Astro","Vascular","OPC","CSF")      

## cell2location
tool =  "cell2location"
Dr = here("processed-data", "spot_deconvo",tool, "IF", "2ndRun_newClass", group)
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
Dr = here("processed-data", "spot_deconvo",tool, "IF", "2ndRun_newClass", group)
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
csv_files = data.frame(files = csv_files[c(33:40),])
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

counts_match <- match(dat$key, temp_df$key)
counts_info <-temp_df[counts_match,]
which(is.na(counts_info), arr.ind=TRUE)

counts_info$tool = tool
counts_info <- subset(counts_info, select = -key)
RCTD <- cbind(dat, counts_info)

column_order <- colnames(cell2location)
RCTD <- RCTD[, column_order]
which(is.na(RCTD), arr.ind=TRUE)

dat1 = rbind(cell2location, tangram, RCTD)
# rmv = c("V11L05-335_B1", "V11U08-081_A1", "V11U08-081_C1", "V11U08-081_D1", "V11U08-084_D1")
# dat1 = dat1[!dat1$sample_id %in% rmv, ]  
which(is.na(dat1), arr.ind=TRUE)

dat2 = dat1
colnames(dat2)[28] <- "Tool"
dat2$Tool=gsub('tangram', 'Tangram', dat2$Tool)

png(here("plots","spot_deconvo","figures","main_figure", "GCboxplot.png"), width = 1150, height = 450, units = "px") 
p =  ggplot(dat2, aes(x = dat2$domain, y = dat2[,"GC"])) + 
    geom_boxplot(aes(fill=Tool), outlier.shape = NA)+theme_bw()+
    #theme_minimal()+
    labs(title = "Predicted % cell types in PRECAST spatial domains", x="", y="")+
    scale_fill_manual(values = c('grey40', 'black','white'))+
    # theme(text = element_text(colour='black'),
    #       legend.position = c(0.9, 0.75),
    #       panel.grid.minor = element_blank(), 
    #       panel.grid.major = element_blank()) +
  
    theme(text = element_text(size = 36, colour = "black"),
      axis.text = element_text(size = 24, colour = "black"),
      axis.text.x = element_text(angle = 90),
      panel.grid.minor = element_blank(), 
      panel.grid.major = element_blank(),
      legend.text = element_text(size = 20, colour = "black"),
      legend.position = c(0.9, 0.75))
print(p)
dev.off()

