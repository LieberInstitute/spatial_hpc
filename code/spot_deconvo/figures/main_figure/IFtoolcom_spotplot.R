setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")
library("here")
library("sessioninfo")
library("SpatialExperiment")
library('ggplot2')
library('dplyr')
library('gridExtra')

spg_in <- here("processed-data", "spot_deconvo", "shared_utilities" , "spg.rds")
spg = readRDS(spg_in)
dat = as.data.frame(colData(spg)) %>% select("key", "array_row", "array_col", "sample_id")
dat$key = sapply(strsplit(dat$key,"_Br"), `[`, 1)

group = "layer"
celltypes = c("Astro","CA1_ProS","CA2.4","Cajal","Choroid","Ependy","GABA.CGE","GABA.LAMP5","GABA.MGE","GC","L2_3.PrS.PaS","L2_3.Prs.Ent",
               "L5","L6_6b","Micro_Macro_T","OPC","Oligo","Sub.1","Sub.2","Thal","Vascular")

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

temp_df$count = rowSums(temp_df[,celltypes])
temp_df[, celltypes] <- lapply(temp_df[, celltypes], function(col) col/temp_df$count)
temp_df$tool = tool
temp_df$key = sapply(strsplit(temp_df$key,"_Br"), `[`, 1)
counts_match <- match(dat$key, temp_df$key)
temp_df <-temp_df[counts_match,]
which(is.na(temp_df), arr.ind=TRUE)

cell2location <- cbind(dat, temp_df)
cell2location <- cell2location[,c(celltypes,"sample_id", "tool", "array_row", "array_col")]
which(is.na(cell2location), arr.ind=TRUE)

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

temp_df$count = rowSums(temp_df[,celltypes])
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

## RCTD
tool =  "RCTD"
Dr = here("processed-data", "spot_deconvo",tool,"2ndRun_newClass_RCTDmarkers", group)
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
temp_df$tool = tool
counts_match <- match(dat$key, temp_df$key)
temp_df <-temp_df[counts_match,]
which(is.na(temp_df), arr.ind=TRUE)

RCTD <- cbind(dat, temp_df)
RCTD <- RCTD[,c(celltypes,"sample_id", "tool", "array_row", "array_col")]
which(is.na(RCTD), arr.ind=TRUE)

dat1 = rbind(cell2location, tangram, RCTD)
which(is.na(dat1), arr.ind=TRUE)

na_color = "#CCCCCC40"
f_names = as_labeller(c('cell2location' = "(iii) cell2location",
                        'RCTD' = "RCTD",
                        'tangram' = "Tangram"))
## proportion plots
#celltypes = colnames(dat1)[6:26]
datb = dat1[which(dat1$sample_id == "V12D07-332_B1"), ]
datb[, celltypes][datb[, celltypes] <= 0.05] <- NA
#datb$tool = gsub('tangram', 'Tangram', datb$tool)
png(here("plots","spot_deconvo","figures","main_figure", "toolcomp.png"), width = 1800, height = 600, units = "px") 
p = ggplot(data = datb, aes(x=array_row, y=array_col, color = GC))+
  geom_point(size = 2.3)+facet_wrap(~tool, labeller = f_names)+theme_void() +
  #scale_color_gradientn(colours = viridisLite::plasma(21), na.value = na_color)+
  scale_color_distiller(type = "seq",palette = rev('Greys'),direction=1, na.value = na_color)+
  labs(color = "min>0.05", x="", y = "")+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        axis.text = element_blank(),
        axis.ticks =element_blank(),
        plot.title = element_text(size = 48, colour = "black", hjust = 0.5),
        legend.text = element_text(size = 12, colour = "black"),
        legend.justification="right",
        legend.box.spacing = unit(-15, "pt"),
        plot.margin = unit(c(0,0,-10,-10), "pt"),
        axis.ticks.length = unit(0, "pt"),
        panel.spacing.x = unit(-25, "pt"),
        strip.text = element_text(size = 52, colour = "black", hjust = 0.5))
print(p)
dev.off()