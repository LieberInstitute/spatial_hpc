setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")
library("here")
library("sessioninfo")
library("SpatialExperiment")
library('ggplot2')
library('dplyr')
library('gridExtra')
library('reshape2')

spe = readRDS(here("processed-data", "spot_deconvo", "shared_utilities", "spe.rds"))
dat = as.data.frame(colData(spe)) %>% select("key", "cluster_collapsed")
count = as.data.frame(colData(spe)) %>% select("key", "count", "cluster_collapsed")

grp = "layer" 
tools = c("cell2location","RCTD","tangram")
df_list <- list()
palette<-c("Choroid"="#4A0404",
           "Vascular"='#A0522D',
           "Ependy"='#D2B48C',
           "Micro_Macro_T"='#71797E',
           'Oligo'='#0000FF',
           'OPC'='#28282B',
           'Astro'='#D3D3D3',
           'GABA.LAMP5'='#40B5AD',
           'GABA.CGE'='#7FFFD4',
           'GABA.MGE'='#AFE1AF',
           'CA2.4'='#F8C8DC',
           'Thal'='#355E3B',
           'GC'='#4CBB17',
           'L5'='#FF7518',
           'L2_3.PrS.PaS'='#FA5F55',
           'CA1_ProS'='#DA70D6',
           'Sub.1'='#FFFF8F',
           'Sub.2'='#FFD700',
           'L6_6b'='#FFAA33',
           'L2_3.Prs.Ent'='#800080',
           'Cajal'='#4682B4'
)

celltypes = names(palette)

  for (tool in tools){
    Dr = here("processed-data", "spot_deconvo", tool, "HE", "2ndRun_newClass", grp)
    if (tool == "RCTD") {Dr = here("processed-data", "spot_deconvo", tool, "2ndRun_newClass_RCTDmarkers", grp)}
    csv_files = data.frame(files = list.files(Dr, pattern = "clusters.csv" , recursive = TRUE, full.names = TRUE))
    if (tool == "RCTD") {csv_files = data.frame(files = csv_files[-c(33:40),])
    csv_files$sample_id = sapply(strsplit(csv_files$files,"/"), `[`, 12)
    }else{csv_files$sample_id = sapply(strsplit(csv_files$files,"/"), `[`, 13)}
    counts_list <- lapply(seq_along(csv_files$files), function(i) 
    {data <- read.csv(csv_files$files[i], row.names = NULL) 
    data$sample_id <- csv_files$sample_id[i]
    return(data)})
    temp_df <- do.call(rbind, counts_list)
    rownames(temp_df) <- NULL
     
    temp_df$count = rowSums(temp_df[,celltypes])
    if (tool == "RCTD") {
      temp_df$key = temp_df$X
      temp_df$key <- gsub("Br2720", "V12F14-051", temp_df$key)
      counts_match = match(count$key, temp_df$key)
      temp_df=temp_df[counts_match,]
      temp_df$count = count$count
      temp_df[, celltypes] <- lapply(temp_df[, celltypes], function(col) col*temp_df$count)
      }
    temp_df[which(temp_df$count == 0),celltypes] = 0
    temp_df = temp_df[,c("key","sample_id",celltypes)]
    temp_df$tool = tool
    if (tool == "tangram"){
      temp_df$key = substr(temp_df$key, 1, nchar(temp_df$key) - 14)
      temp_df$key = gsub("Br2720", "V12F14-051", temp_df$key)}
    if (tool == "RCTD") {
      temp_df$cluster = count$cluster_collapsed
    }else{
      counts_match = match(dat$key, temp_df$key)
      temp_df=temp_df[counts_match,]
      temp_df$cluster = dat$cluster_collapsed
    }
    print(paste0(tool, " ", grp))
    which(is.na(temp_df), arr.ind=TRUE)
    #temp_df$cluster = colData(spe)$cluster_collapsed
    df_list[[length(df_list) + 1]] <- temp_df
  }

final_df <- do.call(rbind, df_list)
which(is.na(final_df), arr.ind=TRUE)
levels(final_df$cluster)[levels(final_df$cluster)=="WM.1"] <- "WM"
levels(final_df$cluster)[levels(final_df$cluster)=="WM.2"] <- "WM"
levels(final_df$cluster)[levels(final_df$cluster)=="WM.3"] <- "WM"

df <- final_df %>%
  group_by(tool, cluster) %>%
  summarize(across(all_of(celltypes), sum)) %>%
  ungroup()

df = as.data.frame(df)
df[, celltypes] <- lapply(df[, celltypes], function(col) col/rowSums(df[,celltypes]))

df = melt(df, id.variable = c("tool", "cluster"))
df$tool = gsub('tangram', 'Tangram', df$tool)


png(here("plots","spot_deconvo","figures","fig_celltypeproportionsINdomains", "sfig_celltypeproportionsINdomains.png"), width = 1200, height = 600, units = "px") 
p = ggplot(df, aes(x = cluster, y = value, fill = variable))+theme_bw() +
  geom_bar(stat = "identity") + facet_grid(~tool) + scale_fill_manual(values = palette) +
  labs(y = "proportion of counts", x="spatial domains", fill = "cell types")+ 
  theme(text = element_text(size = 20, colour = "black"),
        axis.text = element_text(size = 24, colour = "black"),
        axis.text.x = element_text(angle = 90),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        strip.text = element_text(size = 24, color = "black")) +
  guides(fill=guide_legend(ncol =1))
print(p)
dev.off()