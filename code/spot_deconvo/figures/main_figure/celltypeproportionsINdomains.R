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

celltypes = c("Choroid","Vascular","Ependy","Micro_Macro_T",'Oligo',
               'OPC','Astro','GABA.LAMP5','GABA.CGE','GABA.MGE','CA2.4',
               'Thal','GC','L5','L2_3.PrS.PaS','CA1_ProS','Sub.1','Sub.2',
               'L6_6b','L2_3.Prs.Ent','Cajal')

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
levels(final_df$cluster)[levels(final_df$cluster)=="SLM.WM"] <- "SLM.SGZ"

df <- final_df %>%
  group_by(tool, cluster) %>%
  summarize(across(all_of(celltypes), sum)) %>%
  ungroup()

celltypes = c("Choroid", "Vascular","Ependy", "Micro/Macro/T", "Oligo", "OPC",
              "Astro", "GABA.LAMP5", "GABA.CGE", "GABA.MGE", "CA2.4",
              "Thal", "GC", "L5", "L2/3.PrS.PaS", "CA1/ProS", "Sub.1", "Sub.2",
              "L6/6b", "L2/3.Prs.Ent", "Cajal" )
names(df) = c("tool", "cluster", celltypes)
df = as.data.frame(df)
df[, celltypes] <- lapply(df[, celltypes], function(col) col/rowSums(df[,celltypes]))

df = melt(df, id.variable = c("tool", "cluster"))
df$tool = gsub('tangram', 'Tangram', df$tool)
load(here("plots","snRNAseq_palettes.rda"))

df$variable = gsub('CA2.4', 'CA2-4', df$variable)
df$variable = gsub('L5', 'L5/6', df$variable)
df$variable = gsub('L2/3.Prs.Ent', 'L2/3.PrS.Ent', df$variable)

png(here("plots","spot_deconvo","figures","main_figure", "celltypeproportionsINdomains.png"), width = 1200, height = 1000, units = "px") 
p = ggplot(df, aes(x = cluster, y = value, fill = variable ))+theme_bw() +
  geom_bar(stat = "identity") + facet_wrap(~tool, ncol=1, strip.position="right") + 
  scale_fill_manual(values = sn.fine.palette) +
  labs(fill = "fine.cell.class", y = "", x = "")+ 
  theme(text = element_text(size = 36, colour = "black"),
        axis.text = element_text(size = 24, colour = "black"),
        #axis.text.x = element_text(angle = 90),
        axis.text.x = element_blank(),
        #axis.ticks.x=element_blank(),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        strip.text = element_text(size = 32, color = "black"),
        legend.position = "bottom",
        legend.text = element_text(size = 24, colour = "black"),
        legend.title.align=0.5, plot.title = element_text(hjust = 0.5)) +
        guides(fill=guide_legend(nrow =3,title.position="top"))+
  scale_x_discrete(position = "top") 
print(p)
dev.off()