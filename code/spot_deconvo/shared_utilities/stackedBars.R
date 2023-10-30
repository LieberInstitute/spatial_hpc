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

#group = "layer" 
group = c("layer", "broad")
tools = c("cell2location","RCTD","tangram")
df_list <- list()

for (grp in group){
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
    if (grp == "layer"){
      inhb = c("GABA.LAMP5","GABA.MGE","GABA.CGE")
      exct = c("Thal","L5","CA2.4","GC","L2_3.PrS.PaS","CA1_ProS","L6_6b","Sub.2","Sub.1","L2_3.Prs.Ent","Cajal")
      csf = c("Choroid","Ependy")
      celltypes = c(inhb, exct, csf, "Micro_Macro_T", "Astro", "Oligo")
      temp_df$InhN = rowSums(temp_df[,inhb])
      temp_df$CSF = rowSums(temp_df[,csf])
      temp_df$ExcN = rowSums(temp_df[,exct])
      }
    
    celltypes = c('InhN', 'ExcN', 'Astro', 'Micro_Macro_T', 'Oligo', 'Vascular', 'OPC', 'CSF')
    temp_df$count = rowSums(temp_df[,celltypes])
    if (tool == "RCTD") {
    temp_df$key = temp_df$X
    temp_df$key <- gsub("Br2720", "V12F14-051", temp_df$key)
      if (grp == "broad"){
      rmv = which(count$key %in% setdiff(count$key,temp_df$key))
      count = count[-rmv,]}
    counts_match = match(count$key, temp_df$key)
    temp_df=temp_df[counts_match,]
    temp_df$count = count$count
    temp_df[, celltypes] <- lapply(temp_df[, celltypes], function(col) col*temp_df$count)
    }
    #temp_df$sample <-sapply(strsplit(temp_df$key,"-1_"), `[`, 2)
    temp_df$group = grp
    temp_df[which(temp_df$count == 0),celltypes] = 0
    temp_df = temp_df[,c("key","sample_id",celltypes,"group")]
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
}
final_df <- do.call(rbind, df_list)
which(is.na(final_df), arr.ind=TRUE)
levels(final_df$cluster)[levels(final_df$cluster)=="WM.1"] <- "WM"
levels(final_df$cluster)[levels(final_df$cluster)=="WM.2"] <- "WM"
levels(final_df$cluster)[levels(final_df$cluster)=="WM.3"] <- "WM"

df <- final_df %>%
  group_by(tool, group, cluster) %>%
  summarize(across(all_of(celltypes), sum)) %>%
  ungroup()

df = as.data.frame(df)
df[, celltypes] <- lapply(df[, celltypes], function(col) col/rowSums(df[,celltypes]))

df = melt(df, id.variable = c("group", "tool", "cluster"))

df$group=gsub('broad', 'Sub', df$group)
df$group=gsub('layer', 'Fine', df$group)
df$group = factor(df$group, levels = c("Sub", "Fine"))

df$tool = gsub('tangram', 'Tangram', df$tool)

colors = c("green2", "green", "darkorange", "yellow", "magenta", "blue1", "blue2", "blue3")
p = ggplot(df, aes(x = group, y = value, fill = variable))+
  geom_bar(stat = "identity") + facet_grid(tool~cluster) + scale_fill_manual(values = colors) +
  theme(text = element_text(size=24, color="black"))+ theme(legend.position="bottom")+labs(fill = "Sub cell type" )+
  guides(fill = guide_legend(nrow = 1))

ggsave(here("plots","spot_deconvo","shared_utilities","stackedBars_layerVsbroad.png"), plot = p,  width = 18, height = 9)
