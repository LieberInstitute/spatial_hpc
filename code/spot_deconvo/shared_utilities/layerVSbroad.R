setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")
library("here")
library("jaffelab")
library("SpatialExperiment")
library("sessioninfo")
library("tidyverse")
library("reshape")

spg_in <- here("processed-data", "spot_deconvo", "shared_utilities", "spg.rds")
spg = readRDS(spg_in)
# unqiue(spg$brnum)
# Br3942_VSPG Br8325_VSPG
colData(spg)$key = gsub("_Br3942", "", colData(spg)$key)
colData(spg)$key = gsub("_Br8325", "", colData(spg)$key)

group = c("broad", "layer")
tools = c("cell2location", "tangram")
df_list <- list()
for (grp in group){
  
  for (tool in tools){
    Dr = here("processed-data", "spot_deconvo", tool, "IF", "2ndRun_newClass", grp)
    csv_files = data.frame(files = list.files(Dr, pattern = "clusters.csv" , recursive = TRUE, full.names = TRUE))
    counts_list <- lapply(csv_files$files, function(file_path) 
    {data <- read.csv(file_path, row.names = NULL)  
    return(data)})
    temp_df <- do.call(rbind, counts_list)
    rownames(temp_df) <- NULL
    temp_df$key = sapply(strsplit(temp_df$key,"_Br"), `[`, 1)
    
    if (grp == "broad"){
      celltypes = c("Oligo", "Micro_Macro_T", "InhN", "ExcN", "Astro", "Vascular", "OPC", "CSF")
      temp_df$count = rowSums(temp_df[,celltypes])
      temp_df$neuron = temp_df$ExcN + temp_df$InhN
      temp_df$other = temp_df$CSF + temp_df$OPC + temp_df$Vascular
      temp_df$microglia = temp_df$Micro_Macro_T
      temp_df$astrocyte = temp_df$Astro
      temp_df$oligo = temp_df$Oligo
    }else{
      neurons = c("GABA.LAMP5","GABA.MGE","GABA.CGE","Thal","L5","CA2.4","GC","L2_3.PrS.PaS",
                  "CA1_ProS","L6_6b","Sub.2","Sub.1","L2_3.Prs.Ent","Cajal")
      others= c("Choroid","Ependy","Vascular","OPC")
      celltypes = c(neurons, others, "Micro_Macro_T", "Astro", "Oligo")
      temp_df$count = rowSums(temp_df[,celltypes])
      temp_df$neuron = rowSums(temp_df[,neurons])
      temp_df$other = rowSums(temp_df[,others])
      temp_df$microglia = temp_df$Micro_Macro_T
      temp_df$astrocyte = temp_df$Astro
      temp_df$oligo = temp_df$Oligo
    }
    
    celltypes = c('oligo', 'other', 'neuron', 'microglia', 'astrocyte')
    #  temp_df[, celltypes] <- lapply(temp_df[, celltypes], function(col) col / temp_df$count)
    temp_df$sample <-sapply(strsplit(temp_df$key,"-1_"), `[`, 2)
    temp_df$group = grp
    temp_df[which(temp_df$count == 0),celltypes] = 0
    temp_df = temp_df[,c("key","sample",celltypes,"group")]
    temp_df$tool = tool
    df_list[[length(df_list) + 1]] <- temp_df
  }
}
final_df <- do.call(rbind, df_list)
which(is.na(final_df), arr.ind=TRUE)

# Melt the dataframe to long format
# id_vars <- c("key", "sample", "group", "tool")
# final_df <- melt(df, id.vars = id_vars, measure.vars = c("oligo", "other", "neuron", "microglia", "astrocyte"),
#                       variable.name = "celltype", value.name = "proportion")
# colnames(final_df) <- c(id_vars, "celltype", "proportion")
# which(is.na(final_df), arr.ind=TRUE)

#### gather Counts results
Dr = here("processed-data", "spot_deconvo", "groundTruth", "03_CART")
csv_files = data.frame(files = list.files(Dr, pattern = "clusters.csv" , recursive = TRUE, full.names = TRUE))
counts_list <- lapply(csv_files$files, function(file_path) 
{data <- read.csv(file_path, row.names = NULL)  
return(data)})
temp_df <- do.call(rbind, counts_list)
rownames(temp_df) <- NULL
count = temp_df[,c("key","n_cells")]

## RCTD ####
tool =  "RCTD"
df = data.frame()
for (grp in group){
  Dr = here("processed-data", "spot_deconvo", tool, "2ndRun_newClass_RCTDmarkers", grp)
  csv_files = data.frame(files = list.files(Dr, pattern = "clusters.csv" , recursive = TRUE, full.names = TRUE))
  csv_files = csv_files[33:40,]
  counts_list <- lapply(csv_files, function(file_path) 
  {data <- read.csv(file_path, row.names = NULL)  
  return(data)})
  temp_df <- do.call(rbind, counts_list)
  rownames(temp_df) <- NULL
  colnames(temp_df)[colnames(temp_df) == "X"] = "key"
  temp_df$key = sapply(strsplit(temp_df$key,"_Br"), `[`, 1)
  
  if (grp == "broad"){
    temp_df$neuron = temp_df$ExcN + temp_df$InhN
    temp_df$other = temp_df$CSF + temp_df$OPC + temp_df$Vascular
    temp_df$microglia = temp_df$Micro_Macro_T
    temp_df$astrocyte = temp_df$Astro
    temp_df$oligo = temp_df$Oligo
  }else{
    neurons = c("GABA.LAMP5","GABA.MGE","GABA.CGE","Thal","L5","CA2.4","GC","L2_3.PrS.PaS",
                "CA1_ProS","L6_6b","Sub.2","Sub.1","L2_3.Prs.Ent","Cajal")
    others= c("Choroid","Ependy","Vascular","OPC")
    temp_df$neuron = rowSums(temp_df[,neurons])
    temp_df$other = rowSums(temp_df[,others])
    temp_df$microglia = temp_df$Micro_Macro_T
    temp_df$astrocyte = temp_df$Astro
    temp_df$oligo = temp_df$Oligo
  }
  celltypes = c('oligo', 'other', 'neuron', 'microglia', 'astrocyte')
  counts_match = match(count$key, temp_df$key)
  temp_df=temp_df[counts_match,]
  temp_df$count = count$n_cells
  temp_df = temp_df[complete.cases(temp_df), ]
  temp_df[, celltypes] <- lapply(temp_df[, celltypes], function(col) col*temp_df$count)
  temp_df$group = grp
  temp_df$sample <-sapply(strsplit(temp_df$key,"-1_"), `[`, 2)
  temp_df <- temp_df[, c('key','sample',celltypes,'group')]
  which(is.na(temp_df), arr.ind=TRUE)
  df = rbind(df,temp_df)
}
df$tool = tool

# id_vars <- c("key", "sample", "group", "tool")
# RCTD <- melt(df, id.vars = id_vars, measure.vars = c("oligo", "other", "neuron", "microglia", "astrocyte"),
#                  variable.name = "celltype", value.name = "proportion")
# colnames(RCTD) <- c(id_vars, "celltype", "proportion")
# which(is.na(RCTD), arr.ind=TRUE)

final_df = rbind(final_df, df)
df1 <- melt(final_df, id.vars = c("key", "sample", "group", "tool") , measure.vars = c("oligo", "other", "neuron", "microglia", "astrocyte"))
df = inner_join(df1[which(df1$group=="layer"),], df1[which(df1$group=="broad"),], by = c("key", "sample", "tool", "variable"))
# Calculate correlations per celltype, sample, and tool
corr_df <- df %>% group_by(variable, sample, tool) %>%
  summarise(correlation = round(cor(value.x, value.y),2),
            rmse = signif(sqrt(mean((value.x - value.y)^2)), 3))

corr_df <- df %>% group_by(tool,sample) %>%
  summarise(corr = round(cor(value.x, value.y),2),
            rmse = signif(sqrt(mean((value.x - value.y)^2)), 3)
) |>
  group_by(tool) |>
  summarize(
    corr = round(mean(corr), 2), rmse = signif(mean(rmse), 3)
  )

corr_df$corr <- paste("Cor =", corr_df$corr)
corr_df$rmse <- paste("RMSE =", corr_df$rmse)

sum_df <- final_df %>%
  group_by(sample, group, tool) %>%
  summarise(Neuron = sum(neuron), Astrocyte = sum(astrocyte), Microglia = sum(microglia), Oligodendrocyte = sum(oligo))
sum_long <- melt(as.data.frame(sum_df), id.vars = c("sample", "group","tool"), variable.name = "celltype", value.name = "totals")
df1 = inner_join(sum_long[which(sum_long$group=="layer"),], sum_long[which(sum_long$group=="broad"),], by = c("sample", "tool", "variable"))
df1$tool = gsub('tangram', 'Tangram', df1$tool)
corr_df$tool = gsub('tangram', 'Tangram', corr_df$tool)

p = ggplot(df1, aes(x = log(value.y), y = log(value.x), shape = sample, color = variable)) + geom_point(size = 5) +
  facet_wrap(~tool)+ geom_abline(intercept = 0, slope = 1, linetype = 3, color = "black")+
  labs(x = "log(sub counts)", y = "log(fine counts)")+ scale_color_manual(values = c("green", "darkorange", "yellow", "magenta"))+
  scale_shape_manual(values = c(0, 7, 12, 15, 1, 10, 13, 19))+ theme(legend.position="right") +
  guides(shape = guide_legend(ncol = 1, title = "Sample"), color = guide_legend(ncol = 1, title = "Broad cell type")) +
  theme(text = element_text(size=24, color="black"))

ggsave(here("plots","spot_deconvo","shared_utilities","broadVSlayerPoster.png"), plot = p,  width = 18, height = 6)
