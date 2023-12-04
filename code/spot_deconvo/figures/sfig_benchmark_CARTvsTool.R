setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")
library("here")
library("jaffelab")
library("SpatialExperiment")
library("sessioninfo")
library("tidyverse")
library("reshape")

Dr = here("processed-data", "spot_deconvo", "groundTruth", "03_CART")
csv_files = data.frame(files = list.files(Dr, pattern = "clusters.csv" , recursive = TRUE, full.names = TRUE))
counts_list <- lapply(csv_files$files, function(file_path) 
{data <- read.csv(file_path, row.names = NULL)  
return(data)})
CART <- do.call(rbind, counts_list)
CART$sample <-sapply(strsplit(CART$key,"-1_"), `[`, 2)
rownames(CART) <- NULL
CART$tool = "CART"


#grp = "broad"
grp = "layer"
tools <- c("tangram", "cell2location", "RCTD")
df_list <- list()
for (tool in tools){
  Dr = here("processed-data", "spot_deconvo", tool, "IF", "2ndRun_newClass", grp)
  if (tool == "RCTD") {Dr = here("processed-data", "spot_deconvo", tool, "2ndRun_newClass_RCTDmarkers", grp)}
  csv_files = data.frame(files = list.files(Dr, pattern = "clusters.csv" , recursive = TRUE, full.names = TRUE))
  if (tool == "RCTD") {csv_files = data.frame(files = csv_files[33:40,])}
  counts_list <- lapply(csv_files$files, function(file_path) 
  {data <- read.csv(file_path, row.names = NULL)  
  return(data)})
  temp_df <- do.call(rbind, counts_list)
  rownames(temp_df) <- NULL
  if (tool == "RCTD"){temp_df$key = temp_df$X}
  temp_df$key = sapply(strsplit(temp_df$key,"_Br"), `[`, 1)
  if (grp == "broad"){ 
    celltypes = c("Oligo", "Micro_Macro_T", "InhN", "ExcN", "Astro", "Vascular", "OPC", "CSF")
    temp_df$count = rowSums(temp_df[,celltypes])
    temp_df$neuron = temp_df$ExcN + temp_df$InhN
    temp_df$other = temp_df$CSF + temp_df$OPC + temp_df$Vascular
    temp_df$microglia = temp_df$Micro_Macro_T
    temp_df$astrocyte = temp_df$Astro
    temp_df$oligo = temp_df$Oligo
  }else if (grp == "layer"){
    inhb = c("GABA.LAMP5","GABA.MGE","GABA.CGE")
    exct = c("Thal","L5","CA2.4","GC","L2_3.PrS.PaS","CA1_ProS","L6_6b","Sub.2","Sub.1","L2_3.Prs.Ent","Cajal")
    csf = c("Choroid","Ependy")
    celltypes = c(inhb, exct, csf, "Micro_Macro_T", "Astro", "Oligo")
    temp_df$count = rowSums(temp_df[,celltypes])
    temp_df$neuron = rowSums(temp_df[,c(inhb,exct)])
    temp_df$other = rowSums(temp_df[,c(csf)])
    temp_df$microglia = temp_df$Micro_Macro_T
    temp_df$astrocyte = temp_df$Astro
    temp_df$oligo = temp_df$Oligo
  }
  celltypes = c('oligo', 'other', 'neuron', 'microglia', 'astrocyte')
  if (tool == "RCTD"){
    temp_df <- temp_df[temp_df$key %in% CART$key, ]
    cnt = CART$n_cells[CART$key %in%  temp_df$key]
    temp_df[, celltypes] <- lapply(temp_df[, celltypes], function(col) col*cnt)}
  temp_df$sample <-sapply(strsplit(temp_df$key,"-1_"), `[`, 2)
  temp_df[which(temp_df$count == 0),celltypes] = 0
  temp_df = temp_df[,c("key","sample",celltypes)]
  temp_df$tool = tool
  df_list[[length(df_list) + 1]] <- temp_df
}
final_df <- do.call(rbind, df_list)
which(is.na(final_df), arr.ind=TRUE)

RCTD = final_df[which(final_df$tool == "RCTD"),]
rm = setdiff(CART$key, RCTD$key)
CART = CART[!(CART$key %in% rm), ]
final_df = final_df[!(final_df$key %in% rm), ]

CART = CART[,c("key", "sample", "oligo","other", "neuron","microglia","astrocyte","tool")]
df1 = rbind(CART, final_df)
#df1 = melt(df1, id.vars = c("key", "sample", "tool"))

metrics_df = df1 |> group_by(tool, sample) |>
  summarize(
    Oligodendrocyte = sum(oligo),
    Neuron = sum(neuron),
    Microglia = sum(microglia),
    Astrocyte = sum(astrocyte)
  ) 

final_metrics = metrics_df |>
  group_by(tool) |>
  summarize(
    corr = round(mean(corr), 2), rmse = signif(mean(rmse), 3)
  )


# Load the ggplot2 library
library(ggplot2)
library(gridExtra)

metrics_df$tool = gsub('tangram', 'Tangram', metrics_df$tool)
metrics_df = melt(metrics_df, id.vars = c("sample", "tool"))
CART = metrics_df[which(metrics_df$tool == "CART"),]
metrics_df = metrics_df[which(metrics_df$tool != "CART"),]

merged_df <- merge(metrics_df, CART, by = c("sample", "variable"))

png(here("plots","spot_deconvo","figures","fig_benchmark", "sfig_layer_CARTvsTool.png"), width = 1200, height = 400, units = "px") 
p = ggplot(merged_df, aes(x = value.x, y = value.y, color = variable, shape = sample))+geom_point(size = 5)+ theme_bw() +
  facet_wrap(~tool.x, nrow=1)+ scale_color_manual(values = c("green", "darkorange", "yellow", "magenta"))+
  scale_shape_manual(values = c(0, 7, 12, 15, 1, 10, 13, 19)) + 
  geom_abline(intercept = 0, slope = 1, linetype = 3, color = "black") +
  labs( x = "Software estimated", y = "CART estimated", color = "Cell type") + 
  theme(text = element_text(size = 20, colour = "black"),
        axis.text = element_text(size = 24, colour = "black"),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        strip.text = element_text(size = 24, color = "black"))

print(p)
dev.off()
