setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")
library("here")
library("jaffelab")
library("SpatialExperiment")
library("sessioninfo")
library("tidyverse")
library("reshape")

#grp = "broad"
grp = "layer"
tools <- c("CART", "tangram", "cell2location", "RCTD")
df_list <- list()
for (tool in tools){
  Dr = here("processed-data", "spot_deconvo", tool, "IF", "2ndRun_newClass", grp)
  if (tool == "CART") {Dr = here("processed-data", "spot_deconvo", "groundTruth", "03_CART")}
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
  if (tool == "CART"){
    temp_df$count = temp_df$n_cells
  }else if (grp == "broad"){ 
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
  if (tool == "RCTD"){temp_df = temp_df}else{
    temp_df[, celltypes] <- lapply(temp_df[, celltypes], function(col) col / temp_df$count)}
  temp_df$sample <-sapply(strsplit(temp_df$key,"-1_"), `[`, 2)
  temp_df[which(temp_df$count == 0),celltypes] = 0
  temp_df = temp_df[,c("key","sample",celltypes)]
  temp_df$tool = tool
  df_list[[length(df_list) + 1]] <- temp_df
}
final_df <- do.call(rbind, df_list)
which(is.na(final_df), arr.ind=TRUE)

df = melt(final_df, id.variables = c("key", "sample_id", "tool"))

CART = df[which(df$tool == "CART"),]
RCTD = df[which(df$tool == "RCTD"),]
tangram = df[which(df$tool == "tangram"),]
cell2location = df[which(df$tool == "cell2location"),]
rm = setdiff(CART$key, RCTD$key)
CART = CART[!(CART$key %in% rm), ]

sorted_keys <- CART$key
RCTD <- RCTD[RCTD$key %in% sorted_keys, ]
RCTD$actual = CART$value
cell2location <- cell2location[cell2location$key %in% sorted_keys, ]
cell2location$actual = CART$value
tangram <- tangram[tangram$key %in% sorted_keys, ]
tangram$actual = CART$value

# CART = final_df[which(final_df$tool == "CART"),]
# RCTD = final_df[which(final_df$tool == "RCTD"),]
# tangram = final_df[which(final_df$tool == "tangram"),]
# cell2location = final_df[which(final_df$tool == "cell2location"),]
# rm = setdiff(CART$key, RCTD$key)
# CART = CART[!(CART$key %in% rm), ]
# 
# sorted_keys <- CART$key
# RCTD <- RCTD[RCTD$key %in% sorted_keys, ]
# cell2location <- cell2location[cell2location$key %in% sorted_keys, ]
# tangram <- tangram[tangram$key %in% sorted_keys, ]

df1 = rbind(RCTD, cell2location, tangram)

metrics_df = df1 |> filter(variable != "other")|>
  group_by(tool, variable, sample) |>
  summarize(
    corr = round(cor(value, actual), 2),
    rmse = signif(mean((value - actual)**2)**0.5, 3)
  ) 

final_metrics = metrics_df |>
  group_by(tool) |>
  summarize(
    corr = round(mean(corr), 2), rmse = signif(mean(rmse), 3)
  )

final_metrics$corr <- paste("Avg. Cor =", final_metrics$corr)
final_metrics$rmse <- paste("Avg. RMSE =", final_metrics$rmse)
final_metrics$tool = gsub("tangram", "Tangram", final_metrics$tool)

# Load the ggplot2 library
library(ggplot2)
library(gridExtra)

metrics_df$variable = as.character(metrics_df$variable)
metrics_df$variable = gsub('oligo', 'Oligodendrocyte', metrics_df$variable)
metrics_df$variable = gsub('neuron', 'Neuron', metrics_df$variable)
metrics_df$variable = gsub('microglia', 'Microglia', metrics_df$variable)
metrics_df$variable = gsub('astrocyte', 'Astrocyte', metrics_df$variable)

metrics_df$tool = gsub('tangram', 'Tangram', metrics_df$tool)

load(here("plots","snRNAseq_palettes.rda"))
names(sn.broad.palette) = c("Neuron", "Microglia", "Astrocyte", "Oligodendrocyte", "Other")
sn.broad.palette["Other"] = "#0000f4"

png(here("plots","spot_deconvo","figures","fig_benchmark", paste0("sfig_",grp,"_RMSEvsCorr.png")), width = 1200, height = 400, units = "px") 
p = ggplot(metrics_df, aes(x = corr, y = rmse))+
  geom_point(size = 5, aes(color = variable, shape = sample))+ theme_bw() +
  facet_wrap(~tool, nrow=1)+ 
  scale_color_manual(values = sn.broad.palette) +
  scale_shape_manual(values = c(0, 7, 12, 15, 1, 10, 13, 19)) + 
  labs( x = "Correlation", y = "RMSE", color = "Cell type") + 
  theme(text = element_text(size = 30, colour = "black"),
        axis.text = element_text(size = 24, colour = "black"),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        strip.text = element_text(size = 24, color = "black")) + 
  geom_text(data = as.data.frame(final_metrics), label = final_metrics$corr, x= 0.015, y = 0.35, size = 8) +
  geom_text(data = as.data.frame(final_metrics), label = final_metrics$rmse, x= 0.015, y = 0.3, size = 8)

#  geom_text(data = as.data.frame(final_metrics), label = final_metrics$corr, x= 0.01, y = 0.65, size = 8) +
#  geom_text(data = as.data.frame(final_metrics), label = final_metrics$rmse, x= 0.01, y = 0.6, size = 8)


print(p)
dev.off()
