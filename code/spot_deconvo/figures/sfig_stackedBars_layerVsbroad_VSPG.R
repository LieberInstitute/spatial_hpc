setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")
library("here")
library("sessioninfo")
library("SpatialExperiment")
library('ggplot2')
library('dplyr')
library('gridExtra')
library('reshape2')

load(here("processed-data","06_clustering","PRECAST","spe_norm_with_domain.rda"))
dat = as.data.frame(colData(spe_norm)) %>% select("key", "domain", "brnum")
levels(dat$domain)[levels(dat$domain)=="WM.1"] <- "WM"
levels(dat$domain)[levels(dat$domain)=="WM.2"] <- "WM"
levels(dat$domain)[levels(dat$domain)=="WM.3"] <- "WM"
levels(dat$domain)[levels(dat$domain)=="SLM.WM"] <- "SLM.SGZ"
dat <- dat[dat$brnum %in% c("Br3942_VSPG", "Br8325_VSPG"), ]
dat$key <- sub("_Br.*", "", dat$key)

## CART
Dr = here("processed-data", "spot_deconvo", "groundTruth", "03_CART")
csv_files = data.frame(files = list.files(Dr, pattern = "clusters.csv" , recursive = TRUE, full.names = TRUE))
csv_files$sample_id = sapply(strsplit(csv_files$files,"/"), `[`, 11)
counts_list <- lapply(seq_along(csv_files$files), function(i) 
{data <- read.csv(csv_files$files[i], row.names = NULL) 
data$sample_id <- csv_files$sample_id[i]
return(data)})
temp_df <- do.call(rbind, counts_list)
rownames(temp_df) <- NULL
count = temp_df[,c("key", "n_cells")]

merged_df <- merge(dat, count, by = "key")

#group = "layer" 
group = c("layer", "broad")
tools = c("cell2location","RCTD","tangram")
df_list <- list()

for (grp in group){
  for (tool in tools){
    Dr = here("processed-data", "spot_deconvo", tool, "IF", "2ndRun_newClass", grp)
    if (tool == "RCTD") {Dr = here("processed-data", "spot_deconvo", tool, "2ndRun_newClass_RCTDmarkers", grp)}
    csv_files = data.frame(files = list.files(Dr, pattern = "clusters.csv" , recursive = TRUE, full.names = TRUE))
    if (tool == "RCTD") {csv_files = data.frame(files = csv_files[c(33:40),])
    csv_files$sample_id = sapply(strsplit(csv_files$files,"/"), `[`, 12)
    }else{csv_files$sample_id = sapply(strsplit(csv_files$files,"/"), `[`, 13)}
    counts_list <- lapply(seq_along(csv_files$files), function(i) 
    {data <- read.csv(csv_files$files[i], row.names = NULL) 
    data$sample_id <- csv_files$sample_id[i]
    return(data)})
    temp_df <- do.call(rbind, counts_list)
    if (tool == "RCTD") {temp_df$key = temp_df$X}
	temp_df$key <- sub("_Br.*", "", temp_df$key)
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
    counts_match = match(merged_df$key, temp_df$key)
    temp_df=temp_df[counts_match,]
	temp_df$cluster = merged_df$domain
	
	if (tool == "RCTD"){
	temp_df$count = merged_df$n_cells
	}else{
		for (celltype in celltypes) {
		  temp_df[,celltype] = temp_df[,celltype]/temp_df$count
		}
	}
    #temp_df$sample <-sapply(strsplit(temp_df$key,"-1_"), `[`, 2)
    temp_df$group = grp
    temp_df[which(temp_df$count == 0),celltypes] = 0
    temp_df = temp_df[,c("key","sample_id",celltypes,"group", "cluster")]
    temp_df$tool = tool
    print(paste0(tool, " ", grp))
    which(is.na(temp_df), arr.ind=TRUE)
	df_list[[length(df_list) + 1]] <- temp_df
     }
    #temp_df$cluster = colData(spe)$cluster_collapsed
    
  }

final_df <- do.call(rbind, df_list)
which(is.na(final_df), arr.ind=TRUE)

df <- final_df %>%
  group_by(tool, group, cluster) %>%
  summarize(across(all_of(celltypes), sum)) %>%
  ungroup()

df = as.data.frame(df)
df[, celltypes] <- lapply(df[, celltypes], function(col) col/rowSums(df[,celltypes]))

df = melt(df, id.variable = c("group", "tool", "cluster"))

df$group=gsub('broad', 'mid', df$group)
df$group=gsub('layer', 'fine', df$group)
df$group = factor(df$group, levels = c("mid", "fine"))

df$tool = gsub('tangram', 'Tangram', df$tool)
df$variable = gsub('Micro_Macro_T', 'Micro/Macro/T', df$variable)

load(here("plots","snRNAseq_palettes.rda"))
png(here("plots","spot_deconvo","figures","fig_stackedbars_layerVsbroad", "sfig_stackedbars_layerVsbroad.png"), width = 1300, height = 600, units = "px") 
p = ggplot(df, aes(x = group, y = value, fill = variable))+
  geom_bar(stat = "identity") + facet_grid(tool~cluster) + scale_fill_manual(values = sn.mid.palette) +
  labs(y = "proportion of counts", x="spatial domains", fill = "mid.cell.class")+ 
  theme(text = element_text(size = 20, colour = "black"),
        axis.text = element_text(size = 24, colour = "black"),
        axis.text.x = element_text(angle = 90),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        legend.position = "bottom",
        strip.text = element_text(size = 21, color = "black")) +
        guides(fill = guide_legend(nrow = 1))
  print(p)
dev.off()


rctd_df <- df[df$tool == "RCTD", ]

#png(here("plots", "revision_maddy", "mid_RCTDstackedbars.png"), width = 10, height = 20)

# Create the plot with y-axis labels on the right
p <- ggplot(rctd_df, aes(y = cluster, x = value, fill = variable)) +
    geom_bar(stat = "identity") +
    labs(y = "PRECAST domain", x = "") +
    scale_fill_manual(values = sn.mid.palette) +  # Specify color palette
    theme_minimal() +
    theme(
        axis.title.y = element_text(size = 24, hjust = 0.5, vjust = 0.5),     # Center y-axis title
        axis.text.y.right = element_text(size = 24, colour = "black"),        # Enable right y-axis text
        axis.ticks.y.right = element_line(),     # Enable right y-axis ticks
        axis.line.y.right = element_line()        # Add right y-axis line
    ) +
    scale_y_discrete(position = "right")           # Position y-axis on the right

ggsave(here("plots","revision_maddy", "mid_RCTDstackedbars.png"), plot = p,  width = 10, height = 20)

  
