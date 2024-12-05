setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")
suppressPackageStartupMessages(library("SingleCellExperiment"))
suppressPackageStartupMessages(library("here"))
suppressPackageStartupMessages(library("spatialLIBD"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("escheR"))
library(readr)
library(tidyverse)

Dr <- here("processed-data","spot_deconvo","shared_utilities")
spe = readRDS(here(Dr,"spe.rds"))
spe_colData = as.data.frame(colData(spe))
spe_colData$platform = "HE"
## correlation	  
load(here("processed-data", "06_clustering","PRECAST","spe_norm_with_domain.rda"))
spg = spe_norm[,which(spe_norm$slide == "V12D07-332" | spe_norm$slide == "V12D07-335")]
colData_df = as.data.frame(colData(spg))
colData_df$key = sapply(strsplit(colData_df$key,"_Br"), `[`, 1)
nmf_df = read_csv(here("processed-data", "VSPG_image_stitching", "nmf.csv"),col_names = TRUE,show_col_types = FALSE)
spg_colData <- merge(colData_df,nmf_df, by = "key", all.x = TRUE)
spg_colData$platform = "VSPG"

common_cols <- intersect(names(spg_colData), names(spe_colData))
spe_colData <- spe_colData[common_cols]
spg_colData <- spg_colData[common_cols]

combined_df = rbind(spe_colData, spg_colData)
## compare to RCTD
Dr = here("processed-data", "spot_deconvo", "RCTD", "2ndRun_newClass_RCTDmarkers", "layer")
csv_files = data.frame(files = list.files(Dr, pattern = "clusters.csv" , recursive = TRUE, full.names = TRUE))
csv_files = data.frame(files = csv_files[1:40,])
counts_list <- lapply(csv_files$files, function(file_path) 
    {data <- read.csv(file_path, row.names = NULL)  
    return(data)})
temp_df <- do.call(rbind, counts_list)
rownames(temp_df) <- NULL
temp_df$key = temp_df$X
temp_df$key = sapply(strsplit(temp_df$key,"_Br"), `[`, 1)

merged_df <- merge(combined_df,temp_df, by = "key")
RCTD_nmf = merged_df

##Oligo
cor(RCTD_nmf$nmf77, RCTD_nmf$Oligo)
[1] 0.6219841
cor(RCTD_nmf$nmf42, RCTD_nmf$Oligo)
[1] 0.1596806
cor(RCTD_nmf$nmf44, RCTD_nmf$Oligo)
[1] 0.3950035
cor(RCTD_nmf$nmf38, RCTD_nmf$Oligo)
[1] 0.3398897

cor(RCTD_nmf$nmf77[RCTD_nmf$platform == "HE"], RCTD_nmf$Oligo[RCTD_nmf$platform == "HE"], use = "complete.obs")
[1] 0.9441788
cor(RCTD_nmf$nmf42[RCTD_nmf$platform == "HE"], RCTD_nmf$Oligo[RCTD_nmf$platform == "HE"], use = "complete.obs")
[1] 0.2036639
cor(RCTD_nmf$nmf44[RCTD_nmf$platform == "HE"], RCTD_nmf$Oligo[RCTD_nmf$platform == "HE"], use = "complete.obs")
[1] 0.5036317
cor(RCTD_nmf$nmf38[RCTD_nmf$platform == "HE"], RCTD_nmf$Oligo[RCTD_nmf$platform == "HE"], use = "complete.obs")
[1] 0.4283612

cor(RCTD_nmf$nmf77[RCTD_nmf$platform == "VSPG"], RCTD_nmf$Oligo[RCTD_nmf$platform == "VSPG"], use = "complete.obs")
[1] 0.9399053
cor(RCTD_nmf$nmf42[RCTD_nmf$platform == "VSPG"], RCTD_nmf$Oligo[RCTD_nmf$platform == "VSPG"], use = "complete.obs")
[1] 0.1638104
cor(RCTD_nmf$nmf44[RCTD_nmf$platform == "VSPG"], RCTD_nmf$Oligo[RCTD_nmf$platform == "VSPG"], use = "complete.obs")
[1] 0.8763121
cor(RCTD_nmf$nmf38[RCTD_nmf$platform == "VSPG"], RCTD_nmf$Oligo[RCTD_nmf$platform == "VSPG"], use = "complete.obs")
[1] 0.3694046

##Astro
cor(RCTD_nmf$nmf79, RCTD_nmf$Astro)
[1] 0.1241429
cor(RCTD_nmf$nmf81, RCTD_nmf$Astro)
[1] 0.6442813

cor(RCTD_nmf$nmf79[RCTD_nmf$platform == "HE"], RCTD_nmf$Astro[RCTD_nmf$platform == "HE"], use = "complete.obs")
[1] 0.1435998
cor(RCTD_nmf$nmf81[RCTD_nmf$platform == "HE"], RCTD_nmf$Astro[RCTD_nmf$platform == "HE"], use = "complete.obs")
[1] 0.8849263

cor(RCTD_nmf$nmf79[RCTD_nmf$platform == "VSPG"], RCTD_nmf$Astro[RCTD_nmf$platform == "VSPG"], use = "complete.obs")
[1] 0.4254133
cor(RCTD_nmf$nmf81[RCTD_nmf$platform == "VSPG"], RCTD_nmf$Astro[RCTD_nmf$platform == "VSPG"], use = "complete.obs")
[1] 0.9038443

##Ependy
cor(RCTD_nmf$nmf87, RCTD_nmf$Ependy)
[1] 0.7023605
cor(RCTD_nmf$nmf87[RCTD_nmf$platform == "HE"], RCTD_nmf$Ependy[RCTD_nmf$platform == "HE"], use = "complete.obs")
[1] 0.7582062
cor(RCTD_nmf$nmf87[RCTD_nmf$platform == "VSPG"], RCTD_nmf$Ependy[RCTD_nmf$platform == "VSPG"], use = "complete.obs")
[1] 0.8513248

##OPC
cor(RCTD_nmf$nmf36, RCTD_nmf$OPC)
[1] 0.3442669
cor(RCTD_nmf$nmf36[RCTD_nmf$platform == "HE"], RCTD_nmf$OPC[RCTD_nmf$platform == "HE"], use = "complete.obs")
[1] 0.3890879
cor(RCTD_nmf$nmf36[RCTD_nmf$platform == "VSPG"], RCTD_nmf$OPC[RCTD_nmf$platform == "VSPG"], use = "complete.obs")
[1] 0.7379624

##Sub.1
cor(RCTD_nmf$nmf40, RCTD_nmf$Sub.1)
[1] 0.6925582
cor(RCTD_nmf$nmf40[RCTD_nmf$platform == "HE"], RCTD_nmf$Sub.1[RCTD_nmf$platform == "HE"], use = "complete.obs")
[1] 0.7873364
cor(RCTD_nmf$nmf40[RCTD_nmf$platform == "VSPG"], RCTD_nmf$Sub.1[RCTD_nmf$platform == "VSPG"], use = "complete.obs")
[1] 0.8222314

##Sub.2
cor(RCTD_nmf$nmf54, RCTD_nmf$Sub.2)
[1] 0.6688239
cor(RCTD_nmf$nmf54[RCTD_nmf$platform == "HE"], RCTD_nmf$Sub.2[RCTD_nmf$platform == "HE"], use = "complete.obs")
[1] 0.823296
cor(RCTD_nmf$nmf54[RCTD_nmf$platform == "VSPG"], RCTD_nmf$Sub.2[RCTD_nmf$platform == "VSPG"], use = "complete.obs")
[1] 0.924581

##CA1
cor(RCTD_nmf$nmf15, RCTD_nmf$CA1_ProS)
[1] 0.8234017
cor(RCTD_nmf$nmf15[RCTD_nmf$platform == "HE"], RCTD_nmf$CA1_ProS[RCTD_nmf$platform == "HE"], use = "complete.obs")
[1] 0.8928128
cor(RCTD_nmf$nmf15[RCTD_nmf$platform == "VSPG"], RCTD_nmf$CA1_ProS[RCTD_nmf$platform == "VSPG"], use = "complete.obs")
[1] 0.8925183

##CA2.4
cor(RCTD_nmf$nmf61, RCTD_nmf$CA2.4)
[1] 0.4881745
cor(RCTD_nmf$nmf61[RCTD_nmf$platform == "HE"], RCTD_nmf$CA2.4[RCTD_nmf$platform == "HE"], use = "complete.obs")
[1] 0.5589571
cor(RCTD_nmf$nmf61[RCTD_nmf$platform == "VSPG"], RCTD_nmf$CA2.4[RCTD_nmf$platform == "VSPG"], use = "complete.obs")
[1] 0.6903548

##choroid
cor(RCTD_nmf$nmf48, RCTD_nmf$Choroid)
[1] 0.691953
cor(RCTD_nmf$nmf48[RCTD_nmf$platform == "HE"], RCTD_nmf$Choroid[RCTD_nmf$platform == "HE"], use = "complete.obs")
[1] 0.6941381
cor(RCTD_nmf$nmf48[RCTD_nmf$platform == "VSPG"], RCTD_nmf$Choroid[RCTD_nmf$platform == "VSPG"], use = "complete.obs")
[1] 0.5243216

##GABA.LAMP5
cor(RCTD_nmf$nmf35, RCTD_nmf$GABA.LAMP5)
[1] 0.2788109
cor(RCTD_nmf$nmf35[RCTD_nmf$platform == "HE"], RCTD_nmf$GABA.LAMP5[RCTD_nmf$platform == "HE"], use = "complete.obs")
[1] 0.3063416
cor(RCTD_nmf$nmf35[RCTD_nmf$platform == "VSPG"], RCTD_nmf$GABA.LAMP5[RCTD_nmf$platform == "VSPG"], use = "complete.obs")
[1] 0.6612972

cor(RCTD_nmf$nmf47,RCTD_nmf$GABA.LAMP5)
[1] 0.5716698
cor(RCTD_nmf$nmf47[RCTD_nmf$platform == "HE"], RCTD_nmf$GABA.LAMP5[RCTD_nmf$platform == "HE"], use = "complete.obs")
[1] 0.628553
cor(RCTD_nmf$nmf47[RCTD_nmf$platform == "VSPG"], RCTD_nmf$GABA.LAMP5[RCTD_nmf$platform == "VSPG"], use = "complete.obs")
[1] 0.6773249

##GC
cor(RCTD_nmf$nmf5,RCTD_nmf$GC)
[1] 0.7345699
cor(RCTD_nmf$nmf5[RCTD_nmf$platform == "HE"], RCTD_nmf$GC[RCTD_nmf$platform == "HE"], use = "complete.obs")
[1] 0.8640261
cor(RCTD_nmf$nmf5[RCTD_nmf$platform == "VSPG"], RCTD_nmf$GC[RCTD_nmf$platform == "VSPG"], use = "complete.obs")
[1] 0.9607683

##L2_3.Prs.Ent
cor(RCTD_nmf$nmf84,RCTD_nmf$L2_3.Prs.Ent)
[1] 0.225854
cor(RCTD_nmf$nmf84[RCTD_nmf$platform == "HE"], RCTD_nmf$L2_3.Prs.Ent[RCTD_nmf$platform == "HE"], use = "complete.obs")
[1] 0.2368461
cor(RCTD_nmf$nmf84[RCTD_nmf$platform == "VSPG"], RCTD_nmf$L2_3.Prs.Ent[RCTD_nmf$platform == "VSPG"], use = "complete.obs")
[1] 0.4079313

##L6_6b
cor(RCTD_nmf$nmf53,RCTD_nmf$L6_6b)
[1] 0.5072164
cor(RCTD_nmf$nmf53[RCTD_nmf$platform == "HE"], RCTD_nmf$L6_6b[RCTD_nmf$platform == "HE"], use = "complete.obs")
[1] 0.5788333
cor(RCTD_nmf$nmf53[RCTD_nmf$platform == "VSPG"], RCTD_nmf$L6_6b[RCTD_nmf$platform == "VSPG"], use = "complete.obs")
[1] 0.6210496

##L5
cor(RCTD_nmf$nmf51,RCTD_nmf$L5)
[1] 0.4016038
cor(RCTD_nmf$nmf51[RCTD_nmf$platform == "HE"], RCTD_nmf$L5[RCTD_nmf$platform == "HE"], use = "complete.obs")
[1] 0.4455087
cor(RCTD_nmf$nmf51[RCTD_nmf$platform == "VSPG"], RCTD_nmf$L5[RCTD_nmf$platform == "VSPG"], use = "complete.obs")
[1] 0.683681

##Micro_Macro_T
cor(RCTD_nmf$nmf90,RCTD_nmf$Micro_Macro_T)
[1] 0.449008
cor(RCTD_nmf$nmf90[RCTD_nmf$platform == "HE"], RCTD_nmf$Micro_Macro_T[RCTD_nmf$platform == "HE"], use = "complete.obs")
[1] 0.5021027
cor(RCTD_nmf$nmf90[RCTD_nmf$platform == "VSPG"], RCTD_nmf$Micro_Macro_T[RCTD_nmf$platform == "VSPG"], use = "complete.obs")
[1] 0.701673

#colData(spg) <- as(colData(RCTD_nmf), "DataFrame")

man.pal = c("THAL"="#1e1eff","CTX"="#5ffffb", "SUB"="#add294", "PCL-CA1"="#00dc00", "PCL-CA3"="#00a000", "CA4"="#99ff99",
            "GCL"="#005000", "SGZ"="#dfa56e", "ML"="#c1c1c1", "SL"="#444444", "SR"="#828E84", "SLM"="tan4",
            "SO"="#A698AE", "WM"="#7a007a", "CP"="#00006a")
			load(here("processed-data", "06_clustering","PRECAST","spe_norm_with_domain.rda"))

spg = spe_norm[,which(spe_norm$slide == "V12D07-332" | spe_norm$slide == "V12D07-335")]
spg = spe[,which(spe$sample_id == "V11L05-333_A1" | spe$sample_id == "V11L05-333_B1")]

spg$key = sapply(strsplit(spg$key,"_Br"), `[`, 1)
nmf_df = read_csv(here("processed-data", "VSPG_image_stitching", "nmf.csv"),col_names = TRUE,show_col_types = FALSE)
nmf_df <- nmf_df[ , -1]
    
Dr = here("processed-data", "spot_deconvo", "RCTD", "2ndRun_newClass_RCTDmarkers", "layer")
csv_files = data.frame(files = list.files(Dr, pattern = "clusters.csv" , recursive = TRUE, full.names = TRUE))
csv_files = data.frame(files = csv_files[1:40,])
counts_list <- lapply(csv_files$files, function(file_path) 
    {data <- read.csv(file_path, row.names = NULL)  
    return(data)})
temp_df <- do.call(rbind, counts_list)
rownames(temp_df) <- NULL
temp_df$key = temp_df$X
temp_df$key = sapply(strsplit(temp_df$key,"_Br"), `[`, 1)

colData_df <- colData(spg) |>
  data.frame() |>
  left_join(
    nmf_df,
    by = c("key"),
    relationship = "one-to-one"
  )|>
  left_join(
    temp_df,
    by = c("key"),
    relationship = "one-to-one"
  )
  
rownames(colData_df) <- colnames(spg)
colData(spg) <- DataFrame(colData_df)
  
library(patchwork)
sample_ids <- unique(spg$sample_id)
pdf(here("plots", "revision_maddy", "nmf77VsRCTDoligo.pdf"), width = 20, height = 8)  # Adjust dimensions as needed
for (id in sample_ids) {
speb <- spg[, spg$sample_id == id]
speb$domain1 <- "WM"
speb$domain1[!speb$domain %in% c("WM.1", "WM.2", "WM.3")] <- "others"

p = make_escheR(speb) |> add_fill("nmf77") |> add_ground("domain1") + scale_color_manual(values = c("WM" = "#7a007a", "others" = "#CCCCCC40")) +
      scale_fill_distiller(type = "seq",palette = rev('Greys'),direction=1, na.value = "#CCCCCC40")  +  labs(title = id)

p1 = make_escheR(speb) |> add_fill("Oligo") |> add_ground("domain1") + scale_color_manual(values = c("WM" = "#7a007a", "others" = "#CCCCCC40")) +
      scale_fill_distiller(type = "seq",palette = rev('Greys'),direction=1, na.value = "#CCCCCC40")  +  labs(title = id)

combined_plot = p + p1+ plot_layout(ncol = 2)
print(combined_plot)
	  }

dev.off()
  
   
pdf(here("plots", "revision_maddy", "nmf7VsRCTDgc.pdf"), width = 20, height = 8)  # Adjust dimensions as needed
for (id in sample_ids) {
	  
speb <- spg[, spg$sample_id == id]
speb$domain1 <- "GCL"
speb$domain1[!speb$domain %in% c("GCL")] <- "others"

p = make_escheR(speb) |> add_fill("nmf7") |> add_ground("domain1") + scale_color_manual(values = c("GCL" = "#005000", "others"= "#CCCCCC40")) +
      scale_fill_distiller(type = "seq",palette = rev('Greys'),direction=1, na.value = "#CCCCCC40")  +  labs(title = id)

p1 = make_escheR(speb) |> add_fill("GC") |> add_ground("domain1") + scale_color_manual(values = c("GCL" = "#005000", "others"= "#CCCCCC40")) +
      scale_fill_distiller(type = "seq",palette = rev('Greys'),direction=1, na.value = "#CCCCCC40")  +  labs(title = id)

combined_plot = p + p1+ plot_layout(ncol = 2)
print(combined_plot)

}
dev.off()


pdf(here("plots", "revision_maddy", "nmf81VsRCTDastro.pdf"), width = 20, height = 8)  # Adjust dimensions as needed
for (id in sample_ids) {
	  
speb <- spg[, spg$sample_id == id]
speb$domain1 <- "ML/SLM/SL/SR"
speb$domain1[!speb$domain %in% c("ML", "SL.SR", "SR.SLM", "SLM.WM")] <- "others"

p = make_escheR(speb) |> add_fill("nmf81") |> add_ground("domain1") + scale_color_manual(values = c("ML/SLM/SL/SR" = "#dfa56e", "others" = "#CCCCCC40")) +
      scale_fill_distiller(type = "seq",palette = rev('Greys'),direction=1, na.value = "#CCCCCC40")  +  labs(title = id)

p1 = make_escheR(speb) |> add_fill("Astro") |> add_ground("domain1") + scale_color_manual(values = c("ML/SLM/SL/SR" = "#dfa56e", "others" = "#CCCCCC40")) +
      scale_fill_distiller(type = "seq",palette = rev('Greys'),direction=1, na.value = "#CCCCCC40")  +  labs(title = id)

combined_plot = p + p1+ plot_layout(ncol = 2)
print(combined_plot)

}
dev.off()


pdf(here("plots", "revision_maddy", "nmf54VsRCTDsub.2.pdf"), width = 20, height = 8)  # Adjust dimensions as needed
for (id in sample_ids) {
	  
speb <- spg[, spg$sample_id == id]
speb$domain1 <- "SUB"
speb$domain1[!speb$domain %in% c("SUB")] <- "others"

p = make_escheR(speb) |> add_fill("nmf54") |> add_ground("domain1") + scale_color_manual(values = c("SUB" = "#add294", "others" = "#CCCCCC40")) +
      scale_fill_distiller(type = "seq",palette = rev('Greys'),direction=1, na.value = "#CCCCCC40")  +  labs(title = id)

p1 = make_escheR(speb) |> add_fill("Sub.2") |> add_ground("domain1") + scale_color_manual(values = c("SUB" = "#add294", "others" = "#CCCCCC40")) +
      scale_fill_distiller(type = "seq",palette = rev('Greys'),direction=1, na.value = "#CCCCCC40")  +  labs(title = id)

combined_plot = p + p1+ plot_layout(ncol = 2)
print(combined_plot)

}
dev.off()


######## updated 
Dr <- here("processed-data","spot_deconvo","shared_utilities")
spe = readRDS(here(Dr,"spe.rds"))

load(here('processed-data/NMF/spe_nmf_final.rda'))

Dr = here("processed-data", "spot_deconvo", "RCTD", "2ndRun_newClass_RCTDmarkers", "layer")
csv_files = data.frame(files = list.files(Dr, pattern = "clusters.csv" , recursive = TRUE, full.names = TRUE))
csv_files = data.frame(files = csv_files[1:40,])
counts_list <- lapply(csv_files$files, function(file_path) 
    {data <- read.csv(file_path, row.names = NULL)  
    return(data)})
temp_df <- do.call(rbind, counts_list)
rownames(temp_df) <- NULL
temp_df$key = temp_df$X
temp_df$key = sapply(strsplit(temp_df$key,"_Br"), `[`, 1)


colData_df <- colData(spe) |>
  data.frame() |>
  left_join(
    temp_df,
    by = c("key"),
    relationship = "one-to-one"
  )
  
rownames(colData_df) <- colnames(spe)
colData(spe) <- DataFrame(colData_df)

load('/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/plots/spatial_palettes.rda')
srt.palette <- c(srt.palette, others = "#CCCCCC40")
names(srt.palette)[names(srt.palette) == "SLM.WM"] <- "SLM.SGZ"

sample_ids = c("V11L05-333_A1", "V11L05-333_B1")
pdf(here("plots", "revision_maddy", "nmf13VsRCTDL23.pdf"), width = 7, height = 5)  # Adjust dimensions as needed
for (id in sample_ids) {
	  
speb <- spe[, spe$sample_id == id]
#speb$domain1 <- speb$cluster_collapsed
#speb$domain1[!speb$cluster_collapsed %in%  c("ML", "SL.SR" ,"SR.SLM", "SLM.WM")] <- "others"

p = make_escheR(speb) |> add_fill("nmf13") |> add_ground("domain", size = 6) + scale_color_manual(values = srt.palette) +
      scale_fill_distiller(type = "seq",palette = rev('Greys'),direction=1, na.value = "#CCCCCC40", limits = c(0,max(speb$nmf13))) 

p1 = make_escheR(speb) |> add_fill("L2_3.PrS.PaS") |> add_ground("domain", size = 6) + scale_color_manual(values = srt.palette) +
      scale_fill_distiller(type = "seq",palette = rev('Greys'),direction=1, na.value = "#CCCCCC40", limits = c(0,1)) 

print(p)
print(p1)

}
dev.off()


nmf7 , GABA.MGE , c("GABA", "RHP", "SUB.RHP", "SUB"), "#5ffffb"
nmf13, L2_3.PrS.PaS, c("GCL", "CA1", "CA2.4", "RHP"), "#008000"
nmf44, Oligo, c("WM.1", "WM.2", "WM.3"), "#ff3ffc"
nmf81, Astro, c("ML", "SL.SR" ,"SR.SLM", "SLM.WM"), "#dfa56e"


colData_df <- as.data.frame(colData(spe))
nmfs = c("nmf81", "nmf44", "nmf13", "nmf7")
pdf(here("plots", "revision_maddy", "nmfs.pdf"), width = 10, height = 5)  # Adjust dimensions as needed
for (nmf in nmfs) {
	  
p = ggplot(data = colData_df, aes(x = domain,  y = .data[[nmf]], fill = domain))+geom_boxplot()+ scale_fill_manual(values = srt.palette) +
	theme(legend.position = "none", axis.title = element_blank(), axis.text.x = element_text(angle = 90, size = 14), axis.text.y = element_text(size = 14))
print(p)

}
dev.off()

