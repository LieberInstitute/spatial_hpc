setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")
suppressPackageStartupMessages(library("SingleCellExperiment"))
suppressPackageStartupMessages(library("here"))
suppressPackageStartupMessages(library("spatialLIBD"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("escheR"))
library(readr)
library(tidyverse)


load(here('processed-data/NMF/spe_nmf_final.rda'))

Dr = here("processed-data", "spot_deconvo", "RCTD", "2ndRun_newClass_RCTDmarkers", "layer")
csv_files = data.frame(files = list.files(Dr, pattern = "clusters.csv" , recursive = TRUE, full.names = TRUE))
csv_files = data.frame(files = csv_files[-c(33:40),])
counts_list <- lapply(csv_files$files, function(file_path) 
    {data <- read.csv(file_path, row.names = NULL)  
    return(data)})
temp_df <- do.call(rbind, counts_list)
rownames(temp_df) <- NULL
temp_df$X <- gsub("Br2720", "V12F14-051", temp_df$X)
temp_df$key = temp_df$X

colData_df <- colData(spe) |>
  data.frame() |>
  left_join(
    temp_df,
    by = c("key"),
    relationship = "one-to-one"
  )
 
  
rownames(colData_df) <- colnames(spe)
colData(spe) <- DataFrame(colData_df)

cor(colData_df$nmf77, colData_df$Oligo, use = "complete.obs")
[1] 0.9412552
cor(colData_df$nmf44, colData_df$Oligo, use = "complete.obs") 
[1] 0.620673
cor(colData_df$nmf38, colData_df$Oligo, use = "complete.obs") 
[1] 0.5370027
cor(colData_df$nmf42, colData_df$Oligo, use = "complete.obs") 
[1] 0.3680308

cor(colData_df$nmf79, colData_df$Astro, use = "complete.obs") 
[1] 0.4211249
cor(colData_df$nmf81, colData_df$Astro, use = "complete.obs") 
[1] 0.8577056

cor(colData_df$nmf7, colData_df$GABA.MGE, use = "complete.obs") 
[[1] 0.367437
cor(colData_df$nmf13, colData_df$L2_3.PrS.PaS, use = "complete.obs") 
[1] 0.462544


load('/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/plots/spatial_palettes.rda')
srt.palette <- c(srt.palette, others = "#CCCCCC40")
names(srt.palette)[names(srt.palette) == "SLM.WM"] <- "SLM.SGZ"

sample_ids = c("V11L05-333_A1", "V11L05-333_B1")
pdf(here("plots", "revision_maddy", "RCTDs.pdf"), width = 7, height = 5)  # Adjust dimensions as needed
for (id in sample_ids) {
	  
speb <- spe[, spe$sample_id == id]
#speb$domain1 <- speb$cluster_collapsed
#speb$domain1[!speb$cluster_collapsed %in%  c("ML", "SL.SR" ,"SR.SLM", "SLM.WM")] <- "others"

p = make_escheR(speb) |> add_fill("Oligo") |> add_ground("domain", size = 6) + scale_color_manual(values = srt.palette) +
      scale_fill_distiller(type = "seq",palette = rev('Greys'),direction=1, na.value = "#CCCCCC40", limits = c(0,1)) 

p1 = make_escheR(speb) |> add_fill("Astro") |> add_ground("domain", size = 6) + scale_color_manual(values = srt.palette) +
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

