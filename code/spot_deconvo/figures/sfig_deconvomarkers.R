setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")
suppressPackageStartupMessages(library("SingleCellExperiment"))
suppressPackageStartupMessages(library("DeconvoBuddies"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("sessioninfo"))
suppressPackageStartupMessages(library("here"))
suppressPackageStartupMessages(library("spatialLIBD"))
suppressPackageStartupMessages(library("cowplot"))

source(here("code","spot_deconvo","shared_utilities","plottingfunctions.R"))
source(here("code", "spot_deconvo", "shared_utilities", "shared_function.R"))

Dr <- here("processed-data","spot_deconvo","shared_utilities")
plot_dir <- here("plots", "spot_deconvo", "figures", "fig_markergene_identification")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

marker_stats = readRDS(here(Dr,paste0("marker_stats_layer_celltype_class1_noHATAGABAAmy.rds")))
n_markers_per_type  = 25

png(here(plot_dir, "sfig_markergene_identification.png"), width = 1200, height = 1200, units = "px") 

p <- marker_stats |>
  mutate(
    Marker = case_when(
      rank_ratio <= n_markers_per_type ~ paste0(
        "Top-", n_markers_per_type, " Markers"
      ),
      TRUE ~ "Not Marker"
    )
  ) |>
  ggplot(aes(ratio, std.logFC, color = Marker)) +
  geom_point(size = 0.5, alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  facet_wrap(~cellType.target, scales = "free_x", nrow = 5) +
  labs(x = "Mean Ratio") +
  theme_bw() +
  theme(text = element_text(size = 20, colour = "black"),
        axis.text = element_text(size = 24, colour = "black"),
        axis.text.x = element_text(angle = 90),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        strip.text = element_text(size = 24, color = "black"),
        legend.position = c(1, 0),
        legend.justification = c(1, 0))+
        guides(col = guide_legend(override.aes = list(size = 5)))
print(p)
dev.off()
