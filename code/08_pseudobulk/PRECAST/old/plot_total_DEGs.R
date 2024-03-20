###############################
# spatial_HPC project
# Plot DE analysis results
# Anthony Ramnauth, Jan 05 2022
###############################

setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")

suppressPackageStartupMessages({
    library(SpatialExperiment)
    library(here)
    library(SingleCellExperiment)
    library(dplyr)
    library(ggplot2)
    library(sessioninfo)
})

# Load modeling results
modeling_results <- readRDS(file = here::here("processed-data", "08_pseudobulk", "PRECAST",
    "modeling_results.rds"))

# Filter for FDR < 0.05 for each PRECAST cluster
enriched1_filtered_summary <- modeling_results$enrichment %>%
    filter(fdr_1 < 0.05)
enriched2_filtered_summary <- modeling_results$enrichment %>%
    filter(fdr_2 < 0.05)
enriched3_filtered_summary <- modeling_results$enrichment %>%
    filter(fdr_3 < 0.05)
enriched4_filtered_summary <- modeling_results$enrichment %>%
    filter(fdr_4 < 0.05)
enriched5_filtered_summary <- modeling_results$enrichment %>%
    filter(fdr_5 < 0.05)
enriched6_filtered_summary <- modeling_results$enrichment %>%
    filter(fdr_6 < 0.05)
enriched7_filtered_summary <- modeling_results$enrichment %>%
    filter(fdr_7 < 0.05)
enriched8_filtered_summary <- modeling_results$enrichment %>%
    filter(fdr_8 < 0.05)
enriched10_filtered_summary <- modeling_results$enrichment %>%
    filter(fdr_10 < 0.05)
enriched11_filtered_summary <- modeling_results$enrichment %>%
    filter(fdr_11 < 0.05)
enriched12_filtered_summary <- modeling_results$enrichment %>%
    filter(fdr_12 < 0.05)
enriched13_filtered_summary <- modeling_results$enrichment %>%
    filter(fdr_13 < 0.05)
enriched14_filtered_summary <- modeling_results$enrichment %>%
    filter(fdr_14 < 0.05)

df <- data.frame(
    cluster = c(1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14),
    DEGs = c("enriched","enriched","enriched", "enriched","enriched","enriched",
        "enriched","enriched","enriched", "enriched","enriched","enriched","enriched",
         "depleted",  "depleted",  "depleted",  "depleted",  "depleted",  "depleted",
         "depleted",  "depleted",  "depleted",  "depleted",  "depleted",  "depleted",  "depleted"),
    total = c(
        sum(enriched1_filtered_summary$t_stat_1 > 0),
        sum(enriched2_filtered_summary$t_stat_2 > 0),
        sum(enriched3_filtered_summary$t_stat_3 > 0),
        sum(enriched4_filtered_summary$t_stat_4 > 0),
        sum(enriched5_filtered_summary$t_stat_5 > 0),
        sum(enriched6_filtered_summary$t_stat_6 > 0),
        sum(enriched7_filtered_summary$t_stat_7 > 0),
        sum(enriched8_filtered_summary$t_stat_8 > 0),
        sum(enriched10_filtered_summary$t_stat_10 > 0),
        sum(enriched11_filtered_summary$t_stat_11 > 0),
        sum(enriched12_filtered_summary$t_stat_12 > 0),
        sum(enriched13_filtered_summary$t_stat_13 > 0),
        sum(enriched14_filtered_summary$t_stat_14 > 0),
        -sum(enriched1_filtered_summary$t_stat_1 < 0),
        -sum(enriched2_filtered_summary$t_stat_2 < 0),
        -sum(enriched3_filtered_summary$t_stat_3 < 0),
        -sum(enriched4_filtered_summary$t_stat_4 < 0),
        -sum(enriched5_filtered_summary$t_stat_5 < 0),
        -sum(enriched6_filtered_summary$t_stat_6 < 0),
        -sum(enriched7_filtered_summary$t_stat_7 < 0),
        -sum(enriched8_filtered_summary$t_stat_8 < 0),
        -sum(enriched10_filtered_summary$t_stat_10 < 0),
        -sum(enriched11_filtered_summary$t_stat_11 < 0),
        -sum(enriched12_filtered_summary$t_stat_12 < 0),
        -sum(enriched13_filtered_summary$t_stat_13 < 0),
        -sum(enriched14_filtered_summary$t_stat_14 < 0)
    )
)

pdf(file = here::here("plots", "08_pseudobulk", "PRECAST", "PRECAST_total_DEGs_barplot.pdf"),
    width = 6, height = 4)
ggplot(df, aes(x = cluster, y = total, fill = DEGs))+
  geom_col()  +
  scale_fill_manual(values=c('indianred2', 'cornflowerblue'))+
  scale_y_continuous(breaks = seq(-10000, 10000, by = 500)) +
  scale_x_continuous(breaks = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 13),
        legend.text = element_text(size = 10),
        legend.position='right') +
  xlab('PRECAST clusters') + ylab('Total DEGs') +
  geom_text(aes(label = abs(total)), size=2, vjust=c(rep(-.5,13),rep(1.2,13))) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
dev.off()

