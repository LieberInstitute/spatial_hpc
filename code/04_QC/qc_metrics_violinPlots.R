library("SpatialExperiment")
library("scuttle")
library("scran")
library("scater")
library("jaffelab")
library("tidyverse")
library("here")
library("sessioninfo")

load(here("processed-data", "04_QC", "spe_QC.Rdata"), verbose = TRUE)

#### QC plots brain####
pdf(here("plots","04_QC", "QC_outliers_brain.pdf"), width = 21)
## Mito rate
plotColData(spe, x = "brnum", y = "subsets_Mito_percent", colour_by = "high_mito_br") +
  ggtitle("Mito Precent by brain") +
  facet_wrap(~ spe$brnum, scales = "free_x", nrow = 1)

## low sum
plotColData(spe, x = "brnum", y = "sum", colour_by = "low_sum_br") +
  scale_y_log10() +
  ggtitle("Total count by brain") +
  facet_wrap(~ spe$brnum, scales = "free_x", nrow = 1)

## low detected
plotColData(spe, x = "brnum", y = "detected", colour_by = "low_detected_br") +
  scale_y_log10() +
  ggtitle("Detected features by brain") +
  facet_wrap(~ spe$brnum, scales = "free_x", nrow = 1)

# Mito rate vs n detected features
samples <- unique(colData(spe)$sample_id)

for(i in 1:length(samples)){
p = plotColData(spe[,which(spe$sample_id == samples[i])],
            x = "detected", y = "subsets_Mito_percent",
            colour_by = "discard_auto_br", point_size = 2.5, point_alpha = 0.5)
print(p)
}

# Detected features vs total count
pdf(here("plots","04_QC", "QC_outliers_brain_t.pdf"), width = 21)

#for(i in 1:length(samples)){
  for(i in 1:2){
  message(samples[i])
  p  =  plotColData(spe[,spe$sample_id == samples[i]],
              x = "sum", y = "detected",
              colour_by = "discard_auto_br", point_size = 2.5, point_alpha = 0.5) + ggtitle(samples[i])
print(p)
}

dev.off()


#### QC plots sample_id ####
pdf(here("plots", "QC", "QC_outliers_capture_area.pdf"), width = 21)
## Mito rate
plotColData(spe, x = "sample_id", y = "subsets_Mito_percent", colour_by = "high_mito_id") +
  ggtitle("Mito Precent by capture area") +
  facet_wrap(~ spe$brnum, scales = "free_x", nrow = 1)

## low sum
plotColData(spe, x = "sample_id", y = "sum", colour_by = "low_sum_id") +
  scale_y_log10() +
  ggtitle("Total count by capture area") +
  facet_wrap(~ spe$brnum, scales = "free_x", nrow = 1)

## low detected
plotColData(spe, x = "sample_id", y = "detected", colour_by = "low_detected_id") +
  scale_y_log10() +
  ggtitle("Detected features by capture area") +
  facet_wrap(~ spe$brnum, scales = "free_x", nrow = 1)

# Mito rate vs n detected features

plotColData(spe,
            x = "detected", y = "subsets_Mito_percent",
            colour_by = "discard_auto_id", point_size = 2.5, point_alpha = 0.5
)

# Detected features vs total count

plotColData(spe,
            x = "sum", y = "detected",
            colour_by = "discard_auto_id", point_size = 2.5, point_alpha = 0.5
)

dev.off()
