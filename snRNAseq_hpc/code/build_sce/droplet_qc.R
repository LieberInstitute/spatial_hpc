#cd /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/

library("SingleCellExperiment")
library("scran")
library("scater")
library("scry")
library("here")
library("sessioninfo")
library("harmony")
library("BiocSingular")
library("purrr")
library("stringr")
library('tidyverse')
library('dplyr')
library('magrittr')


##Load in QC'd sce object
load(file=here::here("snRNAseq_hpc","processed-data","sce","sce_raw.rda"))

##get droplet score filepaths
droplet_paths <- list.files(here("snRNAseq_hpc","processed-data", "build_sce", "droplet_scores"),
                            full.names = TRUE
)

names(droplet_paths) <- gsub("st", "s", gsub("droplet_scores_round2_|.Rdata", "", basename(droplet_paths)))

e.out <- lapply(droplet_paths, function(x) get(load(x)))

#### Compile drop empty info ####
logs <- list.files(here("snRNAseq_hpc","code", "build_sce", "logs"), pattern = "02_barcode_ranks_empty_drops_qc.[0-9]", full.names = TRUE)
logs <- map(logs, readLines)

extract_knee_lower <- function(log) {
    # Find the index of the line that contains "knee_lower ="
    index <- grep("knee_lower =", log)

    # If the phrase is not found, return NA
    if (length(index) == 0) {
        return(NA)
    } else {
        # Extract the value at the end of the line
        line <- log[index]
        value <- as.numeric(gsub(".*knee_lower =", "", line))
        return(value)
    }
}

# Function to handle duplicates or take max
handle_duplicates <- function(x) {
    if (length(x) > 1) {
        if (x[1] == x[2]) {
            x <- x[1]
        } else {
            x <-x[2]
        }
    }
    return(x)
}

# Apply function to each element of knee_lower_values
knee_lower_values <- sapply(knee_lower_values, handle_duplicates)
print(knee_lower_values)

##make sure e.out is mapped to knee_lower_values correctly
find_eout_name <- function(log) {
    # Extract sample name from the line that starts with "Running Sample:"
    sample_line <- grep("Running Sample:", log, value = TRUE)
    sample_name <- gsub("Running Sample: ", "", sample_line)
    sample_name <- gsub(" \\(.+\\)$", "", sample_name)  # Remove trailing "(1/20)" or similar
    return(sample_name)
}

log_names <- sapply(logs, find_eout_name)
log_names <- sapply(log_names, function(x) x[1])
names(knee_lower_values) <- log_names




##n empty table

FDR_cutoff <- 0.001

drop_summary <- stack(map_int(e.out, nrow)) %>%
    rename(total_n = values) %>%
    left_join(stack(map_int(e.out, ~ sum(.x$FDR < FDR_cutoff, na.rm = TRUE))) %>%
                  rename(non_empty = values)) %>%
    select(Sample = ind, total_n, non_empty) %>%
    left_join(stack(knee_lower) %>% rename(Sample = ind, lower_cutoff = values))

write_csv(drop_summary, file = here("snRNAseq_hpc","processed-data", "build_sce", "drop_summary.csv"))

drop_summary %>%
    arrange(non_empty)

summary(drop_summary$non_empty)
# drop_summary$Sample[which.max(drop_summary$non_empty)]
# Sample total_n non_empty lower_cutoff
# 1  14c-scp  810758      3221           NA
# 2  24c-scp  908623      3602           NA
# 3   2c-scp 1696635      4123           NA
# 4  12c-scp 1128177      4476           NA
# 5  13c-scp 1082692      4477           NA
# 6  38c-scp 1869148      4481           NA
# 7  39c-scp 2112382      4731           NA
# 8  33c-scp 1514385      4735           NA
# 9   1c-scp 1616924      4756           NA
# 10 25c-scp  851106      4851           NA
# 11 36c-scp 1845312      4942           NA
# 12 11c-scp 1466332      5017           NA
# 13 32c-scp 1313783      5049           NA
# 14 27c-scp 1097933      5086           NA
# 15 15c-scp 1184744      5210           NA
# 16 37c-scp 1920223      5312           NA
# 17 23c-scp 1385206      5487           NA
# 18 19c-scp 1285838      5530           NA
# 19 10c-scp 1324261      5666           NA
# 20 22c-scp 1089669      5878           NA
# 21 20c-scp 1177179      6070           NA
# 22 18c-scp 1321550      6094           NA
# 23 21c-scp 1429925      6491           NA
# 24 26c-scp 1095126      6733           NA
# 25 16c-scp 1229773      7377           NA
# 26 17c-scp 1271555     12385           NA
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 3221    4732    5068    5453    5825   12385
# > sum(drop_summary$non_empty)
# [1] 141780
# > drop_summary$Sample[which.max(drop_summary$non_empty)]
# [1] 17c-scp

drop_barplot <- drop_summary %>%
    mutate(empty = total_n - non_empty) %>%
    select(-total_n) %>%
    pivot_longer(!Sample, names_to = "drop_type", values_to = "n_drop") %>%
    ggplot(aes(x = Sample, y = n_drop, fill = drop_type)) +
    geom_col() +
    scale_y_continuous(trans = "log10") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(drop_barplot, filename = here("snRNAseq_hpc","plots", "drop_barplot.png"), width = 9)

## Check empty droplet results
map(e.out, ~ addmargins(table(Signif = .x$FDR <= FDR_cutoff, Limited = .x$Limited, useNA = "ifany")))

#### Eliminate empty droplets ####
e.out.all <- do.call("rbind", e.out)[colnames(sce), ]
sce <- sce[, which(e.out.all$FDR <= 0.001)]

dim(sce)
# [1] 36601 141780

#save drops removed sce
save(sce,file=here("snRNAseq_hpc","processed-data","sce","sce_drops_removed.rda"))

#### Compute QC metrics ####
sce <- scuttle::addPerCellQC(
    sce,
    subsets = list(Mito = which(seqnames(sce) == "chrM")),
    # BPPARAM = BiocParallel::MulticoreParam(4)
)


#### Check for low quality nuc ####
## High mito
# sce$high.mito.sample ## standard name?
sce$high_mito <- isOutlier(sce$subsets_Mito_percent, nmads = 3, type = "higher", batch = sce$Sample)
#FALSE   TRUE
#125821  15959

## low library size
sce$low_lib <- isOutlier(sce$sum, log = TRUE, type = "lower", batch = sce$Sample)
table(sce$low_lib)
# FALSE   TRUE
#139660   2120

## low detected features
# sce$qc.detected
sce$low_genes <- isOutlier(sce$detected, log = TRUE, type = "lower", batch = sce$Sample)
table(sce$low_genes)
# FALSE   TRUE
#137012   4768



## All low sum are also low detected
table(sce$low_lib, sce$low_genes)
#         FALSE   TRUE
#  FALSE 137012   2648
#  TRUE       0   2120


## Annotate nuc to drop
sce$discard_auto <- sce$high_mito | sce$low_lib | sce$low_genes

table(sce$discard_auto)
# FALSE   TRUE
#123421  18359



qc_t <- addmargins(table(sce$Sample, sce$discard_auto))
#           FALSE   TRUE    Sum
#  10c-scp   5313    353   5666
#  11c-scp   3995   1022   5017
#  12c-scp   4370    106   4476
#  13c-scp   3636    841   4477
#  14c-scp   2973    248   3221
#  15c-scp   3919   1291   5210
#  16c-scp   6729    648   7377
#  17c-scp  10981   1404  12385
#  18c-scp   5321    773   6094
#  19c-scp   4618    912   5530
#  1c-scp    4742     14   4756
#  20c-scp   4768   1302   6070
#  21c-scp   5229   1262   6491
#  22c-scp   5847     31   5878
#  23c-scp   5213    274   5487
#  24c-scp   3548     54   3602
#  25c-scp   4289    562   4851
#  26c-scp   6326    407   6733
#  27c-scp   4780    306   5086
#  2c-scp    3145    978   4123
#  32c-scp   4259    790   5049
#  33c-scp   3769    966   4735
#  36c-scp   4082    860   4942
#  37c-scp   4300   1012   5312
#  38c-scp   3673    808   4481
#  39c-scp   3596   1135   4731
#  Sum     123421  18359 141780



round(100 * sweep(qc_t, 1, qc_t[, 3], "/"), 1)
#          FALSE  TRUE   Sum
#  10c-scp  93.8   6.2 100.0
#  11c-scp  79.6  20.4 100.0
#  12c-scp  97.6   2.4 100.0
#  13c-scp  81.2  18.8 100.0
#  14c-scp  92.3   7.7 100.0
#  15c-scp  75.2  24.8 100.0
#  16c-scp  91.2   8.8 100.0
#  17c-scp  88.7  11.3 100.0
#  18c-scp  87.3  12.7 100.0
#  19c-scp  83.5  16.5 100.0
#  1c-scp   99.7   0.3 100.0
#  20c-scp  78.6  21.4 100.0
#  21c-scp  80.6  19.4 100.0
#  22c-scp  99.5   0.5 100.0
#  23c-scp  95.0   5.0 100.0
#  24c-scp  98.5   1.5 100.0
#  25c-scp  88.4  11.6 100.0
#  26c-scp  94.0   6.0 100.0
#  27c-scp  94.0   6.0 100.0
#  2c-scp   76.3  23.7 100.0
#  32c-scp  84.4  15.6 100.0
#  33c-scp  79.6  20.4 100.0
#  36c-scp  82.6  17.4 100.0
#  37c-scp  80.9  19.1 100.0
#  38c-scp  82.0  18.0 100.0
#  39c-scp  76.0  24.0 100.0
#  Sum      87.1  12.9 100.0



#### QC plots ####
pdf(here("snRNAseq_hpc","plots", "QC_auto.pdf"), width = 21)
## Mito rate
plotColData(sce, x = "Sample", y = "subsets_Mito_percent", colour_by = "high_mito") +
    ggtitle("Mito Percent") #+
#   facet_wrap(~ sce$round, scales = "free_x", nrow = 1)

# ## low sum
plotColData(sce, x = "Sample", y = "sum", colour_by = "low_lib") +
    scale_y_log10() +
    ggtitle("Total UMIs") #+
#  facet_wrap(~ sce$round, scales = "free_x", nrow = 1)
# +
#   geom_hline(yintercept = 1000) ## hline doesn't work w/ facet_wrap?

# ## low detected
plotColData(sce, x = "Sample", y = "detected", colour_by = "low_genes") +
    scale_y_log10() +
    ggtitle("Detected genes") #+
# geom_hline(yintercept = 500)+
#   facet_wrap(~ sce$round, scales = "free_x", nrow = 1)

# Mito rate vs n detected features

plotColData(sce,
            x = "detected", y = "subsets_Mito_percent",
            colour_by = "discard_auto", point_size = 2.5, point_alpha = 0.5
)

# Detected features vs total count

plotColData(sce,
            x = "sum", y = "detected",
            colour_by = "discard_auto", point_size = 2.5, point_alpha = 0.5
)

dev.off()

###there's a lot of issues here, especially with mito rate...let's fix that
##Use samples with reasonable QC metric thresholds to compute those that failed. For mito rate, we will do this without differentiating between NeuN and PI sorted guys
good_samples<-colnames(attributes(sce$high_mito)$thresholds[,attributes(sce$high_mito)$thresholds[2,] < 5])


discard_mito <- isOutlier(sce$subsets_Mito_percent,
                          type="higher", batch=sce$Sample,
                          subset=sce$Sample %in% good_samples)

#for detected and sum, we'll do NeuN and PI sorted separately due to differences in distributions for different sort strategies
x<-c('11c-scp','13c-scp','15c-scp','18c-scp','19c-scp','21c-scp','23c-scp','2c-scp','33c-scp','39c-scp')
y<-c('11c-scp','15c-scp','19c-scp','23c-scp','33c-scp','39c-scp')

sce$low_genes <- ifelse(sce$Sample %in% x, sce$low_genes,
                        ifelse(sce$detected <=1000,T,F))

sce$low_lib <- ifelse(sce$Sample %in% y, sce$low_lib,
                      ifelse(sce$sum <=1000,T,F))

sce$high_mito<-discard_mito


## Annotate nuc to drop
sce$discard_semiauto <- sce$high_mito | sce$low_lib | sce$low_genes

table(sce$discard_semiauto)
#FALSE  TRUE
#86095 54875






qc_t <- addmargins(table(sce$Sample, sce$discard_semiauto))
qc_t

# FALSE   TRUE    Sum
# 10c-scp   2600   3066   5666
# 11c-scp   3995   1022   5017
# 12c-scp   2072   2404   4476
# 13c-scp   3586    891   4477
# 14c-scp   1600   1621   3221
# 15c-scp   3919   1291   5210
# 16c-scp   5685   1692   7377
# 17c-scp   8473   3912  12385
# 18c-scp   5271    823   6094
# 19c-scp   4618    912   5530
# 1c-scp    1991   2765   4756
# 20c-scp   4324   1746   6070
# 21c-scp   5140   1351   6491
# 22c-scp   1204   4674   5878
# 23c-scp    982   4505   5487
# 24c-scp    921   2681   3602
# 25c-scp    250   4601   4851
# 26c-scp   2320   4413   6733
# 27c-scp   2042   3044   5086
# 2c-scp    3144    979   4123
# 32c-scp   3721   1328   5049
# 33c-scp   3769    966   4735
# 36c-scp   4001    941   4942
# 37c-scp   4227   1085   5312
# 38c-scp   3454   1027   4481
# 39c-scp   3596   1135   4731
# Sum      86905  54875 141780



round(100 * sweep(qc_t, 1, qc_t[, 3], "/"), 1)
# FALSE  TRUE   Sum
# 10c-scp  45.9  54.1 100.0
# 11c-scp  79.6  20.4 100.0
# 12c-scp  46.3  53.7 100.0
# 13c-scp  80.1  19.9 100.0
# 14c-scp  49.7  50.3 100.0
# 15c-scp  75.2  24.8 100.0
# 16c-scp  77.1  22.9 100.0
# 17c-scp  68.4  31.6 100.0
# 18c-scp  86.5  13.5 100.0
# 19c-scp  83.5  16.5 100.0
# 1c-scp   41.9  58.1 100.0
# 20c-scp  71.2  28.8 100.0
# 21c-scp  79.2  20.8 100.0
# 22c-scp  20.5  79.5 100.0
# 23c-scp  17.9  82.1 100.0
# 24c-scp  25.6  74.4 100.0
# 25c-scp   5.2  94.8 100.0
# 26c-scp  34.5  65.5 100.0
# 27c-scp  40.1  59.9 100.0
# 2c-scp   76.3  23.7 100.0
# 32c-scp  73.7  26.3 100.0
# 33c-scp  79.6  20.4 100.0
# 36c-scp  81.0  19.0 100.0
# 37c-scp  79.6  20.4 100.0
# 38c-scp  77.1  22.9 100.0
# 39c-scp  76.0  24.0 100.0
# Sum      61.3  38.7 100.0


#
pdf(here("snRNAseq_hpc","plots", "QC_semiauto.pdf"), width = 21)
## Mito rate
plotColData(sce, x = "Sample", y = "subsets_Mito_percent", colour_by = "high_mito") +
    ggtitle("Mito Percent") #+
#   facet_wrap(~ sce$round, scales = "free_x", nrow = 1)

# ## low sum
plotColData(sce, x = "Sample", y = "sum", colour_by = "low_lib") +
    scale_y_log10() +
    ggtitle("Total UMIs") #+
#  facet_wrap(~ sce$round, scales = "free_x", nrow = 1)
# +
#   geom_hline(yintercept = 1000) ## hline doesn't work w/ facet_wrap?

# ## low detected
plotColData(sce, x = "Sample", y = "detected", colour_by = "detected") +
    scale_y_log10() +
    ggtitle("Detected genes") +
    facet_wrap(~ sce$sort, scales = "free_x", nrow = 1)

# Mito rate vs n detected features

plotColData(sce,
            x = "detected", y = "subsets_Mito_percent",
            colour_by = "discard_semiauto", point_size = 2.5, point_alpha = 0.5
)

# Detected features vs total count

plotColData(sce,
            x = "sum", y = "detected",
            colour_by = "discard_semiauto", point_size = 2.5, point_alpha = 0.5
)

dev.off()

#### Doublet detection ####
## To speed up, run on sample-level top-HDGs - just take top 2000
set.seed(328)

colData(sce)$doubletScore <- NA

for (i in splitit(sce$Sample)) {
    sce_temp <- sce[, i]
    ## To speed up, run on sample-level top-HVGs - just take top 1000
    normd <- logNormCounts(sce_temp)
    geneVar <- modelGeneVar(normd)
    topHVGs <- getTopHVGs(geneVar, n = 2000)

    dbl_dens <- computeDoubletDensity(normd, subset.row = topHVGs)
    colData(sce)$doubletScore[i] <- dbl_dens
}

summary(sce$doubletScore)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 0.00000  0.08066  0.28000  0.53053  0.66208 20.30150
quantile(sce$doubletScore, probs=seq(0,1,by=0.01),3)
# 0%        1%        2%        3%        4%        5%        6%        7%
#     0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.000000  0.009470
# 8%        9%       10%       11%       12%       13%       14%       15%
#     0.010420  0.012140  0.014754  0.019768  0.020344  0.024376  0.028536  0.031872
# 16%       17%       18%       19%       20%       21%       22%       23%
#     0.036564  0.040000  0.044262  0.049420  0.054870  0.060000  0.064420  0.070536
# 24%       25%       26%       27%       28%       29%       30%       31%
#     0.076096  0.080658  0.087318  0.094048  0.098952  0.104632  0.110600  0.118032
# 32%       33%       34%       35%       36%       37%       38%       39%
#     0.123656  0.132580  0.140000  0.145680  0.152828  0.160854  0.168766  0.177048
# 40%       41%       42%       43%       44%       45%       46%       47%
#     0.184338  0.194040  0.200680  0.210140  0.220000  0.229240  0.238354  0.246876
# 48%       49%       50%       51%       52%       53%       54%       55%
#     0.258632  0.268136  0.280000  0.287560  0.298620  0.310464  0.321708  0.331800
# 56%       57%       58%       59%       60%       61%       62%       63%
#     0.342860  0.354310  0.367404  0.380000  0.393888  0.406866  0.420280  0.435584
# 64%       65%       66%       67%       68%       69%       70%       71%
#     0.449934  0.464548  0.480334  0.498608  0.514206  0.532262  0.552106  0.572928
# 72%       73%       74%       75%       76%       77%       78%       79%
#     0.593040  0.615030  0.636302  0.662082  0.686766  0.712414  0.740000  0.767208
# 80%       81%       82%       83%       84%       85%       86%       87%
#     0.797742  0.829500  0.861940  0.899650  0.937800  0.977408  1.020528  1.073396
# 88%       89%       90%       91%       92%       93%       94%       95%
#     1.129092  1.188122  1.257852  1.333760  1.414416  1.508192  1.614760  1.740000
# 96%       97%       98%       99%      100%
# 1.916148  2.205774  2.774416  4.102993 20.301504


## Visualize doublet scores ##

dbl_df <- colData(sce) %>%
    as.data.frame() %>%
    select(Sample, doubletScore)

dbl_box_plot <- dbl_df %>%
    ggplot(aes(x = reorder(Sample, doubletScore, FUN = median), y = doubletScore)) +
    geom_boxplot() +
    labs(x = "Sample") +
    geom_hline(yintercept = 2.75, color = "red", linetype = "dashed") +
    coord_flip()

ggsave(dbl_box_plot, filename = here(plot_dir, "doublet_scores_boxplot.png"))

dbl_density_plot <- dbl_df %>%
    ggplot(aes(x = doubletScore)) +
    geom_density() +
    labs(x = "doublet score") +
    facet_grid(Sample ~ .) +
    theme_bw()

ggsave(dbl_density_plot, filename = here(plot_dir, "doublet_scores_desnity.png"), height = 17)

## Save
save(sce, file = here::here("snRNAseq_hpc","processed-data", "sce", "sce_no_empty_droplets.Rdata"))

## Save out sample info for easy access
sample_info <- pd %>%
    group_by(Sample, file_id, region, subject, round) %>%
    summarize(
        n = n(),
        n_high_mito = sum(high_mito),
        n_low_sum = sum(low_sum),
        n_low_detect = sum(low_detected),
        n_discard_auto = sum(discard_auto)
    )

write_csv(sample_info, file = here("processed-data", "03_build_sce", "sample_info.csv"))

n_boxplot <- sample_info %>%
    ggplot(aes(x = round, y = n, color = round)) +
    geom_boxplot(outlier.shape = NA) +
    geom_point() +
    ggrepel::geom_text_repel(aes(label = Sample), color = "black") +
    my_theme +
    theme(legend.position = "None")

ggsave(n_boxplot, filename = here(plot_dir, "n_nuclei_boxplot.png"))

#### Save clean data as HDF5 file  ####
load(here("processed-data", "sce", "sce_no_empty_droplets.Rdata"))

sce <- sce[, !sce$discard_semiauto]
dim(sce)
# [1] 36601 86905


save(sce,file=here("snRNAseq_hpc","processed-data","sce","sce_post_qc.rda"))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
