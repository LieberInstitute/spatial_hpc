library("SingleCellExperiment")
library("scuttle")
library("scran")
library("scater")
library("scDblFinder")
library("jaffelab")
library("batchelor")
library("tidyverse")
library("here")
library("sessioninfo")
library("HDF5Array")


## Load raw data
load(here("snRNAseq_hpc","processed-data", "sce", "sce_raw.rda"), verbose = TRUE)

##get droplet score filepaths
droplet_paths <- list.files(here("snRNAseq_hpc","processed-data", "build_sce", "droplet_scores"),
    full.names = TRUE
)

names(droplet_paths) <- gsub("st", "s", gsub("droplet_scores_|.Rdata", "", basename(droplet_paths)))

e.out <- lapply(droplet_paths, function(x) get(load(x)))

#### Compile drop empty info ####
logs <- list.files(here("snRNAseq_hpc","code", "build_sce", "logs"), pattern = "02_barcode_ranks_empty_drops_qc.[0-9]", full.names = TRUE)
logs <- map(logs, readLines)

knee_lower <- map_dbl(logs, ~ parse_number(.x[grepl("knee_lower =", .x)]))
names(knee_lower) <- gsub("st", "s", map_chr(logs, ~str_sub(.x[grepl("Running Sample: ", .x)], " ", 3)))

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
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#   3258    4689    5372    5676    6070   12456
##which sample is problematic?
drop_summary$Sample[which.max(drop_summary$non_empty)]
#[1] 17c-scp

#make Louise-style barplot
drop_barplot <- drop_summary %>%
    mutate(empty = total_n - non_empty) %>%
    select(-total_n) %>%
    pivot_longer(!Sample, names_to = "drop_type", values_to = "n_drop") %>%
    ggplot(aes(x = Sample, y = n_drop, fill = drop_type)) +
    geom_col() +
    scale_y_continuous(trans = "log10") +
    my_theme +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(drop_barplot, filename = here("snRNAseq_hpc","plots", "drop_barplot.png"), width = 9)


#drop_v_cutoff <- drop_summary %>%
#    left_join(round_info) %>%
#    ggplot(aes(x = lower_cutoff, y = non_empty, color = round)) +
#    geom_point() +
#    my_theme

#ggsave(drop_v_cutoff, filename = here(plot_dir, "drop_v_cutoff.png"))

## Check empty droplet results
map(e.out, ~ addmargins(table(Signif = .x$FDR <= FDR_cutoff, Limited = .x$Limited, useNA = "ifany")))

# 
# $`10c-scp`
# Limited
# Signif    FALSE    TRUE    <NA>     Sum
# FALSE   10207       0       0   10207
# TRUE      197    5456       0    5653
# <NA>        0       0 1308401 1308401
# Sum     10404    5456 1308401 1324261
# 
# $`11c-scp`
# Limited
# Signif    FALSE    TRUE    <NA>     Sum
# FALSE   10660       0       0   10660
# TRUE      208    4844       0    5052
# <NA>        0       0 1450620 1450620
# Sum     10868    4844 1450620 1466332
# 
# $`12c-scp`
# Limited
# Signif    FALSE    TRUE    <NA>     Sum
# FALSE    7677       0       0    7677
# TRUE      195    4361       0    4556
# <NA>        0       0 1115944 1115944
# Sum      7872    4361 1115944 1128177
# 
# $`13c-scp`
# Limited
# Signif    FALSE    TRUE    <NA>     Sum
# FALSE    4026       0       0    4026
# TRUE      173    4386       0    4559
# <NA>        0       0 1074107 1074107
# Sum      4199    4386 1074107 1082692
# 
# $`14c-scp`
# Limited
# Signif   FALSE   TRUE   <NA>    Sum
# FALSE   7808      0      0   7808
# TRUE     187   3071      0   3258
# <NA>       0      0 799692 799692
# Sum     7995   3071 799692 810758
# 
# $`15c-scp`
# Limited
# Signif    FALSE    TRUE    <NA>     Sum
# FALSE   12541       0       0   12541
# TRUE       67    5211       0    5278
# <NA>        0       0 1166925 1166925
# Sum     12608    5211 1166925 1184744
# 
# $`16c-scp`
# Limited
# Signif    FALSE    TRUE    <NA>     Sum
# FALSE   18213       0       0   18213
# TRUE      150    7342       0    7492
# <NA>        0       0 1204068 1204068
# Sum     18363    7342 1204068 1229773
# 
# $`17c-scp`
# Limited
# Signif    FALSE    TRUE    <NA>     Sum
# FALSE   15848       0       0   15848
# TRUE      997   11459       0   12456
# <NA>        0       0 1243251 1243251
# Sum     16845   11459 1243251 1271555
# 
# $`18c-scp`
# Limited
# Signif    FALSE    TRUE    <NA>     Sum
# FALSE    8942       0       0    8942
# TRUE      143    5907       0    6050
# <NA>        0       0 1306558 1306558
# Sum      9085    5907 1306558 1321550
# 
# $`19c-scp`
# Limited
# Signif    FALSE    TRUE    <NA>     Sum
# FALSE    9769       0       0    9769
# TRUE      144    5373       0    5517
# <NA>        0       0 1270552 1270552
# Sum      9913    5373 1270552 1285838
# 
# $`1c-scp`
# Limited
# Signif    FALSE    TRUE    <NA>     Sum
# FALSE   10949       0       0   10949
# TRUE      204    4528       0    4732
# <NA>        0       0 1601243 1601243
# Sum     11153    4528 1601243 1616924
# 
# $`20c-scp`
# Limited
# Signif    FALSE    TRUE    <NA>     Sum
# FALSE    5581       0       0    5581
# TRUE      212    5918       0    6130
# <NA>        0       0 1165468 1165468
# Sum      5793    5918 1165468 1177179
# 
# $`21c-scp`
# Limited
# Signif    FALSE    TRUE    <NA>     Sum
# FALSE    7451       0       0    7451
# TRUE      164    6388       0    6552
# <NA>        0       0 1415922 1415922
# Sum      7615    6388 1415922 1429925
# 
# $`22c-scp`
# Limited
# Signif    FALSE    TRUE    <NA>     Sum
# FALSE   20066       0       0   20066
# TRUE      250    5646       0    5896
# <NA>        0       0 1063707 1063707
# Sum     20316    5646 1063707 1089669
# 
# $`23c-scp`
# Limited
# Signif    FALSE    TRUE    <NA>     Sum
# FALSE   26411       0       0   26411
# TRUE      163    5302       0    5465
# <NA>        0       0 1353330 1353330
# Sum     26574    5302 1353330 1385206
# 
# $`24c-scp`
# Limited
# Signif   FALSE   TRUE   <NA>    Sum
# FALSE   6834      0      0   6834
# TRUE     234   3454      0   3688
# <NA>       0      0 898101 898101
# Sum     7068   3454 898101 908623
# 
# $`25c-scp`
# Limited
# Signif   FALSE   TRUE   <NA>    Sum
# FALSE   4516      0      0   4516
# TRUE     786   4321      0   5107
# <NA>       0      0 841483 841483
# Sum     5302   4321 841483 851106
# 
# $`26c-scp`
# Limited
# Signif    FALSE    TRUE    <NA>     Sum
# FALSE    8967       0       0    8967
# TRUE      541    6331       0    6872
# <NA>        0       0 1079287 1079287
# Sum      9508    6331 1079287 1095126
# 
# $`27c-scp`
# Limited
# Signif    FALSE    TRUE    <NA>     Sum
# FALSE    8182       0       0    8182
# TRUE      420    4677       0    5097
# <NA>        0       0 1084654 1084654
# Sum      8602    4677 1084654 1097933
# 
# $`2c-scp`
# Limited
# Signif    FALSE    TRUE    <NA>     Sum
# FALSE   11847       0       0   11847
# TRUE      118    3997       0    4115
# <NA>        0       0 1680673 1680673
# Sum     11965    3997 1680673 1696635


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

table(sce$high_mito, sce$Sample)
#        10c-scp 11c-scp 12c-scp 13c-scp 14c-scp 15c-scp 16c-scp 17c-scp 18c-scp
#  FALSE    5313    4383    4370    3729    2973    4551    6729   10981    5364
#  TRUE      353     634     106     748     248     659     648    1404     730
#
#        19c-scp 1c-scp 20c-scp 21c-scp 22c-scp 23c-scp 24c-scp 25c-scp 26c-scp
#  FALSE    4760   4742    4768    5287    5847    5462    3548    4289    6326
#  TRUE      770     14    1302    1204      31      25      54     562     407
#
#        27c-scp 2c-scp 32c-scp 33c-scp 36c-scp 37c-scp 38c-scp 39c-scp
#  FALSE    4780   3235    4259    4090    4082    4302    3673    3978
#  TRUE      306    888     790     645     860    1010     808     753



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
#87095 54685






qc_t <- addmargins(table(sce$Sample, sce$discard_semiauto))
#           FALSE   TRUE    Sum
#  10c-scp   2600   3066   5666
#  11c-scp   3995   1022   5017
#  12c-scp   2072   2404   4476
#  13c-scp   3636    841   4477
#  14c-scp   1600   1621   3221
#  15c-scp   3919   1291   5210
#  16c-scp   5685   1692   7377
#  17c-scp   8473   3912  12385
#  18c-scp   5321    773   6094
#  19c-scp   4618    912   5530
#  1c-scp    1991   2765   4756
#  20c-scp   4324   1746   6070
#  21c-scp   5229   1262   6491
#  22c-scp   1204   4674   5878
#  23c-scp    982   4505   5487
#  24c-scp    921   2681   3602
#  25c-scp    250   4601   4851
#  26c-scp   2320   4413   6733
#  27c-scp   2042   3044   5086
#  2c-scp    3145    978   4123
#  32c-scp   3721   1328   5049
#  33c-scp   3769    966   4735
#  36c-scp   4001    941   4942
#  37c-scp   4227   1085   5312
#  38c-scp   3454   1027   4481
#  39c-scp   3596   1135   4731
#  Sum      87095  54685 141780


round(100 * sweep(qc_t, 1, qc_t[, 3], "/"), 1)
#          FALSE  TRUE   Sum
#  10c-scp  45.9  54.1 100.0
#  11c-scp  79.6  20.4 100.0
#  12c-scp  46.3  53.7 100.0
#  13c-scp  81.2  18.8 100.0
#  14c-scp  49.7  50.3 100.0
#  15c-scp  75.2  24.8 100.0
#  16c-scp  77.1  22.9 100.0
#  17c-scp  68.4  31.6 100.0
#  18c-scp  87.3  12.7 100.0
#  19c-scp  83.5  16.5 100.0
#  1c-scp   41.9  58.1 100.0
#  20c-scp  71.2  28.8 100.0
#  21c-scp  80.6  19.4 100.0
#  22c-scp  20.5  79.5 100.0
#  23c-scp  17.9  82.1 100.0
#  24c-scp  25.6  74.4 100.0
#  25c-scp   5.2  94.8 100.0
#  26c-scp  34.5  65.5 100.0
#  27c-scp  40.1  59.9 100.0
#  2c-scp   76.3  23.7 100.0
#  32c-scp  73.7  26.3 100.0
#  33c-scp  79.6  20.4 100.0
#  36c-scp  81.0  19.0 100.0
#  37c-scp  79.6  20.4 100.0
#  38c-scp  77.1  22.9 100.0
#  39c-scp  76.0  24.0 100.0
#  Sum      61.4  38.6 100.0


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
plotColData(sce, x = "Sample", y = "detected", colour_by = "low_genes") +
    scale_y_log10() +
    ggtitle("Detected genes") #+
    # geom_hline(yintercept = 500)+
 #   facet_wrap(~ sce$round, scales = "free_x", nrow = 1)

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
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.0000  0.1439  0.4391  0.8215  1.0823 31.4607

quantile(sce$doubletScore, probs=seq(0,1,by=0.01),3)
#        0%         1%         2%         3%         4%         5%         6%
# 0.0000000  0.0000000  0.0000000  0.0000000  0.0000000  0.0000000  0.0000000
#        7%         8%         9%        10%        11%        12%        13%
# 0.0102140  0.0122600  0.0149840  0.0200000  0.0220680  0.0262080  0.0303120
#       14%        15%        16%        17%        18%        19%        20%
# 0.0363000  0.0400000  0.0441360  0.0484000  0.0549760  0.0600000  0.0655200
#       21%        22%        23%        24%        25%        26%        27%
# 0.0713580  0.0772380  0.0820080  0.0882720  0.0943360  0.1000000  0.1055600
#       28%        29%        30%        31%        32%        33%        34%
# 0.1111440  0.1200000  0.1236960  0.1325220  0.1400000  0.1458880  0.1544760
#       35%        36%        37%        38%        39%        40%        41%
# 0.1622720  0.1695900  0.1786720  0.1839000  0.1936860  0.2005640  0.2097760
#       42%        43%        44%        45%        46%        47%        48%
# 0.2200000  0.2278000  0.2366000  0.2452000  0.2547280  0.2642480  0.2748800
#       49%        50%        51%        52%        53%        54%        55%
# 0.2846960  0.2955680  0.3062520  0.3166800  0.3281875  0.3407040  0.3536400
#       56%        57%        58%        59%        60%        61%        62%
# 0.3641220  0.3773440  0.3895840  0.4039920  0.4183220  0.4296280  0.4445760
#       63%        64%        65%        66%        67%        68%        69%
# 0.4598000  0.4744620  0.4904000  0.5070560  0.5241600  0.5402820  0.5605760
#       70%        71%        72%        73%        74%        75%        76%
# 0.5800000  0.6011500  0.6249760  0.6459680  0.6728040  0.6988200  0.7233400
#       77%        78%        79%        80%        81%        82%        83%
# 0.7541700  0.7849380  0.8200000  0.8540880  0.8931740  0.9294740  0.9727337
#       84%        85%        86%        87%        88%        89%        90%
# 1.0178880  1.0621440  1.1167520  1.1695100  1.2263680  1.2933120  1.3600000
#       91%        92%        93%        94%        95%        96%        97%
# 1.4315260  1.5093866  1.5957800  1.6955240  1.8253100  2.0208090  2.3390060
#       98%        99%       100%
# 2.9200000  4.1483280 20.6479520


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

dbl_df %>%
    group_by(Sample) %>%
    summarize(
        median = median(doubletScore),
        q95 = quantile(doubletScore, .95),
        drop = sum(doubletScore >= 2.75),
        drop_percent = 100 * drop / n()
    )

# A tibble: 20 × 5
#   Sample  median   q95  drop drop_percent
#   <chr>    <dbl> <dbl> <int>        <dbl>
# 1 10c-scp  0.317  1.91   109        2.05
# 2 11c-scp  0.384  1.30    56        1.40
# 3 12c-scp  0.301  1.82   103        2.31
# 4 13c-scp  0.401  1.29    29        0.788
# 5 14c-scp  0.300  1.81    36        1.20
# 6 15c-scp  0.369  1.33    84        2.13
# 7 16c-scp  0.270  1.33   130        1.89
# 8 17c-scp  0.36   2.5    445        4.03
# 9 18c-scp  0.266  1.49   113        2.14
#10 19c-scp  0.287  1.25    64        1.40
#11 1c-scp   0.369  1.84    45        1.87
#12 20c-scp  0.380  1.37    71        1.48
#13 21c-scp  0.314  1.27   119        2.26
#14 22c-scp  0.342  2.26    48        3.62
#15 23c-scp  1.21   2.25    21        2
#16 24c-scp  0.686  1.82    10        0.923
#17 25c-scp  1.42   2.02     0        0
#18 26c-scp  0.426  1.79    43        1.32
#19 27c-scp  0.693  1.68    28        1.14
#20 2c-scp   0.387  1.49    44        1.40
 


table(sce$discard_semiauto, sce$doubletScore >= 2.75)
#        FALSE  TRUE
#  FALSE 75653  1598
#  TRUE  35433   841



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
# [1] 36601 77604

save(sce,file=here("snRNAseq_hpc","processed-data","sce","sce_post_qc.rda"))


# sgejobs::job_single('03_droplet_qc', create_shell = TRUE, queue= 'bluejay', memory = '50G', command = "Rscript 03_droplet_qc.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
