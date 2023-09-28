setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")
suppressPackageStartupMessages({
  library("here")
  library("SpatialExperiment")
  library("rtracklayer")
  library("lobstr")
  library("sessioninfo")
})
#-------------------------------------------------------------------------------
#   Add cell counts to HE spatial object
#-------------------------------------------------------------------------------

spe_in <- here("processed-data","02_build_spe","spe_nmf_final.rda")
# sce_in <- "/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/snRNAseq_hpc/processed-data/sce/sce_final.rda"
sce_in <- "/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/snRNAseq_hpc/processed-data/sce/sce_class_final.rda"
spg_in <- here("processed-data","06_clustering", "PRECAST" , "spe_precast_final.rda")
  
out <- here("processed-data", "spot_deconvo", "shared_utilities")
  
load(spe_in)
t = colData(spe)[colData(spe)$slide=="V12F14-051",]
temp = paste0(sapply(strsplit(t$key,"Br"),'[',1),t$sample_id)
colData(spe)$key[spe$slide=="V12F14-051"]=temp

spaceranger_dirs = read.csv(file.path(here::here("code","VistoSeg","code","samples.txt")), header = FALSE, sep = '\t', stringsAsFactors = FALSE, col.names = c('SPpath','sample_id','brain'))
spaceranger_dirs$SPpath = paste0(spaceranger_dirs$SPpath,"outs/spatial/tissue_spot_counts.csv")

segmentations_list <-
  lapply(spaceranger_dirs$sample_id, function(sampleid) {
    file <-spaceranger_dirs$SPpath[spaceranger_dirs$sample_id == sampleid]
    if (!file.exists(file)) {
      return(NULL)
    }
    x <- read.csv(file)
    x$key <- paste0(x$barcode, "_", sampleid)
    return(x)
  })

## Merge them (once the these files are done, this could be replaced by an rbind)
segmentations <-
  Reduce(function(...) {
    merge(..., all = TRUE)
  }, segmentations_list[lengths(segmentations_list) > 0])

## Add the information
segmentation_match <- match(spe$key, segmentations$key)
segmentation_info <-
  segmentations[segmentation_match, -which(
    colnames(segmentations) %in% c("barcode", "tissue", "row", "col", "imagerow", "imagecol", "key")
  )]
colData(spe) <- cbind(colData(spe), segmentation_info)
colData(spe)$sample_id = as.character(colData(spe)$sample_id)
spe$count <- spe$CNmask_dark_blue

#   zellkonverter doesn't know how to convert the 'spatialCoords' slot. We'd
#   ultimately like the spatialCoords in the .obsm['spatial'] slot of the
#   resulting AnnDatas, which corresponds to reducedDims(spe)$spatial in R
load(spg_in)
spg = spe[,which(spe$slide == "V12D07-332" | spe$slide == "V12D07-335")]
reducedDims(spg)$spatial <- spatialCoords(spg)
rownames(spg) <- rowData(spg)$gene_id

reducedDims(spe)$spatial <- spatialCoords(spe)
#reducedDims(spg)$spatial <- spatialCoords(spg)

load(sce_in, verbose = TRUE)
rownames(sce) <- rowData(sce)$gene_id
#colData(sce)$layer.type = gsub("/", "_", colData(sce)$layer.type)
colData(sce)$broad.class = gsub("/", "_", colData(sce)$broad.class)
colData(sce)$cell.class = gsub("/", "_", colData(sce)$cell.class)
colData(sce)$cell.class2 = gsub("/", "_", colData(sce)$cell.class2)
colData(sce)$cell.class3 = gsub("/", "_", colData(sce)$cell.class3)
colData(sce)$cell.class4 = gsub("/", "_", colData(sce)$cell.class4)
sce = sce[,sce$cell.class!="HATA_Amy"]
sce = sce[,sce$cell.class!="GABA.Amy"]

rownames(spe) <- rowData(spe)$gene_id

## EDA on counts ## 
getmode <- function(v) {
  uniqv <- unique(na.omit(v))
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

# mean(spe$CNmask_dark_blue,na.rm=TRUE)
# [1] 5.279843
# median(spe$CNmask_dark_blue,na.rm=TRUE)
# [1] 3
# getmode(spe$CNmask_dark_blue)
# [1] 2
# 
# t(table(spe$CNmask_dark_blue))
# 
# 0     1     2     3     4     5     6     7     8     9    10    11
# [1,] 13171 19026 19881 17491 13528  9995  7013  4705  3271  2253  1694  1170
# 
# 12    13    14    15    16    17    18    19    20    21    22    23
# [1,]   966   799   645   557   496   435   376   329   283   265   228   189
# 
# 24    25    26    27    28    29    30    31    32    33    34    35
# [1,]   196   166   149   147   133   115   108    95    89    81    85    79
# 
# 36    37    38    39    40    41    42    43    44    45    46    47
# [1,]    62    63    67    58    40    60    47    45    57    46    40    34
# 
# 48    49    50    51    52    53    54    55    56    57    58    59
# [1,]    40    53    37    31    34    31    29    27    25    24    20    25
# 
# 60    61    62    63    64    65    66    67    68    69    70    71
# [1,]    17    15    22    23    20    16    21    18    16    22    24    22
# 
# 72    73    74    75    76    77    78    79    80    81    82    83
# [1,]     8    10    20    16    16    18    19    11     8    17     8     9
# 
# 84    85    86    87    88    89    90    91    92    93    94    95
# [1,]    16    17     8    21    10    12    12    16    23    14     8    17
# 
# 96    97    98    99   100   101   102   103   104   105   106   107
# [1,]    16    14    20    16    20    18    11    35    29    31    18    19
# 
# 108   109   110   111   112   113   114   115   116   117   118   119
# [1,]    19    25    26    23    29    25    22    13    24    11    13    16
# 
# 120   121   122   123   124   125   126   127   128   129   131   132
# [1,]     4    10     4     4    11     3     9     3     5     4     2     1
# 
# 133   134   136   140
# [1,]     1     3     1     1

##
saveRDS(sce, here(out,'sce_class1_noHATAGABAAmy.rds'))
saveRDS(spe, here(out,'spe.rds'))
saveRDS(spg, here(out,'spg.rds'))
