setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")
library(here)

load(here::here("processed-data","07_Feature_selection", "nnSVG", "nnSVG.Rdata"))
load(file = here::here("processed-data", "06_Clustering", "spe_modify.Rdata"))

brains = unique(spe$brnum)

df_summary <- as.list(rep(NA, length(brains)))
names(df_summary) <- brains

for (ii in seq_along(brains)){
speb <- spe[, which(spe$brnum == brains[ii])]
sample_part_ids = unique(speb$sample_id)

# number of genes that passed filtering for each sample-part
res_list1 = res_list[which(names(res_list) %in% sample_part_ids)]
sapply(res_list1, nrow)

# match results from each sample-part and store in correct rows
res_ranks <- matrix(NA, nrow = nrow(speb), ncol = length(sample_part_ids))
rownames(res_ranks) <- rownames(speb)
colnames(res_ranks) <- sample_part_ids

for (s in seq_along(sample_part_ids)) {
   
  stopifnot(colnames(res_ranks)[s] == sample_part_ids[s])
  stopifnot(colnames(res_ranks)[s] == names(res_list1)[s])
  
  rownames_s <- rownames(res_list1[[s]])
  res_ranks[rownames_s, s] <- res_list1[[s]][, "rank"]
}

# keep only genes that were not filtered out in all sample-parts
res_ranks <- na.omit(res_ranks)

# calculate average ranks
avg_ranks <- sort(rowMeans(res_ranks))

# summary table
df_summary[[ii]] <- data.frame(
  gene_id = names(avg_ranks), 
  gene_name = rowData(speb)[names(avg_ranks), "gene_name"], 
  gene_type = rowData(speb)[names(avg_ranks), "gene_type"], 
  avg_rank = unname(avg_ranks), 
  row.names = names(avg_ranks)
)

}

temp = lapply(df_summary, head, 20)
result = lapply(temp, "[", ,"gene_name")

Reduce(intersect,result)
# [1] "MBP"     "GFAP"    "PLP1"    "ENC1"    "SLC17A7" "SNAP25"  "UCHL1"  

as.data.frame(result)
# Br2743  Br3942    Br6423  Br6432    Br6471    Br6522    Br8325    Br8492
# 1      MBP     MBP       MBP     MBP       MBP       TTR       MBP       MBP
# 2     GFAP    PLP1      GFAP    GFAP MTRNR2L12       MBP      GFAP      GFAP
# 3     PLP1    GFAP      PLP1    PLP1      GFAP      GFAP      PLP1 MTRNR2L12
# 4     ENC1  SNAP25      NRGN    NRGN       TTR      PLP1    SNAP25      PLP1
# 5  SLC17A7    FTH1      ENC1  SNAP25    SNAP25   SLC17A7      NRGN    CARTPT
# 6     NRGN    NRGN MTRNR2L12   CRYAB      PLP1      ENC1      ENC1    SNAP25
# 7   SNAP25   UCHL1      NNAT   VSNL1      ENC1 MTRNR2L12     UCHL1      NCDN
# 8     CHN1   CRYAB     UCHL1   UCHL1     UCHL1    SNAP25   SLC17A7      CST3
# 9     THY1  TMSB10    SNAP25    NNAT   SLC17A7     UCHL1      THY1     UCHL1
# 10    MOBP SLC17A7   SLC17A7 SLC17A7      CST3      NRGN     S100B       HBB
# 11    CST3   S100B      THY1  TUBB2A     YWHAH      THY1      CHN1       NPY
# 12   UCHL1    PCP4      MOBP   CALM3      THY1      CHN1      NEFL      HBA2
# 13    NEFL    CHN1       CNP    NEFL      HBA2       MT3     CRYAB      ENC1
# 14     NSF    ENC1     YWHAH    ENC1      NCDN      RTN1      MOBP      CHGB
# 15   YWHAH   CALM3     CALM3    RTN1       HBB     NPTXR      FTH1  MTRNR2L1
# 16  TUBB2A    MOBP    MALAT1   S100B    RPL37A     CALM3      SYT1     YWHAH
# 17    SYT1 C9orf16      NEFL  TUBA1B      HPCA      HBA2      HPCA   SLC17A7
# 18   NPTXR  TUBB2A     NPTX1    THY1      CHN1      CST3    TUBB2A     NPTX1
# 19    RTN1    THY1    TUBB2A    FTH1  HSP90AA1     YWHAZ MTRNR2L12    PPFIA2
# 20   YWHAG  DBNDD2     VSNL1    CHN1    ATP1B1      MOBP     CALM3      HBA1
# Br8667
# 1        MBP
# 2       PLP1
# 3    SLC17A7
# 4     SNAP25
# 5       GFAP
# 6      NPTXR
# 7       NRGN
# 8       HPCA
# 9  MTRNR2L12
# 10     UCHL1
# 11      ENC1
# 12      THY1
# 13      CHN1
# 14     CRYAB
# 15     NPTX1
# 16      FTH1
# 17       CCK
# 18     OLFM1
# 19     YWHAH
# 20   C9orf16