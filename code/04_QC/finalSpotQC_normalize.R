# cd /dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/
suppressPackageStartupMessages(library("spatialLIBD"))
suppressPackageStartupMessages(library("SpatialExperiment"))
suppressPackageStartupMessages(library("here"))

load(here("processed-data", "04_QC", "spe_QC.Rdata"), verbose = TRUE)
dim(spe)
# [1]  30359 137442

spe$discard_auto_br <- spe$low_sum_br | spe$low_detected_br
spe$discard_auto_id <- spe$low_sum_id | spe$low_detected_id

table(spe$discard_auto_br)
#  FALSE   TRUE 
# 135215   2227 
table(spe$discard_auto_id)
#  FALSE   TRUE 
# 135640   1802 

# remove spots with overlapping tissue sections
# spebr3942 = spe[,which(spebrnum == "Br3942")]
# spebr3942$tissue_sec_overl0ap <- FALSE
# spebr3942$barcode = rownames(colData(spebr3942))
# spebr3942$barcode = sapply(strsplit(spebr3942$barcode,split="_V"),`[`,1)
# 
# tt = c("TACCGATCCAACACTT-1","GATAAGGGACGATTAG-1","TGTTGGCTGGCGGAAG-1","GCGAGGGACTGCTAGA-1",
#        "GCGCGTTTAAATCGTA-1","GACGACTTTCCAAGAA-1","ATACCCTGGCTCAAAT-1","TAACCGTCCAGTTCAT-1",
#        "CAAGGGAGTGTATTTG-1")
# spebr3942$tissue_sec_overlap[spebr3942$barcode == tt] = TRUE

# remove combined set of low-quality spots
spe <- spe[, colData(spe)$discard_auto_id == FALSE]
dim(spe)
# [1]  30359 135640

# calculate logcounts (log-transformed normalized counts) and store in object
spe <- logNormCounts(spe)
# check
assayNames(spe)

save(spe, file = here::here("processed-data", "04_QC", "spe_QCfinal.Rdata"))
