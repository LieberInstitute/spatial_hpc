setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")

library('here')
library('jaffelab')
library('SingleCellExperiment')
library('ggvenn')
library('UpSetR')

load(file = here::here("processed-data", "08_pseudobulk", "mbkmeans", "DEgenes_brain_cluster15_GCL.Rdata"))
mbkmeans_brain_UP = UP
mbkmeans_brain_DOWN = DOWN
load(file = here::here("processed-data", "08_pseudobulk", "mbkmeans", "DEgenes_captureArea_cluster15_GCL_adjBrnum.Rdata"))
mbkmeans_captureArea_UP = UP
mbkmeans_captureArea_DOWN = DOWN
load(file = here::here("processed-data", "08_pseudobulk", "manual_annotations", "DEgenes_brain_GCL.Rdata"))
manualAnnotation_brain_UP = UP
manualAnnotation_brain_DOWN = DOWN
load(file = here::here("processed-data", "08_pseudobulk", "manual_annotations", "DEgenes_captureArea_GCL.Rdata"))
manualAnnotation_captureArea_UP = UP
manualAnnotation_captureArea_DOWN = DOWN

x = list(
  mbkmeans_brain = mbkmeans_brain_UP$gene,
  mbkmeans_captureArea = mbkmeans_captureArea_UP$gene,
  manualAnnotation_brain = manualAnnotation_brain_UP$gene,
  manualAnnotation_captureArea = manualAnnotation_captureArea_UP$gene
)

pdf(file = here::here("plots", "08_pseudobulk","DE_common_genes.pdf"), width = 10, height = 4)

ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
)
upset(fromList(x), order.by = "freq")
dev.off()

# x
# $mbkmeans_brain
# [1] "PPFIA2" "KCTD4"  "AP1G2"  "SMAD3"  "RFX3"   "THOC3"  "CHRD"   "SEMA5A"
# 
# $mbkmeans_captureArea
# [1] "PPFIA2"   "PKD1"     "TIAM1"    "C1QTNF9B" "KCNK1"    "THOC3"   
# [7] "PPP1R1A"  "GOLGA8A"  "CEBPD"    "TNNT2"    "ITGA7"    "POSTN"   
# [13] "SCART1"  
# 
# $manualAnnotation_brain
# [1] "PPFIA2"     "NECAB3"     "YPEL4"      "ZNF273"     "AP1G2"     
# [6] "AC005837.4"
# 
# $manualAnnotation_captureArea
# [1] "PPFIA2" "ZNF273" "NECAB3" "THOC3"  "ITGA7"  "KCNK1"  "TIAM1"  "CEBPD" 
# [9] "PGA5"  

