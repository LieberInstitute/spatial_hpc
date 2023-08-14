#-------------------------------------------------------------------------------
#   Add cell counts to HE spatial object
#-------------------------------------------------------------------------------

spe_in <- here("processed-data","02_build_spe","spe_nmf_final.rda")
sce_in <- "/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/snRNAseq_hpc/processed-data/sce/sce_final.rda"
out <- here("processed-data", "spot_deconvo", "shared_utilities")
  
load(spe_in)

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
reducedDims(spe)$spatial <- spatialCoords(spe)


load(sce_in, verbose = TRUE)
rownames(sce) <- rowData(sce)$gene_id

saveRDS(sce, here(out,sce.rds))
saveRDS(spe, here(out,spe.rds))