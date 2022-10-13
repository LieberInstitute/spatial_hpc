setwd('/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/')
suppressPackageStartupMessages({
library("here")
library("SpatialExperiment")
library("Seurat")
})

spe_to_seurat <- function(spe){
  ret <- CreateSeuratObject(
    counts=assays(spe)$counts,
    meta.data=data.frame(
      row=spatialCoords(spe)[,1],
      col=spatialCoords(spe)[,2]),
    spot_id = colData(spe)$spot_id
  )
return(ret)
}


# Multiple sample
spe_to_seuratList <- function(spe){
  uniq_sample_id <- colData(spe)$sample_id |> unique()
  
  # Create a seurate object for each unique sample_id
  map(uniq_sample_id,
      .f = function(smp_id, spe){
        # browser()
        ret_spe <- spe[, colData(spe)$sample_id == smp_id]
        ret_seurat <- spe_to_seurat(ret_spe)
        
        return(ret_seurat)
      },
      spe = spe)
}