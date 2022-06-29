setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")

library('basilisk')
library('zellkonverter')
library('spatialLIBD')
library('here')
library('sessioninfo')

load(file = here::here("processed-data", "05_Batch_correction", "spe_harmony.Rdata"),verbose = TRUE)
colData(spe)$dateImg = NULL
anndata_out = here("processed-data", "06_Clustering", "spe_harmony_anndata.h5ad")

#  Append 'spatialCoords' and 'spatialData' slots to 'colData', since in
#  conversion we're treating the spatialExperiment object as if it is a
#  singleCellExperiment, which doesn't have those additional slots.
colData(spe) = cbind(colData(spe), spatialCoords(spe))

write_anndata = function(sce, out_path) {
  invisible(
    basiliskRun(
      fun = function(sce, filename) {
        library('zellkonverter')
        library('reticulate')
        
        # Convert SCE to AnnData:
        adata <- SCE2AnnData(sce)
        
        #  Write AnnData object to disk
        adata$write(filename = filename)
        
        return()
      },
      env = zellkonverterAnnDataEnv(),
      sce = sce,
      filename = out_path
    )
  )
}

write_anndata(spe, anndata_out)

session_info()