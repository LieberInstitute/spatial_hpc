setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")
suppressPackageStartupMessages(library("basilisk"))
suppressPackageStartupMessages(library("SpatialExperiment"))
suppressPackageStartupMessages(library("SingleCellExperiment"))
suppressPackageStartupMessages(library("zellkonverter"))
suppressPackageStartupMessages(library("sessioninfo"))
suppressPackageStartupMessages(library("here"))

spe_in <- here("processed-data", "NMF", "spe_nmf_final.rda")
spe_out <- here("processed-data", "10-image_stitching", "spe.h5ad")

write_anndata <- function(sce, out_path) {
  invisible(
    basiliskRun(
      fun = function(sce, filename) {
        library("zellkonverter")
        library("reticulate")
        
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

load(spe_in, verbose = TRUE)
#   zellkonverter doesn't know how to convert the 'spatialCoords' slot. We'd
#   ultimately like the spatialCoords in the .obsm['spatial'] slot of the
#   resulting AnnData, which corresponds to reducedDims(spe)$spatial in R
reducedDims(spe)$spatial <- spatialCoords(spe)
colData(spe)$dateImg = as.character(colData(spe)$dateImg)

write_anndata(spe, spe_out)

session_info()
