###############################################################################
#  Functions
###############################################################################

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

