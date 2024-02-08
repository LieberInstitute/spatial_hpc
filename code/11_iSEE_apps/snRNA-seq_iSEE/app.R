library("SingleCellExperiment")
library("iSEE")
library("shiny")
library("paletteer")
library("scuttle")

load("sce_iSEE.rda", verbose = TRUE)
load("snrna_palettes_isee.rda")

## Don't run this on app.R since we don't want to run this every single time
#lobstr::obj_size(sce)




source("initial.R", print.eval = TRUE)


#rse_gene <- registerAppOptions(rse_gene, color.maxlevels = length(Sample_ID)
iSEE(
    sce,
    appTitle = "HPC snRNAseq data",
    initial = initial,
    colormap = ExperimentColorMap(colData = list(
        fine.cell.class = function(n) {
            cols <- sn.fine.palette
        },
        mid.cell.class = function(n) {
            cols <- sn.mid.palette
        },
        broad.cell.class = function(n) {
            cols <- sn.broad.palette
        }

    ))
)
