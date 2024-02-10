library("rsconnect")

# source("token.R")

options(repos = BiocManager::repositories())
rsconnect::deployApp(
    appFiles = c("app.R", "sce_iSEE.rda", "snrna_palettes_isee.rda", "initial.R"),
    appName = "HPC_snRNAseq_data",
    account = "libd",
    server = "shinyapps.io"
)
