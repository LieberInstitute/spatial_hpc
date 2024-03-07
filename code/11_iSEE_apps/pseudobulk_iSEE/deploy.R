library("rsconnect")

#source("token.R")

options(repos = BiocManager::repositories())
rsconnect::deployApp(
    appFiles = c("app.R", "spe_pseudo.rda", "spatial_palettes_isee.rda", "initial.R"),
    appName = "HPC_bulkseq_data",
    account = "libd",
    server = "shinyapps.io"
)
