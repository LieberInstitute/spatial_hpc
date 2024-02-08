library("rsconnect")

# source("token.R")

options(repos = BiocManager::repositories())
rsconnect::deployApp(
    appFiles = c("app.R", "rse_gene_TrkB_KO_LS_n8_wm.Rdata", "initial.R"),
    appName = "bulkseq_lateral_septum",
    account = "libd",
    server = "shinyapps.io"
)
