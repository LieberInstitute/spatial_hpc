library("rsconnect")
library("here")

app_dir <- here::here("code", "03_shiny_app_basic")

## Or you can go to your shinyapps.io account and copy this
## Here we do this to keep our information hidden.
##
# load(here("code", "03_shiny_app_basic", ".deploy_info.Rdata"), verbose = TRUE)
# rsconnect::setAccountInfo(
#   name = deploy_info$name,
#   token = deploy_info$token,
#   secret = deploy_info$secret
# )

## You need this to enable shinyapps to install Bioconductor packages
options(repos = BiocManager::repositories())

## Deploy the app, that is, upload it to shinyapps.io
rsconnect::deployApp(
    appDir = here("code", "03_shiny_app_basic", "01_spatialLIBD_app_round9"),
    appFiles = c(
        "app.R",
        "spe_round9.rds"
    ),
    appName = "Visium_HPC_round9",
    account = "libd",
    server = "shinyapps.io"
)
