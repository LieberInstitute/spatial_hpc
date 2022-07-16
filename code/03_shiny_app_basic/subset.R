library("spatialLIBD")
library("lobstr")
library("here")

## Load the main spe object
load(here::here("processed-data", "02_build_spe", "spe_basic.Rdata"),
    verbose = TRUE
)
## Check how big it is in memory
lobstr::obj_size(spe) 
# 5.075627 B
## That's too big for shinyapps.io. Aim to have an object near 2GB.

## Subset the spe object outside of shinyapps.io. Otherwise, the peak memory is
## still affected by loading the object.
## Also, running lobstr::obj_size() takes a while to run, which we don't need
## to run every time someone accesses the shiny app.
imgData(spe) <-
    imgData(spe)[!imgData(spe)$image_id %in% c("hires", "detected", "aligned"), ]
lobstr::obj_size(spe) 
# 2.351571 B
## Ok, this seems reasonable.

## Save the reduced version of the spe object in the shiny app directory
## instead of using soft links.
save(spe, file = here::here("code", "03_shiny_app_basic", "spe.Rdata"))
