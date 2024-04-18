library("here")
library("withr")

with_dir(here(), {
    rmarkdown::render("README.Rmd", "html_document")
    system(paste0("mv README.html ", file.path(dirname(here::here()), "spatial_hpc_website", "index.html")))
})

with_dir(
    file.path(dirname(here::here()), "spatial_hpc_website"),
    system("git ci -am -'Updated website with code/update_website.R'; git push origin gh-pages")
)
