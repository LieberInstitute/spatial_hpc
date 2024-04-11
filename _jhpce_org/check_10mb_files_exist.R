## Check that files exist at JHPCE
gitignore <- readLines(".gitignore")
gitignore <- gitignore[!gitignore %in% c(".Rproj.user", ".Rhistory", ".RData", ".Ruserdata")]
jhpce_file_exists <- file.exists(gitignore)
table(jhpce_file_exists)
# jhpce_file_exists
# FALSE  TRUE 
#     1   882
gitignore[!jhpce_file_exists]
# [1] "plots/02_build_spe/referenceMapping.pdf"

## At JHPCE I ran:
# chmod g+w plots/02_build_spe/
## On my laptop I then ran:
# scp plots/02_build_spe/referenceMapping.pdf et:/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/plots/02_build_spe/

table(jhpce_file_exists)
# jhpce_file_exists
# TRUE 
#  883 
