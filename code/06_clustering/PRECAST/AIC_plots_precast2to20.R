
setwd('/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/')
suppressPackageStartupMessages({
    library(dplyr)
    library(purrr)
    library(Seurat)
    library(SpatialExperiment)
    library(PRECAST)
    library(spatialLIBD)
    library(ggplot2)
    library(gridExtra)
    library("here")
})

##write loop to get all the ICs
kvals=2:20
###load final PRECAST clusters
aicList<-list()
for(k in kvals){
print('loading object ',k)
load(file = here("processed-data", "06_clustering", "PRECAST",
                 paste0("nnSVG_2000_",k,"_HE.Rdata")))
PRECASTObj <- selectModel(PRECASTObj,criteria='AIC')
aicList[k]<-PRECASTObj@resList$icMat[2]
print(aicList[k])
rm(PRECASTObj)
gc()
}

aic<-as.numeric(unlist(aicList))

# Adjust the indices to start from 2
aic_df <- data.frame(Index = 2:20, AIC = aic)

# Plot using ggplot2 with y-axis labels in scientific notation
ggplot(aic_df, aes(x = Index, y = AIC)) +
  geom_line() + # Add line
  geom_point() + # Add points for all
  geom_point(data = aic_df[18:19,], aes(x = Index, y = AIC), color = "red") + # Make the last 3 points red, adjusted for new indexing
  geom_vline(xintercept = 4.9, linetype = "dotted", color = "red", size = 1) + # Add dotted vertical line at 4.9
  geom_vline(xintercept = 17.1, linetype = "dotted", color = "red", size = 1) + # Add another dotted vertical line at 17.1
  scale_y_continuous(labels = scales::scientific) + # Format y-axis labels in scientific notation
  theme_minimal() + # Use a minimal theme
  labs(title = "Line Plot of AIC Values",
       x = "Index",
       y = "AIC Value")


