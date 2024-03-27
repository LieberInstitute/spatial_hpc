setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")

library(ggplot2)
library(here)

df = read.csv(file = here::here("processed-data", "06_Clustering", "BayesSpace", "1st_run", "fasthplus_results.csv"))
pdf(file = here::here("plots", "06_Clustering","BayesSpace_fasthplus.pdf"))
ggplot(data = df, aes(x = k, y = fasthplus, group = 1)) +
  geom_line() +
  geom_point()
dev.off()
