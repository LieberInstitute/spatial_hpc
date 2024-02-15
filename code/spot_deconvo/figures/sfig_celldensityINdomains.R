setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")
suppressPackageStartupMessages(library("SingleCellExperiment"))
suppressPackageStartupMessages(library("here"))
suppressPackageStartupMessages(library("spatialLIBD"))
suppressPackageStartupMessages(library("ggplot2"))

Dr <- here("processed-data","spot_deconvo","shared_utilities")
plot_dir <- here("plots", "spot_deconvo", "figures", "fig_celldensity")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

spe = readRDS(here(Dr,"spe.rds"))
dat = as.data.frame(colData(spe))
dat = dat[dat$count<50, ]
#levels(dat$cluster_collapsed)[levels(dat$cluster_collapsed)=="WM.1"] <- "WM"
#levels(dat$cluster_collapsed)[levels(dat$cluster_collapsed)=="WM.2"] <- "WM"
#levels(dat$cluster_collapsed)[levels(dat$cluster_collapsed)=="WM.3"] <- "WM"

load('/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/plots/spatial_palettes.rda')
colors = srt.palette
colors = unname(colors)
png(here("plots","spot_deconvo","figures","fig_celldensity", "sfig_celldensityINdomains.png"), width = 1200, height = 600, units = "px") 
p = ggplot(data = dat, aes(x = dat$cluster_collapsed, y=count, fill = cluster_collapsed))+
  geom_boxplot(outlier.shape = NA, show.legend = FALSE)+theme_bw()+scale_fill_manual(values = colors)+
  labs(y = "nuclei count per spot", x="spatial domains")+ 
  theme(text = element_text(size = 20, colour = "black"),
        axis.text = element_text(size = 24, colour = "black"),
        axis.text.x = element_text(angle = 90),
        axis.line = element_line(size=2, colour = "black"),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        panel.border = element_blank())
print(p)
dev.off()

getmode <- function(v) {
  uniqv <- unique(na.omit(v))
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

mean(dat$count[dat$cluster_collapsed == "ML"])
#[1] 2.755694
getmode(dat$count[dat$cluster_collapsed == "ML"])
#[1] 1
median(dat$count[dat$cluster_collapsed == "ML"])
#[1] 2

mean(dat$count[dat$cluster_collapsed == "GCL"])
#[1] 9.716166
getmode(dat$count[dat$cluster_collapsed == "GCL"])
#[1] 5
median(dat$count[dat$cluster_collapsed == "GCL"])
#[1] 9