library(RcppML)
library(singlet)
library(here)
library(sessioninfo)
setwd('/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/')

load(here::here('snRNAseq_hpc','processed-data','NMF','logcounts.rda'))
cvnmf<-cross_validate_nmf(
    mtx,
    ranks=c(5,10,50,100,125,150,200),
    n_replicates = 3,
    tol = 1e-03,
    maxit = 100,
    verbose = 3,
    L1 = 0.1,
    L2 = 0,
    threads = 0,
    test_density = 0.2
)

pdf(file=here::here('plots','figures','supp_figures','cvnmf_lineplot.pdf'),h=3,w=8)
df<-cvnmf
filtered_df <- df %>%
     group_by(k, rep) %>%
     filter(iter == max(iter)) %>%
     ungroup()

ggplot(filtered_df, aes(x = k, y = test_error, color = as.factor(rep))) +
    geom_line() +
    labs(title = "Test Error vs. k for Different Reps",
         x = "k",
         y = "Test Error",
         color = "Rep") +
     theme_minimal()
dev.off()


save(cvnmf,file=here::here('snRNAseq_hpc','processed-data','NMF','cvnmf_3reps.rda'))
