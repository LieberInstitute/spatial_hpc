library(SpatialExperiment)
library(scater)
library(RcppML)
library(ggspavis)
library(here)
print('Loading spe!')
load(file=here("processed-data","02_build_spe","spe_bayes_clus.Rdata"))
counts<-assay(spe,'logcounts')
library(Matrix)
counts_sparse <- as(counts, "dgCMatrix")
print('Running NMF!')
x<-nmf(counts,
50,
tol = 1e-06,
maxit = 1000,
verbose = TRUE,
seed = 1512,
L1 = c(0, 0),
mask_zeros = FALSE,
diag = TRUE,
nonneg = TRUE
)

save(x,file=here("processed-data","NMF","RcppML_nmf_firstTry.rda"))

##bind patterns
pats<-t(x$h)
colnames(pats)<-paste("Pattern", 1:50, sep = "_")
colData(spe)<-cbind(colData(spe),pats)

patterns <- colnames(colData(spe))[75:124]

brains <- unique(spe$brnum)

p<-list()
for (i in seq_along(patterns)) {
    p[[i]]<-list()
    for (j in seq_along(brains)) {
        speb <- spe[, (colData(spe)$brnum == brains[j])]
        speb$sample_id <- droplevels(speb$sample_id)
        speb$sample_id <- as.character(speb$sample_id)
        samples <- unique(speb$sample_id)
        speb$sample_id <- factor(speb$sample_id, levels = samples)
        samples
        speb$brnum <- droplevels(speb$brnum)

        print(paste0("printing plot",i,"_",j))
        p[[i]][[j]]<-plotVisium(
            speb,
            spots = TRUE,
            fill = patterns[i],
            highlight = NULL,
            facets = "sample_id",
            image = TRUE,
            assay = "logcounts",
            trans = "identity",
            x_coord = NULL,
            y_coord = NULL,
            y_reverse = TRUE,
            sample_ids = NULL,
            image_ids = NULL,
            palette = "viridis"
        )+ ggplot2::theme(legend.text = ggplot2::element_text(size = 12),
                           plot.title = ggplot2::element_text(size = 18))
    }}
    
for (i in seq_along(patterns)) {
    pdf(
        file = here::here("plots", "NMF",
                          paste0(patterns[i], ".pdf")))
    print("printing next one")

    for (j in seq_along(p[[i]])) {
        print(p[[i]][[j]])
    }

    dev.off()
}

save(spe,file=here("processed-data","02_build_spe,"spe_nmf.Rdata")