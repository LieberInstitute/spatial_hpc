bayesSpace_name <- colnames(colData(spe))[72]

brains <- unique(spe$brnum)

k<-12
p<-list()
for (i in seq_along(bayesSpace_name)) {
    p[[i]]<-list()
    for (j in seq_along(brains)) {
        speb <- spe[, (colData(spe)$brnum == brains[j])]
        speb$sample_id <- droplevels(speb$sample_id)
        speb$sample_id <- as.character(speb$sample_id)
        samples <- unique(speb$sample_id)
        speb$sample_id <- factor(speb$sample_id, levels = samples)
        samples
        speb$brnum <- droplevels(speb$brnum)

        cols <- Polychrome::palette36.colors(k[i])
        names(cols) <- levels(speb[[bayesSpace_name[i]]])
        print(paste0("printing plot",i,"_",j))
        p[[i]][[j]]<-plotVisium(
            speb,
            spots = TRUE,
            fill = bayesSpace_name[i],
            highlight = NULL,
            facets = "sample_id",
            image = FALSE,
            assay = "logcounts",
            trans = "identity",
            x_coord = NULL,
            y_coord = NULL,
            y_reverse = TRUE,
            sample_ids = NULL,
            image_ids = NULL,
            palette = cols
        )+ ggplot2::theme(legend.text = ggplot2::element_text(size = 12),
                           plot.title = ggplot2::element_text(size = 18))
    }}

for (i in seq_along(bayesSpace_name)) {
    pdf(
        file = here::here("plots", "06_clustering", "PRECAST",
                          paste0(bayesSpace_name[i], ".pdf")))
    print("printing next one")

    for (j in seq_along(p[[i]])) {
        print(p[[i]][[j]])
    }

    dev.off()
}



clus<-as.list(colData(spe)[,c(20,50:74)])
heat<-compareClusterings(clus,T)
heat<-ifelse(heat==0,1,heat)
pdf(here::here('plots','06_clustering','ARI.pdf'))
pheatmap::pheatmap(heat)
dev.off()

library(tidyr)
library(ggplot2)

# Create a data frame to store ARI values for each brnum and clustering
ari_df <- data.frame()

for (br in levels(speb$brnum)) {
    speb_subset <- speb[, speb$brnum == br]
    ref <- colData(speb_subset)[, 20]
    clus <- as.list(colData(speb_subset)[, c(50:74)])

    # Calculate ARI values for each clustering against the reference
    ari_values <- sapply(clus, function(x) pairwiseRand(ref, x, mode = "index", adjusted = TRUE))

    # Convert the ARI values to a data frame and add a brnum column
    ari_br <- data.frame(t(ari_values))
    colnames(ari_br) <- colnames(colData(speb_subset)[, 50:74])
    ari_br$brnum <- br

    # Combine the ARI data frame with the current ari_br
    ari_df <- rbind(ari_df, ari_br)
}

# Convert the data frame to a long format
ari_long <- gather(ari_df, clustering, ARI, -brnum)

# Add a column for the algorithm
ari_long$algorithm <- ifelse(grepl("bayes", ari_long$clustering), "BayesSpace", "PRECAST")

# Create a boxplot
pdf(here::here('plots','ARI_boxplot.pdf'))
ggplot(ari_long, aes(x = clustering, y = ARI, fill = algorithm)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    labs(title = "ARI Comparison for Different Clusterings vs. Reference", x = "Clustering", y = "ARI") +
    facet_wrap(~ algorithm, scales = "free_x", strip.position = "bottom") +
    scale_fill_manual(values = c("BayesSpace" = "#6da7de", "PRECAST" = "#de6d85"))
dev.off()


# Create a boxplot

pdf(here::here('plots','ARI_boxplot.pdf'))
ggplot(ari_long, aes(x = clustering, y = ARI)) +
    geom_boxplot() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    labs(title = "ARI Comparison for Different Clusterings vs. Reference", x = "Clustering", y = "ARI")
dev.off()