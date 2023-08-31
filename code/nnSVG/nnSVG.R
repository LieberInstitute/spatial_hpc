library(SpatialExperiment)
library(here)
library(nnSVG)
library(dplyr)
library(scater)
library(scran)

load(file=here::here("processed-data","05_preprocess_batchCorrection","spe_norm.rda"))
dim(spe)
##filter out mito genes...just trying something here
spe[!which(seqnames(spe)=='chrM'),]
dim(spe)
##filter genes for lowly expressed guys
n_umis <- 100
ix_low_genes <- rowSums(counts(spe)) < n_umis
table(ix_low_genes)

spe <- spe[!ix_low_genes, ]
dim(spe)

ix_zero_genes <- rowSums(counts(spe)) == 0
table(ix_zero_genes)

if (sum(ix_zero_genes) > 0) {
    spe <- spe[!ix_zero_genes, ]
}

dim(spe)
colData(spe)$sample_id<-factor(colData(spe)$sample_id)

sample_ids <- levels(colData(spe)$sample_id)

res_list <- as.list(rep(NA, length(sample_ids)))
names(res_list) <- sample_ids

for (s in seq_along(sample_ids)) {

    # select sample
    ix <- colData(spe)$sample_id == sample_ids[s]
    spe_sub <- spe[, ix]

    dim(spe_sub)

    # run nnSVG filtering for mitochondrial genes and low-expressed genes
    # (note: set 'filter_mito = TRUE' in most datasets)
    spe_sub <- filter_genes(
        spe_sub,
        filter_genes_ncounts = 3,
        filter_genes_pcspots = 0.5,
        filter_mito=F
    )

    # remove any zeros introduced by filtering
    ix_zeros <- colSums(counts(spe_sub)) == 0
    if (sum(ix_zeros) > 0) {
        spe_sub <- spe_sub[, !ix_zeros]
    }

    dim(spe_sub)

    # re-calculate logcounts after filtering
    spe_sub <- computeLibraryFactors(spe_sub)
    spe_sub <- logNormCounts(spe_sub)

    # run nnSVG
    set.seed(123)
    spe_sub <- nnSVG(spe_sub,n_threads=12,verbose=F)

    # store results for this sample
    res_list[[s]] <- rowData(spe_sub)
    rm(spe_sub)
}

sapply(res_list, nrow)

# match results from each sample and store in matching rows
res_ranks <- matrix(NA, nrow = nrow(spe), ncol = length(sample_ids))
rownames(res_ranks) <- rownames(spe)
colnames(res_ranks) <- sample_ids

for (s in seq_along(sample_ids)) {
    stopifnot(colnames(res_ranks)[s] == sample_ids[s])
    stopifnot(colnames(res_ranks)[s] == names(res_list)[s])

    rownames_s <- rownames(res_list[[s]])
    res_ranks[rownames_s, s] <- res_list[[s]][, "rank"]
}


# remove genes that were filtered out in all samples
ix_allna <- apply(res_ranks, 1, function(r) all(is.na(r)))
res_ranks <- res_ranks[!ix_allna, ]

dim(res_ranks)

avg_ranks <- rowMeans(res_ranks, na.rm = TRUE)


# calculate number of samples where each gene is within top 100 ranked SVGs
# for that sample
n_withinTop1000 <- apply(res_ranks, 1, function(r) sum(r <= 1000, na.rm = TRUE))


# summary table
df_summary <- data.frame(
    gene_id = names(avg_ranks),
    gene_name = rowData(spe)[names(avg_ranks), "gene_name"],
    gene_type = rowData(spe)[names(avg_ranks), "gene_type"],
    overall_rank = rank(avg_ranks),
    average_rank = unname(avg_ranks),
    n_withinTop1000 = unname(n_withinTop1000),
    row.names = names(avg_ranks)
)

# sort by average rank
df_summary <- df_summary[order(df_summary$average_rank), ]
head(df_summary)

# calculate number of samples where each gene is within top 1000 ranked SVGs
# for that sample
n_withinTop1000 <- apply(res_ranks, 1, function(r) sum(r <= 1000, na.rm = TRUE))


# top n genes
# (note: NAs in this example due to subsampling genes for faster runtime)
top1000genes <- df_summary$gene_name[1:1000]

# summary table of "replicated" SVGs (i.e. genes that are highly ranked in at
# least x samples)
df_summaryReplicated <- df_summary[df_summary$n_withinTop1000 >= 2, ]

# re-calculate rank within this set
df_summaryReplicated$overall_rank <- rank(df_summaryReplicated$average_rank)

dim(df_summaryReplicated)

head(df_summaryReplicated)

# top "replicated" SVGs
topSVGsReplicated <- df_summaryReplicated$gene_name

save(res_ranks,topSVGsReplicated,df_summaryReplicated,df_summary,top1000genes,
     file=here::here('processed-data','nnSVG','nnSVG_outs_mitoFilter.rda'))


