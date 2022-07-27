

load(here::here("processed-data","07_Feature_selection", "nnSVG", "nnSVG.Rdata"))
load(file = here::here("processed-data", "06_Clustering", "spe_modify.Rdata"))

sample_part_ids = unique(spe$sample_id)

# number of genes that passed filtering for each sample-part
sapply(res_list, nrow)


# match results from each sample-part and store in correct rows
res_ranks <- matrix(NA, nrow = nrow(spe), ncol = length(sample_part_ids))
rownames(res_ranks) <- rownames(spe)
colnames(res_ranks) <- sample_part_ids

for (s in seq_along(sample_part_ids)) {
  stopifnot(colnames(res_ranks)[s] == sample_part_ids[s])
  stopifnot(colnames(res_ranks)[s] == names(res_list)[s])
  
  rownames_s <- rownames(res_list[[s]])
  res_ranks[rownames_s, s] <- res_list[[s]][, "rank"]
}

# keep only genes that were not filtered out in all sample-parts
res_ranks <- na.omit(res_ranks)

# calculate average ranks
avg_ranks <- sort(rowMeans(res_ranks))

# summary table
df_summary <- data.frame(
  gene_id = names(avg_ranks), 
  gene_name = rowData(spe)[names(avg_ranks), "gene_name"], 
  gene_type = rowData(spe)[names(avg_ranks), "gene_type"], 
  avg_rank = unname(avg_ranks), 
  row.names = names(avg_ranks)
)

head(df_summary, 20)
mySVGs = head(df_summary$gene_name, 20)
tonySVGs <- c(
  "MBP",
  "PLP1",
  "GFAP",
  "MTRNR2L12",
  "TTR",
  "SNAP25",
  "SLC17A7",
  "UCHL1",
  "NPTXR",
  "NPTX2",
  "THY1",
  "HPCA",
  "ENC1",
  "NPTX1",
  "CHN1",
  "NRGN",
  "YWHAH",
  "OLFM1",
  "NNAT",
  "NCDN",
  "CRYAB",
  "CST3"
)

SVG = setdiff(mySVGs, tonySVGs)

# Locate the marker genes
SVG_search <- rowData(spe)$gene_search[match(SVG, rowData(spe)$gene_name)]
# brains <- unique(spe$brnum)
brains = c("Br6423","Br6432","Br2743","Br8325","Br3942","Br6471","Br8667","Br8492","Br6522")
for (i in SVG_search) {
  gene <- i
  gene_name <- strsplit(i, ";")[[1]][1]
  pdf(here("plots", "07_Feature_selection", "nnSVG", paste0(gene_name, ".pdf")), width = 21, height = 20)
  
  ii <- 1
  speb <- spe[, which(spe$brnum == brains[ii])]
  samples <- unique(speb$sample_id)
  samples
  p1 <- vis_gene(spe = speb, sampleid = samples[1], geneid = gene, spatial = TRUE, assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2, ... = paste0("_", brains[ii]))
  p2 <- vis_gene(spe = speb, sampleid = samples[2], geneid = gene, spatial = TRUE, assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2, ... = paste0("_", brains[ii]))
  p3 <- vis_gene(spe = speb, sampleid = samples[3], geneid = gene, spatial = TRUE, assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2, ... = paste0("_", brains[ii]))
  p4 <- vis_gene(spe = speb, sampleid = samples[4], geneid = gene, spatial = TRUE, assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2, ... = paste0("_", brains[ii]))
  
  grid.arrange(p1, p2, p3, p4, nrow = 2)
  
  ##
  ii <- 2
  speb <- spe[, which(spe$brnum == brains[ii])]
  samples <- unique(speb$sample_id)
  samples
  
  p1 <- vis_gene(spe = speb, sampleid = samples[1], geneid = gene, spatial = TRUE, assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2, ... = paste0("_", brains[ii]))
  p2 <- vis_gene(spe = speb, sampleid = samples[2], geneid = gene, spatial = TRUE, assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2, ... = paste0("_", brains[ii]))
  p3 <- vis_gene(spe = speb, sampleid = samples[1], geneid = gene, spatial = TRUE, assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2, ... = paste0("_", brains[ii]))
  p4 <- vis_gene(spe = speb, sampleid = samples[2], geneid = gene, spatial = TRUE, assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2, ... = paste0("_", brains[ii]))
  
  grid.arrange(p1, p3, p2, p4, nrow = 2)
  
  ##
  ii <- 3
  speb <- spe[, which(spe$brnum == brains[ii])]
  samples <- unique(speb$sample_id)
  samples
  
  p1 <- vis_gene(spe = speb, sampleid = samples[3], geneid = gene, spatial = TRUE, assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2, ... = paste0("_", brains[ii]))
  p2 <- vis_gene(spe = speb, sampleid = samples[1], geneid = gene, spatial = TRUE, assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2, ... = paste0("_", brains[ii]))
  p3 <- vis_gene(spe = speb, sampleid = samples[3], geneid = gene, spatial = TRUE, assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2, ... = paste0("_", brains[ii]))
  p4 <- vis_gene(spe = speb, sampleid = samples[1], geneid = gene, spatial = TRUE, assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2, ... = paste0("_", brains[ii]))
  
  grid.arrange(p1, p2, p3, p4, nrow = 2)
  
  p1 <- vis_gene(spe = speb, sampleid = samples[2], geneid = gene, spatial = TRUE, assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2, ... = paste0("_", brains[ii]))
  p2 <- vis_gene(spe = speb, sampleid = samples[4], geneid = gene, spatial = TRUE, assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2, ... = paste0("_", brains[ii]))
  p3 <- vis_gene(spe = speb, sampleid = samples[2], geneid = gene, spatial = TRUE, assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2, ... = paste0("_", brains[ii]))
  p4 <- vis_gene(spe = speb, sampleid = samples[4], geneid = gene, spatial = TRUE, assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2, ... = paste0("_", brains[ii]))
  
  grid.arrange(p1, p3, p2, p4, nrow = 2)
  
  ##
  ii <- 4
  speb <- spe[, which(spe$brnum == brains[ii])]
  samples <- unique(speb$sample_id)
  samples
  
  p1 <- vis_gene(spe = speb, sampleid = samples[1], geneid = gene, spatial = TRUE, assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2, ... = paste0("_", brains[ii]))
  p2 <- vis_gene(spe = speb, sampleid = samples[2], geneid = gene, spatial = TRUE, assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2, ... = paste0("_", brains[ii]))
  p3 <- vis_gene(spe = speb, sampleid = samples[5], geneid = gene, spatial = TRUE, assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2, ... = paste0("_", brains[ii]))
  p4 <- vis_gene(spe = speb, sampleid = samples[3], geneid = gene, spatial = TRUE, assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2, ... = paste0("_", brains[ii]))
  p5 <- vis_gene(spe = speb, sampleid = samples[4], geneid = gene, spatial = TRUE, assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2, ... = paste0("_", brains[ii]))
  grid.arrange(p1, p2, p3, p4, p5, nrow = 2)
  
  ##
  ii <- 5
  speb <- spe[, which(spe$brnum == brains[ii])]
  samples <- unique(speb$sample_id)
  samples
  
  p1 <- vis_gene(spe = speb, sampleid = samples[1], geneid = gene, spatial = TRUE, assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2, ... = paste0("_", brains[ii]))
  p2 <- vis_gene(spe = speb, sampleid = samples[2], geneid = gene, spatial = TRUE, assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2, ... = paste0("_", brains[ii]))
  p3 <- vis_gene(spe = speb, sampleid = samples[3], geneid = gene, spatial = TRUE, assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2, ... = paste0("_", brains[ii]))
  p4 <- vis_gene(spe = speb, sampleid = samples[4], geneid = gene, spatial = TRUE, assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2, ... = paste0("_", brains[ii]))
  
  grid.arrange(p1, p2, p3, p4, nrow = 2)
  
  ##
  ii <- 6
  speb <- spe[, which(spe$brnum == brains[ii])]
  samples <- unique(speb$sample_id)
  samples
  
  p1 <- vis_gene(spe = speb, sampleid = samples[1], geneid = gene, spatial = TRUE, assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2, ... = paste0("_", brains[ii]))
  p2 <- vis_gene(spe = speb, sampleid = samples[2], geneid = gene, spatial = TRUE, assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2, ... = paste0("_", brains[ii]))
  p3 <- vis_gene(spe = speb, sampleid = samples[3], geneid = gene, spatial = TRUE, assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2, ... = paste0("_", brains[ii]))
  
  grid.arrange(p1, p2, p3, nrow = 2)
  
  ##
  ii <- 7
  speb <- spe[, which(spe$brnum == brains[ii])]
  samples <- unique(speb$sample_id)
  samples
  
  p1 <- vis_gene(spe = speb, sampleid = samples[2], geneid = gene, spatial = TRUE, assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2, ... = paste0("_", brains[ii]))
  p2 <- vis_gene(spe = speb, sampleid = samples[1], geneid = gene, spatial = TRUE, assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2, ... = paste0("_", brains[ii]))
  p3 <- vis_gene(spe = speb, sampleid = samples[4], geneid = gene, spatial = TRUE, assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2, ... = paste0("_", brains[ii]))
  p4 <- vis_gene(spe = speb, sampleid = samples[3], geneid = gene, spatial = TRUE, assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2, ... = paste0("_", brains[ii]))
  
  grid.arrange(p1, p2, p3, p4, nrow = 2)
  
  ##
  ii <- 8
  speb <- spe[, which(spe$brnum == brains[ii])]
  samples <- unique(speb$sample_id)
  samples
  
  p1 <- vis_gene(spe = speb, sampleid = samples[1], geneid = gene, spatial = TRUE, assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2, ... = paste0("_", brains[ii]))
  p2 <- vis_gene(spe = speb, sampleid = samples[2], geneid = gene, spatial = TRUE, assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2, ... = paste0("_", brains[ii]))
  p3 <- vis_gene(spe = speb, sampleid = samples[1], geneid = gene, spatial = TRUE, assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2, ... = paste0("_", brains[ii]))
  p4 <- vis_gene(spe = speb, sampleid = samples[2], geneid = gene, spatial = TRUE, assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2, ... = paste0("_", brains[ii]))
  
  grid.arrange(p1, p2, p3, p4, nrow = 2)
  
  ##
  ii <- 9
  speb <- spe[, which(spe$brnum == brains[ii])]
  samples <- unique(speb$sample_id)
  samples
  
  p1 <- vis_gene(spe = speb, sampleid = samples[1], geneid = gene, spatial = TRUE, assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2, ... = paste0("_", brains[ii]))
  p2 <- vis_gene(spe = speb, sampleid = samples[2], geneid = gene, spatial = TRUE, assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2, ... = paste0("_", brains[ii]))
  p3 <- vis_gene(spe = speb, sampleid = samples[3], geneid = gene, spatial = TRUE, assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2, ... = paste0("_", brains[ii]))
  p4 <- vis_gene(spe = speb, sampleid = samples[4], geneid = gene, spatial = TRUE, assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2, ... = paste0("_", brains[ii]))
  
  grid.arrange(p1, p2, p3, p4, nrow = 2)
  
  dev.off()
  message(i)
}

