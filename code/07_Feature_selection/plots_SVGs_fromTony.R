gene = "NPTX2; ENSG00000106236"

human_markers_search <- rowData(spe)$gene_search[match("NPTX2", rowData(spe)$gene_name)]

ii = 1
speb = spe[,which(spe$brnum == brains[ii])]
samples = unique(speb$sample_id)
samples
p1 = vis_gene(spe = speb,sampleid = samples[1],geneid = gene, spatial = TRUE,assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2,... = paste0("_",brains[ii]))
p2 = vis_gene(spe = speb,sampleid = samples[2],geneid = gene, spatial = TRUE,assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2,... = paste0("_",brains[ii]))
p3 = vis_gene(spe = speb,sampleid = samples[3],geneid = gene, spatial = TRUE,assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2,... = paste0("_",brains[ii]))
p4 = vis_gene(spe = speb,sampleid = samples[4],geneid = gene, spatial = TRUE,assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2,... = paste0("_",brains[ii]))

grid.arrange(p1,p2,p3,p4,nrow = 2)

##
ii = 2
speb = spe[,which(spe$brnum == brains[ii])]
samples = unique(speb$sample_id)
samples

p1 = vis_gene(spe = speb,sampleid = samples[1],geneid = gene, spatial = TRUE,assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2,... = paste0("_",brains[ii]))
p2 = vis_gene(spe = speb,sampleid = samples[2],geneid = gene, spatial = TRUE,assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2,... = paste0("_",brains[ii]))
p3 = vis_gene(spe = speb,sampleid = samples[1],geneid = gene, spatial = TRUE,assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2,... = paste0("_",brains[ii]))
p4 = vis_gene(spe = speb,sampleid = samples[2],geneid = gene, spatial = TRUE,assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2,... = paste0("_",brains[ii]))

grid.arrange(p1,p3,p2,p4,nrow = 2) 

##
ii = 3
speb = spe[,which(spe$brnum == brains[ii])]
samples = unique(speb$sample_id)
samples

p1 = vis_gene(spe = speb,sampleid = samples[3],geneid = gene, spatial = TRUE,assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2,... = paste0("_",brains[ii]))
p2 = vis_gene(spe = speb,sampleid = samples[1],geneid = gene, spatial = TRUE,assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2,... = paste0("_",brains[ii]))
p3 = vis_gene(spe = speb,sampleid = samples[3],geneid = gene, spatial = TRUE,assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2,... = paste0("_",brains[ii]))
p4 = vis_gene(spe = speb,sampleid = samples[1],geneid = gene, spatial = TRUE,assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2,... = paste0("_",brains[ii]))

grid.arrange(p1,p2,p3,p4,nrow = 2)

p1 = vis_gene(spe = speb,sampleid = samples[2],geneid = gene, spatial = TRUE,assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2,... = paste0("_",brains[ii]))
p2 = vis_gene(spe = speb,sampleid = samples[4],geneid = gene, spatial = TRUE,assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2,... = paste0("_",brains[ii]))
p3 = vis_gene(spe = speb,sampleid = samples[2],geneid = gene, spatial = TRUE,assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2,... = paste0("_",brains[ii]))
p4 = vis_gene(spe = speb,sampleid = samples[4],geneid = gene, spatial = TRUE,assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2,... = paste0("_",brains[ii]))

grid.arrange(p1,p3,p2,p4,nrow = 2)

##
ii = 4
speb = spe[,which(spe$brnum == brains[ii])]
samples = unique(speb$sample_id)
samples

p1 = vis_gene(spe = speb,sampleid = samples[1],geneid = gene, spatial = TRUE,assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2,... = paste0("_",brains[ii]))
p2 = vis_gene(spe = speb,sampleid = samples[2],geneid = gene, spatial = TRUE,assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2,... = paste0("_",brains[ii]))
p3 = vis_gene(spe = speb,sampleid = samples[5],geneid = gene, spatial = TRUE,assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2,... = paste0("_",brains[ii]))
p4 = vis_gene(spe = speb,sampleid = samples[3],geneid = gene, spatial = TRUE,assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2,... = paste0("_",brains[ii]))
p5 = vis_gene(spe = speb,sampleid = samples[4],geneid = gene, spatial = TRUE,assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2,... = paste0("_",brains[ii]))
grid.arrange(p1,p2,p3,p4,p5,nrow = 2)

##
ii = 5
speb = spe[,which(spe$brnum == brains[ii])]
samples = unique(speb$sample_id)
samples

p1 = vis_gene(spe = speb,sampleid = samples[1],geneid = gene, spatial = TRUE,assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2,... = paste0("_",brains[ii]))
p2 = vis_gene(spe = speb,sampleid = samples[2],geneid = gene, spatial = TRUE,assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2,... = paste0("_",brains[ii]))
p3 = vis_gene(spe = speb,sampleid = samples[3],geneid = gene, spatial = TRUE,assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2,... = paste0("_",brains[ii]))
p4 = vis_gene(spe = speb,sampleid = samples[4],geneid = gene, spatial = TRUE,assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2,... = paste0("_",brains[ii]))

grid.arrange(p1,p2,p3,p4,nrow = 2)

##
ii = 6
speb = spe[,which(spe$brnum == brains[ii])]
samples = unique(speb$sample_id)
samples

p1 = vis_gene(spe = speb,sampleid = samples[1],geneid = gene, spatial = TRUE,assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2,... = paste0("_",brains[ii]))
p2 = vis_gene(spe = speb,sampleid = samples[2],geneid = gene, spatial = TRUE,assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2,... = paste0("_",brains[ii]))
p3 = vis_gene(spe = speb,sampleid = samples[3],geneid = gene, spatial = TRUE,assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2,... = paste0("_",brains[ii]))

grid.arrange(p1,p2,p3,nrow = 2)

##
ii = 7
speb = spe[,which(spe$brnum == brains[ii])]
samples = unique(speb$sample_id)
samples

p1 = vis_gene(spe = speb,sampleid = samples[2],geneid = gene, spatial = TRUE,assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2,... = paste0("_",brains[ii]))
p2 = vis_gene(spe = speb,sampleid = samples[1],geneid = gene, spatial = TRUE,assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2,... = paste0("_",brains[ii]))
p3 = vis_gene(spe = speb,sampleid = samples[4],geneid = gene, spatial = TRUE,assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2,... = paste0("_",brains[ii]))
p4 = vis_gene(spe = speb,sampleid = samples[3],geneid = gene, spatial = TRUE,assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2,... = paste0("_",brains[ii]))

grid.arrange(p1,p2,p3,p4,nrow = 2)

##
ii = 8
speb = spe[,which(spe$brnum == brains[ii])]
samples = unique(speb$sample_id)
samples

p1 = vis_gene(spe = speb,sampleid = samples[1],geneid = gene, spatial = TRUE,assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2,... = paste0("_",brains[ii]))
p2 = vis_gene(spe = speb,sampleid = samples[2],geneid = gene, spatial = TRUE,assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2,... = paste0("_",brains[ii]))
p3 = vis_gene(spe = speb,sampleid = samples[1],geneid = gene, spatial = TRUE,assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2,... = paste0("_",brains[ii]))
p4 = vis_gene(spe = speb,sampleid = samples[2],geneid = gene, spatial = TRUE,assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2,... = paste0("_",brains[ii]))

grid.arrange(p1,p2,p3,p4,nrow = 2)

##
ii = 9
speb = spe[,which(spe$brnum == brains[ii])]
samples = unique(speb$sample_id)
samples

p1 = vis_gene(spe = speb,sampleid = samples[1],geneid = gene, spatial = TRUE,assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2,... = paste0("_",brains[ii]))
p2 = vis_gene(spe = speb,sampleid = samples[2],geneid = gene, spatial = TRUE,assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2,... = paste0("_",brains[ii]))
p3 = vis_gene(spe = speb,sampleid = samples[3],geneid = gene, spatial = TRUE,assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2,... = paste0("_",brains[ii]))
p4 = vis_gene(spe = speb,sampleid = samples[4],geneid = gene, spatial = TRUE,assayname = "logcounts", minCount = 0, alpha = 0.7, viridis = TRUE, point_size = 2,... = paste0("_",brains[ii]))

grid.arrange(p1,p2,p3,p4,nrow = 2)

dev.off()
