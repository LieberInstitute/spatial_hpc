
layer_stat_cor <-
    function(stats1, stats2,
             reverse = FALSE,
             top_n = NULL) {

        start_time <- Sys.time()

        ## Extract statistics from first dataset
        tstats1 <- stats1$enrichment
        tstats1 <- tstats1[ , grepl("^t_stat_", colnames(tstats1))]
        colnames(tstats1) <- gsub("^t_stat_", "", colnames(tstats1))
        rownames(tstats1) <- stats1$anova$ensembl

        if (reverse) {
            tstats1 <- tstats1 * -1
            colnames(tstats1) <- vapply(strsplit(colnames(tstats1), "-"), function(x) {
                paste(rev(x), collapse = "-")
            }, character(1))
        }

        cat("Time after processing first dataset: ", Sys.time() - start_time, "\n")

        ## Extract statistics from second dataset
        start_time <- Sys.time()

        tstats2 <- stats2$enrichment
        tstats2 <- tstats2[ , grepl("^t_stat_", colnames(tstats2))]
        colnames(tstats2) <- gsub("^t_stat_", "", colnames(tstats2))
        rownames(tstats2) <- stats2$anova$ensembl

        if (reverse) {
            tstats2 <- tstats2 * -1
            colnames(tstats2) <- vapply(strsplit(colnames(tstats2), "-"), function(x) {
                paste(rev(x), collapse = "-")
            }, character(1))
        }

        cat("Time after processing second dataset: ", Sys.time() - start_time, "\n")

        ## Subset to the 'top_n' marker genes if specified
        if (!is.null(top_n)) {
            start_time <- Sys.time()

            stopifnot(top_n < nrow(tstats1) & top_n > 0)
            top_n_index <- unique(as.vector(apply(tstats1, 2, function(t) {
                order(t, decreasing = TRUE)[seq_len(top_n)]
            })))
            tstats1 <- tstats1[top_n_index, , drop = FALSE]
            tstats2 <- tstats2[top_n_index, , drop = FALSE]

            cat("Time after subsetting to top n marker genes: ", Sys.time() - start_time, "\n")
        }

        start_time <- Sys.time()
        common_genes <- intersect(rownames(tstats1), rownames(tstats2))
        if (length(common_genes) == 0) {
            stop("It looks like 'stats1' and 'stats2' do not share ENSEMBL gene names on the rownames or none of the genes are matching.", call. = FALSE)
        }

        tstats1 <- tstats1[common_genes, , drop = FALSE]
        tstats2 <- tstats2[common_genes, , drop = FALSE]

        cat("Time after selecting common genes: ", Sys.time() - start_time, "\n")

        ## Compute correlation
        start_time <- Sys.time()
        cor_res <- cor(tstats1, tstats2)

        cat("Time after computing correlation: ", Sys.time() - start_time, "\n")

        return(cor_res)
    }

# Load necessary packages
library(dplyr)
library(tidyr)




order<-c("DG.GC.PROX1.PDLIM5","DG.GC.PROX1.SGCZ","DG.MC.ARHGAP24.DLC1","CA3.CFAP299.SYN3",
  "CA2.CFAP299.HGF","CA1.dorsal.GRIK1.GRM3","CA1.ventral.ACVR1C.SYT13","SUB.proximal.ROBO1.COL5A2",
  "SUB.distal.FN1.NTNG1","EC.L5.BCL11B.ADRA1A","SUB.proximal.ROBO1.SEMA3E","EC.L56.TLE4.NXPH2","EC.L6b.TLE4.CCN2","EC.L6.TLE4.SULF1",
  "EC.L6.THEMIS.CDH13","EC.L6.THEMIS.RGS12","EC.L5.RORB.TLL1","EC.L5.RORB.TPBG","EC.L2.CUX2.CALB1","EC.L2.CUX2.PDGFD","EC.L3.PCP4.ADARB2",
  "EC.L2.CUX2.IL1RAPL2","EC.L2.RELN.BMPR1B","EC.L2.RELN.BCL11B","EC.L2.CUX2.LAMA3","CR.RELN.NDNF",
  "InN.LHX6.AC008415.1","InN.PVALB.MEPE","InN.PVALB.PLCL1","InN.PVALB.PLEKHH2","InN.SST.ADAMTS12",
  "InN.SST.EPB41L4A","InN.SST.OTOF","InN.SST.NPY","InN.LAMP5.CHST9",
  "InN.LAMP5.KIT","InN.LAMP5.NMBR","InN.NR2F2.DDR2","InN.NR2F2.MIR4300HG","InN.NR2F2.PTPRK",
  "InN.NR2F2.SLC17A8","InN.NR2F2.ANO2","InN.VIP.ABI3BP","InN.VIP.NOX4","InN.VIP.CHRNA2","InN.VIP.SCML4",
  "InN.VIP.SCTR","InN.VIP.PENK","InN.MEIS2.SHISAL2B",
  "Astro.AQP4.GFAP","Astro.AQP4.CHRDL1","Oligo.CPXM2.KANK4","Oligo.OPALIN.LAMA2","Oligo.OPALIN.LINC01098",
  "Oligo.OPALIN.SLC5A11","OPC.PDGFRA.EGR1","OPC.PDGFRA.GRIA4","COP.GPR17.ADAM33","Micro.C1QB.CD83",
  "Micro.C1QB.P2RY12","Macro.F13A1.COLEC12","Myeloid.LSP1.LYZ","T.SKAP1.CD247","aEndo.DKK2.FBLN5",
  "Endo.CLDN5.VWF","PC.CLDN5.ABCC9","vSMC.ABCC9.P2RY14","aSMC.ACTA2.TAGLN","VLMC.COL1A1.COL1A2")

pdf(file=here::here('plots','figures','figure_2','heatmap_sestan.pdf'),w=12,h=8)
pheatmap(x2,clustering_distance_cols = 'correlation',clustering_distance_rows = 'correlation',treeheight_row = 0,treeheight_col = 0,cluster_row=F,cluster_cols=F,color = viridis(20))
dev.off()




x<-layer_stat_cor(sestan,stats)
x<-x[order,levels(sce$fine.type)]
pdf(file=here::here('plots','figures','figure_3','heatmap_sestan_presentation.pdf'),w=18,h=12)
pheatmap(x,cluster_cols=F,cluster_rows=F,color=magma(25),legend=T,
         annotation_col=as.data.frame(annotation),
         annotation_row=annotation2,
         annotation_colors = palettes,annotation_legend=F,
         show_rownames = T)
dev.off()

pdf(file=here::here('plots','figures','figure_3','heatmap_sestan_legend.pdf'),w=16,h=12)
pheatmap(heat_cor2,cluster_cols=F,cluster_rows=F,color=magma(25),legend=T)
         #annotation_col=as.data.frame(annotation),
         #annotation_row=annotation2,
         #annotation_colors = palettes,annotation_legend=F,
         #show_rownames = T)
dev.off()


pdf(file=here::here('plots','figures','figure_2','heatmap_sestan.pdf'),h=7.5,w=8.5)
pheatmap(x,clustering_distance_cols = 'correlation',clustering_distance_rows = 'correlation',color=viridis(15))
dev.off()
