setwd("/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/")
suppressPackageStartupMessages({
  library(here)
  library(SpatialExperiment)
  library(tidyverse)
  library(sessioninfo)
  library(spatialLIBD)
})

#spe = readRDS('/dcs04/lieber/lcolladotor/spatialHPC_LIBD4035/spatial_hpc/processed-data/spot_deconvo/shared_utilities/spe.rds')
#spe$neuropil_pos <- ifelse(spe$Pmask_dark_blue < 0.005, TRUE, FALSE)

#############SLM.SGZ-SL.SR#############		
#spe_sub <- spe[, colData(spe)$cluster_collapsed %in% c("SLM.WM", "SL.SR") &  colData(spe)$neuropil_pos == TRUE]
#spe_sub$cluster_collapsed <- droplevels(spe_sub$cluster_collapsed)
#neuropil_pseudo <-
#  registration_pseudobulk(
#    spe_sub,
#    var_registration = "cluster_collapsed",
#    var_sample_id = "sample_id",
#    covars = c("age", "sex"),
#    min_ncells = 50,
#    pseudobulk_rds_file = here("processed-data", "revision_maddy", "neuropilDE", "SLM.SGZ-SL.SR_neuropil_pseudo.rds")
#    )
#	
#dim(neuropil_pseudo)
	#[1] 3077   16	

neuropil_pseudo = readRDS(here("processed-data", "revision_maddy", "neuropilDE", "SLM.SGZ-SL.SR_neuropil_pseudo.rds"))
dx_res <- registration_stats_enrichment(
	  neuropil_pseudo,
	  block_cor = "sample_id",
	  covars = c("age", "sex"),
	  var_registration = "cluster_collapsed",
	  gene_ensembl = "gene_id",
	  gene_name = "gene_name"
	)
	

	dx_res$sig <- with(dx_res, 
	                   (logFC_SLM.WM > log2(3) & fdr_SLM.WM < 0.05) | 
	                   (logFC_SLM.WM < -log2(3) & fdr_SLM.WM < 0.05)
	)
 
	top_genes_up <- dx_res %>% 
	  filter(sig == TRUE & logFC_SLM.WM > log2(3)) %>% 
	  arrange(desc(abs(logFC_SLM.WM))) %>% 
	  slice(1:5)
  	top_genes_down <- dx_res %>% 
  	  filter(sig == TRUE & logFC_SLM.WM < -log2(3)) %>% 
  	  arrange(desc(abs(logFC_SLM.WM))) %>% 
  	  slice(1:5)
	top_genes = rbind(top_genes_up,top_genes_down)
			
	p <- ggplot(data = dx_res, aes(x = logFC_SLM.WM, y = -log10(p_value_SLM.WM), color = sig)) +
	  geom_point() + 
	  scale_color_manual(values = c("FALSE" = "grey", "TRUE" = "red"), breaks = c("TRUE", "FALSE"), name = "FDR<0.05") + 
	  theme_minimal() +
	  geom_text(data = top_genes, aes(label = gene), color = "black", vjust = -0.5, hjust = 0.5) +
	  #geom_hline(yintercept = -log10(0.03), lty = "dashed", color = "royalblue") + 
	  #geom_vline(xintercept = -log2(2), lty = "dashed", color = "royalblue") + 
	  #geom_vline(xintercept = log2(2), lty = "dashed", color = "royalblue") + 
	  labs(
	    x = expression(log[2]~"fold change"), 
	    y = expression(-log[10]~FDR), 
	    title = "Differential expression in neuropil spots of SLM.SGZ vs. SL.SR"
	  )
	
#p <- ggplot(data = dx_res1, aes(x = logFC_SLM.WM, y = -log10(p_value_SLM.WM), color = sig_03)) +
#     geom_point() + 
#     scale_color_manual(values = c("FALSE" = "grey", "TRUE" = "red"), breaks = c("TRUE", "FALSE"), name = "FDR<0.03") + 
#     theme_minimal() +
#     geom_text(data = top_genes, aes(label = gene), vjust = -0.5, hjust = 0.5)
	  
ggsave(here("plots", "revision_maddy", "neuropilDE", "SLM.SGZ-SL.SR.png"), plot = p, width = 8, height = 6, dpi = 300)
				
 unique(spe$cluster_collapsed)

 [1] CA1      WM.2     WM.1     Choroid  SL.SR    SLM.WM   WM.3     SR.SLM  
 [9] SUB.RHP  GABA     Vascular RHP      SUB      CA2-4    GCL      ML      
16 Levels: GCL CA2-4 CA1 SUB SUB.RHP RHP GABA SL.SR ML SR.SLM SLM.WM ... Choroid

#############SLM.SGZ-E#############			

#spe_sub <- spe[, colData(spe)$cluster_collapsed %in% c("SLM.WM", "SL.SR", "SR.SLM", "ML") &  colData(spe)$neuropil_pos == TRUE]
#spe_sub$neuropil_cluster = 2
#spe_sub$neuropil_cluster[spe_sub$cluster_collapsed == "SLM.WM"] = 1
#
#neuropil_pseudo <-
#  registration_pseudobulk(
#    spe_sub,
#    var_registration = "neuropil_cluster",
#    var_sample_id = "sample_id",
#    covars = c("age", "sex"),
#    min_ncells = 50,
#    pseudobulk_rds_file = here("processed-data", "revision_maddy", "neuropilDE", "SLM.SGZ-E_neuropil_pseudo.rds")
#    )
#	
#dim(neuropil_pseudo)
#	#[1] 15719    62	

neuropil_pseudo = readRDS(here("processed-data", "revision_maddy", "neuropilDE", "SLM.SGZ-E_neuropil_pseudo.rds"))
dx_res <- registration_stats_enrichment(
  neuropil_pseudo,
  block_cor = "sample_id",
  covars = c("age", "sex"),
  var_registration = "neuropil_cluster",
  gene_ensembl = "gene_id",
  gene_name = "gene_name"
)
	

dx_res$sig <- with(dx_res, 
                   (logFC_1> log2(3) & fdr_1 < 0.05) | 
                   (logFC_1< -log2(3) & fdr_1 < 0.05)
) 
top_genes_up <- dx_res %>% 
  filter(sig == TRUE & logFC_1 > log2(3)) %>% 
  arrange(desc(abs(logFC_1))) %>% 
  slice(1:5)
top_genes_down <- dx_res %>% 
  filter(sig == TRUE & logFC_1 < -log2(3)) %>% 
  arrange(desc(abs(logFC_1))) %>% 
  slice(1:5)
top_genes = rbind(top_genes_up,top_genes_down)
			
p <- ggplot(data = dx_res, aes(x = logFC_1, y = -log10(p_value_1), color = sig)) +
  geom_point() + 
  scale_color_manual(values = c("FALSE" = "grey", "TRUE" = "red"), breaks = c("TRUE", "FALSE"), name = "FDR<0.05") + 
  theme_minimal() +
  geom_text(data = top_genes, aes(label = gene), color = "black", vjust = -0.5, hjust = 0.5) +
  #geom_hline(yintercept = -log10(0.03), lty = "dashed", color = "royalblue") + 
  #geom_vline(xintercept = -log2(2), lty = "dashed", color = "royalblue") + 
  #geom_vline(xintercept = log2(2), lty = "dashed", color = "royalblue") + 
  labs(
    x = expression(log[2]~"fold change"), 
    y = expression(-log[10]~FDR), 
    title = "Differential expression in neuropil spots of SLM.SGZ vs. ML/SL.SR/SR.SLM"
  )

ggsave(here("plots", "revision_maddy", "neuropilDE", "SLM.SGZ-E.png"), plot = p, width = 8, height = 6, dpi = 300)
				
#############SLM.SGZ-MLM#############				
				
#spe_sub <- spe[, colData(spe)$cluster_collapsed %in% c("SLM.WM", "ML") &  colData(spe)$neuropil_pos == TRUE]
#spe_sub$cluster_collapsed <- droplevels(spe_sub$cluster_collapsed)
#
#neuropil_pseudo <-
#  registration_pseudobulk(
#    spe_sub,
#    var_registration = "cluster_collapsed",
#    var_sample_id = "sample_id",
#    covars = c("age", "sex"),
#    min_ncells = 50,
#    pseudobulk_rds_file = here("processed-data", "revision_maddy", "neuropilDE", "SLM.SGZ-ML_neuropil_pseudo.rds")
#    )
#
#dim(neuropil_pseudo)
#	#[1] 15719    62	
neuropil_pseudo = readRDS(here("processed-data", "revision_maddy", "neuropilDE", "SLM.SGZ-ML_neuropil_pseudo.rds"))

dx_res <- registration_stats_enrichment(
  neuropil_pseudo,
  block_cor = "sample_id",
  covars = c("age", "sex"),
  var_registration = "cluster_collapsed",
  gene_ensembl = "gene_id",
  gene_name = "gene_name"
)

dx_res$sig <- with(dx_res, 
                   (logFC_SLM.WM > log2(3) & fdr_SLM.WM < 0.05) | 
                   (logFC_SLM.WM < -log2(3) & fdr_SLM.WM < 0.05)
) 
top_genes_up <- dx_res %>% 
  filter(sig == TRUE & logFC_SLM.WM > log2(3)) %>% 
  arrange(desc(abs(logFC_SLM.WM))) %>% 
  slice(1:5)
top_genes_down <- dx_res %>% 
  filter(sig == TRUE & logFC_SLM.WM < -log2(3)) %>% 
  arrange(desc(abs(logFC_SLM.WM))) %>% 
  slice(1:5)
top_genes = rbind(top_genes_up,top_genes_down)

			
p <- ggplot(data = dx_res, aes(x = logFC_SLM.WM, y = -log10(p_value_SLM.WM), color = sig)) +
  geom_point() + 
  scale_color_manual(values = c("FALSE" = "grey", "TRUE" = "red"), breaks = c("TRUE", "FALSE"), name = "FDR<0.05") + 
  theme_minimal() +
  geom_text(data = top_genes, aes(label = gene), color = "black", vjust = -0.5, hjust = 0.5) +
  #geom_hline(yintercept = -log10(0.03), lty = "dashed", color = "royalblue") + 
  #geom_vline(xintercept = -log2(2), lty = "dashed", color = "royalblue") + 
  #geom_vline(xintercept = log2(2), lty = "dashed", color = "royalblue") + 
  labs(
    x = expression(log[2]~"fold change"), 
    y = expression(-log[10]~FDR), 
    title = "Differential expression in neuropil spots of SLM.SGZ vs. ML"
  )


ggsave(here("plots", "revision_maddy", "neuropilDE", "SLM.SGZ-ML.png"), plot = p, width = 8, height = 6, dpi = 300)


#############SLM.SGZ-SR.SLM#############

#spe_sub <- spe[, colData(spe)$cluster_collapsed %in% c("SLM.WM", "SR.SLM") &  colData(spe)$neuropil_pos == TRUE]
#spe_sub$cluster_collapsed <- droplevels(spe_sub$cluster_collapsed)
#neuropil_pseudo <-
#  registration_pseudobulk(
#    spe_sub,
#    var_registration = "cluster_collapsed",
#    var_sample_id = "sample_id",
#    covars = c("age", "sex"),
#    min_ncells = 50,
#    pseudobulk_rds_file = here("processed-data", "revision_maddy", "neuropilDE", "SLM.SGZ-SR.SLM_neuropil_pseudo.rds")
#    )
#
#dim(neuropil_pseudo)
	#[1] 15719    62	
neuropil_pseudo = readRDS(here("processed-data", "revision_maddy", "neuropilDE", "SLM.SGZ-SR.SLM_neuropil_pseudo.rds"))
dx_res <- registration_stats_enrichment(
  neuropil_pseudo,
  block_cor = "sample_id",
  covars = c("age", "sex"),
  var_registration = "cluster_collapsed",
  gene_ensembl = "gene_id",
  gene_name = "gene_name"
)

dx_res$sig <- with(dx_res, 
                   (logFC_SLM.WM > log2(3) & fdr_SLM.WM < 0.05) | 
                   (logFC_SLM.WM < -log2(3) & fdr_SLM.WM < 0.05)
)
 
top_genes_up <- dx_res %>% 
  filter(sig == TRUE & logFC_SLM.WM > log2(3)) %>% 
  arrange(desc(abs(logFC_SLM.WM))) %>% 
  slice(1:5)
top_genes_down <- dx_res %>% 
  filter(sig == TRUE & logFC_SLM.WM < -log2(3)) %>% 
  arrange(desc(abs(logFC_SLM.WM))) %>% 
  slice(1:5)
top_genes = rbind(top_genes_up,top_genes_down)

			
p <- ggplot(data = dx_res, aes(x = logFC_SLM.WM, y = -log10(p_value_SLM.WM), color = sig)) +
  geom_point() + 
  scale_color_manual(values = c("FALSE" = "grey", "TRUE" = "red"), breaks = c("TRUE", "FALSE"), name = "FDR<0.05") + 
  theme_minimal() +
  geom_text(data = top_genes, aes(label = gene), color = "black", vjust = -0.5, hjust = 0.5) +
  #geom_hline(yintercept = -log10(0.03), lty = "dashed", color = "royalblue") + 
  #geom_vline(xintercept = -log2(2), lty = "dashed", color = "royalblue") + 
  #geom_vline(xintercept = log2(2), lty = "dashed", color = "royalblue") + 
  labs(
    x = expression(log[2]~"fold change"), 
    y = expression(-log[10]~FDR), 
    title = "Differential expression in neuropil spots of SLM.SGZ vs. SR.SLM"
  )


ggsave(here("plots", "revision_maddy", "neuropilDE", "SLM.SGZ-SR.SLM.png"), plot = p, width = 8, height = 6, dpi = 300)


