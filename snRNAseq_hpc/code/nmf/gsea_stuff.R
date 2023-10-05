
sign(de.dge[[i]]$logFC) * sqrt(de.dge[[i]]$F)

pats<-list()
for(i in 1:100){
    pats[[i]]<-x@w[,i]
   # names(pats[[i]])<-rownames(x@w)
    pats[[i]]<-pats[[i]][order(pats[[i]],decreasing=T)]
}
names(pats)<-colnames(x@w)

pats<-list()
for(i in 1:100){
    pats[[i]]<-z[,i]
    # names(pats[[i]])<-rownames(pa
    pats[[i]]<-pats[[i]][order(pats[[i]],decreasing=T)]
    pats[[i]]<-pats[[i]][!is.na(pats[[i]])]
}
names(pats)<-colnames(x@w)

gseaRes<-list()
for(i in 1:length(celltype)){
    print(Sys.time())
    print(paste0('running GSEA for pattern ',i))
    gseaRes[[i]]<-fgsea(
        trimmedList,
        celltype[[i]],
        minSize = 10,
        maxSize = 1000,
        gseaParam = 1
    )}

fdr<-list()
for(i in 1:length(gseaRes)){
    fdr[[i]]<-gseaRes[[i]]$padj
    names(fdr[[i]])<-gseaRes[[i]]$pathway
}


nes<-list()
for(i in 1:length(gseaRes)){
    nes[[i]]<-gseaRes[[i]]$NES
    names(nes[[i]])<-gseaRes[[i]]$pathway
}

celltype<-list()
for(i in levels(sce$fine.type)){
    celltype[[i]]<-sign(statsstats$enrichment[[paste0("t_stat_",i)]]
    names(celltype[[i]])<-rownames(stats$enrichment)
    celltype[[i]]<-celltype[[i]][order(celltype[[i]],decreasing=T)]
}

visclust<-list()
for(i in levels(spe$cluster)){
   visclust[[i]]<-spe_stats$enrichment[[paste0("t_stat_",i)]]*sign(spe_stats$enrichment[[paste0("logFC_",i)]])
    names(visclust[[i]])<-rownames(spe_stats$enrichment)
    visclust[[i]]<-visclust[[i]][order(visclust[[i]],decreasing=T)]
}

twas<-list()
for(i in unique(dis$twas)){
    twas[[i]]<-dis$gene_symbol[dis$twas==i]
    twas[[i]]<-unique(twas[[i]])
}
type_scz_twas<-list()
for(i in 1:length(celltype)){
type_scz_twas[[i]]<-fgsea(
    scz_gene,
    celltype[[i]],
    minSize = 10,
    maxSize = 1000,
    gseaParam = 1
)}

vis_scz_gene<-list()
for(i in 1:length(visclust)){
    vis_scz_gene[[i]]<-fgsea(
        scz_gene,
        visclust[[i]],
        minSize = 10,
        maxSize = 1000,
        gseaParam = 1
    )}

type_heat<-list()
for(i in names(type_scz_gene)){
    type_heat[[i]]<-type_scz_gene[[i]]$padj
    names(type_heat[[i]])<-type_scz_gene[[i]]$pathway}
type_heat<-do.call(cbind,type_heat)
    if(length(type_heat[is.na(type_heat)])>0){
        type_heat[is.na(type_heat)]<-1
        type_heat=-log10(type_heat)
    }else{
        type_heat=-log10(type_heat)
    }


vis_heat<-list()
for(i in names(vis_dx)){
    vis_heat[[i]]<-vis_dx[[i]]$padj
    names(vis_heat[[i]])<-vis_dx[[i]]$pathway
}

twas_dx<-list()
for(i in unique(twas$dx)){
    twas_dx[[i]]<-twas$symbol[twas$dx==i]
    twas_dx[[i]]<-unique(twas_dx[[i]])
    twas_dx[[i]]<-twas_dx[[i]][twas_dx[[i]] %in% rownames(sce)]
}

gwas_dx<-list()
for(i in unique(gwas$dx)){
    gwas_dx[[i]]<-gwas$symbol[gwas$dx==i]
    gwas_dx[[i]]<-unique(gwas_dx[[i]])
    gwas_dx[[i]]<-gwas_dx[[i]][gwas_dx[[i]] %in% rownames(sce)]
}

dx<-list()
for(i in unique(dx_dx$dx)){
    dx[[i]]<-dx_dx$symbol[dx_dx$dx==i]
    dx[[i]]<-unique(dx[[i]])
    #dx[[i]]<-dx[[i]][dx[[i]] %in% rownames(sce)]
}

mdd<-twas[twas$twas=='MDD',]
mdd<-mdd[!is.na(mdd$fdr),]

twas_dx$MDD<-mdd$symbol

type_twas<-list()
for(i in 1:length(celltype)){
    type_twas[[i]]<-fgsea(
        twas_dx,
        celltype[[i]],
        minSize = 5,
        maxSize = 200,
        gseaParam = 1
    )}

vis_twas<-list()
for(i in 1:length(visclust)){
    vis_twas[[i]]<-fgsea(
        twas_dx,
        visclust[[i]],
        minSize = 5,
        maxSize = 200,
        gseaParam = 1
    )}



z<-for(i in 1:length(pats)){
all(is.finite(pats[[i]]))}


type_twas<-list()
for(i in 1:length(celltype)){
    type_twas[[i]]<-fgsea(
        twas_dx,
        celltype[[i]],
        minSize = 10,
        maxSize = 2000,
        gseaParam = 1
    )}

vis_twas<-list()
for(i in 1:length(visclust)){
    vis_twas[[i]]<-fgsea(
        twas_dx,
        visclust[[i]],
        minSize = 10,
        maxSize = 2000,
        gseaParam = 1
    )}

vis_heat<-list()
for(i in names(vis_scz_gene)){
    vis_heat[[i]]<-vis_scz_gene[[i]]$padj
    names(vis_heat[[i]])<-vis_scz_gene[[i]]$pathway}
vis_heat<-do.call(cbind,vis_heat)
if(length(vis_heat[is.na(vis_heat)])>0){
    vis_heat[is.na(vis_heat)]<-1
    vis_heat=-log10(vis_heat)
}else{
    vis_heat=-log10(vis_heat)
}

set.seed(1)
gesecaRes <- geseca(genesets, x@w, minSize = 5, maxSize = 5000, eps=1e-100,center=F)


list_of_lists <- list(twas_dx = twas_dx, gwas_dx = gwas_dx, dx = dx)

# Loop through each list in the list of lists
for (list_name in names(list_of_lists)) {
    current_list <- list_of_lists[[list_name]]

    # Filter list based on the names of its elements
    filtered_list <- current_list[names(current_list) %in% c('MDD','ASD','ADHD','AUD','AD','Anorexia','SCZ','GAD','BIP','NDD','OUD')]

    # Replace the original list with the filtered list
    list_of_lists[[list_name]] <- filtered_list

    # Or do something else with filtered_list, e.g., print or store it
    print(filtered_list)
}

# Initialize an empty list to store the consensus data
consensusList <- list()

# Loop through each list in genelists
for (name in names(genelists)) {
    current_list <- genelists[[name]]

    # Loop through each named vector in the current list
    for (vector_name in names(current_list)) {
        current_vector <- current_list[[vector_name]]

        # If this name has not been seen before, initialize it in consensusList
        if (is.null(consensusList[[vector_name]])) {
            consensusList[[vector_name]] <- character(0)
        }

        # Append the current vector to the existing data in consensusList
        consensusList[[vector_name]] <- c(consensusList[[vector_name]], current_vector)
    }
}

for(i in 1:length(consensusList)){
    consensusList[[i]]<-unique(consensusList[[i]])
    consensusList[[i]]<-consensusList[[i]][consensusList[[i]] %in% rownames(sce)]
}

# Count the occurrence of each gene in each list
count_gene_occurrence <- function(gene, list_of_lists) {
    count <- 0
    for (lst in list_of_lists) {
        if (gene %in% lst) {
            count <- count + 1
        }
    }
    return(count)
}

# Initialize an empty list to store the trimmed data
trimmedList <- list()

# Loop through each named vector in the consensusList
for (vector_name in names(consensusList)) {
    current_vector <- consensusList[[vector_name]]

    # Check if the length exceeds 500
    if (length(current_vector) > 500) {
        original_lists <- list(
            genelists$twas_dx[[vector_name]],
            genelists$gwas_dx[[vector_name]],
            genelists$dx[[vector_name]]
        )

        # Filter genes that occur in at least 2/3 lists
        filtered_vector <- current_vector[sapply(current_vector, function(g) {
            count_gene_occurrence(g, original_lists)
        }) >= 2]

        # If still too long, take only those in 3/3 geneList lists
        if (length(filtered_vector) > 500) {
            filtered_vector <- filtered_vector[sapply(filtered_vector, function(g) {
                count_gene_occurrence(g, original_lists)
            }) == 3]
        }

        # Store the filtered or trimmed vector
        trimmedList[[vector_name]] <- filtered_vector
    } else {
        # If the length is already under the limit, just store it as is
        trimmedList[[vector_name]] <- current_vector
    }
}


relevant_tissues <- c(
    "Brain", "Brain Cortex", "Amygdala", "Anterior Cingulate Cortex", "Brain Hippocampus",
    "Cerebellum", "Cortex", "Frontal Cortex", "Nucleus Accumbens Basal Ganglia", "Substantia Nigra",
    "Brain Hypothalamus", "Caudate Basal Ganglia", "Cerebellar Hemisphere", "Putamen Basal Ganglia",
    "Dorsal Lateral Prefrontal Cortex", "Brain Substantia Nigra", "Brain Dorsolateral Prefrontal Cortex",
    "Hypothalamus", "Brain Caudate Basal Ganglia", "Brain Putamen Basal Ganglia",
    "Brain Cerebellar Hemisphere", "Brain Cerebellum", "Brain Frontal Cortex",
    "Brain Nucleus Accumbens", "Brain Prefrontal Cortex", "Brain Hippocampal", "Brain Spinal Cord",
    "Brain Caudate", "Brain Cerebrum Hypothalamus", "Adult Cerebral Cortex", "Hippocampus", "Nucleus Accumbens",
    "Nucleus accumbens basal ganglia", "Brain Anterior cingulate cortex BA24", "Brain Caudate basal ganglia",
    "Brain Frontal Cortex BA9", "Brain Nucleus accumbens basal ganglia", "Brain Putamen basal ganglia",
    "Brain Amygdala", "Brain Spinal cord cervical c-1", "Brain Substantia nigra",
    "Aortic Artery", "Artery Tibia", "Artery Coronary", "Artery Aorta", "Artery Tibial",
    "CD14", "Monocyte", "Lymphoblastoid Cell Lines", "Lymphocyte Cell", "Immune Cell"
)




# Now trimmedList contains vectors with length at most 500, prioritizing genes that occur in multiple original lists
set.seed(1)
aucRes<-AUCell_run(counts(sce),trimmedList)
