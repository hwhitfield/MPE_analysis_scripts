# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
#   PERFORM AIMS SUBTYPING FOR FIGURE 3 DATA -- Malignant & Mesothelial cells
#   
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


##### ----- Get directories
setwd("./MPE_paper_analysis/")
data_dir <- paste0(getwd(),"/data/")
matrix_dir <- paste0(data_dir,"raw_matrices/")
fig_dir <- paste0(getwd(),"/figures/")

source(paste0(getwd(),"/scripts/helper_functions.R"))




##### ----- Load & process data ----------------------------------------------------------

ALL_COLDATA  <- read.table(file=paste0(data_dir,"ALL_COLDATA.tsv"),sep="\t")


### --- Load count matrices into SingleCellExperiment objects
BCB_names <- c("BCB112","BCB139",
               "BCB66","BCB20","BCB21",
               "BCB66_E","BCB20_E","BCB21_E")
Non_immune_sce_lst <- lapply(BCB_names, 
                             function(bcb){
                               loadSCE(paste0(bcb,"_"), 
                                       paste0(matrix_dir,bcb),
                                       ALL_COLDATA)})
names(Non_immune_sce_lst) <- BCB_names


### --- Subset to non-immune
Non_immune_cells <- ALL_COLDATA[ALL_COLDATA$CellType_Fig3 %in% c("Malignant","Mesothelial"),]$CellID
Non_immune_sce_lst <- lapply(Non_immune_sce_lst, function(x){x[,x$CellID %in% Non_immune_cells]})


### --- Scran normalise
Non_immune_sce_lst <- lapply(Non_immune_sce_lst, function(x) {scranNorm(x)})



### --- CLUSTERING FOR PSEUDOBULK -------------------------------

Cluster_KNN <- function(sce_object, K){
  g <- buildSNNGraph(sce_object[as.vector(rowSums(logcounts(sce_object)>0)>5),], 
                     k=K,  assay.type = "logcounts")
  clust <- igraph::cluster_walktrap(g)$membership
  sce_object[[paste0("Clust_Knn",K)]] <- as.vector(clust)
  return(sce_object)
}

Non_immune_sce_lst <- lapply(Non_immune_sce_lst, function(sce_x){Cluster_KNN(sce_x,K=30)})
Non_immune_sce_lst <- lapply(Non_immune_sce_lst, function(sce_x){scater::runTSNE(sce_x)})


PatientClusters <- unlist(lapply(Non_immune_sce_lst, function(x){x$Clust_Knn30}))
names(PatientClusters) <- unlist(lapply(Non_immune_sce_lst, function(x){x$CellID}))

ALL_COLDATA$ClusterKNN_patient <- as.vector(PatientClusters[ALL_COLDATA$CellID])
write.table(ALL_COLDATA , file=paste0(data_dir,"ALL_COLDATA.tsv"),sep="\t")

#save(Non_immune_sce_lst, file=paste0(data_dir,"Non_immune_sce_lst.Rdata"))



### --- GET GENE IDS FOR GENEFU -------------------------------

GetGeneIDs <- function(sce_object, ensdb){
  require(BiocFileCache)
  require(AnnotationHub)
  ## Change name for genefu
  rowData(sce_object)$EntrezGene.ID <- mapIds(ensdb,
                                              keys=rownames(sce_object),
                                              keytype="SYMBOL", column="ENTREZID")
  
  dge_x <- dge_x[!is.na(dge_x$genes$EntrezGene.ID),]
  dge_x <- dge_x[!(duplicated(dge_x$genes$EntrezGene.ID)),]
  rownames(dge_x) <- dge_x$genes$EntrezGene.ID
  
  return(sce_object)
}

ensdb <- AnnotationHub()[["AH73881"]] #GRCh38, ensembl97 hsa 
Non_immune_sce_lst <- lapply(Non_immune_sce_lst, function(sce_x){GetGeneIDs(sce_x,ensdb)})




### --- PSEUDOBULK KNN CLUSTERS WITHIN EACH PATIENTS -------------------------------

Pseudobulk_sce <- function(sce_object, cluster_labels,cpm_threshold=0.1, sample_threshold=0.1){
  require(edgeR)
  Summed_Sce <- scater::aggregateAcrossCells(sce_object, id=cluster_labels)
  DGE_Summed_sce <- edgeR::DGEList(counts(Summed_Sce), 
                                   genes = rowData(sce_object))
  
  ## Filter samples with low library size
  DGE_Summed_sce <- DGE_Summed_sce[,DGE_Summed_sce$samples$lib.size > 1e4]
  
  return(DGE_Summed_sce)
}

DGE_lst <- lapply(Non_immune_sce_lst, function(x){Pseudobulk_sce(x, x$Clust_Knn30)})




### --- CLASSIFICATION WITH GENEFU -------------------------------

library(genefu)


Pseudobulk_classify <- function(dge_x){
  require(genefu)
  
  # Don't need CPM because classification is independent between samples (single sample)
  # Should use a variance stabilising transformation like log
  annot_df <- dge_x$genes
  log_data <- t(log2(dge_x$counts +1))
  
  ### AIMS model -- might be inappropriate because many zeros
  AIMS_subtypes <- genefu::molecular.subtyping(sbt.model="AIMS", data=log_data, 
                                               annot=annot_df, verbose=TRUE, do.mapping=TRUE)
  
  AMIS_subtypes_dict <- as.vector(AIMS_subtypes$subtype)
  names(AMIS_subtypes_dict) <- rownames(AIMS_subtypes$subtype)
  dge_x$samples$AIMS <- as.vector(AMIS_subtypes_dict[rownames(dge_x$samples)])

  return(dge_x)
}


DGE_lst <- lapply(DGE_lst, function(x){Pseudobulk_classify(x)})



### --- Append subtype
cellID_to_subtype <- c()
for (i_name in names(Non_immune_sce_lst)){
  sce_x <- Non_immune_sce_lst[[i_name]]
  dge_x <- DGE_lst[[i_name]]
  
  subtype_dict <- setNames(dge_x$samples$AIMS,as.character(rownames(dge_x$samples)))
  
  sce_x$AIMS_subtype <- as.vector(subtype_dict[as.character(sce_x$Clust_Knn30)])
  cellID_to_subtype <- c(cellID_to_subtype,setNames(sce_x$AIMS_subtype, sce_x$CellID))
  Non_immune_sce_lst[[i_name]] <- sce_x
}

cellID_to_subtype[setdiff(ALL_COLDATA$CellID,names(cellID_to_subtype))] <- NA

ALL_COLDATA$AIMS_subtype <- as.vector(cellID_to_subtype[ALL_COLDATA$CellID])
write.table(ALL_COLDATA , file=paste0(data_dir,"ALL_COLDATA.tsv"),sep="\t")





##### ----- SUPP FIGURE 8 TSNE PER PATIENT -----------------------------------


library(ggRNA)
Colour_values_AIMS <- c("#CABEE9", "#7C7189", "#FAE093", "#D04E59", "#2F3D70")
names(Colour_values_AIMS) <- c( "LumA", "LumB", "Her2","Normal", "Basal")

plot_lst <- lapply(Non_immune_sce_lst, function(x){
  ggDIMRED(x, dimred="TSNE", colour_by = "AIMS_subtype",col_pal = Colour_values_AIMS,point_size=1,point_alpha = 0.85)+
    theme_blank() + theme(legend.position = "none")
  })
names(plot_lst) <- names(Non_immune_sce_lst)

plot_lst <- plot_lst[c("BCB20","BCB20_E",
                       "BCB21","BCB21_E",
                       "BCB66","BCB66_E",
                       "BCB139","BCB112")]

plot_leg <- ggRNA::GetLegend(Colour_values_AIMS, "Subtype")

pdf(paste0(fig_dir,"SuppFig8_SubtypetSNEs.pdf"), width=6,height=10)
cowplot::plot_grid(cowplot::plot_grid(plotlist = plot_lst,ncol=2, 
                            labels=names(plot_lst), label_fontface ="plain" ),
                   plot_leg,ncol=1, rel_heights = c(6,1))
dev.off()



##### ----- LEHMAN BASAL SUBTYPING -----------------------------------

library(VISION)

### --- Get TNBC samples
load(file=paste0(data_dir,"Non_immune_sce_lst.Rdata"))

TNBC_sce_lst <- Non_immune_sce_lst[names(Non_immune_sce_lst) %in% c("BCB112", "BCB66","BCB66_E")]
rm(Non_immune_sce_lst)
gc()


### --- Load lehmann genes
L_dict <- c("LAR", "BS1", "BS2", "Mes", "Mes-Stem", "Immune")
names(L_dict) <- c("LuminalAR","BasalLike2","BasalLike1","Mesenchymal","MesStem", "Immuno")

LehmannGenes <- read.csv(paste0(data_dir,"LehmannTNBCSubtypes_processed.csv"), 
                         sep=",",header=TRUE)
LehmannGenes$SUBTYPE <- as.vector(L_dict[LehmannGenes$SUBTYPE])
LehmannGenes <- LehmannGenes[LehmannGenes$REGULATION == "UP",]
Gene_lst <- split(LehmannGenes$GENE, LehmannGenes$SUBTYPE)

LehmannGenes <- LehmannGenes[!(duplicated(LehmannGenes$GENE)),]
UniqueGene_lst <- split(LehmannGenes$GENE, LehmannGenes$SUBTYPE)


### --- Convert signatures to VISION objects

MakeVisionObj <- function(gene_vec, sig_str, row_genes){
  require(VISION)
  gene_vec <- gene_vec[!is.na(gene_vec)]
  gene_vec <- unique(gene_vec)
  gene_vec <- gene_vec[gene_vec %in% row_genes]
  
  genes <- rep(1,length(gene_vec))
  names(genes) <- gene_vec
  return(VISION::createGeneSignature(name = sig_str, sigData = genes))
}

data_genes <- unique(unlist(lapply(TNBC_sce_lst, function(x){rownames(x)})))
obj_lst <- lapply(names(Gene_lst), function(x){MakeVisionObj(Gene_lst[[x]], x, data_genes)})

obj_lst_2 <- lapply(names(UniqueGene_lst), function(x){MakeVisionObj(Gene_lst[[x]], x, data_genes)})


### --- Run Vision 

runVISION <- function(sce_x, sig_lst){
  
  ## Use normalised expression data, but not log-transformed
  counts(sce_x) <- as(counts(sce_x), "dgCMatrix")
  norm_counts <- scuttle::normalizeCounts(sce_x, transform="none", exprs_values = "counts")
  
  VISION_Obj <- Vision(data =  norm_counts,pool=FALSE,
                       signatures = sig_lst)
  options(mc.cores = 1)
  VISION_RESULTS <- analyze(VISION_Obj)
  
  ## Append signature scores
  for (iSig in colnames(getSignatureScores(VISION_RESULTS))){
    sce_x[[iSig]] <- as.vector(getSignatureScores(VISION_RESULTS)[, iSig])
  }
  
  return(sce_x)
}

TNBC_sce_lst <- lapply(TNBC_sce_lst, function(x){runVISION(x, obj_lst)})

#save(TNBC_sce_lst, file=paste0(data_dir,"TNBC_sce_lst.Rdata"))


### --- Select subtype with max score
maxScore <- function(sce_x){
  coldata <- as.data.frame(colData(sce_x)[,colnames(colData(sce_x)) 
                                          %in% c("LAR", "BS1", "BS2", "Mes", "Mes-Stem", "Immune")])
  sce_x$LehmannSubtype <- colnames(coldata)[apply(coldata, 1, which.max)]
  return(sce_x)
}

TNBC_sce_lst <- lapply(TNBC_sce_lst, function(x){maxScore(x)})




##### ----- SUPP FIGURE 9 --------------------------------------------------


library(Manu)
Colours_Lehmann <- get_pal("Hoiho")
names(Colours_Lehmann) <-  c("LAR", "BS1", "BS2", "Mes", "Mes-Stem", "Immune")

plot_lst <- lapply(TNBC_sce_lst, function(x){
  ggDIMRED(x, dimred="TSNE", colour_by = "LehmannSubtype",
           col_pal = Colours_Lehmann,point_size=1,point_alpha = 0.85)+
              theme_blank() + theme(legend.position = "none")
})
names(plot_lst) <- names(TNBC_sce_lst)

plot_lst <- plot_lst[c("BCB112", "BCB66","BCB66_E","BCB139")]


plot_leg <- ggRNA::GetLegend(Colours_Lehmann, "Subtype")

pdf(paste0(fig_dir,"SUPP9_PatientTSNEs_byLehmann.pdf"), width=8,height=3)
cowplot::plot_grid(cowplot::plot_grid(plotlist = plot_lst[c("BCB112", "BCB66","BCB66_E")],
                                  nrow=1, labels=names(plot_lst), label_fontface ="plain" ),
                   plot_leg,ncol=1, rel_heights = c(6,1))
dev.off()




