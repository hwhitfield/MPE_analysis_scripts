# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
#   MPE PROCESSING & ANALYSIS HELPER FUNCTIONS
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

##### ----- Load libraries
requiredPackages <- c("SingleCellExperiment", "scater","scran","DropletUtils","Matrix",
                      "BiocFileCache","AnnotationHub","SingleR","edgeR","GSEABase",
                      "cowplot", "ggplot2", "gridExtra", "grid","scMerge", "viridis")

for (pkg in requiredPackages){
  suppressWarnings(suppressMessages(library(pkg, character.only = T)))
}




##### ----- Loading & processing data -----------------------------------

loadSCE <- function(prefix,dir, metadat_x){
  xx <- edgeR::read10X(mtx = paste0(prefix,"matrix.mtx.gz"), 
                       genes = "features.tsv.gz", 
                       barcodes = paste0(prefix,"barcodes.tsv.gz"),
                       path =dir,
                       DGEList = FALSE)
  dat_x <- xx$counts
  colnames(dat_x) <- paste0(prefix, as.vector(xx$samples$Barcode))
  
  metadat_x <- metadat_x[metadat_x$CellID %in% colnames(dat_x),]
  rownames(metadat_x) <- metadat_x$CellID
  dat_x <- dat_x[,metadat_x$CellID]
  
  
  sce_x <- SingleCellExperiment(list(counts = dat_x), 
                                rowData = xx$genes,
                                colData=metadat_x) 
  sce_x <- sce_x[!(duplicated(rowData(sce_x)$Symbol)),]
  rownames(sce_x) <- rowData(sce_x)$Symbol
  sce_x$CellID <- colnames(dat_x)
  return(sce_x)
}

sceBind <- function(sce_lst,colDat_to_keep=NULL,rowsumThresh=5){
  require(scMerge)
  
  if(is.null(names(sce_lst))){
    sce_names <- paste0("batch_", 1:length(sce_lst))
  } else {
    sce_names <- names(sce_lst)
  }
  
  sce_lst <- lapply(sce_lst,function(x){asCountMat(x)})
  
  sce_obj <- sce_cbind(sce_lst,method ="union",
                       cut_off_batch =0,
                       cut_off_overall =0,
                       exprs ="counts",
                       batch_names = sce_names,
                       colData_names=colDat_to_keep)
  
  counts(sce_obj) <- as(counts(sce_obj), "dgCMatrix")
  
  ## Normalise
  sce_obj <- sce_obj[rowSums(counts(sce_obj)>0)>rowsumThresh]
  set.seed(1000)
  clusters <- quickCluster(sce_obj)
  sce_obj <- scran::computeSumFactors(sce_obj, cluster=clusters, min.mean=0.1) 
  sce_obj <- logNormCounts(sce_obj)
  
  ## Run dimension reduction
  sce_obj <- scater::runTSNE(sce_obj, name="TSNE", exprs_values="logcounts")
  
  return(sce_obj)
}

asCountMat <- function(sce_obj, mat_type="matrix"){
  if (mat_type=="matrix"){
    counts(sce_obj) <- as.matrix(counts(sce_obj))
  } else if (mat_type=="dgCMatrix"){
    counts(sce_obj) <- as(counts(sce_obj), "dgCMatrix")
  }
  return(sce_obj)
}

scranNorm <- function(sce_object){
  set.seed(1000)
  
  clusters <- quickCluster(sce_object)
  sce_object <- scran::computeSumFactors(sce_object, cluster=clusters, min.mean=0.1) 
  sce_object <- logNormCounts(sce_object)
  
  sce_object$ConvSizeFact <-  sizeFactors(sce_object)
  sce_object$log10Sum <-  log10(sce_object$sum)
  return(sce_object)
}


GetPosCellIDs <- function(sce_list, gene){
  return(as.vector(unlist(lapply(sce_list, 
                                 function(x){
                                   as.vector(x[,as.vector(counts(x[rownames(x)==gene,])>0)]$CellID)
                                 }))))
}



GetPseudoBulk <- function(sce_x, cluster_str, min_cell_numb=50,
                          cpm_thresh=0.5,sample_thresh=0.1){  
  
  ## Drop irrelevant clusters
  clust_vec <- sce_x[[cluster_str]]
  clust_to_keep <- names(table(clust_vec)[table(clust_vec)>min_cell_numb])
  sce_x <- sce_x[,clust_vec %in% clust_to_keep]
  
  Summed_Sce <- aggregateAcrossCells(sce_x, id=sce_x[[cluster_str]])
  DGE_Summed_sce <- DGEList(counts(Summed_Sce))
  
  ## Discared labels with very small library sizes
  discarded <- isOutlier(DGE_Summed_sce$samples$lib.size, log=TRUE,
                         type="lower")
  DGE_Summed_sce <- DGE_Summed_sce[,!discarded]
  
  ## Remove genes that are lowly expressed 
  keep = rowMeans(edgeR::cpm(y=DGE_Summed_sce, log = TRUE) >= cpm_thresh) >= sample_thresh
  DGE_Summed_sce <- DGE_Summed_sce[keep, ]
  
  ## Normalise by library size
  DGE_Summed_sce <- calcNormFactors(DGE_Summed_sce)  
  
  DGE_Summed_sce$logCPM <- edgeR::cpm(DGE_Summed_sce, log=TRUE, prior.count = 1) 
  
  return(DGE_Summed_sce)
}



ConvertSeurat <- function(sce_obj){
  
  if (!(is.null(SingleCellExperiment::altExpNames(sce_obj)))){
    sce_obj <- SingleCellExperiment::removeAltExps(sce_obj)
  }
  return(Seurat::as.Seurat(sce_obj))
  
}

SeuratIntegrate <- function(sce_obj,split_col=NULL, n_features=2000){
  require(Seurat)
  
  if (is.list(sce_obj)){
    seurat_lst <- lapply(sce_obj,function(x){ConvertSeurat(x)})
  } else if (is.null(split_col)){
    message("Please provide `split_col`")
    return(NULL)
  } else {
    seurat_obj <- ConvertSeurat(sce_obj)
    seurat_lst <- Seurat::SplitObject(seurat_obj,split.by =split_col)
  }
  
  features <- Seurat::SelectIntegrationFeatures(object.list = seurat_lst,nfeatures=n_features)
  enrich.anchors <- Seurat::FindIntegrationAnchors(object.list = seurat_lst, anchor.features = features)
  obj.combined <- Seurat::IntegrateData(anchorset = enrich.anchors)
  sce_obj <- as.SingleCellExperiment(obj.combined)
  sce_obj <- RunDimRed(sce_obj)
  return(sce_obj)
}




##### ----- Plotting themes -----------------------------------
rl=1
PlainBar_theme <- theme(legend.text = element_text(size = rel(rl)*1.3),
                        legend.text.align = 0.5,legend.title.align = 0.5,
                        legend.spacing.x = unit(0.2, 'cm'),
                        legend.spacing.y = unit(0.1, 'cm'),
                        legend.margin=margin(0,0,0,0),
                        legend.box.margin=margin(-5,-5,-5,-5),
                        legend.box.spacing= unit(0.1, 'cm'),
                        legend.title = element_text(size = rel(rl)*1.6, face = 'bold'),
                        legend.position = 'bottom',
                        panel.border = element_rect(colour = NA),
                        panel.grid.minor = element_blank(),
                        panel.grid.major = element_blank(),
                        axis.text.y  = element_text(angle=90, hjust=0.5, size = rel(rl)), 
                        axis.text.x = element_blank(), 
                        axis.line.y = element_blank(),axis.line.x = element_blank(),
                        axis.ticks.y=element_blank(), axis.ticks.x=element_blank(),
                        axis.title.y =element_blank(),  axis.title.x=element_blank())



MARKER_THEME <- theme(axis.ticks = element_blank(), panel.border=element_blank(),
                      axis.title = element_blank(), axis.text=element_blank(), 
                      axis.line = element_blank(),legend.title=element_blank())



