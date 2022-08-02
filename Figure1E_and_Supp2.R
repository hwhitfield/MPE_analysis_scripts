
##### ----- Get directories


data_dir <- paste0(getwd(),"/MPE_paper_analysis/data/")
matrix_dir <- paste0(getwd(),"/MPE_paper_analysis/raw_matrices/")
pal_etal_dir <- "~/data/_datasets/sc/Smyth_HumanBrCa/"
pal_save_dir <- paste0(data_dir,"pal_etal_data/")
fig_dir <- paste0(getwd(),"/MPE_paper_analysis/figures/")

source(paste0(getwd(),"/MPE_paper_analysis/scripts/helper_functions.R"))


### COLOURS ------------
Colour_values <- c("#9e5476",
                   "#d18b79","#CA7B83",
                   "#B2DBA0", "#70A18F", "#637c8f",
                   "#886598",
                   "#e7ebbc",
                   "#6074ab", 
                   "#C7BBDD",
                   "grey", "grey","grey")

names(Colour_values) <- c("Epithelial",
                          "CAFs","Endothelial",
                          "TAMs","Myeloid", "DCs",  
                          "Pericytes",
                          "B cells",
                          "T cells",
                          "Plasma cells",
                          "12","13","5")

Colour_values2 <- c("#9e5476",#"#9e5476", 
                    "#d18b79",# "#C7BBDD",
                    "#CA7B83","#C7BBDD",
                   "#b2dba0","#70a18f", "#637c8f",#"#303E70",
                   "#e7ebbc", "#6074ab","#8bbde6")

names(Colour_values2) <- c("Malignant",#"TNBC", 
                           "Mesothelial",
                           "CAFs","Endothelial",
                          "TAMs", "Myeloid", "DCs",#"Pericytes",
                          "B cells","T cells", "NK cells")

Colour_values_3 <- c("#9e5476","#2F3D70")
names(Colour_values_3) <- c("ER", "TNBC")


##### ----- Load primary brca data & process --------------------------

library(Seurat)

### --- All
ERTotal <- readRDS(file=paste0(pal_etal_dir, "SeuratObject_ERTotal.rds"))
TNBC <- readRDS(file=paste0(pal_etal_dir, "SeuratObject_TNBC.rds"))
pal_metadata_all <- rbind(ERTotal@meta.data[,colnames(ERTotal@meta.data) %in% c("orig.ident","nCount_RNA","nFeature_RNA","group","seurat_clusters")], 
                          TNBC@meta.data[,colnames(TNBC@meta.data) %in% c("orig.ident","nCount_RNA","nFeature_RNA","group","seurat_clusters")])
pal_metadata_all$CellID <- rownames(pal_metadata_all)
rm(ERTotal, TNBC)
gc()

### --- Get celltypes
TNBCSub_cellTypes <- c("T cells","TAMs","Plasma cells",
                       "CAFs","T cells", "B cells","DCs",
                       "Endothelial", "Pericytes", "Myeloid")
names(TNBCSub_cellTypes) <- as.character(0:9)

ERSub_cellTypes <- c("T cells","TAMs","CAFs",
                       "Pericytes","5", "Endothelial",
                     "TAMs","B cells", "Myeloid", "CAFs",
                     "Plasma cells", "12", "13")
names(ERSub_cellTypes) <- as.character(0:12)

### --- Sub
ERTotalSub <- readRDS(file=paste0(pal_etal_dir, "SeuratObject_ERTotalSub.rds"))
ERTotalSub@meta.data$celltype <- as.vector(ERSub_cellTypes[as.vector(ERTotalSub$seurat_clusters)])

TNBCSub <- readRDS(file=paste0(pal_etal_dir, "SeuratObject_TNBCSub.rds"))
TNBCSub@meta.data$celltype <- as.vector(TNBCSub_cellTypes[as.vector(TNBCSub$seurat_clusters)])

pal_metadata_sub <- rbind(ERTotalSub@meta.data[,colnames(ERTotalSub@meta.data) %in% c("celltype","nCount_RNA","nFeature_RNA","group","seurat_clusters")], 
                          TNBCSub@meta.data[,colnames(TNBCSub@meta.data) %in% c("celltype","nCount_RNA","nFeature_RNA","group","seurat_clusters")])
pal_metadata_sub$CellID <- rownames(pal_metadata_sub)
rm(ERTotalSub, TNBCSub)
gc()


### --- Add cancer cell annot
celltype_dict <- setNames(pal_metadata_sub$celltype,pal_metadata_sub$CellID)
celltype_dict[setdiff(pal_metadata_all$CellID, pal_metadata_sub$CellID)] <- "Epithelial"
pal_metadata_all$celltype <- as.vector(celltype_dict[pal_metadata_all$CellID])

pal_metadata_all$PatientID <- pal_metadata_all$group
write.table(pal_metadata_all,paste0(pal_etal_dir, "pal_metadata_all.tsv"),sep="\t")




##### ----- SUPP FIGURE 2 PROP BARS ---------------------------------------------------


pal_metadata_all <- read.table(paste0(pal_etal_dir, "pal_metadata_all.tsv"),sep="\t")

plot_dat <- pal_metadata_all[,colnames(pal_metadata_all) %in% c("celltype","PatientID", "CellID")]
plot_dat <- ggRNA::GetProps_perPatient(plot_dat, "celltype","PatientID")
plot_dat[is.na(plot_dat)] <- 0
plot_dat$celltype <- rownames(plot_dat)
dat_melt = reshape2::melt(plot_dat, variable.name="PatientID")

ordering <- c("Epithelial","Endothelial", "CAFs", "Pericytes",
              "T cells","TAMs", "Plasma cells", "B cells",  "DCs","Myeloid",
              "12", "13", "5")
dat_melt$celltype <- factor(dat_melt$celltype, levels=rev(ordering))

pdf(paste0(fig_dir,"SUPP2_cellNumbs.pdf"),width=8,height=4)
ggplot(dat_melt, aes(fill=celltype, y=value, x=PatientID)) + 
  geom_bar(position="fill", stat="identity",na.rm=T)+labs(y="CellNumber",fill="")+theme_bw()+
  scale_fill_manual(values=Colour_values[!(names(Colour_values) %in% c("12","13","5"))])+
  scale_y_continuous(expand=c(0,0),labels = scales::percent_format()) +
  theme(axis.text.x=element_text(angle=45,hjust=1),axis.ticks.x=element_blank(),
        axis.ticks.y = element_line(size = 0.3))
  
dev.off()





##### ----- LOAD DATA ---------------------------------------------------


### --- Load primary brca data
pal_files <- list.files(pal_dir)
pal_files <- pal_files[grepl("GSM",pal_files)]
pal_files <- gsub("barcodes.tsv.gz","",pal_files)
pal_files <- gsub("matrix.mtx.gz","",pal_files)
pal_files<- unique(pal_files)

## Drop lymph node/tumour pairs
pal_files <- pal_files[!(grepl("LN-",pal_files))]
pal_files <- pal_files[!(grepl("0056|0064|0173|0029",pal_files))]

## ER/TN only
pal_files <- pal_files[grepl("_ER-|_TN-",pal_files)]

pal_names <- gsub(".*_", "", pal_files)
pal_names <- gsub("MH|AH|PM|SH|Tum","",pal_names)
pal_names <- trimws(pal_names,whitespace = "[-]")
pal_names <- gsub("-","_",pal_names)

pal_dict <- setNames(pal_files,pal_names)
pal_dict["ER_0040_T"] <- pal_dict["ER_0040"]

pal_metadata_all$CellID <- paste0(as.vector(pal_dict[pal_metadata_all$PatientID]), 
                                  gsub(".*_", "",rownames(pal_metadata_all)))

PAL_lst <- lapply(pal_files, 
                  function(bcb){loadSCE(bcb,
                                        pal_dir,
                                        pal_metadata_all)})
names(PAL_lst) <- pal_names


save(PAL_lst, file=paste0(pal_save_dir, "PAL_lst.Rdata"))



##### ----- Get data

### --- Merge
PAL_lst <- lapply(PAL_lst,function(x){asCountMat(x)})
PAL_sce <- sce_cbind(PAL_lst,method ="union",
                     cut_off_batch =0,
                     cut_off_overall =0,
                     exprs ="counts",
                     batch_names = names(PAL_lst),
                     colData_names=c("celltype", "PatientID", "CellID","seurat_clusters"))
counts(PAL_sce) <- as(counts(PAL_sce), "dgCMatrix")

rm(PAL_lst)
gc()

### --- Normalise
PAL_sce <- scater::addPerCellQC(PAL_sce)
PAL_sce <- scranNorm(PAL_sce)

save(PAL_sce, file=paste0(pal_save_dir, "PAL_sce.Rdata"))




##### ----- FIGURE 1E CELLS DETECTED HEATMAP -----------------------------------

load(file=paste0(pal_save_dir, "PAL_sce.Rdata")) # -- 112,685 cells
load(paste0(data_dir,"Fig1_sce.Rdata"))


### --- Fix up cell types 
PAL_sce <- PAL_sce[,!(PAL_sce$celltype %in% c("13", "12", "5", "Pericytes", "Plasma cells","TAMs","DCs"))]
PAL_sce$celltype <- replace(PAL_sce$celltype, PAL_sce$celltype=='Epithelial', "Malignant")

Fig1_sce <- Fig1_sce[,!(Fig1_sce$CellType_Fig1 == "Unassigned")]
Fig1_sce$celltype <- sub("_"," ",Fig1_sce$CellType_Fig1)


### --- Get prop data
cell_ordering <- c("Malignant",
                   "Mesothelial","CAFs","Endothelial",
                   "Myeloid", "NK cells","T cells", "B cells")
genes = c("KLRF1", "DES",
          "CD19",
          "CD79A","CD3D","CDH2", 
          "ESR1","CDH1","EPCAM") 

getPropData <- function(sce_x,genes_to_plot,group_str="celltype"){
  prop_lst <- list()
  for (i_ct in unique(sce_x[[group_str]])){
    prop_lst[[i_ct]] <- rowMeans(logcounts(sce_x[rownames(sce_x) %in% genes_to_plot,sce_x[[group_str]]==i_ct])>0)
  }
  plot_dat <- as.data.frame(prop_lst)
  plot_dat$gene <- factor(rownames(plot_dat), levels=genes_to_plot)
  
  plot_dat <- reshape2::melt(plot_dat, value.name="Detection", variable.name ="CellType")
  plot_dat$CellType <- gsub("[.]", " ",as.vector(plot_dat$CellType))
  plot_dat$CellType <- factor(plot_dat$CellType, levels=levels(sce_x[[group_str]]))
  plot_dat$Detection <- plot_dat$Detection*100
  return(plot_dat)
}


MPE_props <- getPropData(Fig1_sce,genes)
PAL_props <- getPropData(PAL_sce,genes)
MPE_props$Dataset <- rep("MPE", nrow(MPE_props))
PAL_props$Dataset <- rep("Primary", nrow(PAL_props))

prop_dat <- rbind(MPE_props,PAL_props)

prop_dat$CellType <- factor(prop_dat$CellType, levels=cell_ordering)

pdf(paste0(fig_dir,"FIG1E_cellsDetected.pdf"),width=8,height=3)
ggplot(data=prop_dat, aes(x=CellType, y=gene, fill=Detection))+
  geom_tile()+facet_grid(~Dataset, scales="free")+ labs(fill="% Cells \nDetected")+
  scale_x_discrete(expand=c(0,0))+scale_y_discrete(expand=c(0,0))+
  scale_fill_distiller(palette ="PuRd",direction=1)+
  theme(strip.text = element_text(size = 16),
        strip.background = element_blank(),strip.placement = "outside",
        axis.ticks = element_blank(), axis.text.x=element_text(angle=45,hjust=1,size=12),
        axis.title.y=element_blank(),axis.title.x=element_blank(),
        axis.text.y=element_text(size=12))
                       
dev.off()


