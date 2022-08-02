# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
#   PRODUCE FIGURE 3 -- Malignant & Mesothelial cells
#   
#   To recreate Figure 3
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


library(ggRNA)
library(Seurat)
library(reshape2)

##### ----- Get directories
setwd("./MPE_paper_analysis/")
data_dir <- paste0(getwd(),"/data/")
matrix_dir <- paste0(data_dir,"raw_matrices/")
fig_dir <- paste0(getwd(),"/figures/")

source(paste0(getwd(),"/scripts/helper_functions.R"))




### COLOURS ------------

Colour_values_PatientID <- c("#077187", "#6AAAB7", "#8E2043", "#BC7A8F", 
                             "#FEA090","#FECFC7", "#3E5496", "#0A9086", 
                             "#E0607E", "#E1D164")
names(Colour_values_PatientID) <- c("BCB66","BCB66_E", "BCB20","BCB20_E", "BCB21", 
                                    "BCB21_E", "BCB112",  "BCB139",  "BCB90",  "BCB114")

Colour_values <- c("#9e5476", "#d18b79")
names(Colour_values) <- c("Malignant", "Mesothelial")

Colour_values_AIMS <- c("#CABEE9", "#7C7189", "#FAE093", "#D04E59", "#2F3D70")
names(Colour_values_AIMS) <- c( "LumA", "LumB", "Her2","Normal", "Basal")


##### ----- Load data ----------------------------------------------------------

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


### --- Combine

### Use refined celltype labels
Non_immune_cells <- ALL_COLDATA[ALL_COLDATA$CellType_Fig3 %in% c("Malignant","Mesothelial"),]$CellID
Non_immune_sce_lst <- lapply(Non_immune_sce_lst, function(x){x[,x$CellID %in% Non_immune_cells]})

Fig3_sce <- sceBind(Non_immune_sce_lst, colDat_to_keep =colnames(ALL_COLDATA),rowsumThresh=0)
Fig3_sce$PatientID <- Fig3_sce$batch

#save(Fig3_sce, file=paste0(data_dir, "Fig3_sce.Rdata"))




##### ----- INTEGRATE BATCHES -----------------------------------


## -- Add chemistry version
SeqVers_dict <- c("ChemistryV3","ChemistryV3",  "ChemistryV3",  "ChemistryV3", "ChemistryV3", "ChemistryV2", "ChemistryV2", "ChemistryV2")
names(SeqVers_dict) <- c("BCB139",  "BCB112","BCB21",  "BCB66", "BCB20", "BCB21_E",  "BCB66_E", "BCB20_E")
Fig3_sce$ChemV <- as.vector(SeqVers_dict[Fig3_sce$PatientID])

## -- Integrate by chemistry version for plotting tSNEs
seurat_obj <- Seurat::as.Seurat(Fig3_sce)
seurat_lst <- Seurat::SplitObject(seurat_obj,split.by ="ChemV")
features <- Seurat::SelectIntegrationFeatures(object.list = seurat_lst,nfeatures=2000)
enrich.anchors <- Seurat::FindIntegrationAnchors(object.list = seurat_lst, anchor.features = features)
obj.combined <- Seurat::IntegrateData(anchorset = enrich.anchors)
Fig3_merged_sce <- Seurat::as.SingleCellExperiment(obj.combined)
Fig3_merged_sce <- scater::runTSNE(Fig3_merged_sce)

#save(Fig3_merged_sce, file=paste0(data_dir, "Fig3_merged_sce.Rdata"))



##### ----- PLOT FIG 3A TSNES -----------------------------------

# load(file=paste0(data_dir, "Fig3_merged_sce.Rdata"))

alph =0.8
siz=5
A1 <- ggDIMRED(Fig3_merged_sce, dimred="TSNE", colour_by = "PatientID",
               point_size=siz, point_alpha =alph, col_pal=Colour_values_PatientID)
A2 <- ggDIMRED(Fig3_merged_sce, dimred="TSNE", colour_by = "CellType_Fig3",
               point_size=siz,point_alpha =alph, col_pal=Colour_values)
A3 <- ggDIMRED(Fig3_merged_sce, dimred="TSNE", colour_by = "AIMS_subtype",
               point_size=siz,point_alpha =alph,col_pal = Colour_values_AIMS)


png(paste0(fig_dir,"Fig3A_byPatientID.png"),width=2000, height=1800)
A1+theme_blank()
dev.off()

pdf(paste0(fig_dir,"Fig3A_byPatientID.pdf"),width=2000, height=1800)
A1+theme_blank()
dev.off()

pdf(paste0(fig_dir,"Fig3A_byPatientID.pdf"),  width=7, height=6)
ggDIMRED(Fig3_merged_sce, dimred="TSNE", colour_by = "PatientID",
         point_size=1, point_alpha =1, col_pal=Colour_values_PatientID)+theme_blank()
dev.off()


reso <- 500
length <- 3.25*reso/72
png(paste0(fig_dir,"Fig3A_byCellType.png"), units="in", width=length, height=length*0.9, res=reso)
A2+theme_blank()
dev.off()

reso <- 500
length <- 3.25*reso/72
png(paste0(fig_dir,"Fig3A_bySubtype.png"), units="in", width=length, height=length*0.9, res=reso)
A3+theme_blank()
dev.off()


##### ----- PLOT FIG 3B PROPORTIONAL BARS -----------------------------------

# load(file=paste0(data_dir, "Fig3_sce.Rdata"))

CellType_Ordering <- c("Malignant", "Mesothelial", "Unassigned")
SampleCols <- c("BCB139",  "BCB112", "BCB66","BCB20", "BCB21", "BCB66_E", "BCB20_E","BCB21_E" )

## Get proportions
CellType_Props <- ggRNA::GetProps_perPatient(Fig3_sce, label="CellType_Fig3", group = "PatientID")
CellType_Props$CellType <- rownames(CellType_Props)
CellType_Props <- CellType_Props[order(factor(CellType_Props$CellType, levels=rev(unique(CellType_Ordering)))),]
CellType_Props <- reshape2::melt(CellType_Props, id.vars="CellType", value.name="Proportion", variable.name="PatientID")

## Add cell number info
SampleNumbs_dict <- unlist(lapply(SampleCols, function(x)  {ncol(Fig3_sce[,Fig3_sce$PatientID == x])}))
names(SampleNumbs_dict) <- SampleCols
CellType_Props$CellNumb <- as.vector(SampleNumbs_dict[as.vector(CellType_Props$PatientID)])

## Reorder
CellType_Props$CellType <- factor(CellType_Props$CellType, levels=rev(CellType_Ordering))
CellType_Props$PatientID <- factor(CellType_Props$PatientID, levels=rev(SampleCols))

## Add Enrichment info
EnrichmentDict <- c("NonEnriched", "NonEnriched", "NonEnriched", "NonEnriched", "NonEnriched",
                    "Enriched", "Enriched", "Enriched")
names(EnrichmentDict) <- SampleCols
CellType_Props$Enrich <- as.vector(EnrichmentDict[as.vector(CellType_Props$PatientID)])
CellType_Props$Enrich <- factor(CellType_Props$Enrich, levels=c("NonEnriched","Enriched"))


PropBar_pp <- ggplot(CellType_Props, aes(x=PatientID, y=Proportion, fill=CellType)) + 
  geom_bar(stat = "identity", position = "stack", width=0.95) +
  xlab("\nPatientID") +ylab("\nLabel Proportion") +
  scale_fill_manual(values=Colour_values, breaks=CellType_Ordering)+theme_bw()+  coord_flip()+
  scale_y_discrete(limits =CellType_Ordering)+
  PlainBar_theme+ labs(fill = "Cell Type")+
  facet_grid(rows=vars(Enrich), scales="free", space="free_y")+
  theme(strip.background = element_blank(), strip.text = element_blank())

cellNumb_bar <- ggplot(data=CellType_Props, aes(x=PatientID, y=CellNumb)) + 
  geom_bar(stat="identity", position="dodge", width=0.93,  fill="dimgrey")+
  coord_flip()+theme_bw()+PlainBar_theme+
  theme(axis.text.x=element_text(size = rel(rl)*1.2, hjust=0.8),
        axis.text.y=element_blank(),axis.line.x=element_line(colour="black",size=0.2))+
  facet_grid(rows=vars(Enrich), scales="free", space="free_y")+
  scale_y_continuous(labels=c("0","2k", "4k", "6k", "8k"), expand=c(0,0),
                    limits=c(0, signif(max(CellType_Props$CellNumb),2)+500),breaks = c(0, 2000, 4000, 6000,8000))


### Combine plots
pdf(paste0(fig_dir,"Fig3_PropBar_ByEnrichment.pdf"),  width=7, height=8)
grid.arrange(grobs=list(PropBar_pp+theme(legend.position = "none",plot.margin=unit(c(0,0,0.5,0), "cm")), 
                        cellNumb_bar+theme(plot.margin=unit(c(0,0,0,0), "cm"),panel.grid.major.x=element_line(size=0.2),
                                           strip.background = element_blank(), strip.text = element_blank())),
             ncol=2, nrow=1, 
             widths = c(5,1))
dev.off()




##### ----- PLOT FIG 3C MARKER GENE BOXPLOT -----------------------------------

# load(file=paste0(data_dir, "Fig3_sce.Rdata"))
markers <- list(Epithelial=c("EPCAM","CLDN4", "IMP3", "MUC1", "CDH1", "CEACAM5", "SCGB2A2", "CLDN3", "CLDN7"),
                Mesothelial=c("DES", "CALB2",  "WT1", "UPK3B",  "CDH2", "COL1A2", "S100A4", "MSLN"),
                Joint=c("KRT8", "KRT19", "VIM",  "CD44", "ACTB"))

library(scfunc)

### --- Calculate logTPM
sce_tmp <- GetGeneLengths(Fig3_sce)  #data_dir
TPM_dat <- scuttle::calculateTPM(sce_tmp)
TPM_dat <- TPM_dat[rownames(TPM_dat) %in% as.vector(unlist(markers)),] 
TPM_dat <- as.data.frame(log2(TPM_dat+1))

### --- Average because sparsity
E <- as.vector(colMeans(TPM_dat[rownames(TPM_dat) %in% markers$Epithelial,]))
M <- as.vector(colMeans(TPM_dat[rownames(TPM_dat) %in% markers$Mesothelial,]))
J <- as.vector(colMeans(TPM_dat[rownames(TPM_dat) %in% markers$Joint,]))

plot_df <- data.frame(AverageExpression=c(E, M, J), 
                      CellType=rep(Fig3_sce$CellType_Fig3, 3),
                      MarkerType=c(rep("Epithelial Genes", length(E)),
                                   rep("Mesothelial Genes", length(M)),
                                   rep("Overlapping", length(J))))


pp_box <- ggplot(plot_df, aes(x=CellType , y=AverageExpression, fill=CellType))+
  geom_boxplot(alpha=1,outlier.shape=1,lwd=0.3, outlier.size = 0.5) + facet_wrap(~MarkerType, nrow=1) +
  scale_fill_manual(values=Colour_values)+xlab("")+ylab("Average Marker Expression")+
  scale_y_continuous(expand=c(0,0))+
  theme(legend.position="none",
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),axis.ticks.x=element_blank(),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 14, margin = margin(t = 0, r = 10, b = 0, l = 0)),
        axis.title.x = element_blank(),
        strip.text.x = element_text(size = rel(rl)*1.4))+
      coord_cartesian(ylim = c(0, round(max(plot_df$AverageExpression))))

pdf(paste0(fig_dir,"Fig3_markersBoxplot.pdf"), width=8, height=4)
print(pp_box)
dev.off()

png(paste0(fig_dir,"Fig3_markersBoxplot.png"), width=2000, height=1000, res=280)
print(pp_box)
dev.off()




##### ----- SUPP FIGURE 7 UPSET PLOT -----------------------------------


library(ComplexUpset)
marker_genes <- c("EPCAM", "DES")
dat_x <- as.data.frame(t(as.matrix(logcounts(Fig3_sce[rownames(Fig3_sce) %in% marker_genes,]) > 0)))
dat_x$Celltype <- as.vector(setNames(Fig3_sce$CellType_Fig3, Fig3_sce$CellID)[rownames(dat_x)])

table(dat_x[dat_x$Celltype=="Malignant",]$EPCAM)
table(dat_x[dat_x$Celltype=="Malignant",]$DES)

pdf(paste0(fig_dir, "SUPP_Upset.pdf"), width = 7, height=6)
ComplexUpset::upset(data=dat_x, intersect=marker_genes, name= "Markers",
                    width_ratio=0.27,height_ratio =0.45,
                    annotations = list(
                      'Cell Type \n (% cells)'=list(
                        aes=aes(x=intersection, fill=Celltype),
                        geom=list(
                          geom_bar(stat='count', position='fill'),
                          theme(legend.position = "top", 
                                axis.text.y =element_text(size = rel(rl)*0.9), 
                                axis.text.x=element_blank(),axis.title.x=element_blank(),
                              #  axis.text.x = element_text(size = rel(rl)*0.3),
                                axis.title.y =element_text(size = rel(rl)*1.1)), 
                            #    axis.title.x =element_text(size = rel(rl)*1.1)),
                          scale_y_continuous(labels=scales::percent_format(), expand=c(0,0)),
                          scale_fill_manual(values = Colour_values[as.character(unique(dat_x$Celltype))]) 
                        )
                      )
                    ),
                    themes=upset_modify_themes(list(overall_sizes=theme(axis.text=element_text(size=6))))
                    )
dev.off()



