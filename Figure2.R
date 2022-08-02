# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
#   PRODUCE FIGURE 2 -- Immune cells, non-enriched samples
#   
#   To recreate Figure 2
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


##### ----- Load libraries
requiredPackages <- c("reshape2","ProjecTILs","Seurat","ggRNA","viridis")

for (pkg in requiredPackages){
  suppressWarnings(suppressMessages(library(pkg, character.only = T)))
}


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

Immune_Colour_values <- c("#b2dba0", "#637c8f", "#70a18f",
                          "#e7ebbc", "#8bbde6", "#6074ab", 
                          "grey")

names(Immune_Colour_values) <- c("Macrophage", "DC", "Myeloid",  
                                 "B_cells","NK_cells", "T_cells",
                                 "Unassigned")


##### ----- PREP DATA FOR PLOT -------------------------------------------------


### --- Load colData

ALL_COLDATA  <- read.table(file=paste0(data_dir,"ALL_COLDATA.tsv"),sep="\t")


### --- Load count matrices into SingleCellExperiment objects
BCB_names <- c("BCB90","BCB114","BCB112","BCB139",
               "BCB66","BCB20","BCB21")
non_enriched_lst <- lapply(BCB_names, 
                           function(bcb){
                             loadSCE(paste0(bcb,"_"), 
                                     paste0(matrix_dir,bcb),
                                     ALL_COLDATA)})
names(non_enriched_lst) <- BCB_names


### --- Remove non-immune cells
immune_cells <- ALL_COLDATA[!(ALL_COLDATA$CellType %in% c("Malignant","Mesothelial", "Unassigned")),]$CellID

## Keep potentially immune 'unassigned' cells
Epithelial_marker_pos <- unlist(lapply(c("EPCAM","CLDN4",  "CLDN3", "CLDN7"), 
                                       function(x){GetPosCellIDs(non_enriched_lst,x)}))
Mesothelial_marker_pos <- unlist(lapply(c("WT1", "MSLN","DES"), 
                                        function(x){GetPosCellIDs(non_enriched_lst,x)}))
unassigned_Immune_cells <- setdiff(ALL_COLDATA[ALL_COLDATA$CellType =="Unassigned",]$CellID,
                                   c(Epithelial_marker_pos,Mesothelial_marker_pos))
unassigned_Immune_cells <- setdiff(unassigned_Immune_cells,
                                   ALL_COLDATA[(ALL_COLDATA$Izar_label %in% c("Malignant","Fibroblast")) |
                                                 (ALL_COLDATA$LUAD_label %in% c("Epithelial cells")) |
                                                 ALL_COLDATA$CellAtlas_labels %in% c("Epithelial_cells"),]$CellID)

non_enriched_lst <- lapply(non_enriched_lst, function(x){x[,x$CellID %in% c(immune_cells, unassigned_Immune_cells)]})

### --- Combine patients
Fig2_sce <- sceBind(non_enriched_lst, colDat_to_keep =colnames(ALL_COLDATA),rowsumThresh=0)
Fig2_sce$PatientID <- Fig2_sce$batch
save(Fig2_sce, file=paste0(data_dir, "Fig2_sce.Rdata"))





##### ----- SUBCLUSTER T, NK & B CELLS -----------------------------------------


SubsetCells <- function(sce_x,celltypes){

  subset_sce_x <- sce_x[,sce_x$CellType %in% celltypes]
  subset_sce_x <- subset_sce_x[!(grepl("^MT-",rownames(subset_sce_x))),]
  
  ## Process
  subset_sce_x <- subset_sce_x[rowSums(counts(subset_sce_x)>0)>5,]
  subset_sce_x <- scranNorm(subset_sce_x)
  
  ## Dimension reduction
  subset_sce_x <- scater::runTSNE(subset_sce_x)

  return(subset_sce_x)
}

#load(file=paste0(data_dir, "Fig2_sce.Rdata"))


Tcell_sce <- SubsetCells(Fig2_sce,c("T_cells"))
NKTcell_sce <- SubsetCells(Fig2_sce, c("T_cells","NK_cells"))

Bcell_sce <- SubsetCells(Fig2_sce,c("B_cells", "Unassigned"))
#save(Bcell_sce, file=paste0(data_dir, "Fig2_Bcell_sce.Rdata"))


### ---- SUPP Myeloid tSNE
Myeloid_sce <- SubsetCells(Fig2_sce,c("Myeloid","Macrophage", "DC"))

myeloid_leg <- GetLegend(Immune_Colour_values[c("Myeloid","Macrophage", "DC")], "Cell Type")
pdf(paste0(fig_dir,"SUPP_MyeloidtSNE.pdf"),  width=7, height=6)
cowplot::plot_grid(ggRNA::ggDIMRED(Myeloid_sce, dimred="TSNE", colour_by="CellType",col_pal = Immune_Colour_values,
                                   point_alpha=0.9, point_size=1.6, ordering=c("Myeloid","Macrophage", "DC"))+theme_blank(),
                   myeloid_leg, ncol=1, rel_heights = c(6,1))
dev.off()




##### ----- T CELL SUBTYPING ANALYSIS ------------------------------------------

### --- Predict T-cell subtypes usnig ProjecTILs & Seurat
require(ProjecTILs)
require(Seurat)
ref <- ProjecTILs::load.reference.map()

T.seurat <- CreateSeuratObject(logcounts(Tcell_sce))
query.projected <- make.projection(T.seurat, ref=ref, filter.cells = F,skip.normalize=TRUE)
query.predict <- cellstate.predict(ref=ref, query=query.projected)
table(query.predict$functional.cluster)

## Append to sce objects
prediction_dict <- query.predict$functional.cluster
prediction_dict[setdiff(NKTcell_sce$CellID,names(prediction_dict))] <- "NK_cells"
NKTcell_sce$ProjecTILs_pred <- as.vector(prediction_dict[NKTcell_sce$CellID])
# save(NKTcell_sce, file=paste0(data_dir, "Fig2_NKTcell_sce.Rdata"))





##### ----- PLOT FIG 2A TSNE ------------------------------------------


# load(file=paste0(data_dir, "Fig2_sce.Rdata"))



A <- ggRNA::ggDIMRED(Fig2_sce, dimred="TSNE", colour_by="CellType",col_pal = Immune_Colour_values,
                          point_alpha=1, point_size=1.1)+theme_blank()
B <- ggRNA::ggDIMRED(Fig2_sce, dimred="TSNE", colour_by="PatientID",col_pal = Colour_values_PatientID,
              point_alpha=1, point_size=1.1)+theme_blank()
B_leg <- GetLegend(Colour_values_PatientID[names(Colour_values_PatientID) %in% Fig2_sce$PatientID], "PatientID")

pdf(paste0(fig_dir,"Fig2A_byCellType.pdf"), width=7, height=6)
A
dev.off()


png(paste0(fig_dir,"Fig2A_byCellType.png"),  width=2000, height=1800, res=250)
A
dev.off()

pdf(paste0(fig_dir,"SuppFig2A_byPatientID.pdf"),  width=7, height=6)
cowplot::plot_grid(B, B_leg, ncol=1, rel_heights = c(6,1))
dev.off()




##### ----- PLOT FIG 2B HEATMAP ------------------------------------------


gene_lst <- c("CD1C","ITGAX","HLA-DRB5", 
              "CD68", "CD64", "CCR5", "ITGAM", "FUT4", "CD14", 
              "KLRF1", "NCAM1","FCGR3A", "CD3D", "CD3G","CD8A","CD8B", 
              "IL5RA",  "IL3RA","SDC1", "JCHAIN", "CD19","MS4A1","CD79A")
common_name_dict <- setNames(c("CD56","CD16","CD64","CD20",
                               "CD11b","CD11c","CD15","CD125","CD123","CD138"),
                             c("NCAM1","FCGR3A","FCGR1A","MS4A1",
                               "ITGAM","ITGAX","FUT4","IL5RA","IL3RA","SDC1" ))
common_name_dict[setdiff(gene_lst,names(common_name_dict))] <- setdiff(gene_lst,names(common_name_dict))


celltypes <- c("DC","Myeloid","Macrophage",  "NK_cells","T_cells","Unassigned","B_cells")

ggheatmap <- ggRNA::scHeatmap(Fig2_sce, "CellType", gene_lst,
                       exprs_val="TPM",log=TRUE,scale_str = "genes",
                       cluster_genes =TRUE, cluster_samples=FALSE,
                       groups_to_plot=celltypes, min_cells=100)

row_names <- unique(as.vector(ggheatmap$data$variable))
alt_names <- as.vector(common_name_dict[row_names])
new_names <- paste0(row_names, "(",alt_names,")")
rows <- unlist(lapply(1:length(row_names),function(i){ifelse(row_names[[i]]==alt_names[[i]],row_names[[i]],new_names[[i]])}))

pdf(paste0(fig_dir,"FIG2_markerHeatmap.pdf"),  width=8, height=6)
ggheatmap+ scale_y_discrete(labels= rows,expand=c(0,0))+
  scale_fill_gradientn(colours= hcl.colors(100, palette = "Mako"))+
  theme(axis.text.x=element_blank())+
  guides(fill=guide_colorbar(title="",title.position ="bottom",direction="vertical",
                             barwidth =unit(0.025, units="npc"),
                             barheight = unit(0.5, units="npc")))
                              
dev.off()




##### ----- PLOT FIG 2C PROPORTIONAL BAR ----------------------------------------

### Get proportions
CellType_Props <- GetProps_perPatient(Fig2_sce, label="CellType", "PatientID")
CellType_Props$CellType <- rownames(CellType_Props)
CellType_Props <- reshape2::melt(CellType_Props, id.vars="CellType", value.name="Proportion", variable.name="PatientID")
CellType_Props <- CellType_Props[order(-CellType_Props$Proportion),]

## Add cell type
CellType_Ordering <- c("T_cells","NK_cells", "Myeloid", "Macrophage","DC", "B_cells","Unassigned")
CellType_Props$CellType <- factor(CellType_Props$CellType, levels=rev(CellType_Ordering))

## Add patientID
SampleCols <- c("BCB90", "BCB114","BCB21", "BCB20", "BCB139",  "BCB112", "BCB66")
CellType_Props$PatientID <- factor(CellType_Props$PatientID, levels=SampleCols)

## Add cell number info
SampleNumbs_dict <- unlist(lapply(SampleCols, function(x)  {ncol(Fig2_sce[,Fig2_sce$PatientID == x])}))
names(SampleNumbs_dict) <- SampleCols
CellType_Props$CellNumb <- as.vector(SampleNumbs_dict[as.vector(CellType_Props$PatientID)])


### Plot  
PropBar_pp <- ggplot(CellType_Props, aes(x = PatientID, y = Proportion, fill=CellType)) +  
  geom_bar(stat = "identity", position = "stack", width=0.95) +
  scale_fill_manual(values=Immune_Colour_values, breaks=CellType_Props$CellType)+theme_bw()+  coord_flip()+
  scale_y_discrete(limits =CellType_Props$CellType)+
  PlainBar_theme+theme(strip.background = element_blank(), strip.text = element_blank())


## Get cell numbs 
CellType_Numb <- CellType_Props[!(duplicated(CellType_Props$PatientID)),colnames(CellType_Props) %in% c("PatientID", "CellNumb")]
CellType_Numb$CellNumb <- as.numeric(CellType_Numb$CellNumb)

cellNumb_bar <- ggplot(data=CellType_Props, aes(x=PatientID, y=CellNumb)) + 
  geom_bar(stat="identity", position="dodge", width=0.93,  fill="dimgrey")+
  coord_flip()+theme_bw()+PlainBar_theme+
  theme(axis.text.x=element_text(size = rel(1)*1.2, hjust=0.8),axis.ticks.x=element_line(colour="black"),
        axis.text.y=element_blank(),axis.line.x=element_line(colour="black",size=0.2))+
  scale_y_continuous(labels=c("0","2k", "4k", "6k"), limits=c(0, 6000),
                     breaks = c(0, 2000, 4000, 6000),expand=c(0,0))


pdf(paste0(fig_dir, "FIG2_PropBar.pdf", ""), width=8, height=6)
grid.arrange(grobs=list(PropBar_pp+theme(legend.position = "none",plot.margin=unit(c(0,0,0.5,0), "cm")), 
                        cellNumb_bar+theme(plot.margin=unit(c(0,0,0,0), "cm"))),
             ncol=2, nrow=1, 
             widths = c(5,1))
dev.off()






##### ----- PLOT SUPP FIG B CELL TSNES ------------------------------------------


load(file=paste0(data_dir, "Fig2_Bcell_sce.Rdata"))


markers_to_plot <- c("CD52", "IGHM", "JCHAIN", "IRF4")
alph=0.65
siz=1.8
bar_h=0.05
bar_w=0.015

line_size=0.2


plot_lst <- lapply(markers_to_plot, function(x){
  ggDIMRED(Bcell_sce, dimred="TSNE",colour_by=x,point_size=siz,point_alpha=alph,npc_units=0.05,exprs_quant=1)+
    MARKER_THEME+scale_color_viridis(name=x,limits=c(0,3.5), na.value="#FDE825")+
    theme(legend.key.height=unit(bar_h, units="npc"),legend.key.width=unit(bar_w, units="npc"),
          axis.line.y.left=element_line(size=line_size), axis.line.x.bottom=element_line(size=line_size)) 
      })
names(plot_lst) <- markers_to_plot


Bcell_tsne <- ggDIMRED(Bcell_sce, dimred="TSNE",colour_by="CellType", 
                       col_pal=setNames(c("#E9F3A3","grey"),c("B_cells","Unassigned")),
                        ordering=NULL,point_size=siz,point_alpha=0.7)+
                          MARKER_THEME + 
                          theme(legend.position = "none",
                                 axis.line.y.left=element_line(size=line_size), 
                                  axis.line.x.bottom=element_line(size=line_size))+
                           guides(colour=guide_legend(ncol=1,override.aes = list(size=3))) 

Patient_tsne <- ggDIMRED(Bcell_sce, dimred="TSNE",colour_by="PatientID", 
                       col_pal=Colour_values_PatientID,ordering=NULL,point_size=siz,point_alpha=0.7)+
                    MARKER_THEME + theme(legend.position = "none",
                                         axis.line.y.left=element_line(size=line_size), 
                                         axis.line.x.bottom=element_line(size=line_size))+
                          guides(colour=guide_legend(ncol=1,override.aes = list(size=3))) 

CellType_leg <- GetLegend(setNames(c("#E9F3A3","grey"),c("B cells","Unassigned")), "Cell Type",n_row = 2,leg_size =3)
PatientID_leg <- GetLegend(Colour_values_PatientID, "PatientID", n_row = 5,leg_size =3)

pdf(paste0(fig_dir,"FIG2_BcellPlasmaMarkers.pdf"),  width=12, height=5)
cowplot::plot_grid(
              cowplot::plot_grid(CellType_leg,PatientID_leg,ncol=1)+theme(plot.margin = margin(0,-20,0,0)),
              cowplot::plot_grid(Bcell_tsne+theme(plot.margin = margin(0,0,5,25)),
                   plot_lst[[1]], 
                   plot_lst[[2]],
                   Patient_tsne+theme(plot.margin = margin(0,0,5,25)),
                   plot_lst[[3]], 
                   plot_lst[[4]], 
                   nrow=2,
                   labels=c("",markers_to_plot[1:2], "", markers_to_plot[3:4]),
                   label_size=14, label_fontface ="plain",
                   label_x=c(0,0.6,0.6,0,0.53,0.6), label_y=0.95),
                ncol=2, rel_widths  = c(1,4))
dev.off()




##### ----- PLOT FG 2D BAR PLOTS ------------------------------------------

MPE_T_sce <- Fig2_sce[,Fig2_sce$CellType == "T_cells"]

MPE_T_sce$CD4pos <- as.vector(counts(MPE_T_sce[rownames(MPE_T_sce) == "CD4",])>0)
MPE_T_sce$CD8pos <- as.vector(colSums(counts(MPE_T_sce[rownames(MPE_T_sce) %in% c("CD8A", "CD8B"),])>0)>0)

MPE_T_sce$Marker <- MPE_T_sce$CellType
MPE_T_sce$Marker <- replace(MPE_T_sce$Marker,
                            MPE_T_sce$CD4pos & MPE_T_sce$CD8pos,  "CD4+/CD8+")
MPE_T_sce$Marker <- replace(MPE_T_sce$Marker,
                            MPE_T_sce$CD4pos & !(MPE_T_sce$CD8pos),  "CD4+")
MPE_T_sce$Marker <- replace(MPE_T_sce$Marker,
                            !(MPE_T_sce$CD4pos) & MPE_T_sce$CD8pos,  "CD8+")
MPE_T_sce$Marker <- replace(MPE_T_sce$Marker,
                            !(MPE_T_sce$CD4pos) & !(MPE_T_sce$CD8pos),  NA)

plot_dat <- as.data.frame(colData(MPE_T_sce)[,colnames(colData(MPE_T_sce)) %in% c("Marker","PatientID")])

table(!(is.na(plot_dat$Marker)))/nrow(plot_dat)*100
plot_dat <- plot_dat[!(is.na(plot_dat$Marker)),]

plot_dat$Marker <- factor(plot_dat$Marker, levels=rev(c("CD8+","CD4+","CD4+/CD8+")))
plot_dat$PatientID <- factor(plot_dat$PatientID, levels=c("BCB90", "BCB114","BCB21", "BCB20", "BCB139",  "BCB112", "BCB66"))


library(dplyr)

ggbar <- plot_dat %>% 
  count(PatientID, Marker) %>%       
  group_by(PatientID) %>%
  mutate(pct= prop.table(n) * 100) %>%
  ggplot() + aes(PatientID, pct, fill=Marker) + 
  geom_bar(stat="identity")+
  theme_bw()+
  theme(axis.ticks.y=element_blank(),axis.text.y=element_text(size=16),
        axis.title.x=element_text(size=14,margin = unit(c(5, 0, 0, 0), "mm")),
        strip.text=element_text(size=20))+
  scale_fill_manual(values=(c("grey","#CABEE9","#7C7189")))+
  scale_y_continuous(expand=c(0,0))+coord_flip()+ 
  labs(y="T cell Percent (%)",x="")

pdf(paste0(fig_dir,"FIG2_TcellMarkers.pdf"), width=8, height=6)
ggbar
dev.off()



##### ----- PLOT FIG 2E BARPLOT ------------------------------------------


### --- Refine predictions based on CD4/CD8 expression
## CD4+ cells predicted to be CD8-types and
## CD8+ cells predicted to be CD8 Th1/Tfh/CD4 were excluded
NKTcell_sce$CD4pos <- as.vector(counts(NKTcell_sce[rownames(NKTcell_sce) == "CD4",])>0)
NKTcell_sce$CD8pos <- as.vector(colSums(counts(NKTcell_sce[rownames(NKTcell_sce) %in% c("CD8A", "CD8B"),])>0)>0)

NKTcell_sce$ProjecTILs_pred <- replace(NKTcell_sce$ProjecTILs_pred, 
                                       NKTcell_sce$ProjecTILs_pred %in% c("Tfh", "Th1", "CD4_NaiveLike") & NKTcell_sce$CD8pos,
                                        "Excluded")
NKTcell_sce$ProjecTILs_pred <- replace(NKTcell_sce$ProjecTILs_pred, 
                                       grepl("CD8",NKTcell_sce$ProjecTILs_pred) & NKTcell_sce$CD4pos,
                                       "Excluded")

### --- Prep plot data
library(dplyr)
count_pct <- function(df) {
  return(
    df %>%
      tally %>% 
      mutate(n_pct = 100*n/sum(n))
  )
}

plot_dat <- as.data.frame(colData(NKTcell_sce)[,colnames(colData(NKTcell_sce)) 
                                           %in% c("ProjecTILs_pred", "PatientID", "CellID")])
plot_dat <- plot_dat %>% group_by(PatientID,ProjecTILs_pred) %>% count_pct
plot_dat$PatientID <- factor(plot_dat$PatientID, levels=c("BCB66","BCB112","BCB139","BCB20","BCB21","BCB114","BCB90"))
plot_dat$ProjecTILs_pred <- factor(plot_dat$ProjecTILs_pred,
                                   levels=c("CD8_NaiveLike", "CD8_EarlyActiv", "CD8_EffectorMemory",
                                            "Th1","CD4_NaiveLike","Tfh", "NK_cells", "Excluded"))

pdf(paste0(fig_dir, "FIG2_NKTcellSubtypes_byPatient.pdf"),width=10, height=6)
ggplot(data=plot_dat, aes(x=ProjecTILs_pred, y=n_pct, fill=PatientID))+
  geom_bar(position="dodge", stat="identity",width=0.9,colour="black",size=0.1)+
  labs(y="% of Patient Cells", x="")+
  theme(axis.text.x = element_text(angle=45, hjust=1, size=12), axis.ticks.x = element_blank())+
  scale_fill_manual(values=Colour_values_PatientID[c("BCB66","BCB112","BCB139","BCB20","BCB21","BCB114","BCB90")])+
  scale_y_continuous(expand = c(0, 0), limits = c(0,85))+scale_x_discrete(expand=c(0,0))
dev.off()



