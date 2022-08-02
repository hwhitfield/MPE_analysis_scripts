# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
#   PRODUCE FIGURE 1 -- All cells non-enriched samples
#   
#   To recreate Figure 1
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


library(ggRNA)

##### ----- Get directories
#wd="/stornext/Home/data/allstaff/w/whitfield.h/PhD_Project_MPE/MPE_paper_analysis/"

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

Colour_values <- c("#9e5476", "#d18b79",  "#70a18f",
                   "#e7ebbc",
                   "#8bbde6", "#6074ab", 
                   "grey")

names(Colour_values) <- c("Malignant", "Mesothelial",  "Myeloid",  
                          "B_cells",
                          "NK_cells", "T_cells",
                          "Unassigned")



##### ----- Load data ----------------------------------------------------------

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


### --- Combine
Fig1_sce <- sceBind(non_enriched_lst, colDat_to_keep = colnames(ALL_COLDATA),rowsumThresh=0)
Fig1_sce$PatientID <- Fig1_sce$batch
#save(Fig1_sce, file=paste0(data_dir, "Fig1_sce.Rdata"))





##### ----- PLOT FIG 1A -----------------------------------

# load(file=paste0(data_dir, "Fig1_sce.Rdata"))
Fig1_sce$CellType_Fig1 <- replace(Fig1_sce$CellType,
                                  Fig1_sce$CellType %in% c("DC","Macrophage"),
                                  "Myeloid")

A1 <- ggRNA::ggDIMRED(Fig1_sce, dimred="TSNE", colour_by = "CellType_Fig1",
                      point_size=4, point_alpha =0.8, col_pal=Colour_values)
A2 <- ggRNA::ggDIMRED(Fig1_sce, dimred="TSNE", colour_by = "PatientID",point_alpha =0.7, col_pal=Colour_values_PatientID)


reso <- 500
length <- 3.25*reso/72
png(paste0(fig_dir,"Fig1A_byCellType.png"), units="in", width=length, height=length*0.9, res=reso)
A1+theme_blank()
dev.off()

pdf(paste0(fig_dir,"Fig1A_byCellType.pdf"),  width=7, height=6)
ggRNA::ggDIMRED(Fig1_sce, dimred="TSNE", colour_by = "CellType_Fig1",
                point_size=1, point_alpha =1, col_pal=Colour_values)+theme_blank()
dev.off()

patient_cols <- GetLegend(Colour_values_PatientID, "PatientID",n_row=2)
pdf(paste0(fig_dir,"Fig1A_byPatientID.pdf"),  width=7, height=7)
cowplot::plot_grid(A2+theme_blank(), patient_cols, ncol=1, rel_heights = c(7,1))
dev.off()


##### ----- PLOT FIG 1A PROPORTIONAL BAR ----------------------------------------


CellType_Props <- ggRNA::GetProps(Fig1_sce, label="CellType_Fig1")
CellType_Props$CellType <- rownames(CellType_Props)
CellType_Props$Sample <- rep("AllPatients", length(rownames(CellType_Props)))
CellType_Props <- CellType_Props[,!(colnames(CellType_Props) == "Freq")]
CellType_Props <- CellType_Props[order(-CellType_Props$Proportion),]
CellType_Props$CellType <- factor(CellType_Props$CellType, levels=rev(unique(CellType_Props$CellType)))

PropBar <- ggplot(CellType_Props, aes(x = Sample, y=Proportion, fill=CellType)) + geom_bar(position="dodge", stat="identity") + 
  geom_bar(stat = "identity")+
  scale_fill_manual(values=Colour_values, breaks=CellType_Props$CellType)+ 
  coord_flip()+theme_bw()+PlainBar_theme+
  scale_y_continuous(expand = c(0, 0), limits = c(0,100))+
  theme(axis.text.y=element_blank(), legend.position = "none")


pdf(paste0(fig_dir, "Fig1A_prop.pdf", ""), width=6, height=1)
PropBar
dev.off()




##### ----- PLOT FIG 1B PROPORTIONAL BAR ----------------------------------------

### Get proportions
CellType_Props <- GetProps_perPatient(Fig1_sce, label="CellType_Fig1", "PatientID")
CellType_Props$CellType <- rownames(CellType_Props)
CellType_Props <- reshape2::melt(CellType_Props, id.vars="CellType", value.name="Proportion", variable.name="PatientID")

## Add cell type
CellType_Props <- CellType_Props[order(-CellType_Props$Proportion),]
CellType_Props$CellType <- factor(CellType_Props$CellType, levels=rev(unique(CellType_Props$CellType)))

## Add patientID
SampleCols <- c("BCB90", "BCB114","BCB21", "BCB20", "BCB139",  "BCB112", "BCB66")
CellType_Props$PatientID <- factor(CellType_Props$PatientID, levels=SampleCols)

## Add cell number info
SampleNumbs_dict <- unlist(lapply(SampleCols, function(x)  {ncol(Fig1_sce[,Fig1_sce$PatientID == x])}))
names(SampleNumbs_dict) <- SampleCols
CellType_Props$CellNumb <- as.vector(SampleNumbs_dict[as.vector(CellType_Props$PatientID)])


### Plot  
PropBar_pp <- ggplot(CellType_Props, aes(x = PatientID, y = Proportion, fill=CellType)) +  
  geom_bar(stat = "identity", position = "stack", width=0.95) +
  scale_fill_manual(values=Colour_values, breaks=CellType_Props$CellType)+theme_bw()+  coord_flip()+
  scale_y_discrete(limits =CellType_Props$CellType)+PlainBar_theme+
  theme(strip.background = element_blank(), strip.text = element_blank(), legend.position = "none")

## Get cell numbs 
CellType_Numb <- CellType_Props[!(duplicated(CellType_Props$PatientID)),colnames(CellType_Props) %in% c("PatientID", "CellNumb")]
CellType_Numb$CellNumb <- as.numeric(CellType_Numb$CellNumb)

cellNumb_bar <- ggplot(data=CellType_Props, aes(x=PatientID, y=CellNumb)) + 
  geom_bar(stat="identity", position="dodge", width=0.93,  fill="dimgrey")+
  coord_flip()+theme_bw()+PlainBar_theme+
  theme(axis.text.x=element_text(size = rel(1)*1.2, hjust=0.8),axis.ticks.x=element_line(colour="black"),
        axis.text.y=element_blank(),axis.line.x=element_line(colour="black",size=0.2))+
  scale_y_continuous(labels=c("0","4k", "8k", "12k"), limits=c(0, 13000),
                     breaks = c(0, 4000, 8000, 12000),expand=c(0,0))

pdf(paste0(fig_dir, "FIG1B_PropBar.pdf"), width=8, height=6)
grid.arrange(grobs=list(PropBar_pp+theme(plot.margin=unit(c(0,0,0.5,0), "cm")), 
                        cellNumb_bar+theme(plot.margin=unit(c(0,0,0,0), "cm"))),
             ncol=2, nrow=1, 
             widths = c(5,1))
dev.off()



##### ----- PLOT FIG 1C CORRELATION HEATMAP ----------------------------------------

### --- Calculate TPM for Pseudobulked data


tmp_sce <- Fig1_sce[!(grepl("RPL|MT|RP", rownames(Fig1_sce))),]
tmp_sce$PseudobulkClust <- paste0(tmp_sce$PatientID,"_",tmp_sce$CellType_Fig1)
tmp_sce <- tmp_sce[,!(tmp_sce$CellType_Fig1 == "Unassigned")]
Fig1_dge <- GetPseudoBulk(tmp_sce, "PseudobulkClust", 
              min_cell_numb=80, cpm_thresh=0.5,sample_thresh=0.1)

tpm_dat <- scater::calculateTPM(Fig1_dge$counts)

pdf(paste0(fig_dir, "FIG1B_corrHeatmap.pdf"), width=8, height=6)
corrHeatmap(tpm_dat,
            setNames(sub(".*?_", "", colnames(Fig1_dge)),colnames(Fig1_dge)),
            col_pal=Colour_values, annot_str="",scale="none")
dev.off()







##### ----- PLOT FIG 1D GENE EXPRESSION HEATMAP ----------------------------------------

## Get genes with SingleR
library(SingleR)
tmp_sce <- Fig1_sce[!(grepl("RPL|MT|RP", rownames(Fig1_sce))),]
tmp_sce$PseudobulkClust <- paste0(tmp_sce$PatientID,"_",tmp_sce$CellType_Fig1)
tmp_sce <- tmp_sce[,!(tmp_sce$CellType_Fig1 == "Unassigned")]
sce_trained <- SingleR::trainSingleR(ref=tmp_sce,labels=tmp_sce$CellType_Fig1, assay.type = "logcounts", de.method="wilcox")

jList <- list()
for (jLabel in unique(names(sce_trained$search$extra))){
  jList[[jLabel]] <- unique(unlist(sce_trained$search$extra[[jLabel]]))
}

## Get celltype-unique genes 
Heatmap_genes <- names(table(unlist(jList))[table(unlist(jList))<2])
Heatmap_genes <- c(Heatmap_genes, 
                   intersect(jList$Mesothelial, jList$Malignant),
                   "DES","MSLN","WT1", "CDH2","VIM", "COL1A1", "COL6A2", "KRT19","CDH1","CLDN4", "CLDN7")
## Drop some ubiquitous or uninteresting genes
drop_genes <- c("ZFP36L2","HLA-A", "HLA-C","HLA-E","CD52","BTG1","CD48","CD37","TXNIP","HCST",
                "CDC42","FTL","YBX1","TPT1","GNB1","EIF1","FAU","GNB2L1","NEAT1","TOMM7","JUN")
Heatmap_genes <- Heatmap_genes[!(Heatmap_genes %in% drop_genes)]

### --- Pseudobulk sample ordering in heatmap, subgroups with more than 80 cells
ordering <- c( unique(tmp_sce$PseudobulkClust)[grepl("Mesothelial",unique(tmp_sce$PseudobulkClust))],
               unique(tmp_sce$PseudobulkClust)[grepl("Malignant",unique(tmp_sce$PseudobulkClust))],
               unique(tmp_sce$PseudobulkClust)[grepl("Myeloid",unique(tmp_sce$PseudobulkClust))],
               unique(tmp_sce$PseudobulkClust)[grepl("B_cells",unique(tmp_sce$PseudobulkClust))],
               unique(tmp_sce$PseudobulkClust)[grepl("T_cells",unique(tmp_sce$PseudobulkClust))],
               unique(tmp_sce$PseudobulkClust)[grepl("NK_cells",unique(tmp_sce$PseudobulkClust))])
ordering <- ordering[ordering %in% names(table(tmp_sce$PseudobulkClust)[table(tmp_sce$PseudobulkClust)>80])]

ggheatmap <- ggRNA::scHeatmap(tmp_sce, "PseudobulkClust", unique(Heatmap_genes),
                              exprs_val="TPM",log=TRUE,scale_str = "genes",
                              cluster_genes =TRUE, cluster_samples=FALSE,
                              lower_thresh=0.05, upper_thresh=0.999,rl=0.7,
                              groups_to_plot=rev(ordering), label_side="bottom")

#png(save_path,  width=2000, height=2000,  res=200)
pdf(paste0(fig_dir, "FIG1B_geneExpressionHeatmap.pdf"), width=6, height=12)
ggheatmap+theme(axis.text.y=element_text(margin=margin(r=0)))+
          guides(fill=guide_colorbar(title="Average logTPM Expression",
                                     title.position ="bottom",direction="vertical",
                                     barwidth =unit(0.04, units="npc"),
                                     barheight = unit(0.2, units="npc")))
dev.off()
