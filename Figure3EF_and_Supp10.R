
##### ----- Get directories


setwd("./MPE_paper_analysis/")
data_dir <- paste0(getwd(),"/data/")
matrix_dir <- paste0(data_dir,"raw_matrices/")
pal_etal_dir <- "~/data/_datasets/sc/Smyth_HumanBrCa/"
pal_save_dir <- paste0(data_dir,"pal_etal_data/")
fig_dir <- paste0(getwd(),"/figures/")

source(paste0(getwd(),"/scripts/helper_functions.R"))


### COLOURS ------------
Colour_values_subtype <- setNames(c("#CABEE9","#2F3D70","#BC8E7D"),
                                  c("ER","TN", "TNB1"))


Colour_values_PatientID <- c("#077187", "#6AAAB7", "#8E2043", 
                             "#FEA090","#FECFC7", "#3E5496", "#0A9086")
names(Colour_values_PatientID) <- c("BCB66","BCB66_E", "BCB20","BCB20_E", 
                                    "BCB21_E", "BCB112",  "BCB139")

Colour_values_AIMS <- c("#CABEE9", "#7C7189", "#FAE093", "#D04E59", "#2F3D70")
names(Colour_values_AIMS) <- c( "LumA", "LumB", "Her2","Normal", "Basal")


### --- LOAD & PROCESS DATA ----------------------------------------------------


### --- Get malignant data
load(file=paste0(data_dir, "Fig3_sce.Rdata"))
Fig3_sce <- Fig3_sce[,Fig3_sce$CellType_Fig3 == "Malignant"]
Fig3_sce <- Fig3_sce[,!(Fig3_sce$PatientID == "BCB21")]



### --- Merge MPE data to select MKI67 clusters 

library(scfunc)
malignantMPE_sce <- SeuratIntegrate(Fig3_sce[as.vector(rowSums(counts(Fig3_sce)>0)>0), ],
                                            split_col="Patient", n_features = 3000)

malignantMPE_sce <- runTSNE(malignantMPE_sce)
malignantMPE_sce <- runUMAP(malignantMPE_sce)
malignantMPE_sce <- runPCA(malignantMPE_sce)
malignantMPE_sce <- runTSNE(malignantMPE_sce,dimred="PCA",name="PCA_TSNE")

g <- buildSNNGraph(malignantMPE_sce,k=50, assay.type = "logcounts")
clust <- igraph::cluster_walktrap(g)$membership
malignantMPE_sce$Clust_Knn50 <- as.character(as.vector(clust))

#save(malignantMPE_sce, file=paste0(data_dir, "merged_malignantMPE_sce.Rdata"))



### --- FIGURE 3E --------------------------------------------------------------


### --- Select MPE clusters for pseudobulk
library(RColorBrewer)
library(pals)
MKI67_cols <- setNames(c(stepped()[13:16],stepped2()[17:20]),
                       c("3","5","7", "8",
                         "1", "2", "4", "6"))
MKI67_cols <- MKI67_cols[as.character(1:8)]

clust_cols <- setNames(brewer.pal(n = 8, name = "Dark2"),
                       as.character(1:8))
A <- ggDIMRED(malignantMPE_sce , dimred="PCA_TSNE",colour_by = "MKI67",point_size =4,npc_units=0.05)+
                theme(legend.key.height=unit(0.17, units="npc"),
                      legend.key.width=unit(0.03, units="npc"),
                      axis.ticks = element_blank(), panel.border=element_blank(),
                      axis.title = element_blank(), axis.text=element_blank(), 
                      axis.line = element_blank(),legend.title=element_blank(),
                      legend.position="left", legend.text = element_text(size=18)) 
B <- ggDIMRED(malignantMPE_sce , dimred="PCA_TSNE",colour_by = "Clust_Knn50",
              col_pal = MKI67_cols,point_size =4)+theme_blank()

plot_dat <- as.data.frame(colData(malignantMPE_sce)[colnames(colData(malignantMPE_sce)) %in% c("Clust_Knn50", "CellID")])
plot_dat$MKI67 <- as.vector(logcounts(malignantMPE_sce[rownames(malignantMPE_sce)=="MKI67",]))
plot_dat$Clust_Knn50 <- as.character(plot_dat$Clust_Knn50)
C <- ggplot(data=plot_dat, aes(x=Clust_Knn50, y= MKI67, fill=Clust_Knn50))+
  geom_boxplot(lwd=0.2) + scale_fill_manual(values=MKI67_cols)+ labs(y="MKI67 Expression") +
  guides(fill = guide_legend(ncol=2,title="Cluster")) +
  theme_bw() + 
  theme(axis.ticks.x=element_blank(), axis.text.x=element_blank(),
        axis.title.x=element_blank(), 
        axis.title.y=element_text(size=18),axis.text.y=element_text(size=16),
        legend.position = "none")+scale_y_continuous(expand=c(0,0))

gg_leg <- GetLegend(MKI67_cols, "KNN Cluster", n_row=1,leg_size =6)

pdf(paste0(fig_dir,"FIG3E_Malignant_MKI67clust_boxes.pdf"),width=7,height=3)
cowplot::plot_grid(C,gg_leg, ncol=1, rel_heights = c(1,0.3))
dev.off()

reso <- 350
length <- 3.25*reso/72
png(paste0(fig_dir,"FIG3E_MalignantTSNE_MKI67.png"), units="in", width=length, height=length*0.9, res=reso)
A
dev.off()

png(paste0(fig_dir,"FIG3E_MalignantTSNE_Clust.png"), units="in", width=length, height=length*0.9, res=reso)
B
dev.off()



### --- SUPP FIGURE 10 --------------------------------------------------------------

patient_tsne <- cowplot::plot_grid(ggDIMRED(malignantMPE_sce, dimred="PCA_TSNE",colour_by = "PatientID",
                                            col_pal = Colour_values_PatientID,point_size =0.8)+theme_blank(),
                                   GetLegend(Colour_values_PatientID, "PatientID",n_row=2,leg_size =1,rl=0.6),
                                   ncol=1, rel_heights = c(7,1))


pdf(paste0(fig_dir,"SUPPFIG10_MKI67clusts.pdf"),width=7,height=8)
patient_tsne
dev.off()







### --- DIFFERENTIAL EXPRESSION ANALYSIS --------------------------------------


### --- Pseudo-bulk within MPE data 

BuildDGE <- function(sce_x, group_str="Subtype",
                     libsize_filt = 1e+6, cpm_filt=0.1){
  summed_x <- aggregateAcrossCells(sce_x, id=sce_x$PB_clust)
  
  dge_x <- DGEList(counts(summed_x),
                     group=summed_x[[group_str]],
                     samples=as.data.frame(colData(summed_x)))
  
  dge_x$samples$MKI67_clust <- ifelse(grepl("MKI67_clust",dge_x$samples$PB_clust), 
                                        "MKI67_clust", "other_clust")
  
  ## Discard samples with small library sizes
  print("Retained samples:")
  print(table(dge_x$samples$lib.size > libsize_filt))
  dge_x <- dge_x[,dge_x$samples$lib.size > libsize_filt]
  
  ## Filter genes
  keep = rowMeans(edgeR::cpm(y=dge_x, log = TRUE) >= cpm_filt) >= 0.1
  dge_x <- dge_x[keep, ]
  print("Retained genes:")
  print(table(keep))
  
  ## Normalise by library size
  dge_x <- calcNormFactors(dge_x)  
  dge_x$logCPM <- edgeR::cpm(dge_x, log=TRUE, prior.count = 1) 
  #ggRLE(MPE_dge$logCPM, MPE_dge$samples, "MKI67_clust",isLog = TRUE)
  return(dge_x)
}

plot3PCA <- function(dge_x){
  require(ggRNA)
  patient_cols <- Colour_values_PatientID[names(Colour_values_PatientID) %in% dge_x$samples$PatientID]
  A <- MDS2gg(dge_x, "PatientID",col_pal=patient_cols,shape_str = "MKI67_clust")
  B <- ggPCA(dge_x$logCPM, dge_x$samples,
             "PatientID",  col_pal=patient_cols,shape_str = "MKI67_clust")
  C <- ggPCA(dge_x$batch_corrected, dge_x$samples,
             "PatientID",  col_pal=patient_cols, shape_str = "MKI67_clust")
  
  plot_leg <- cowplot::get_legend(A+theme(legend.direction = "horizontal",
                                          legend.box="horizontal")+
                                    guides(colour=guide_legend(nrow=1),
                                           shape=guide_legend(nrow=1)))
  plot_lst <- list(A+theme(legend.position="none"),
                   B+theme(legend.position="none"),
                   C+theme(legend.position="none"))
  return(cowplot::plot_grid(cowplot::plot_grid(
    plotlist =plot_lst, nrow=1, 
    labels = c("MDS", "PCA", "PCA - corrected"),
    label_fontface = "plain"),
    plot_leg, ncol=1, rel_heights = c(6,1)))
  
}

library(ggRNA)

### --- Get pseudobulk clusts
# load(file=paste0(data_dir, "merged_malignantMPE_sce.Rdata"))
malignantMPE_sce$MKI67_clust <- ifelse(malignantMPE_sce$Clust_Knn50 %in% c("3","5","7","8"), 
                                       "MKI67_clust","other_clust")
malignantMPE_sce$PB_clust <- paste0(malignantMPE_sce$PatientID, "_",malignantMPE_sce$MKI67_clust)
PB_clust_dict <- setNames(malignantMPE_sce$PB_clust, malignantMPE_sce$CellID)


### --- Get non-integrated malignant count data
load(file=paste0(data_dir, "Fig3_sce.Rdata"))
Fig3_sce <- Fig3_sce[,Fig3_sce$CellType_Fig3 == "Malignant"]
Fig3_sce <- Fig3_sce[,!(Fig3_sce$PatientID == "BCB21")]
Fig3_sce$Subtype <- ifelse(grepl("139|20|21",Fig3_sce$PatientID), "ER", "TN")
Fig3_sce$PB_clust <- as.vector(PB_clust_dict[Fig3_sce$CellID])



##### ----- TN patients
TN_dge <- BuildDGE(Fig3_sce[,!(Fig3_sce$Subtype =="ER")], group_str="Subtype",
                    libsize_filt = 1e+6, cpm_filt=0.1)

### Design matrix
patient_factor <- factor(TN_dge$samples$PatientID)
MKI67_factor <- factor(TN_dge$samples$MKI67_clust, levels=c("other_clust", "MKI67_clust"))
design_mat_TN <- model.matrix(~0 + patient_factor + MKI67_factor)

TN_dge <- estimateDisp(TN_dge, design_mat_TN, robust = TRUE)

### Correct for patient
design_to_preserve <- model.matrix(~0 + MKI67_factor)
TN_dge$batch_corrected <- limma::removeBatchEffect(cpm(TN_dge, log=TRUE, prior.count = 1),
                                                    design = design_to_preserve,
                                                    batch=patient_factor)

plot3PCA(TN_dge)


##### ----- Luminal patients
ER_dge <- BuildDGE(Fig3_sce[,Fig3_sce$Subtype =="ER"], group_str="Subtype",
                   libsize_filt = 1e+5, cpm_filt=0.1)

### Design matrix
patient_factor <- factor(ER_dge$samples$PatientID)
MKI67_factor <- factor(ER_dge$samples$MKI67_clust, levels=c("other_clust", "MKI67_clust"))
design_mat_ER <- model.matrix(~0 + patient_factor + MKI67_factor)

ER_dge <- estimateDisp(ER_dge, design_mat_ER, robust = TRUE)

### Correct for patient
design_to_preserve <- model.matrix(~0 + MKI67_factor)
ER_dge$batch_corrected <- limma::removeBatchEffect(cpm(ER_dge, log=TRUE, prior.count = 1),
                                                   design = design_to_preserve,
                                                   batch=patient_factor)

plot3PCA(ER_dge)




##### ----- Perform Differential expression

qfit_TN <- glmQLFit(TN_dge, design_mat_TN)
qfit_ER <- glmQLFit(ER_dge, design_mat_ER)


## Test DE 
TN_de_obj <- glmQLFTest(qfit_TN, coef=4) 
summary(decideTestsDGE(TN_de_obj, adjust.method="BH"))

ER_de_obj <- glmQLFTest(qfit_ER, coef=5) 
summary(decideTestsDGE(ER_de_obj, adjust.method="BH"))




### --- KEGG ENRICHMENT: FIGURE 3F ----------------------------------------------

## -- Get gene ids
library(biomaRt)
ensembl <- useEnsembl(biomart='ENSEMBL_MART_ENSEMBL', dataset="hsapiens_gene_ensembl", host="https://oct2018.archive.ensembl.org")
annot <- biomaRt::getBM(attributes=c('ensembl_gene_id','version', 'external_gene_name',  'gene_biotype','entrezgene', 'hgnc_symbol'),
                        mart = ensembl)
annot <- annot[!(annot$hgnc_symbol==""),]

idx <- match(rownames(TN_de_obj), annot$hgnc_symbol)
TN_de_obj <- TN_de_obj[which(!is.na(idx)),]
idx <- idx[!is.na(idx)] 
TN_de_obj$genes <- annot[idx,] 

idx <- match(rownames(ER_de_obj), annot$hgnc_symbol)
ER_de_obj <- ER_de_obj[which(!is.na(idx)),]
idx <- idx[!is.na(idx)] 
ER_de_obj$genes <- annot[idx,] 


## -- KEGG
KEGGEnrich_TN <- kegga(TN_de_obj, species="Hs", 
                    geneid=TN_de_obj$genes$entrezgene, FDR=0.05)
KEGGEnrich_ER <- kegga(ER_de_obj, species="Hs", 
                    geneid=ER_de_obj$genes$entrezgene, FDR=0.05)


##### ----- Check overlap with primary

PAL_TN_KEGG_pathways <- c("DNA replication","Spliceosome","Base excision repair",
                          "Nucleotide excision repair","Mismatch repair","Cellular senescence",
                          "Cell cycle","Homologous recombination","Fanconi anemia pathway",
                          "RNA transport")

PAL_ER_KEGG_pathways <- c("Cell cycle", "DNA replication",
                          "Mismatch repair", "Homologous recombination",
                          "Oxidative phosphorylation",
                          "Cellular senescence", "Metabolic pathways",
                          "Nucleotide excision repair",
                          "Progesterone-mediated oocyte maturation",
                          "Base excision repair")

topTN <- topKEGG(KEGGEnrich_TN, number=20)
TN_KEGG_pathways <- topTN[topTN$P.Up < 0.05,]$Pathway

topER <- topKEGG(KEGGEnrich_ER, number=20)
ER_KEGG_pathways <- topER[topER$P.Up < 0.05,]$Pathway

## Notably the only one missing is Oxidative phosphorylation & Metabolic pathways & "RNA transport"
## 10/13 KEGG pathways discussed in paper

TN_pathways <- topTN[topTN$Pathway %in% unique(c(PAL_TN_KEGG_pathways, PAL_ER_KEGG_pathways)),]
ER_pathways <- topER[topER$Pathway %in% unique(c(PAL_TN_KEGG_pathways, PAL_ER_KEGG_pathways)),]

TN_pathways$Subtype <- rep('TN',nrow(TN_pathways))
ER_pathways$Subtype <- rep('ER',nrow(ER_pathways))

plot_dat <- rbind(TN_pathways, ER_pathways)
plot_dat <- plot_dat[rev(order(plot_dat$P.Up)),]

plot_dat <- data.frame(P=-log10(plot_dat$P.Up),
                    Pathway=plot_dat$Pathway,
                    Subtype=plot_dat$Subtype)

setdiff(TN_pathways$Pathway,ER_pathways$Pathway)
plot_dat[nrow(plot_dat)+1,] = c(0, "Spliceosome", "ER")
plot_dat$P <- as.numeric(plot_dat$P)

plot_dat$Pathway <- factor(plot_dat$Pathway, levels=unique(plot_dat$Pathway))
plot_dat <- plot_dat[order(plot_dat$Pathway),]


plot_dat$PAL <-  ifelse((plot_dat$Pathway %in% PAL_TN_KEGG_pathways) &
                          plot_dat$Pathway %in% PAL_ER_KEGG_pathways, "Both", NA)
plot_dat$PAL  <-  replace(plot_dat$PAL,
                          (plot_dat$Pathway %in% PAL_ER_KEGG_pathways) &
                            is.na(plot_dat$PAL),"ER")
plot_dat$PAL  <-  replace(plot_dat$PAL,
                          (plot_dat$Pathway %in% PAL_TN_KEGG_pathways) &
                            is.na(plot_dat$PAL),"TN")
KEGG_cols <- setNames(c("#2F3D70","#CABEE9","#8AAAD0"),
                      c("TN","ER","Both"))

### --- PLOT FIGURE 3F
A <- ggplot(data=plot_dat[plot_dat$Subtype=="ER",], aes(x=P, y=Pathway, fill=PAL))+
  geom_bar(stat="identity")+theme_bw() + 
  theme(axis.ticks.y=element_blank(), axis.text.y=element_blank(),
        legend.position = "none",plot.margin=margin(0,5,0,0))+
  labs(x="-Log10(P-value)", y="")+scale_x_reverse(expand=c(0,0)) + 
  scale_fill_manual(values=KEGG_cols)

B <- ggplot(data=plot_dat[plot_dat$Subtype=="TN",], aes(x=P, y=Pathway, fill=PAL))+
  geom_bar(stat="identity")+theme_bw() + scale_x_continuous(expand=c(0,0))+
  theme(axis.ticks.y=element_blank(), axis.text.y=element_text(size=13,hjust=0.5),
        legend.position = "none", plot.margin=margin(0,0,0,0))+
  labs(x="-Log10(P-value)", y="") + scale_fill_manual(values=KEGG_cols)


gg_leg <- ggRNA::GetLegend(KEGG_cols, "Enriched in Primary", n_row = 3)

pdf(paste0(fig_dir,"FIG3F_KEGG.pdf"),width=13,height=3)
cowplot::plot_grid(cowplot::plot_grid(A,B,rel_widths =c(1,1.8)),
                   gg_leg,nrow=1,rel_widths = c(6,1))
dev.off()








