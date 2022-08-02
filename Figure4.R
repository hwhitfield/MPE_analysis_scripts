# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
#   PRODUCE FIGURE 4 -- Malignant & Mesothelial cells
#   
#   To recreate Figure 4
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


##### ----- Get directories
data_dir <- paste0(getwd(),"/data/")
matrix_dir <- paste0(data_dir,"raw_matrices/")
fig_dir <- paste0(getwd(),"/figures/")

source(paste0(getwd(),"/scripts/helper_functions.R"))


requiredPackages <- c("scales", "dplyr", "pals","GSEABase",
                      "VISION","gplots","org.Hs.eg.db","ggRNA")

for (pkg in requiredPackages){
  suppressWarnings(suppressMessages(library(pkg, character.only = T)))
}




### COLOURS ------------

Colour_values_PatientID <- c("#077187", "#6AAAB7", "#8E2043", "#BC7A8F", 
                             "#FEA090","#FECFC7", "#3E5496", "#0A9086", 
                             "#E0607E", "#E1D164")
names(Colour_values_PatientID) <- c("BCB66","BCB66_E", "BCB20","BCB20_E", "BCB21", 
                                    "BCB21_E", "BCB112",  "BCB139",  "BCB90",  "BCB114")

Colour_values <- c("#9e5476", "#d18b79")
names(Colour_values) <- c("Malignant", "Mesothelial")

cols <- c("#FFC0CB", "#94376C")
names(cols) <- c("FGFR2-", "FGFR2+")


##### ----- Load data ----------------------------------------------------------



### --- Load data from figure 3

load(file=paste0(data_dir, "Fig3_sce.Rdata"))



FGF_sigs <- c("PID_FGF_PATHWAY","REACTOME_PI_3K_CASCADE_FGFR2",
              "REACTOME_DOWNSTREAM_SIGNALING_OF_ACTIVATED_FGFR2")

facet_labels <- c("PID FGF Pathway", "Downstream Signaling of FGFR2", "PI-3K Cascade: FGFR2")
names(facet_labels) <- c("PID_FGF_PATHWAY","REACTOME_DOWNSTREAM_SIGNALING_OF_ACTIVATED_FGFR2","REACTOME_PI_3K_CASCADE_FGFR2")



##### ----- Pseudo-bulk celltypes within each patient --------------------------



### --- Only patients with enough of both celltype
Fig3_sce$PseudobulkClust <- paste0(Fig3_sce$PatientID,"_",Fig3_sce$CellType_Fig3)
tmp_sce <- Fig3_sce[,Fig3_sce$PatientID %in% c("BCB66", "BCB139", "BCB112", "BCB21_E")]

### --- Pseudobulk (includes gene filtering & calcNormFactors
Fig4_dge <- GetPseudoBulk(tmp_sce, "PseudobulkClust", cpm_thresh=0.5,sample_thresh=0.1)
Fig4_dge$samples$Patient <- sub('_[^_]*$', '', rownames(Fig4_dge$samples))
Fig4_dge$samples$CellType <- sub(".*_", "", rownames(Fig4_dge$samples))

## Check normalisation
#ggRLE(Fig4_dge$counts, Fig4_dge$samples, "CellType",isLog = FALSE)
#ggRLE(Fig4_dge$logCPM, Fig4_dge$samples, "Patient",isLog = TRUE)

### --- FIGURE 4A PCA

## Check for transcriptomic signal of celltype differences
pdf(paste0(fig_dir, "Fig4PCA.pdf"), width=5, height=3)
ggPCA(Fig4_dge$logCPM,Fig4_dge$samples,col_str="CellType", col_pal=Colour_values,
      shape_str = "Patient",isLog=TRUE,point_size =6)
dev.off()



##### ----- Differential expression analysis ------------------------------------


### --- Build design matrix
Patient <- as.vector(Fig4_dge$samples$Patient)
CellType <- as.vector(Fig4_dge$samples$CellType)
design_mat <- model.matrix( ~0 + Patient + CellType)

### --- Estimate dispersions
Fig4_dge <- edgeR::estimateDisp(Fig4_dge, design_mat)
Fig4_qfit <- edgeR::glmQLFit(Fig4_dge, design_mat)


### --- Test differential expression
Fig4_qtest <- edgeR::glmQLFTest(Fig4_qfit, coef=5) 
summary(decideTestsDGE(Fig4_qtest, adjust.method="BH"))


### --- Output results
Fig4_toptags <- topTags(Fig4_qtest, n=Inf, p.value=0.05, adjust.method = "BH", sort.by = "PValue")
Fig4_toptags <- Fig4_toptags[Fig4_toptags$table$FDR < 0.05 ,]
write.table(Fig4_toptags$table, paste0(fig_dir, "SuppTable3.csv", ""), sep=",", row.names = FALSE)




##### ----- GENE SET ENRICHMENT ANALYSES: SUPP TABLE 5 -----------------------

### --- Get ENTREZ gene id for limma
entrez_geneid <- mapIds(org.Hs.eg.db, keys=rownames(Fig4_qtest), 
                 column="ENTREZID", 
                 keytype="SYMBOL")

### --- GENE ONTOLOGY
GOEnrich <- limma::goana(Fig4_qtest, species="Hs", geneid=entrez_geneid, FDR=0.05)
Table_GO = limma:topGO(GOEnrich, n=30, ontology=c("BP", "CC", "MF")) #  truncate.term=30,
write.table(Table_GO, paste0(fig_dir, "SuppTable5_1.csv"), sep=",", row.names = TRUE)

### --- KEGG
KEGGEnrich <- limma::kegga(Fig4_qtest, species="Hs", entrez_geneid)
Table_KEGG =  limma:topKEGG(KEGGEnrich, n=30)
write.table(Table_KEGG, paste0(fig_dir, "SuppTable5_2.csv"), sep=",", row.names = TRUE)

### --- Hallmarks

Hs.hallmark <- readRDS(url("http://bioinf.wehi.edu.au/MSigDB/v7.1/Hs.h.all.v7.1.entrez.rds"))
idx_cam <- ids2indices(Hs.hallmark,id=entrez_geneid)

cam <- limma::camera(Fig4_dge$counts, index=idx_cam, design=design_mat, contrast=5) 
cam = cam[order(cam$PValue),]
options(digits=2)
Table_HM <- head(cam,30)
write.table(Table_HM, paste0(fig_dir, "SuppTable5_3.csv"), sep=",", row.names = TRUE)




##### ----- FIG 4B UPPER HEATMAP, DE GENES ------------------------------------


### --- Get top up- and down- regulated genes
top_tags <- topTags(Fig4_qtest, n=Inf, p.value=0.05, adjust.method = "BH", sort.by = "PValue")
top_tags <- top_tags[top_tags$table$FDR < 0.05 ,]
top_tags <- top_tags[!(grepl("RPL|MT|RP", rownames(top_tags))),]

TopGenes_DN <- rownames(top_tags[top_tags$table$logFC < 0 ,])[1:20]
TopGenes_UP <- rownames(top_tags[top_tags$table$logFC > 0 ,])[1:20]

### --- Remove patient effect for heatmap
design_to_preserve <- model.matrix( ~0 + CellType)
Fig4_dge$patient_corrected <- removeBatchEffect(edgeR::cpm(Fig4_dge, log=TRUE, prior.count = 1), 
                                                design = design_to_preserve, batch=Patient)

topgenes_log_expression <- Fig4_dge$patient_corrected[rownames(Fig4_dge$patient_corrected) %in% c(TopGenes_DN, TopGenes_UP),]
topgenes_log_expression <- t(scale(t(topgenes_log_expression)))


par(mfrow=c(1,1))
pdf(paste0(fig_dir, "FIG4B_DEheatmap.pdf", ""),width = 6, height = 12)
heatmap_dat <- heatmap.2(topgenes_log_expression, 
                         labRow=  rownames(topgenes_log_expression),  
                         labCol=colnames(topgenes_log_expression),
                         main = "",
                         cexRow=1,margin=c(1,1),
                         trace="none", 
                         lhei = c(1,7,1),lwid = c(5,1), lmat=rbind( c(3, 0), c(1,2),c(4,0)),
                         dendrogram = "column",Rowv=TRUE,
                         col=viridis,
                         scale = "none",
                         key.xlab="Scaled logCPM", keysize=2,key.title=NA)  
dev.off()


## Get column ordering for lower heatmap
heatmap_colnames <- colnames(topgenes_log_expression)[heatmap_dat$colInd]




##### ----- Calculate gene signature VISION scores -----------------------------


MakeVisionObj <- function(gene_vec, sig_str, row_genes){
  require(VISION)
  gene_vec <- gene_vec[!is.na(gene_vec)]
  gene_vec <- unique(gene_vec)
  gene_vec <- gene_vec[gene_vec %in% row_genes]
  
  genes <- rep(1,length(gene_vec))
  names(genes) <- gene_vec
  return(VISION::createGeneSignature(name = sig_str, sigData = genes))
}


### --- Get gene signatures

## Tirosh et al.
## Download & save as csv : https://science.sciencemag.org/content/352/6282/189
tirosh_dat <-  read.csv2(paste0(data_dir, "aad0501_table_s3.csv"), sep=",", 
                         skip=4, header=TRUE, na.strings=c(""))
CAF_genes <- gsub("'", "", as.vector(tirosh_dat$CAFs))

## Tan et al.
## Downloaded from: https://www.embopress.org/doi/full/10.15252/emmm.201404208
tan_emt <- read.csv(paste0(data_dir,"Thiery_EMTsignature_both_tumour_cellLine.txt"),sep='\t')
tan_cellLine_genes  <- as.character(tan_emt$HGNC.symbol)[tan_emt$epiMes_cellLine %in% 'epi']
tan_tumour_genes  <- as.character(tan_emt$HGNC.symbol)[tan_emt$epiMes_tumor %in% 'epi']


Tirosh_obj <- MakeVisionObj(CAF_genes, "TIROSH_CAF", rownames(Fig3_sce))
Tan_CL_obj <- MakeVisionObj(tan_cellLine_genes, "TAN_EPI_CL", rownames(Fig3_sce))
Tan_T_obj <- MakeVisionObj(tan_tumour_genes, "TAN_EPI_TUM", rownames(Fig3_sce))

## C2 MSigDB
## Download from: http://www.gsea-msigdb.org/gsea/msigdb/collections.jsp
conn <- file(paste0(data_dir, "c2.all.v7.5.1.symbols.gmt"),open="r")
MSigDB <- GSEABase::getGmt(conn)



MSigDB_lst <- c()
for (iSig in FGF_sigs){
  iGeneIDs <- geneIds(MSigDB[[iSig]])
  genes_vec <- rep(1,length(iGeneIDs))
  names(genes_vec) <- iGeneIDs
  SigObject <- createGeneSignature(name = iSig, sigData = genes_vec)
  MSigDB_lst <- c(MSigDB_lst, createGeneSignature(name = iSig, sigData = genes_vec))
}


### --- Run Vision 

# Use normalised expression data, but not log-transformed
Fig3_sce <- scater::logNormCounts(Fig3_sce, transform="none", exprs_values = "counts")

VISION_Obj <- Vision(data =  normcounts(Fig3_sce),pool=FALSE,
                     signatures = c(Tirosh_obj, Tan_CL_obj,Tan_T_obj, MSigDB_lst))
options(mc.cores = 1)
VISION_RESULTS <- analyze(VISION_Obj)


### --- Append scores

for (iSig in colnames(getSignatureScores(VISION_RESULTS))){
  Fig3_sce[[iSig]] <- as.vector(getSignatureScores(VISION_RESULTS)[, iSig])
}

#save(Fig3_sce, file=paste0(data_dir, "Fig3_sce.Rdata"))






##### ----- FIG 4B LOWER HEATMAP -----------------------------------------------



## Subset to cells used in differential expression analysis
coldata <- as.data.frame(colData(Fig3_sce[,Fig3_sce$PatientID %in% c("BCB66", "BCB139", "BCB112", "BCB21_E")]))
coldata <- coldata[, colnames(coldata) %in% c("PseudobulkClust","TIROSH_CAF","TAN_EPI_CL","TAN_EPI_TUM")]

plot_df <- reshape2::melt(coldata, id.vars="PseudobulkClust", variable.name="GeneSet", value.name="Score")
plot_df$Group <- paste0(plot_df$PseudobulkClust, "_", plot_df$GeneSet)
plot_df_mean <- aggregate(Score ~ Group, data=plot_df, mean)
plot_df_mean$CellType <- trimws(regmatches(plot_df_mean$Group,regexpr("_M\\s*(.*?)\\s*_",plot_df_mean$Group)), whitespace = "_")
plot_df_mean$GeneSet <- sub(".*?_T", "T", plot_df_mean$Group)

plot_df_mean$PB <- sub('_T.*', "", plot_df_mean$Group)
plot_df_mean$PB <- factor(plot_df_mean$PB, levels=heatmap_colnames)


score_names_dict <-  c("Tan Epi (CellLine)", "Tan Epi (Tumour)","Tirosh CAF")
names(score_names_dict) <-  c("TAN_EPI_CL", "TAN_EPI_TUM", "TIROSH_CAF")
rl=1
colour_vec <- ocean.matter(200)
pp_box <- ggplot(plot_df_mean , aes(x=GeneSet , y=PB, fill=Score))+geom_tile() + 
  ylab("Pseudo-bulk Cluster")+labs(fill = "Average Score")+theme_bw()+
  theme(axis.text.x=element_blank(),
        axis.title.x = element_blank(), 
        axis.ticks=element_blank(),
        axis.text.y=element_text(size=rel(rl)*1.3),
        axis.title.y = element_text(size = rel(rl)*1.5, margin = margin(t = 0, r = 0, b = 0, l = 0)),
        legend.position="bottom", 
        legend.title = element_text(size = rel(rl)*1.2),
        legend.text = element_text(size = rel(rl)*1.2, vjust=1),
        panel.grid.minor = element_blank(),    panel.border = element_rect(colour = NA),
        panel.grid.major = element_blank())+ 
  scale_fill_gradientn(colours=colour_vec, 
                       limits = c(min(plot_df_mean$Score), floor(max(plot_df_mean$Score))), 
                       na.value=colour_vec[length(colour_vec)])+
  guides(fill=guide_colorbar(barwidth =unit(2, "inch"), barheight = unit(0.2, "inch")), label="Average Score")+
  coord_equal()+
  scale_x_discrete(labels=as.vector(score_names_dict[unique(plot_df_mean$GeneSet)]), position="top")+
  coord_flip()


pdf(paste0(fig_dir, "FIG4B_DEheatmap_lower.pdf"),width = 6, height = 3)
pp_box
dev.off()




##### ----- FIG 4D VIOLIN PLOTS ------------------------------------------------



### --- Split epithelial cells up based on receptor expression

tmp_sce <- Fig3_sce[,Fig3_sce$CellType_Fig3=="Malignant"]
tmp_sce$FGFR_status <- ifelse(as.vector(counts(tmp_sce[rownames(tmp_sce) == "FGFR2",])>0), "FGFR2+", "FGFR2-")


coldata <- as.data.frame(colData(tmp_sce))
coldata <- coldata[,colnames(coldata) %in% c("FGFR_status",FGF_sigs)]

plot_data <- reshape2::melt(coldata, value.name="Score", variable.name="GeneSet", id.vars="FGFR_status") 
plot_data$Score <- scales::rescale(plot_data$Score, to=c(0,1))

rl=1
FIG4C <- ggplot(plot_data, aes(x=FGFR_status, y=Score, fill=FGFR_status))+
  geom_violin(alpha=1,lwd=0.3)+#theme_bw()+
  facet_wrap(~GeneSet, nrow=1, labeller = as_labeller(facet_labels))+
  xlab("FGFR2 Status")+scale_fill_manual(values=cols)+
  theme(strip.text.x = element_text(size = rel(rl)*1.1), panel.grid.minor = element_blank(), 
        axis.text.x = element_blank(),axis.ticks.x=element_blank(), axis.title.x=element_blank(),
        axis.title.y=element_text(size = rel(rl)*1.3), legend.position = "bottom")

pdf(paste0(fig_dir, "FIG4C_FGFRviolins.pdf"),width = 7, height = 3)
FIG4C+scale_y_continuous(limits=c(layer_scales(FIG4C)$y$range$range[1],1))
dev.off()


coldata$PatientID <- as.vector(setNames(Fig3_sce$PatientID, Fig3_sce$CellID)[rownames(coldata)])
plot_data <- reshape2::melt(coldata, value.name="Score", variable.name="GeneSet", id.vars=c("FGFR_status", "PatientID"))

## Exclude BCB21 for low cell numbers
table(paste0(tmp_sce$PatientID,"_",tmp_sce$FGFR_status))

FIG4C_byPatient <- ggplot(plot_data[!(grepl("BCB21",plot_data$PatientID)),],
                aes(x=FGFR_status, y=Score, fill=PatientID))+
  geom_violin(alpha=1,lwd=0.3,scale='width')+
  facet_wrap(~GeneSet, nrow=1, labeller = as_labeller(facet_labels),scales="free")+
  xlab("FGFR2 Status")+scale_fill_manual(values=cols)+
  theme(strip.text.x = element_text(size = rel(rl)*1.1), panel.grid.minor = element_blank(), 
        axis.title.x=element_blank(),axis.title.y=element_text(size = rel(rl)*1.3), legend.position = "bottom")+
  scale_fill_manual(values=Colour_values_PatientID)

FIG4C_byPatient+scale_y_continuous(limits=c(layer_scales(FIG4C_byPatient)$y$range$range[1],1))


ggplot(plot_data[!(grepl("BCB21",plot_data$PatientID)),],
       aes(x=PatientID, y=Score, fill=FGFR_status))+
  geom_violin(alpha=0.8,lwd=0.3,scale='width')+
  facet_wrap(~GeneSet, nrow=1, labeller = as_labeller(facet_labels),scales="free")+
  xlab("FGFR2 Status")+scale_fill_manual(values=cols)+
  theme(strip.text.x = element_text(size = rel(rl)*1.1), panel.grid.minor = element_blank(), 
        axis.title.x=element_blank(),axis.title.y=element_text(size = rel(rl)*1.3), legend.position = "bottom")

Fig4C_meanScores <- plot_data %>% group_by(PatientID, FGFR_status) %>%
    dplyr::summarize(MeanScore = mean(Score, na.rm=F)) %>%
  ggplot() + aes(x=PatientID, y=FGFR_status, fill=MeanScore) + geom_tile() +
  scale_x_discrete(expand=c(0,0))+scale_y_discrete(expand=c(0,0))+
  scale_fill_distiller(palette ="PuRd",direction=1)

Fig4C_medScores <- plot_data %>% group_by(PatientID, FGFR_status) %>%
  dplyr::summarize(MedScore = median(Score, na.rm=F)) %>%
  ggplot() + aes(x=PatientID, y=FGFR_status, fill=MedScore) + geom_tile() +
  scale_x_discrete(expand=c(0,0))+scale_y_discrete(expand=c(0,0))+scale_fill_viridis()
  scale_fill_distiller(palette ="PuRd",direction=1)



pdf(paste0(fig_dir, "FIG4Supp_FGFRstatus_perPatient.pdf"),width = 7, height = 4)
ggRNA::PropBar(tmp_sce, "FGFR_status", "PatientID")
dev.off()




##### ----- PLOT FIG 4E FGFR2 EXPRS -----------------------------------

#load(file=paste0(data_dir, "Fig3_sce.Rdata"))
#tmp_sce <- Fig3_sce[,Fig3_sce$CellType_Fig3=="Malignant"]
#tmp_sce$FGFR_status <- ifelse(as.vector(counts(tmp_sce[rownames(tmp_sce) == "FGFR2",])>0), "FGFR2+", "FGFR2-")
#coldata <- as.data.frame(colData(tmp_sce))
#coldata <- coldata[,colnames(coldata) %in% c("FGFR_status",FGF_sigs)]

tmp_sce <- tmp_sce[,tmp_sce$FGFR_status=="FGFR2+"]
tmp_sce <- tmp_sce[,!(tmp_sce$PatientID %in% c("BCB21", "BCB20","BCB20_E"))] ## Too few malignant/mesothelial cells

#SampleOrdering <- c("BCB139",  "BCB112", "BCB66","BCB20", "BCB66_E", "BCB20_E","BCB21_E")
SampleOrdering <- c("BCB21_E","BCB66_E","BCB66","BCB112","BCB139")
rl=0.8

### --- Get FGFR2+ cell numbers
numb_dat <- as.data.frame(colData(tmp_sce)) %>% 
      dplyr::count(PatientID, FGFR_status) 
numb_dat$PatientID <- factor(numb_dat$PatientID, levels=SampleOrdering)

cellNumb_bar <- ggplot(data=numb_dat, aes(x=n, y=PatientID)) + 
  geom_bar(stat="identity", position="dodge", width=0.93,  fill="dimgrey")+
  theme_bw()+labs(x="Number FGFR2+ Cells")+
  theme(axis.text.x=element_text(size = rel(rl)*1.2),
        axis.text.y=element_blank(),
        axis.line.x=element_line(colour="black",size=0.2),
        panel.border = element_rect(colour = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.title.y =element_blank(), axis.ticks.y=element_blank(),
        axis.title.x=element_text(size = rel(rl)*1.3))+ # margin=margin(10,0,0,0)
        scale_x_reverse(expand=c(0,0)) 


### --- Get FGFR2 log counts in FGFR2+ cells
vioin_dat <- data.frame(logCount=as.vector(logcounts(tmp_sce[rownames(tmp_sce)=="FGFR2",])),
                       PatientID=tmp_sce$PatientID)
vioin_dat$PatientID <- factor(vioin_dat$PatientID, levels=SampleOrdering)
patientViolins <- ggplot(data=vioin_dat, aes(x=logCount, y=PatientID, fill=PatientID)) + 
                  geom_violin(alpha=1,lwd=0.3,scale='width',trim=F)+ 
                  theme_bw()+ labs(x="FGFR2 Log-Expression") +
                  scale_y_discrete(expand=c(0,0))+ 
                  scale_fill_manual(values=Colour_values_PatientID)+
                    theme(axis.ticks.y=element_blank(),
                          axis.text.y=element_blank(),
                          axis.title.y=element_blank(),
                          axis.text.x=element_text(size = rel(rl)*1.2),
                          axis.title.x=element_text(size= rel(rl)*1.3),
                          legend.position = "none")
patientViolins <- patientViolins +scale_x_continuous(expand=c(0,0), limits=c(0, layer_scales(patientViolins)$x$range$range[2]))
  


### --- Get FGFR2 signalling scores
coldata$PatientID <- as.vector(setNames(Fig3_sce$PatientID, Fig3_sce$CellID)[rownames(coldata)])
plot_data <- reshape2::melt(coldata, value.name="Score", variable.name="GeneSet", id.vars=c("FGFR_status", "PatientID"))
plot_data <- plot_data[!(plot_data$PatientID == "BCB21"),]
plot_data <- plot_data[plot_data$GeneSet == "REACTOME_PI_3K_CASCADE_FGFR2",]
plot_data$PatientID <- factor(plot_data$PatientID, levels=rev(SampleOrdering))
plot_data$FGFR_status <- factor(plot_data$FGFR_status, levels=c("FGFR2-","FGFR2+"))

plot_data$Score <- scales::rescale(plot_data$Score, to=c(0,1))

df_median <- plot_data %>% 
  group_by(PatientID,FGFR_status) %>%
  summarise(medianScore = median(Score))
df_median$FGFR2col<- as.vector(cols[df_median$FGFR_status])

## Plot PI-3K pathway
patientBoxes <- ggplot()+ 
                geom_hline(data = df_median,size = 0.4,  linetype="solid",
                      aes(yintercept = medianScore, colour =FGFR_status)) +
                geom_boxplot(data=plot_data,
                      aes(x=FGFR_status, y=Score, fill=FGFR_status),
                      alpha=1,lwd=0.3, width=0.95,outlier.shape=NA,na.rm=T)+ 
                facet_wrap(~PatientID, ncol=1,scales="free")+
                ylab("Scaled Score - PI-3K Cascade: FGFR2")+theme_bw() +
                scale_fill_manual(values=cols)+  
                scale_colour_manual(values=cols) +
                theme(strip.text.x = element_blank(),
                      axis.text.x=element_blank(),
                      axis.text.y=element_text(size=rel(rl)*0.9),
                      axis.ticks.x=element_blank(),
                      panel.grid.minor = element_blank(), 
                      panel.grid.major = element_blank(), 
                      legend.position = "none",
                      axis.title.x=element_blank(),
                      axis.title.y=element_text(size =rel(rl)*1.3),
                      panel.spacing = unit(0.02, "lines"))+
                scale_y_continuous(position = "right",expand=c(0,0))+
                coord_cartesian(ylim =quantile(plot_data$Score, c(0.05, 0.95)))

### --- Combine plots
pdf(paste0(fig_dir,"Fig4_FGFR2_perPatient.pdf"),width=8, height=6)
grid.arrange(grobs=list(cellNumb_bar+theme(plot.margin=unit(c(0,-0.1,0.01,0), "cm")), 
                        patientViolins+theme(plot.margin=unit(c(0.07,0,0,0), "cm")),
                        patientBoxes+theme(plot.margin=unit(c(0.08,0,0.8,0), "cm"))),
             nrow=1, widths = c(0.7,1.5,0.8))
dev.off()







