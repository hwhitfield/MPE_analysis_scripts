

##### ----- Get directories


data_dir <- paste0(getwd(),"/data/")
matrix_dir <- paste0(data_dir,"raw_matrices/")
pal_etal_dir <- "~/data/_datasets/sc/Smyth_HumanBrCa/"
pal_save_dir <- paste0(data_dir,"pal_etal_data/")
fig_dir <- paste0(getwd(),"figures/")

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



##### ----- LOAD & PROCESS DATA ---------------------------------------------------

load(file=paste0(pal_save_dir, "/PAL_sce.Rdata")) # -- 112,685 cells from Figure1E_and_Supp2.R

## -- Subset & process immune
PAL_immune_sce <- PAL_sce[,!(PAL_sce$celltype %in% c("12", "13", "5", "Epithelial", "Endothelial", "CAFs", "Pericytes"))]

PAL_immune_sce <- PAL_immune_sce[as.vector(rowSums(logcounts(PAL_immune_sce)>0)>5),]

g <- buildSNNGraph(PAL_immune_sce, 
                   k=30,  assay.type = "logcounts")
clust <- igraph::cluster_walktrap(g)$membership
PAL_immune_sce$Clust_Knn30 <- as.vector(clust)

#save(PAL_immune_sce,file=paste0(pal_save_dir, "/PAL_immune_sce.Rdata"))




##### ----- FIGURE 2D & SUPP FIG 6 ---------------------------------------------------

load(file=paste0(pal_save_dir, "/PAL_immune_sce.Rdata"))

### --- Get T cell marker %
PAL_T_sce <- PAL_immune_sce[, PAL_immune_sce$celltype == "T cells"]
rm(PAL_immune_sce)
gc()

PAL_T_sce$CD4pos <- as.vector(counts(PAL_T_sce[rownames(PAL_T_sce) == "CD4",])>0)
PAL_T_sce$CD8pos <- as.vector(colSums(counts(PAL_T_sce[rownames(PAL_T_sce) %in% c("CD8A", "CD8B"),])>0)>0)

PAL_T_sce$Marker <- PAL_T_sce$celltype
PAL_T_sce$Marker <- replace(PAL_T_sce$Marker,
                            PAL_T_sce$CD4pos & PAL_T_sce$CD8pos,  "CD4+/CD8+")
PAL_T_sce$Marker <- replace(PAL_T_sce$Marker,
                            PAL_T_sce$CD4pos & !(PAL_T_sce$CD8pos),  "CD4+")
PAL_T_sce$Marker <- replace(PAL_T_sce$Marker,
                            !(PAL_T_sce$CD4pos) & PAL_T_sce$CD8pos,  "CD8+")
PAL_T_sce$Marker <- replace(PAL_T_sce$Marker,
                            !(PAL_T_sce$CD4pos) & !(PAL_T_sce$CD8pos),  NA)


### --- Prep for plotting
PAL_T_sce$Subtype <- ifelse(grepl("ER",PAL_T_sce$PatientID), "ER","TNBC")

plot_dat <- rbind(as.data.frame(colData(PAL_T_sce)[,colnames(colData(PAL_T_sce)) %in% c("Marker","PatientID","Subtype")]),
      as.data.frame(colData(MPE_T_sce)[,colnames(colData(MPE_T_sce)) %in% c("Marker","PatientID","Subtype")]))

plot_dat <- plot_dat[!(is.na(plot_dat$Marker)),]
plot_dat$Marker <- factor(plot_dat$Marker, levels=rev(c("CD8+","CD4+","CD4+/CD8+")))


gg_dat <- plot_dat %>% 
  count(Dataset, PatientID, Marker) %>%       
  group_by(Dataset,PatientID) %>%
  mutate(pct= prop.table(n) * 100) 
gg_dat$Subtype <- ifelse(grepl("TN",gg_dat$PatientID), "TNBC Patients","ER+ Patients")

### --- Plot
gg_bar <- ggplot(data=gg_dat) + aes(PatientID, pct, fill=Marker) + 
  geom_bar(stat="identity")+#facet_wrap(~Subtype,ncol=1, scales="free")+
  theme_bw()+
  theme(axis.ticks.y=element_blank(),axis.text.y=element_text(size=16),
        axis.title.x=element_text(size=14,margin = unit(c(5, 0, 0, 0), "mm")),
        strip.text=element_text(size=20))+
  scale_fill_manual(values=(c("grey","#CABEE9","#7C7189")))+
  scale_y_continuous(expand=c(0,0))+coord_flip()+ 
  labs(y="T cell Percent (%)",x="")

pdf(paste0(fig_dir,"SUPPFIG2_TcellMarkers_Primary.pdf"), width=7, height=5)
gg_bar
dev.off()




