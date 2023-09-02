######################################################################################################
##
##  FIGURE 4 -- IZAR DATA
##
## 



saveFolder = ""


## -- Load packages & functions
source("./helper_functions.R")



## -- Colours
IZAR_cols <- c("#9e5476", "#d18b79",  "#70a18f",
               "#e7ebbc", "pink",
               "#8bbde6", "#6074ab")

names(IZAR_cols) <- c("Malignant", "Fibroblasts",  "Macrophage",  
                      "Bcells", "Erythrocytes",
                      "DC", "Tcells")

patient_cols <- c("#FA8345", "#CC281B", "#DE007F",
                  "#FFE184", "#F2F28B", "#1F87BD",
                  "#3749A4", "#A8E686")
names(patient_cols) <- c("Patient 1", "Patient 2", "Patient 2.1",
                         "Patient 3", "Patient 4", "Patient 5",
                         "Patient 5.1", "Patient 6")

clust_cols <- c(blend_colours("#BC829C", "#67374D", 5), # Cancer clusters
                blend_colours("#DFAEA1", "#B7593F", 4), # Fibroblast clusters
                blend_colours("#94BAAC", "#44685B", 4), # macrophages
                "#A3CAEB", "#61A5DD", # DC
                "#e7ebbc", # B cell
                "#6074ab", # T cell
                "pink")
names(clust_cols) <- as.character(1:18)


tsneTheme <- theme(panel.grid.minor = element_blank(),
                   panel.grid.major = element_blank(),
                   axis.text.y = element_blank(), axis.text.x = element_blank(), 
                   axis.ticks.y=element_blank(), axis.ticks.x=element_blank(),
                   axis.title.y=element_blank(), axis.title.x=element_blank(),
                   axis.line.y = element_blank(),axis.line.x = element_blank())


##### ----- LOAD DATA --------------------------------------------------------

## -- Scores
load(file=paste(saveFolder,"Izar_scores.Rdata", sep="")) # melted_df_IZAR

## -- SCE object
load(file=paste0(saveFolder, "Izar_sce.Rdata"))
Izar_sce <- Izar_sce[,!(Izar_sce$CellType %in% c("-1", "NotAssigned"))]
Izar_sce$patient_id <- replace(Izar_sce$patient, grepl(".1",Izar_sce$sample_ID) & grepl("2",Izar_sce$patient), "2.1")
Izar_sce$patient_id <- replace(Izar_sce$patient_id, grepl(".1",Izar_sce$sample_ID) & grepl("5",Izar_sce$patient_id), "5.1")
Izar_sce$patient_id <- paste0("Patient ", Izar_sce$patient_id)

tsne_df <- data.frame(tsne_x=as.numeric(Izar_sce$TSNE_x), tsne_y=as.numeric(Izar_sce$TSNE_y), 
                      patient= Izar_sce$patient_id, celltype=Izar_sce$CellType, 
                      cluster=Izar_sce$clst, cellID = colnames(Izar_sce))


## -- Original tSNE & metadata from paper can be found at:
## https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE146026

## -- Fix custers
tsne_df$clust2 <- tsne_df$cluster == "2"
tsne_df$clust5 <- tsne_df$cluster == "5"
tsne_df$clust20 <- tsne_df$cluster == "20" 
tsne_df$clust7<- tsne_df$cluster == "7" 

clust_dict <- as.character(1:18)
names(clust_dict) <- c("21","20", "17", "4",
                       "7","9", "3", "2",
                       "13","11", "10", "5",
                       "1","18", "19", "8",
                       "15", "16")
tsne_df$correct_clust <- as.vector(clust_dict[tsne_df$cluster])
melted_df_IZAR$Clust <- as.vector(clust_dict[melted_df_IZAR$Clust])



##### ----- CURATE GENE SET THEMES ---------------------------------------------
proliferation_related <- as.vector(unique(melted_df_IZAR$geneset))[grepl("PROLIF|STEM_CELL_UP",as.vector(unique(melted_df_IZAR$geneset)))] 
metabolism_related <- as.vector(unique(melted_df_IZAR$geneset))[grepl("METABOL|EFFLUX|TRANSPORT|RIBOS|TRANSLATION", as.vector(unique(melted_df_IZAR$geneset)))] 

cancer_related <- as.vector(unique(melted_df_IZAR$geneset))[
  grepl("CANCER|METASTA|AML|MYEL|GBM|BRCA|TUMOR_FIELD|TOMA|NOJIMA_SFRP2|FLUOROURACIL|TNF|PAX8|RASSF1", 
        as.vector(unique(melted_df_IZAR$geneset)))]


EMT_related <- as.vector(unique(melted_df_IZAR$geneset))[grepl("WNT|ECM|JUNCTION|EMT|ZEB", as.vector(unique(melted_df_IZAR$geneset)))]
CAF_related <- as.vector(unique(melted_df_IZAR$geneset))[grepl("COLLAGEN|MUSCLE|FIBRO|VEGFA", as.vector(unique(melted_df_IZAR$geneset)))]
IMMUNE_related <- as.vector(unique(melted_df_IZAR$geneset))[grepl("PHOCYTE|_B_|_T_|NEUTROPHIL|MACROPHAG|ASTHMA|COMPLEMENT|ANTIGEN|INFECTION", as.vector(unique(melted_df_IZAR$geneset)))]

IMMUNE_related <- IMMUNE_related[!(IMMUNE_related %in% cancer_related)]
cancer_related <- cancer_related[!(cancer_related %in% c(CAF_related, EMT_related, proliferation_related))]

proliferation_related <- c("YAMASHITA_LIVER_CANCER_STEM_CELL_UP", "ZHANG_PROLIFERATING_VS_QUIESCENT")
cancer_related <- c("LIU_PROSTATE_CANCER_UP", "NOJIMA_SFRP2_TARGETS_DN","LANDIS_BREAST_CANCER_PROGRESSION_UP",
                    
                    "KANG_FLUOROURACIL_RESISTANCE_UP",
                    "SATO_SILENCED_EPIGENETICALLY_IN_PANCREATIC_CANCER",
                    "HOFMANN_MYELODYSPLASTIC_SYNDROM_LOW_RISK_DN",
                    "NOUSHMEHR_GBM_SILENCED_BY_METHYLATION",
                    
                    "MACLACHLAN_BRCA1_TARGETS_UP", "CHEN_LUNG_CANCER_SURVIVAL", 
                    "HWANG_PROSTATE_CANCER_MARKERS",
                    "YANG_BREAST_CANCER_ESR1_BULK_UP", 
                    "ZHOU_TNF_SIGNALING_30MIN",  "LI_LUNG_CANCER")

EMT_related <- c("WNT_SIGNALING", "REACTOME_TIGHT_JUNCTION_INTERACTIONS", "AIGNER_ZEB1_TARGETS") 
CAF_related <- c("WESTON_VEGFA_TARGETS_12HR",
                 "REACTOME_ECM_PROTEOGLYCANS","REACTOME_COLLAGEN_DEGRADATION","REACTOME_COLLAGEN_BIOSYNTHESIS_AND_MODIFYING_ENZYMES",
                 "MISHRA_CARCINOMA_ASSOCIATED_FIBROBLAST_DN", "REACTOME_SMOOTH_MUSCLE_CONTRACTION")
IMMUNE_related <- c("AKL_HTLV1_INFECTION_UP",
                    "REACTOME_INITIAL_TRIGGERING_OF_COMPLEMENT",
                    "LIAN_NEUTROPHIL_GRANULE_CONSTITUENTS","HADDAD_T_LYMPHOCYTE_AND_NK_PROGENITOR_DN",
                    "MORI_LARGE_PRE_BII_LYMPHOCYTE_DN",
                    "KEGG_ASTHMA")
B_plasma <- c("TARTE_PLASMA_CELL_VS_B_LYMPHOCYTE_UP","SHIN_B_CELL_LYMPHOMA_CLUSTER_8",
              "PELLICCIOTTA_HDAC_IN_ANTIGEN_PRESENTATION_UP", "PELLICCIOTTA_HDAC_IN_ANTIGEN_PRESENTATION_DN")

geneset_themes <- c(rep("Metabolism-related", length(metabolism_related)),
                    rep("EMT-related", length(EMT_related)),
                    rep("Proliferation", length(proliferation_related)),
                    rep("Cancer", length(cancer_related)),
                    
                    rep("CAF/Fibroblast", length(CAF_related)),
                    rep("Immune-related", length(IMMUNE_related)),
                    rep("B/Plasma", length(B_plasma)))
names(geneset_themes) <- c(metabolism_related, EMT_related, proliferation_related, cancer_related, CAF_related, IMMUNE_related,B_plasma)


bio_theme_cols <- c(usecol(pal_unikn_pref)[1:4], usecol(pal_unikn_pref)[6:7])
names(bio_theme_cols) <- c("Immune-related", "EMT-related", "Proliferation","CAF/Fibroblast",
                           "B/Plasma", "Cancer")

geneset_cols <- as.vector(bio_theme_cols[as.vector(geneset_themes)])
names(geneset_cols) <- names(geneset_themes)
geneset_cols[setdiff(as.vector(unique(melted_df_IZAR$geneset)), names(geneset_themes))] <- "grey"

geneset_themes[setdiff(as.vector(unique(melted_df_IZAR$geneset)), names(geneset_themes))] <- "none"

genesets_of_interest <- names(geneset_themes[!(geneset_themes %in% c("none", "Metabolism-related"))])



##### ----- PLOT FIG 4B -----------------------------------------------------

## -- Set row ordering
celltype_ordering <- c("Bcells","DC","Erythrocytes", "Tcells", "Macrophage",
                       "Malignant", "Fibroblasts")
clust_to_cellType_df <-  tsne_df[!(duplicated(tsne_df$correct_clust)),colnames(tsne_df) %in% c("celltype", "correct_clust")]
clust_to_cellType_df <- clust_to_cellType_df[order(match(clust_to_cellType_df$celltype,celltype_ordering)),]

clust_to_cellType_dict <- clust_to_cellType_df$celltype
names(clust_to_cellType_dict) <- clust_to_cellType_df$correct_clust

## -- Mean per cluster ------------------------------------------
plot_df_clust <- melted_df_IZAR  %>%
  group_by(Clust, geneset) %>%
  summarize(
    n = n(),
    mean_score = mean(score, na.rm=T)
  ) 
short_df_clust <- reshape2::dcast(plot_df_clust, Clust ~ geneset)
rownames(short_df_clust) <- short_df_clust$Clust
short_df_clust <- short_df_clust[,!(colnames(short_df_clust)=="Clust")]

## -- Get col ordering 
short_df_clust <- short_df_clust[names(clust_to_cellType_dict),]
hcr <- hclust(as.dist(1-cor.pairs(short_df_clust[,genesets_of_interest]+1,cor.method="spearman")))
ddr <- as.dendrogram(hcr)
colInd <- order.dendrogram(ddr)
#short_df_clust <- short_df_clust[,colInd]
short_df_clust <- short_df_clust[,match(names(geneset_cols), colnames(short_df_clust))]

## Get colour row bar
heatmap_row_cols_clust <- IZAR_cols[as.vector(clust_to_cellType_dict[rownames(short_df_clust)])]
names(heatmap_row_cols_clust) <- rownames(short_df_clust)

## Get colour col bar
heatmap_col_cols_clust <- geneset_cols[colnames(short_df_clust)]
names(heatmap_col_cols_clust) <- colnames(short_df_clust)

heatmap_dat <- t(short_df_clust[,genesets_of_interest][,colInd])

svg(paste0(SAVE_DIR,"FIG3_scoreHeatmap_genesOfInterest.svg"), width=5, height=7)
heatmap3(heatmap_dat, scale = "row", main="",RowSideLabs="", ColSideLabs="",margins=c(6,1),
         RowSideColors=heatmap_col_cols_clust[rownames(heatmap_dat)],
         ColSideColors = heatmap_row_cols_clust,
         xlab="", ylab="",	Rowv =TRUE, Colv=TRUE, showColDendro=TRUE,showRowDendro=TRUE,
         col=rev(brewer.pal(name = "YlGnBu", n=9)), cexRow=0.6,) 
dev.off()


svg(paste0(SAVE_DIR,"FIG3_scoreHeatmap_legends.svg"), width=7, height=5)
cowplot::plot_grid(GetLegend(bio_theme_cols, "Biological Theme",n_row=3, by_row=TRUE), 
                   GetLegend(IZAR_cols, "CellType",n_row=3, by_row=TRUE), 
                   ncol=1)
dev.off()


##### ----- PLOT FIG 4A -----------------------------------------------------


### --- Original tSNE by cell type & patient

rl=1.1
originalTSNE_gg <- ggplot(aes(x=tsne_x, y=tsne_y, colour=patient), data = tsne_df)+geom_point(size=1, alpha=0.9)+scale_colour_manual(values=patient_cols)+
  theme_bw()+
  theme( panel.grid.minor = element_blank(),
         panel.grid.major = element_blank(),
         axis.line.y = element_blank(),axis.line.x = element_blank())+ labs(colour = "Sample")+
  guides(color = guide_legend(override.aes = list(size = 4)))

celltypeTSNE_gg <- ggplot(aes(x=tsne_x, y=tsne_y, colour=celltype), data = tsne_df)+geom_point(size=1, alpha=0.9)+scale_colour_manual(values=IZAR_cols)+
  theme_bw()+
  theme( panel.grid.minor = element_blank(),
         panel.grid.major = element_blank(),
         axis.line.y = element_blank(),axis.line.x = element_blank())+ labs(colour = "Cell Type")+
  guides(color = guide_legend(override.aes = list(size = 4)))

svg(paste0(SAVE_DIR,"tsnes/FIG3_originalTSNE.svg"), width=6, height=5)
print(originalTSNE_gg)
dev.off()


svg(paste0(SAVE_DIR,"tsnes/FIG3_originalTSNE_byCellType.svg"), width=6, height=5)
print(celltypeTSNE_gg)
dev.off()

### --- Original tSNE by scDECAF score

genesets_of_interest <- c("ZHOU_TNF_SIGNALING_30MIN", 
                          "ZHANG_PROLIFERATING_VS_QUIESCENT", 
                          "AIGNER_ZEB1_TARGETS",
                          "YAMASHITA_LIVER_CANCER_STEM_CELL_UP", 
                          "WESTON_VEGFA_TARGETS_12HR", "PELLICCIOTTA_HDAC_IN_ANTIGEN_PRESENTATION_UP")

## -- Plots
plot_lst <- lapply(genesets_of_interest,function(gs_i){
  tsnePanel(melted_df_IZAR, tsne_df, gs_i, bar_width = 7, bar_height = 0.7)})
names(plot_lst) <- genesets_of_interest

## -- Fix legends
plot_lst$CellType <- celltypeTSNE_gg+tsneTheme+theme(legend.position = "bottom")+
  guides(color = guide_legend(nrow = 2, override.aes = list(size = 4)))
plot_lst$Patient <- originalTSNE_gg+tsneTheme+theme(legend.position = "bottom")+
  guides(color = guide_legend(nrow = 2, override.aes = list(size = 4)))

## -- Save
for (i_plot in names(plot_lst)){
  png(paste0(SAVE_DIR,"tsnes/FIG3_tSNE_", i_plot,".png"), width=800, height=800, res=160)
  print(plot_lst[[i_plot]]+labs(x="tsne_x", y="tsne_y"))
  dev.off()
  
}



##### ----- PLOT FIG 4C -----------------------------------------------------

genesets_to_plot <- c("TARTE_PLASMA_CELL_VS_B_LYMPHOCYTE_UP", "MARTIN_NFKB_TARGETS_UP", "WNT_SIGNALING",
                      "LANDIS_BREAST_CANCER_PROGRESSION_UP")

ggBoxplot <- function(i_gs){
  plot_df <- melted_df_IZAR[melted_df_IZAR$geneset == i_gs,]
  gg <- ggplot(plot_df, aes(x=reorder(Clust,score,na.rm = TRUE), y=score, fill=CellType)) + geom_boxplot() +
    facet_wrap(~geneset, scales = "free_x") + theme_bw()+labs(x="Cluster", y="Score", fill="Cell type")+
    scale_fill_manual(values=IZAR_cols) + 
    theme(strip.text.x=element_text(size=10), axis.ticks.x= element_blank(),
          legend.title.align = 0.5)
  return(gg)
}

## -- Plot
plot_lst <- lapply(genesets_to_plot, function(x){ggBoxplot(x)})
names(plot_lst) <- genesets_to_plot

## -- Only one legend
celltype_leg <- cowplot::get_legend(plot_lst$MARTIN_NFKB_TARGETS_UP+guides(fill=guide_legend(nrow=4)))
plot_lst <- lapply(plot_lst, function(x){x+theme(legend.position = "none")})
plot_lst$Legend <- celltype_leg

## -- Save
for (i_plot in names(plot_lst)){
  png(paste0(SAVE_DIR,"boxplots/FIG3_boxplots_", i_plot,".png"), width=1200, height=600, res=160)
  print(plot_lst[[i_plot]]+theme(strip.background = element_blank(),
                                 strip.text.x = element_blank()))
  dev.off()
  
}




