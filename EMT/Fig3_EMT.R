######################################################################################################
##
##  FIGURE 3 -- EMT DATA
##
## 


DATA_DIR = ""
SAVE_DIR = ""


## -- Load packages & functions
source("./helper_functions.R")

## -- Colours
Spec_Div <- brewer.pal(11, "Spectral")
time_cols <- c("grey", Spec_Div[1],   # 0 & 8h
               Spec_Div[4],  # 1day
               Spec_Div[8],  # 3-7 day
               Spec_Div[10],  # 8h_rm
               Spec_Div[11])   #1-3d_rm
names(time_cols) <- factor(x=c("0d","8h", "1d", "3-7d", "8h_rm", "1-3d_rm"), 
                           levels=c("0d", "8h", "1d", "3-7d", "8h_rm", "1-3d_rm"))

EMT_cols <- c("#D04E59", "#FAE093", "#7C7189")
names(EMT_cols) <- c("Control", "EMT", "MET")



##### ----- DIFFERENTIAL EXPRESSION ANALYSIS -----------------------------------

load(file=paste0(DATA_DIR, "CookedEMT_DGE_all.Rdata"))
CookEMT_DGE <- DGE_Summed_sce

CellLine_vec <- as.vector(CookEMT_DGE$samples$CellLine)
Time_vec <- as.vector(CookEMT_DGE$samples$Time)
design_mat <- model.matrix( ~0 + Time_vec + CellLine_vec)
design_to_preserve <- model.matrix( ~0 + Time_vec)

### --- Perform DE
CookEMT_DGE <- estimateDisp(CookEMT_DGE, design_mat, robust = TRUE)
qfit <- glmQLFit(CookEMT_DGE, design_mat)
EMT_DE <- glmQLFTest(qfit, contrast = c(-1,0,0,1,0,0,0,0)) 
MET_DE <- glmQLFTest(qfit, contrast = c(-1,1,0,0,0,0,0,0)) 

DE_obj <- topTags(EMT_DE, n=Inf, p.value=0.05, adjust.method = "BH", sort.by = "PValue")
table(DE_obj$table$logFC>0)
write.table(DE_obj$table[DE_obj$table$FDR < 0.05,],
            file=paste0(SAVE_DIR, "deGenes.txt", ""),sep="\t", col.names=TRUE, row.names=FALSE)

### --- Get Signatures
C2_gmt <- getGmt("/stornext/Home/data/allstaff/w/whitfield.h/data_load/genesets/c2.all.v7.1.symbols.gmt")
C2_genesets <- lapply(C2_gmt, function(x){geneIds(x)})
names(C2_genesets) <- names(C2_gmt)
C2_genesets <- C2_genesets[lapply(C2_genesets, length) > 10]
idx_C2 <- ids2indices(C2_genesets, id=rownames(CookEMT_DGE))

### --- Control Vs EMT
idx_cam <- ids2indices(C2_genesets ,id=rownames(CookEMT_DGE$counts))
cam <- camera(CookEMT_DGE$counts, index=idx_cam, design=design_mat, c(-1,0,0,1,0,0,0,0)) 
cam <- cam[cam$PValue < 0.05,]
cam <- cam[cam$FDR < 0.05,]
cam$Geneset <- rownames(cam)
write.table(cam, paste0(SAVE_DIR, "EMT_de_enriched.csv", ""), sep=",", row.names = FALSE)

### --- Batch correction
CookEMT_DGE$batch_corrected <- removeBatchEffect(cpm(CookEMT_DGE, log=TRUE, prior.count = 1), design = design_to_preserve, batch=CellLine_vec)


##### ----- PLOT FIG 3A -----------------------------------------------------

### --- Plot legend separately
SideBar_dict <- setNames(CookEMT_DGE$samples$Time,rownames(CookEMT_DGE$samples))
SideBarColours_dict <- as.vector(time_cols[SideBar_dict])
names(SideBarColours_dict) <- names(SideBar_dict)
SideBar_leg <- GetLegend(time_cols, "Time")


svg(paste0(SAVE_DIR,"FIG2_pbEMT_heatmap_legend.svg"), width=4, height=2)
plot(SideBar_leg)
dev.off()

### --- TOP 50 UP & DN 
n_genes = Inf
featureTab <- GetTT(EMT_DE, n_genes)
DNgenes <- featureTab[featureTab$table$logFC < 0,]
UPgenes <- featureTab[featureTab$table$logFC > 0,]
featureSet <- c(rownames(UPgenes[order(UPgenes$table$PValue),][1:50,]), 
                rownames(DNgenes[order(DNgenes$table$PValue),][1:50,]))

topgenes_log_expression.full <- CookedEMT_DGE$batch_corrected[rownames(CookedEMT_DGE$batch_corrected) %in% featureSet,]
topgenes_log_expression.full <- t(scale(t(topgenes_log_expression.full)))

### --- Heatmap
col_order <- rownames(CookedEMT_DGE$samples[order(match(CookedEMT_DGE$samples$Time,c("0d", "8h", "1d", "3-7d", "8h_rm", "1-3d_rm"))),])

svg(paste0(SAVE_DIR,"FIG2_pbEMT_heatmap_50UpDn_inOrder.svg"), width=5, height=7,antialias ="none")
heatmap.2(topgenes_log_expression.full[,col_order], trace="none",
          scale = "column", dendrogram = "row",Rowv=TRUE, Colv=FALSE,cexRow=0.2,
          labRow=rep("", nrow(topgenes_log_expression.full)),
          labCol=unlist(lapply(strsplit(colnames(topgenes_log_expression.full[,col_order]), "_"), function(x){x[[1]]})),
          #margin=c(1,1),trace="none", 
          #lhei = c(4,10),lwid = c(6,6), lmat=rbind( c(3, 4), c(1,2)),
          col=rev(brewer.pal(name = "YlGnBu", n=9)), ColSideColors=SideBarColours_dict[col_order],
          key=FALSE) #   lhei = c(2,10),lwid = c(2,6)
dev.off()




##### ----- PLOT FIG 3B -----------------------------------------------------


final_genesets <- c("VERRECCHIA_RESPONSE_TO_TGFB1_C1", 
                    "HOLLERN_EMT_BREAST_TUMOR_UP",
                    "MCBRYAN_PUBERTAL_TGFB1_TARGETS_DN")

### PDFs
for (i_gs in final_genesets){
  pdf(paste0(SAVE_DIR, "barcodePlots/FIG2_pbEMT_barcode_", i_gs,".pdf"), width=7, height=5)
  barcodeplot(EMT_DE$table$logFC, index=idx_C2[[i_gs]],
              labels=c("CONTROL", "EMT"), xlab="Log-fold Change")
  dev.off()
}


##### ----- PLOT FIG 3C -----------------------------------------------------

ggBox <- function(emt_df, x_gs){
  plot_df <- emt_df[emt_df$geneset==x_gs,]
  
  gg <- ggplot(plot_df, aes(x=status, y=score, fill=status)) + geom_boxplot() +
    facet_wrap(~method, scales = "free_x") + 
    scale_fill_manual(values=EMT_cols) + theme_bw()+
    theme(strip.text.x=element_text(size=14), axis.ticks.x=element_blank(),
          axis.title.x=element_blank())+labs(y="Score", fill="Status")+ 
    guides(fill=guide_legend(nrow=1))
  return(gg)
}

load(file=paste(saveFolder,"CookEMT_TGFB1scores_ALLcellLines.Rdata", sep=""))

for (x_gs in final_genesets){
  pdf(paste0(SAVE_DIR,"boxPlots/FIG2_pbEMT_boxplots_", x_gs,".pdf"), width=10, height=5)
  print(ggBox(CookedEMT_TGFB1_melt, x_gs)+theme(axis.text.x=element_blank(), legend.position = "none"))
  dev.off()
}


