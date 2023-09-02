######################################################################################################
##
##  HELPER FUNCTIONS
##
## 


## LIBRARIES ----------------- 

requiredPackages <-
  c("DropletUtils","SingleCellExperiment", "scater","scran", "edgeR",
    "tidyverse", "dplyr", "data.table",
    "PMA","reshape2", 
    
    ## Scoring
    "GSEABase", "VISION", "scDECAF",  "org.Hs.eg.db",
    
    ## Plotting
    "ggplot2", "gplots", "pals", "gridExtra", "ggpubr", "unikn", "RColorBrewer",
    "heatmap3", "Manu", "heatmap3")

for(pkg in requiredPackages){
  suppressWarnings(suppressMessages(library(pkg, character.only = T)))
}

##### ----- GENESET FUNCTIONS ---------------------------------------

GetMSigDB <- function(database="Hallmarks"){
  
  if (substr(getwd(), 2, 9) == "stornext"){
    Hallmark_Geneset_DIR="/stornext/Home/data/allstaff/w/whitfield.h/data_load/genesets/h.all.v7.0.symbols.gmt"
    C2_DIR="/stornext/Home/data/allstaff/w/whitfield.h/data_load/genesets/c2.all.v7.1.symbols.gmt"
    SC_DIR="/stornext/Home/data/allstaff/w/whitfield.h/data_load/genesets/scsig.all.v1.0.1.symbols.gmt"
  } else {
    Hallmark_Geneset_DIR="C:\\Users\\whitfield.h\\Desktop\\Projects\\B_I_G_Data\\GeneSets\\MSigDB\\msigdb.v7.2.symbols.gmt"
    C2_DIR="C:\\Users\\whitfield.h\\Desktop\\Projects\\B_I_G_Data\\GeneSets\\MSigDB\\c2.all.v7.1.symbols.gmt"
    SC_DIR="C:\\Users\\whitfield.h\\Desktop\\Projects\\B_I_G_Data\\GeneSets\\MSigDB\\scsig.all.v1.0.1.symbols.gmt"
  }
  
  if (database == "Hallmarks") {
    Hallmark_Genesets <- list()
    conn <- file(Hallmark_Geneset_DIR,open="r")
    linn <-readLines(conn)
    for (i in 1:length(linn)){
      iLine <- linn[i]
      iLine <- strsplit(iLine, "\t")
      iName <- iLine[[1]][1]
      iList <- iLine[[1]]
      iList <- iList[3:length(iList)]
      Hallmark_Genesets[[iName]] <- iList
      # Hallmark_Genesets <- c(Hallmark_Genesets, list(GeneSet(setName=iName, geneIds=iList)))
    }
    close(conn)
    return(Hallmark_Genesets)
    
  } else if (database == "C2") {
    C2_Genesets <- list()
    conn <- file(C2_DIR,open="r")
    linn <-readLines(conn)
    for (i in 1:length(linn)){
      iLine <- linn[i]
      iLine <- strsplit(iLine, "\t")
      iName <- iLine[[1]][1]
      iList <- iLine[[1]]
      iList <- iList[3:length(iList)]
      C2_Genesets[[iName]] <- iList
      #C2_Genesets <- c(C2_Genesets, list(GeneSet(setName=iName, geneIds=iList)))
    }
    close(conn)
    return(C2_Genesets)
  } else if (database == "Both"){
    Genesets <- list()
    conn <- file(Hallmark_Geneset_DIR,open="r")
    linn <-readLines(conn)
    for (i in 1:length(linn)){
      iLine <- linn[i]
      iLine <- strsplit(iLine, "\t")
      iName <- iLine[[1]][1]
      iList <- iLine[[1]]
      iList <- iList[3:length(iList)]
      Genesets[[iName]] <- iList
      
    }
    close(conn)
    
    conn <- file(C2_DIR,open="r")
    linn <-readLines(conn)
    for (i in 1:length(linn)){
      iLine <- linn[i]
      iLine <- strsplit(iLine, "\t")
      iName <- iLine[[1]][1]
      iList <- iLine[[1]]
      iList <- iList[3:length(iList)]
      Genesets[[iName]] <- iList
      
    }
    close(conn)
    
    return(Genesets)
  } else if (database == "SC"){
    SC_Genesets <- list()
    conn <- file(SC_DIR,open="r")
    linn <-readLines(conn)
    for (i in 1:length(linn)){
      iLine <- linn[i]
      iLine <- strsplit(iLine, "\t")
      iName <- iLine[[1]][1]
      iList <- iLine[[1]]
      iList <- iList[3:length(iList)]
      SC_Genesets[[iName]] <- iList
      # Hallmark_Genesets <- c(Hallmark_Genesets, list(GeneSet(setName=iName, geneIds=iList)))
    }
    close(conn)
    return(SC_Genesets)
  }
  
}




MakeVisionObject <- function(genesetList, biDirection=TRUE){
  MSigDB_signatures <- c()
  message("# Converting genesets to Vision objects")
  
  if (biDirection){
    if (length(genesetList[grepl("UP", names(genesetList))]) > 0){
      
      x <- genesetList[grepl("UP", names(genesetList))]
      y <- unlist(strsplit(names(x), "_UP"))
      
      
      for (iDirGeneSet in y){
        
        # print(paste0(" -> Extracting ", iDirGeneSet, " UP and DN..."))
        
        UP_str <- paste0(iDirGeneSet, "_UP", "")
        DN_str <- paste0(iDirGeneSet, "_DN", "")
        
        UP_genes <- genesetList[names(genesetList) == UP_str][[UP_str]]
        DN_genes <- genesetList[names(genesetList) == DN_str][[DN_str]]
        
        ### Make signature object
        UPDN_genes <- c(rep(1,length(UP_genes)), rep(-1, length(DN_genes)))
        names(UPDN_genes) <- c(UP_genes, DN_genes)
        iDirSet <- VISION::createGeneSignature(name = iDirGeneSet, sigData = UPDN_genes)
        
        MSigDB_signatures <- c(MSigDB_signatures, iDirSet)
      }
      
      genesetList <- genesetList[!(names(genesetList) %in% c(paste0(y, "_UP"), paste0(y, "_DN")))]
    }
  }
  
  for (iSig in names(genesetList)){
    # print(paste0("-> Extracting ",iSig,""))
    iGeneIDs <- genesetList[names(genesetList) %in% c(iSig)][[iSig]]
    
    ### Make signature object
    iSig_genes <- rep(1,length(iGeneIDs))
    names(iSig_genes) <- iGeneIDs
    SigObject <- createGeneSignature(name = iSig, sigData = iSig_genes)
    
    MSigDB_signatures <- c(MSigDB_signatures, SigObject)
    
  }
  
  return(MSigDB_signatures)
}


GetVisionObject_MSigDB <- function(){
  
  if (substr(getwd(), 2, 9) == "stornext"){
    Hallmark_Geneset_DIR="/stornext/Home/data/allstaff/w/whitfield.h/BIG_DATA/GeneSets/h.all.v7.0.symbols.gmt"
    C2_DIR="/stornext/Home/data/allstaff/w/whitfield.h/BIG_DATA/GeneSets/c2.all.v7.1.symbols.gmt"
  } else {
    Hallmark_Geneset_DIR="C:\\Users\\whitfield.h\\Desktop\\Projects\\B_I_G_Data\\GeneSets\\MSigDB\\msigdb.v7.2.symbols.gmt"
    C2_DIR="C:\\Users\\whitfield.h\\Desktop\\Projects\\B_I_G_Data\\GeneSets\\MSigDB\\c2.all.v7.1.symbols.gmt"
  }
  
  
  message("# Loading C2 genesets from MSigDB")
  conn <- file(C2_DIR,open="r")
  C2_db <- getGmt(conn)
  close(conn)
  
  MSigDB_signatures <- c()
  message("# Converting genesets to Vision objects")
  if (length(Sig_names_lst[grepl("DN", Sig_names_lst)]) > 0){
    
    x <- Sig_names_lst[grepl("DN", Sig_names_lst)]
    y <- unlist(strsplit(x, "_DN"))
    
    
    for (iDirGeneSet in y){
      
      print(paste0(" -> Extracting ", iDirGeneSet, " UP and DN..."))
      
      UP_str <- paste0(iDirGeneSet, "_UP", "")
      UP_genes <- geneIds(C2_db[names(C2_db) == UP_str])[[UP_str]]
      DN_str <- paste0(iDirGeneSet, "_DN", "")
      DN_genes <- geneIds(C2_db[names(C2_db) == DN_str])[[DN_str]]
      
      ### Make signature object
      UPDN_genes <- c(rep(1,length(UP_genes)), rep(-1, length(DN_genes)))
      names(UPDN_genes) <- c(UP_genes, DN_genes)
      iDirSet <- createGeneSignature(name = iDirGeneSet, sigData = UPDN_genes)
      
      MSigDB_signatures <- c(MSigDB_signatures, iDirSet)
    }
    
    Sig_names_lst <- Sig_names_lst[!(grepl("DN|UP", Sig_names_lst))]
  }
  
  for (iSig in names(C2_db)){
    print(paste0("-> Extracting ",iSig,""))
    iGeneIDs <- geneIds(C2_db[names(C2_db) %in% c(iSig)])[[iSig]]
    
    ### Make signature object
    iSig_genes <- rep(1,length(iGeneIDs))
    names(iSig_genes) <- iGeneIDs
    SigObject <- createGeneSignature(name = iSig, sigData = iSig_genes)
    
    MSigDB_signatures <- c(MSigDB_signatures, SigObject)
    
  }
  
  return(MSigDB_signatures)
}


##### ----- SCORING FUNCTIONS ------------------------------------------------


RunCookEMT_scDECAF <- function(x_dat, genesetList, genesetThresh=5, nVar=3000,standardise = TRUE){
  
  ### --- Get HVGs
  features_lst <- VariableFeatures(object = x_dat)
  
  ### --- Subset data
  data_df <- as.matrix(GetAssayData(object = x_dat, slot = 'scale.data'))
  data_df <- data_df[features_lst[features_lst %in% rownames(data_df)],]
  x_dat <- x_dat[features_lst[features_lst %in% rownames(x_dat)],]
  features_lst <- features_lst[features_lst %in% rownames(data_df)]
  
  ### --- Get Embedding
  x_dat <- RunUMAP(object = x_dat, features=features_lst, reduction.key = 'sctUMAP_', reduction.name = 'sctumap')
  cell_embed <- as.data.frame(x_dat@reductions$sctumap@cell.embeddings)
  
  ### --- Overlapping genes
  target <- scDECAF::genesets2ids(data_df[match(features_lst, rownames(data_df)),], genesetList) # genesetList=pb_markers
  genesetDrop <- unlist(lapply(as.vector(colSums(target)), function(x){ifelse(x<genesetThresh, FALSE, TRUE)}))
  target <- target[,genesetDrop]
  
  ### --- Run scDECAF
  message(paste0("#####  Running scDECAF "))
  message(paste0(" Standardise : ", standardise))
  message(paste0(" HVGs : ", as.character(length(features_lst))))
  message(paste0(" # of Genesets : ", dim(target)[2]))
  output <- scDECAF::scDECAF(data = data_df, gs = target, standardize = standardise, 
                             hvg = features_lst, k = 10, embedding = cell_embed, cca.k=10,
                             n_components = 10, max_iter = 2, thresh = 0.5)
  
  return(as.data.frame(attr(output,"raw_scores")))
}

RunCookEMT_VISION <- function(x_dat, genesetList){
  ### --- Remove missing genes in genesets
  genesetList_filt <- lapply(genesetList, function(x){x[x %in% rownames(x_dat)]})
  x_dat <- x_dat[intersect(rownames(x_dat), unique(as.vector(unlist(genesetList_filt)))),]
  
  ### --- Get non-log norm counts
  dat_obj <- SingleCellExperiment(assays=list(counts=as.matrix(GetAssayData(x_dat, "counts")), 
                                              logcounts=as.matrix(GetAssayData(x_dat, "data"))))
  dat_obj <- scater::logNormCounts(dat_obj, log=FALSE, exprs_values = "counts")
  VISION_counts <- normcounts(dat_obj)
  VISION_counts <- VISION_counts[intersect(rownames(VISION_counts), unique(as.vector(unlist(genesetList_filt)))),]
  genesetList_filt <- lapply(genesetList_filt, function(x){x[x %in% rownames(VISION_counts)]})
  
  ### --- VISION objs
  genesigs_obj <- MakeVisionObject(genesetList_filt, biDirection=FALSE)
  genesetNames <- unlist(lapply(genesigs_obj, function(x){attributes(x)$name}))
  VISION_Obj <- Vision(data = VISION_counts, signatures = genesigs_obj, pool=FALSE, min_signature_genes=10)
  
  ### --- Run VISION
  options(mc.cores = 2)
  VISION_RESULTS <- analyze(VISION_Obj)
  
  return(as.data.frame(getSignatureScores(VISION_RESULTS)))
  
}


##### ----- PLOTTING ------------------------------------------------------


GetLegend <- function(col_dict, lab_str, title_pos="top",
                      n_row=1, by_row=TRUE,rl=1){
  require(cowplot)
  legend_dat_themes <- data.frame(Sample=factor(names(col_dict), levels=names(time_cols_orig)), 
                                  prop=rep(1, length(col_dict)), colours = as.vector(col_dict))
  legend_bar <- ggplot(legend_dat_themes, aes(x=Sample, y=prop, fill=Sample)) + 
    geom_bar(position="dodge", stat="identity")+
    scale_fill_manual(values=col_dict)+theme_bw()+
    guides(fill=guide_legend(nrow=n_row,byrow=by_row, override.aes = list(size=4),title.position=title_pos), 
           color = guide_legend(override.aes = list(size = 4)),title.position=title_pos)+labs(fill=lab_str)
  legend_themes <- cowplot::get_legend(legend_bar)
  return(legend_themes)
}


tsnePanel <- function(melted_df_IZAR, tsne_df, gs_x, bar_width=5, bar_height=1, alpha_x=0.7, rl=1){
  dat_x <- melted_df_IZAR[melted_df_IZAR$geneset==gs_x,]
  tsne_df$gs_score <- dat_x$score[match(dat_x$cellID, tsne_df$cellID)]
  
  gg <- ggplot(aes(x=tsne_x, y=tsne_y, colour=gs_score), data = tsne_df)+geom_point(size=1, alpha=alpha_x)+
    theme_bw()+
    theme( legend.position = "bottom",legend.title = element_text(size=rel(rl)*0.6),
           legend.text  = element_text(size=rel(rl)*0.4), legend.margin = margin(-10,0,0,1),
           panel.grid.minor = element_blank(),panel.grid.major = element_blank(),
           axis.title.y = element_blank(), axis.title.x = element_blank(), 
           axis.text.y = element_blank(), axis.text.x = element_blank(), 
           axis.ticks.y=element_blank(), axis.ticks.x=element_blank(),
           axis.line.y = element_blank(),axis.line.x = element_blank())+scale_colour_viridis()+
    labs(colour="Score", title="")+ guides(colour = guide_colourbar(barwidth=bar_width,barheight=bar_height))
  return(gg)
}
