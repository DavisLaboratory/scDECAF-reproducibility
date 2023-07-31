######################################################################################################
##
##  PROCESS EMT DATA 
##
## 


## -- Load packages & functions
source("./helper_functions.R")


##### ----- LOAD DATASET -----------------------------------------------

## Download data --------------------------------------
EMT_SAVE_DIR <- "/stornext/Home/data/allstaff/w/whitfield.h/data/_datasets/sc/scEMT_3cellLines/raw_data/"
SCE_SAVE_DIR <- "/stornext/Home/data/allstaff/w/whitfield.h/data/_datasets/sc/scEMT_3cellLines/"

ftp_address = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE147nnn/GSE147405/suppl/"
cellLine_files = c("GSE147405_A549_EGF_TimeCourse","GSE147405_A549_TGFB1_TimeCourse","GSE147405_A549_TNF_TimeCourse",
                   "GSE147405_DU145_EGF_TimeCourse", "GSE147405_DU145_TGFB1_TimeCourse","GSE147405_DU145_TNF_TimeCourse",
                   "GSE147405_MCF7_EGF_TimeCourse", "GSE147405_MCF7_TGFB1_TimeCourse","GSE147405_MCF7_TNF_TimeCourse",
                   "GSE147405_OVCA420_EGF_TimeCourse", "GSE147405_OVCA420_TGFB1_TimeCourse","GSE147405_OVCA420_TNF_TimeCourse")
cellLine_file_types = c("_UMI_matrix.csv.gz","_metadata.csv.gz")

count_files = c("GSE147405_TimeCourse_Mix1", "GSE147405_TimeCourse_Mix2", 
                "GSE147405_TimeCourse_Mix3a", "GSE147405_TimeCourse_Mix3b",
                "GSE147405_TimeCourse_Mix4a", "GSE147405_TimeCourse_Mix4b")
count_file_types = c("_barcode_annotations.csv.gz", "_barcode_counts.csv.gz")

## Get cell-line files
for (i_file in cellLine_files){
  download.file(url=paste0(ftp_address,i_file,cellLine_file_types[[1]]),
                destfile=paste0(EMT_SAVE_DIR,i_file,cellLine_file_types[[1]]))
  download.file(url=paste0(ftp_address,i_file,cellLine_file_types[[2]]),
                destfile=paste0(EMT_SAVE_DIR,i_file, cellLine_file_types[[2]]))
}

## Get count files
for (i_file in count_files){
  download.file(url=paste0(ftp_address,i_file,count_file_types[[1]]),
                destfile=paste0(EMT_SAVE_DIR,i_file,count_file_types[[1]]))
  download.file(url=paste0(ftp_address,i_file,count_file_types[[2]]),
                destfile=paste0(EMT_SAVE_DIR,i_file, count_file_types[[2]]))
}

cellLine_names = unlist(lapply(strsplit(cellLine_files, "_"), function(x){paste(x[2:3], collapse="_")}))
names(cellLine_names) <- cellLine_files


##### ----- PROCESS DATASET -----------------------------------------------

## Chosen for their high number of cells & clearer trajectories in original paper
trimmed_cellLine_lst <- cellLine_files[grepl("TGFB1", cellLine_files)]
cellLine_names = unlist(lapply(strsplit(trimmed_cellLine_lst, "_"), function(x){paste(x[2:3], collapse="_")}))
names(cellLine_names) <- trimmed_cellLine_lst

features_n=3000
EMT_sce_lst <- list()
EMT_sce_lst_scTransform <- list()
for (iSample in trimmed_cellLine_lst){
  
  ## READ DATA -------------------------
  message(paste0("Reading in ",iSample,"..."))
  dt = fread(paste0(EMT_SAVE_DIR, iSample, cellLine_file_types[[1]]))
  dt_meta = fread(paste0(EMT_SAVE_DIR, iSample, cellLine_file_types[[2]]))
  
  dt_rownames <- dt$V1
  dt <- dt[,2:ncol(dt)]
  dt <- as.matrix(sapply(dt, as.numeric))  
  rownames(dt) <- dt_rownames
  
  ## MAKE SCE & FILTER -------------------------
  EMT_obj <- SingleCellExperiment(list(counts=dt)) # as(dt, "dgCMatrix")
  colData(EMT_obj) <- DataFrame(dt_meta)
  EMT_obj$CellID <- EMT_obj$V1
  colnames(EMT_obj) <- EMT_obj$CellID
  EMT_obj <- EMT_obj[as.vector(rowSums(counts(EMT_obj))>0),]
  EMT_obj <- EMT_obj[!(duplicated(rownames(EMT_obj))),]
  
  EMT_sce_lst[[cellLine_names[[iSample]]]] <- EMT_obj
  
  ## MAKE SEURAT OBJECT -------------------------
  EMT_obj <- as.Seurat(EMT_obj, counts="counts", data=NULL)
  
  ## NORMALISE -------------------------
  EMT_obj <- SCTransform(object = EMT_obj, assay="originalexp", verbose = FALSE, variable.features.n = features_n)
  
  ## SAVE -------------------------
  EMT_sce_lst_scTransform[[cellLine_names[[iSample]]]] <- EMT_obj
  save(EMT_sce_lst_scTransform, file=paste0(saveFolder, "CookEMT_TGFB1List.Rdata"))
  save(EMT_sce_lst, file=paste0(saveFolder, "CookEMT_TGFB1List_sce.Rdata"))
}

EMT_sce_lst_scTransform <- lapply(EMT_sce_lst_scTransform, function(x){RenameAssays(object = x, originalexp = 'RNA')})

save(EMT_sce_lst_scTransform, file=paste0(saveFolder, "CookEMT_TGFB1List.Rdata"))


##### ----- SCORE DATASET w/ C2 -----------------------------------------------


### --- Load data
load(file=paste0(saveFolder, "CookEMT_TGFB1List.Rdata")) # EMT_sce_lst_scTransform

C2_genesets <- GetMSigDB(database="C2")
C2_genesets_filt <- C2_genesets[lapply(C2_genesets, length) > 10]


CookEMT_TGFB1_sctransform_C2filt_score_lst <- list()

for (x_name in names(EMT_sce_lst_scTransform)){
  x_dat <- EMT_sce_lst_scTransform[[x_name]]
  
  message(paste0("Running VISION for : ", x_name))
  CookEMT_TGFB1_sctransform_C2filt_score_lst[[paste0("VISION_","sctransform","_allgenes_", x_name)]] <- RunCookEMT_VISION(x_dat,C2_genesets_filt)
  message(paste0("Running scDECAF for : ", x_name))
  CookEMT_TGFB1_sctransform_C2filt_score_lst[[paste0("scDECAF_","sctransform","_","3k_",x_name)]] <- RunCookEMT_scDECAF(x_dat, C2_genesets_filt)
  message("Save intermediate...")
  save(CookEMT_TGFB1_sctransform_C2filt_score_lst, 
       file=paste(saveFolder,"CookEMT_TGFB1scores_ALLcellLines.Rdata", sep=""))
}

## Merge scores
melted_score_lst <- lapply(names(CookEMT_TGFB1_sctransform_C2filt_score_lst), function(x_name){
  message(paste0("Processing ", x_name))
  info_vec <- strsplit(x_name,"_")[[1]]
  df_x <- as.data.frame(CookEMT_TGFB1_sctransform_C2filt_score_lst[[x_name]])
  df_x$cellID <- rownames(df_x)
  
  melted_df_x <- melt(df_x, value.name="score", id.vars="cellID", variable.name = "geneset")
  n_rows <- dim(melted_df_x)[1]
  melted_df_x$method <- rep(info_vec[1], n_rows)
  melted_df_x$norm_method <- rep(info_vec[2], n_rows)
  melted_df_x$n_genes <- rep(info_vec[3], n_rows)
  melted_df_x$cellLine <- rep(info_vec[4], n_rows)
  melted_df_x$stim <- rep(info_vec[5], n_rows)
  return(melted_df_x)
})

names(melted_score_lst) <- names(CookEMT_TGFB1_sctransform_C2filt_score_lst)
CookedEMT_TGFB1_melt <- rbindlist(melted_score_lst, use.names = TRUE)

## ADD COLDATA ----------------------
load(file=paste0("/stornext/Home/data/allstaff/w/whitfield.h/data/_datasets/sc/scEMT_3cellLines/",
                 "CookEMT_colData.Rdata")) # colData_comb

time_dict <- colData_comb$time
names(time_dict) <- colData_comb$CellID
pseudotime_dict <- colData_comb$pseudotime
names(pseudotime_dict) <- colData_comb$CellID

CookedEMT_TGFB1_melt$time <- as.vector(time_dict[CookedEMT_TGFB1_melt$cellID])
CookedEMT_TGFB1_melt$pseudotime <- as.vector(pseudotime_dict[CookedEMT_TGFB1_melt$cellID])

save(CookedEMT_TGFB1_melt, 
     file=paste(saveFolder,"CookEMT_TGFB1scores_ALLcellLines.Rdata", sep=""))



##### ----- PSEUDOBULK DATASET -----------------------------------------------

library(scMerge)

load(file=paste0(saveFolder, "CookEMT_TGFB1List_sce.Rdata")) # EMT_sce_lst

### -- Merge
cols_to_keep <- colnames(colData(sce_lst[[1]]))[2:ncol(colData(sce_lst[[1]]))]
EMT_sce <- scMerge::sce_cbind(EMT_sce_lst, 
                     method="union",
                     cut_off_batch=0,
                     cut_off_overall=0,
                     exprs = c("counts"),
                     batch_names = as.vector(cellLine_names), 
                     colData_names=cols_to_keep)

EMT_sce$Sample <- paste0(EMT_sce$CellLine, "_TGFB1_", EMT_sce$time)

## Drop OCTA420 due to low cell numbers
sort(table(EMT_sce$Sample))

EMT_sce <- EMT_sce[,!(EMT_sce$CellLine=="OVCA420")]
EMT_sce$time <- replace(EMT_sce$time, EMT_sce$time %in% c("1d_rm", "3d_rm"), "1-3d_rm")
EMT_sce$time <- replace(EMT_sce$time, EMT_sce$time %in% c("3d", "7d"), "3-7d")

### --- Aggregate
Summed_Sce <- aggregateAcrossCells(EMT_sce, id=EMT_sce$Sample)
DGE_Summed_sce <- DGEList(counts(Summed_Sce))
discarded <- isOutlier(DGE_Summed_sce$samples$lib.size, log=TRUE,
                       type="lower")
DGE_Summed_sce <- DGE_Summed_sce[,!discarded]
keep = rowMeans(edgeR::cpm(y=DGE_Summed_sce, log = TRUE) >= 0.5) >= 0.1 
DGE_Summed_sce <- DGE_Summed_sce[keep, ]
DGE_Summed_sce <- calcNormFactors(DGE_Summed_sce)

DATA_logCPM <- edgeR::cpm(DGE_Summed_sce, log=TRUE, prior.count = 1, normalized.lib.sizes=TRUE, lib.size=DGE_Summed_sce$samples$lib.size) 
DGE_Summed_sce$logCPM <- DATA_logCPM

### --- Save
DGE_Summed_sce$samples$CellLine <- unlist(lapply(strsplit(rownames(DGE_Summed_sce$samples), "_"), function(x){x[[1]]}))
DGE_Summed_sce$samples$Time <- unlist(lapply(strsplit(rownames(DGE_Summed_sce$samples), "_"), function(x){paste0(x[3:length(x)], collapse="_")}))
save(DGE_Summed_sce, file=paste0(DATA_DIR, "CookedEMT_DGE_all.Rdata"))




