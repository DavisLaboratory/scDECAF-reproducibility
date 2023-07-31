######################################################################################################
##
##  PROCESS IZAR DATA 
##
## 

## -- Load packages & functions
source("./helper_functions.R")

saveFolder=""

##### ----- LOAD DATASET -----------------------------------------------

## From: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE146026
IZAR_PATH="/stornext/Home/data/allstaff/w/whitfield.h/data/_datasets/sc/GSE146026_Izar_HGSOC_ascites_10x_log.tsv"

### Load log counts
Izar <- read.table(IZAR_PATH, sep="\t", row.names=1, header=TRUE, stringsAsFactors = FALSE)

## Get cell meta data
ColData <- Izar[1:7, ]
ColData <- t(ColData)

## Get log counts
LogCounts <- Izar[8:dim(Izar)[1], ]

## Make SCE object
Izar_sce <- SingleCellExperiment(list(logcounts=as(as.matrix(LogCounts), "dgCMatrix")), colData=DataFrame(ColData))

### From Supp Table -- Set cell names
Cluster_numb <-  c("21", "20", "17", "4", "7",
                   "9", "3", "2", "13",
                   "11", "10", "5", "1", 
                   "18", "19",
                   "8", "15", "16", 
                   "6", "12", "14")
CellType <- c(rep("Malignant", 5), 
              rep("Fibroblasts",4), 
              rep("Macrophage", 4), 
              rep("DC",2), 
              "Bcells", "Tcells", "Erythrocytes", 
              rep("NotAssigned", 3))
ClusterNames <- setNames(c(CellType, "-1"), c(Cluster_numb, "-1"))
Izar_sce$CellType <- ClusterNames[as.vector(Izar_sce$clst)]

### Save data
save(Izar_sce, file=paste0(saveFolder, "Izar_sce.Rdata"))


##### ----- SCORE DATASET w/ PRUNING -------------------------------------------


### --- Load Genesets
C2_genesets <- GetMSigDB(database="C2")
C2_genesets_filt <- C2_genesets[lapply(C2_genesets, length) > 10]

### --- Get gene overlap
dat_obj <- Izar_sce[intersect(rownames(Izar_sce), unique(as.vector(unlist(C2_genesets)))),]

### --- Get HVGs
dat_varModel <-  scran::modelGeneVar(dat_obj, block=dat_obj$CellType, density.weights=FALSE)
features_lst <- scran::getTopHVGs(dat_varModel , n=3000)

### --- Subset data
dat_obj <- dat_obj[features_lst[features_lst %in% rownames(dat_obj)],]
features_lst <- features_lst[features_lst %in% rownames(dat_obj)]

### --- Cell embeddings
dat_obj <- scater::runUMAP(dat_obj, scale=FALSE,  exprs_values = "logcounts")
cell_embed = as.data.frame(reducedDim(dat_obj, "UMAP"))

### --- Remove missing genes in genesets
genes_to_keep <- intersect(rownames(dat_obj), unique(as.vector(unlist(C2_genesets))))
C2_genesets <- C2_genesets[lapply(C2_genesets, length) > 10]
C2_genesets_filt <- lapply(C2_genesets, function(x){x[x %in% rownames(dat_obj)]})


### --- Prune genesets
pruned_genesets <- scDECAF::pruneGenesets(data=as.matrix(logcounts(Izar_sce)), 
                                          genesetlist = C2_genesets_filt,
                                          hvg=features_lst,
                                          embedding =cell_embed,
                                          min_gs_size =10,lambda = exp(-2),
                                          suppress_plot=TRUE)

pruned_genesets <- as.character(pruned_genesets)
# From 5529 genesets down to 91

### --- Run scDECAF
data_df <- as.matrix(logcounts(dat_obj))
target <- scDECAF::genesets2ids(data_df[match(features_lst, rownames(data_df)),], C2_genesets_filt[pruned_genesets]) # genesetList=pb_markers

genesetThresh=5
standardise=TRUE
genesetDrop <- unlist(lapply(as.vector(colSums(target)), function(x){ifelse(x<genesetThresh, FALSE, TRUE)}))
target <- target[,genesetDrop]
output <- scDECAF::scDECAF(data = data_df, gs = target, standardize = standardise, 
                           hvg = features_lst, k = 10, embedding = cell_embed, cca.k=10,
                           n_components = 10, max_iter = 2, thresh = 0.5)
IZAR_scDECAFscores_pruned <- attr(output,"raw_scores")


### --- Tidy score dataframe
IZAR_scDECAFscores_pruned <- as.data.frame(IZAR_scDECAFscores_pruned)
IZAR_scDECAFscores_pruned$cellID <- rownames(IZAR_scDECAFscores_pruned)
melted_df_IZAR <- melt(IZAR_scDECAFscores_pruned, value.name="score", variable.name = "geneset",id.vars="cellID")

CellTypeID <- Izar_sce$CellType
names(CellTypeID) <- colnames(Izar_sce)
ClustID <- Izar_sce$clst
names(ClustID) <- colnames(Izar_sce)

melted_df_IZAR$CellType <- as.vector(CellTypeID[melted_df_IZAR$cellID])
melted_df_IZAR$Clust <- as.vector(ClustID[melted_df_IZAR$cellID])

save(melted_df_IZAR, file=paste(saveFolder,"Izar_scores.Rdata", sep=""))
