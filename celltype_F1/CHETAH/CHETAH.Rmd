---
  title: "CHETAH"
output: pdf_document
---
  
  ```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

# Pancreas using Seger as reference and baron as query

```{r}
library(scRNAseq)
library(CHETAH)
library(scran)
library(scater)

# Preprocessing Baron and Seger
baron <- BaronPancreasData('human')
seger <- SegerstolpePancreasData()

clusters <- quickCluster(baron)
baron <- computeSumFactors(baron, clusters=clusters)
summary(sizeFactors(baron))
baron <- logNormCounts(baron)
baron <- runTSNE(baron)

clusters <- quickCluster(seger)
seger <- computeSumFactors(seger, clusters=clusters)
summary(sizeFactors(seger))
seger <- logNormCounts(seger)

label_isnot_na <- which(!is.na(seger$`cell type`))
seger <- seger[, label_isnot_na]
seger$celltypes <- seger$`cell type`

# Running CHETAH
input <- CHETAH::CHETAHclassifier(input = baron, ref_cells = seger)

```
