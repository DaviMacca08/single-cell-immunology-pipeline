# =========================================================
# Environment & packages
# =========================================================

suppressPackageStartupMessages({
  
  library(tidyverse)
  library(Seurat)
  library(SingleCellExperiment)
  library(harmony)
  library(DESeq2)
  
  library(celldex)
  library(SingleR)
  library(Azimuth)
  
  library(ggplot2)
  library(patchwork)
  library(pheatmap)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(clusterProfiler)
  library(msigdbr)
  library(ReactomePA)
  library(fgsea)
  
  library(CellChat)
  library(liana)
  
  library(matrixStats)
  
  library(slingshot)
  library(condiments)
  library(scater)
  library(scran)

  # GEO data acquisition & annotation
  library(GEOquery)
  library(hgu133plus2.db)
})