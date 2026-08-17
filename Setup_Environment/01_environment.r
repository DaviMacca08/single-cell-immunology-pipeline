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
  
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(clusterProfiler)
  library(msigdbr)
  library(ReactomePA)
  library(fgsea)
  
  library(CellChat)
  library(liana)
  
})