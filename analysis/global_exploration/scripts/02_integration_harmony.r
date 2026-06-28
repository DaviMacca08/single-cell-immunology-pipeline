# =========================================================
# Bulk integration workflow (Merge + PCA + Harmony)
# Input: filtered single-sample Seurat objects
# Output: integrated Seurat object + QC plots
# =========================================================


# =========================================================
#                  Libraries & Setup
# =========================================================

source("Setup_Environment/00_paths.R")
source("Setup_Environment/01_environment.R")
source("Setup_Environment/02_io_helpers.R")
source("Setup_Environment/03_checks.R")
source("Setup_Environment/04_seed.R")

set_seed(1234)

message("=== Starting Integration pipeline ===")


# =========================================================
#                  Load filtered samples
# =========================================================

message("Loading filtered Seurat objects...")

setwd(paths$filtered)

rds_files <- list.files(pattern = "\\.rds$")

if (length(rds_files) == 0) {
  stop("ERROR: No .rds files found in filtered directory")
}

sample_list <- list()

for (f in rds_files) {
  
  sample_name <- tools::file_path_sans_ext(f)
  
  message("Loading sample: ", sample_name)
  
  sample_list[[sample_name]] <- readRDS(f)
}

message("Loaded ", length(sample_list), " samples")


# =========================================================
#                  Sanity checks
# =========================================================

if (length(sample_list) < 2) {
  stop("ERROR: At least 2 samples required for integration")
}

cell_count <- tibble(
  sample = names(sample_list), 
  n_cells = sapply(sample_list, ncol)
  )

print(cell_count)


# =========================================================
#                  Merge samples
# =========================================================

message("Merging Seurat objects...")

seurat_merge <- merge(
  x = sample_list[[1]],
  y = sample_list[-1],
  add.cell.ids = names(sample_list),
  project = "IBD_gut_immune_scRNA"
)

check_seurat(seurat_merge)

message("Merge completed: ", ncol(seurat_merge), " total cells")


# =========================================================
#                  Normalization + HVGs
# =========================================================

message("Normalizing data...")

seurat_merge <- NormalizeData(seurat_merge, normalization.method = "LogNormalize", scale.factor = 10000)

message("Finding HVGs...")

seurat_merge <- FindVariableFeatures(seurat_merge, selection.method = "vst", nfeatures = 2000)

top10 <- head(VariableFeatures(seurat_merge), 10)


# =========================================================
#                  HVG plot
# =========================================================

hvg_plot <- VariableFeaturePlot(seurat_merge) + ggtitle("Highly Variable Genes (Top 10 Labeled)")
hvg_plot <- LabelPoints(plot = hvg_plot, points = top10, repel = TRUE)

open_pdf(filename = "HVG_feature_plot_labeled.pdf", dir = paths$plots_integration_hvg)

print(hvg_plot)

close_pdf()

# =========================================================
#                  Scaling
# =========================================================

message("Scaling data...")

seurat_merge <- ScaleData(
  seurat_merge,
  features = VariableFeatures(seurat_merge),
  vars.to.regress = "nCount_RNA"
)


# =========================================================
#                  PCA
# =========================================================

message("Running PCA...")

seurat_merge <- RunPCA(
  seurat_merge,
  features = VariableFeatures(seurat_merge)
)

print(seurat_merge[["pca"]], dims = 1:15, nfeatures = 10)


# ---- PCA diagnostics ----
# Loadings plot
pca_loadings <- VizDimLoadings(seurat_merge, dims = 1:2, reduction = "pca") + 
  ggtitle("Top Gene Loadings: PC1 - PC2")

save_plot(pca_loadings, filename = "PCA_LoadingsPlot.png", dir = paths$plots_integration_hvg)

#Embedding plot
pca_embedding <- DimPlot(seurat_merge, reduction = "pca") + NoLegend() +
  ggtitle("Cell Distribution in PCA Space")

save_plot(pca_embedding, filename = "PCA_EmbeddingPlot.png", dir = paths$plots_integration_hvg)


# Elbow plot
open_pdf(filename = "PCA_ElbowPlot.pdf", dir = paths$plots_integration_hvg, width = 9, height = 9)
print(ElbowPlot(seurat_merge, ndims = 50) + ggtitle("Elbow Plot"))
close_pdf()

# Heatmap 
message("Generating PC heatmaps (diagnostic mode)...")
for (i in 1:15) {
  
  DimHeatmap(
    seurat_merge, 
    dims = i, 
    balanced = TRUE, 
    fast = TRUE
  )
  
  message(paste("Showing PC", i))
  
  gc()
  
}


# =========================================================
#                  UMAP (pre-Harmony)
# =========================================================

message("Running UMAP (pre-Harmony)...")

seurat_merge <- RunUMAP(seurat_merge, dims = 1:30, reduction = "pca")

meta_vars <- c("sample_id", "donor", "disease", "condition", "compartment")

DimPlot_Before <- list()

for (v in meta_vars){
  
  # Check metadata exists
  if (!v %in% colnames(seurat_merge@meta.data)) {
    
    warning(paste(v, "not found in metadata"))
    next
  }
  
  # Check metadata has more than 1 level
  n_groups <- length(unique(seurat_merge@meta.data[[v]]))
  
  if (n_groups < 2) {
    
    warning(paste(v, "has less than 2 groups"))
    next
  }
  
  # doing plot
  DimPlot_Before[[v]] <- DimPlot(seurat_merge, group.by = v, reduction = "umap") +
    ggtitle(paste("Before Harmony -", v))
  
  # save plot
  save_plot(DimPlot_Before[[v]], filename = paste0("BeforeHarmony_", v, ".png"), dir = paths$plots_integration_hvg,
            width = 8, height = 6)
  
}


# =========================================================
#                  Harmony integration
# =========================================================

message("Running Harmony integration...")

cd_harmony <- RunHarmony(
  seurat_merge,
  group.by.vars = "sample_id",
  plot_convergence = FALSE
)

message("Harmony completed")


# =========================================================
#                  UMAP + clustering (Harmony space)
# =========================================================

cd_harmony <- RunUMAP(cd_harmony, reduction = "harmony", dims = 1:30)

cd_harmony <- FindNeighbors(cd_harmony, reduction = "harmony", dims = 1:30)

cd_harmony <- FindClusters(cd_harmony, resolution = 0.6)

message("Clustering completed")

clusters <- DimPlot(cd_harmony, reduction = "umap", group.by = "seurat_clusters", label  = TRUE, repel = TRUE) + NoLegend()

save_plot(clusters, filename = "DimPlot_clusters.png", dir = paths$plots_integration_hvg)

# =========================================================
#                  Post-Harmony plots
# =========================================================

DimPlot_After <- list()
  
for (v in meta_vars){
  
  # Check metadata exists
  if (!v %in% colnames(cd_harmony@meta.data)) {
    
    warning(paste(v, "not found"))
    next
  }
  
  # Check metadata has more than 1 level
  if (length(unique(cd_harmony@meta.data[[v]])) < 2) {
    
    warning(paste(v, "has only one group"))
    next
  }
  
  # doing plot
  DimPlot_After[[v]] <- DimPlot(cd_harmony, group.by = v, reduction = "umap") +
    ggtitle(paste("After Harmony -", v))
  
  # save plot
  save_plot(DimPlot_After[[v]], filename = paste0("AfterHarmony_", v, ".png"), dir = paths$plots_integration_hvg,
            width = 8, height = 6)
  
}


# =========================================================
#                  Before vs After comparison
# =========================================================

for (v in meta_vars) {
  
  if (!v %in% names(DimPlot_Before) || !v %in% names(DimPlot_After)) {
    warning("Missing plot for: ", v)
    next
  }
  
  p <- DimPlot_Before[[v]] + DimPlot_After[[v]]
  
  save_plot(p, filename = paste0(v, "_before_vs_after_harmony.png"), dir = paths$plots_integration_hvg,
            width = 20, height = 12)
}


# =========================================================
#                  Save outputs
# =========================================================

save_rds(cd_harmony, filename = "srt_obj_merge_harmony.rds", dir = paths$objects_global)
save_rds(seurat_merge, filename = "srt_obj_merge_pca.rds", dir = paths$objects_global)

message("Integration pipeline completed successfully")


# =========================================================
#                  Session information
# =========================================================

message("Saving session information (Harmony integration stage)...")

session_file <- file.path(paths$logs,"sessionInfo_harmony.txt")

writeLines(capture.output(sessionInfo()),session_file)

save_session_info(dir = paths$logs, filename = "sessionInfo_harmony.txt", label = "Harmony integration stage")

message("Session information saved at: ", paths$logs)


# =========================================================
#                  Final pipeline message
# =========================================================

message("=================================================")
message("PIPELINE STEP COMPLETED: Integration (Harmony)")
message("=================================================")