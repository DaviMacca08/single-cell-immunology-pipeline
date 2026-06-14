# =========================================================
# Consensus Cell Type Annotation Layer
# Input:
#   - Seurat object (post-Harmony clustering)
#   - Cluster vs SingleR table (CellAtlas and MonacoImmune)
#   - Cluster vs Azimuth table (PBMC reference)
#
# Output:
#   - Seurat object with consensus annotation
#   - Cluster-level annotation table
#   - Conflict report table
#   - Diagnostic plots
# =========================================================


# =========================================================
#                  Libraries & Setup
# =========================================================

source("scripts/00_setup/00_paths.R")
source("scripts/00_setup/01_environment.R")
source("scripts/00_setup/02_io_helpers.R")
source("scripts/00_setup/04_seed.R")

set_seed(1234)

message("Starting consensus cell type annotation...")


# =========================================================
#                   Load all input data
# =========================================================

message("Loading Seurat object...")

cd_harmony <- readRDS(
  file.path(paths$objects, "srt_obj_annotated_SingleR_Azimuth.rds")
)

message("Loading cluster-level annotation tables...")

cluster_singler_Monaco <- read.csv(file.path(paths$tables_singler, "MonacoImmune_MajorityVote.csv"))

cluster_singler_CellAtlas <- read.csv(file.path(paths$tables_singler, "CellAtlas_MajorityVote.csv"))

cluster_azimuth_PBMC <- read.csv(file.path(paths$tables_azimuth, "PBMC_Azimuth_MajorityVote.csv"))


# =========================================================
#                  Initial sanity checks
# =========================================================

message(" Running sanity checks...")

stopifnot(inherits(cd_harmony, "Seurat"))
stopifnot("seurat_clusters" %in% colnames(cd_harmony@meta.data))

stopifnot(!is.null(cluster_singler_CellAtlas))
stopifnot(!is.null(cluster_singler_Monaco))
stopifnot(!is.null(cluster_azimuth_PBMC))


# =========================================================
#                  Standardize cluster IDs
# =========================================================

message("Standardize cluster identifiers...")

cluster_singler_CellAtlas$cluster <- as.character(cluster_singler_CellAtlas$cluster)
cluster_singler_Monaco$cluster <- as.character(cluster_singler_Monaco$cluster)
cluster_azimuth_PBMC$cluster <- as.character(cluster_azimuth_PBMC$cluster)

cluster_ids <- as.character(cd_harmony$seurat_clusters)


# =========================================================
#             Build cluster-level annotation table
# =========================================================

message("Building annotation comparison table...")

cellatlas_tbl <- cluster_singler_CellAtlas |> 
  dplyr::select(cluster, label) |> 
  dplyr::rename(cellatlas = label)

monaco_tbl <- cluster_singler_Monaco |> 
  dplyr::select(cluster, label) |> 
  dplyr::rename(monaco = label)

azimuth_tbl <- cluster_azimuth_PBMC |> 
  dplyr::select(cluster, label) |> 
  dplyr::rename(azimuth = label)

cluster_annotations <- cellatlas_tbl |> 
  left_join(monaco_tbl, by = "cluster") |> 
  left_join(azimuth_tbl, by = "cluster")


# =========================================================
#            Add cluster-level annotations
# =========================================================

message(" Mapping CellAtlas majority labels...")

cellatlas_map <- setNames(
  cluster_singler_CellAtlas$label,
  cluster_singler_CellAtlas$cluster
)

cd_harmony$cellatlas_majority <- unname(cellatlas_map[cluster_ids])

message("Mapping MonacoImmune majority labels...")

monacoimmune_map <- setNames(
  cluster_singler_Monaco$label,
  cluster_singler_Monaco$cluster
)

cd_harmony$monacoimmune_majority <- unname(monacoimmune_map[cluster_ids])

message("Mapping Azimuth PBMC majority labels...")

pbmc_map <- setNames(
  cluster_azimuth_PBMC$label,
  cluster_azimuth_PBMC$cluster
)

cd_harmony$pbmc_majority <- unname(pbmc_map[cluster_ids])


# =========================================================
#                Diagnostic UMAP Plots 
# =========================================================

message("Generating diagnostics UMAP plots...")

cluster_plots <- DimPlot(cd_harmony, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE) + NoLegend() +
  ggtitle("Seurat clusters")

cellatlas_plot <- DimPlot(cd_harmony, reduction = "umap", group.by = "cellatlas_majority", label = TRUE, repel = TRUE) + NoLegend() + 
  ggtitle("CellAtlas Mapping") 

monaco_plot <- DimPlot(cd_harmony, reduction = "umap", group.by = "monacoimmune_majority", label = TRUE, repel = TRUE) + NoLegend() +
  ggtitle("MonacoImmune Mapping")


pbmc_plot <- DimPlot(cd_harmony, reduction = "umap", group.by = "pbmc_majority", label = TRUE, repel = TRUE) + NoLegend() + 
  ggtitle("Azimuth PBMC Mapping")

# -----------------------
#  Save diagnostic plots
# -----------------------

message("Saving diagnostic UMAP comparisons...")

save_plot(cellatlas_plot + cluster_plots, filename = "Cluster_CellAtlas.pdf", dir = paths$plot_cellannotation, 
          width = 14, height = 6)

save_plot(monaco_plot + cluster_plots, filename = "Cluster_MonacoImmune.pdf", dir = paths$plot_cellannotation, 
          width = 14, height = 6)

save_plot(pbmc_plot + cluster_plots, filename = "Cluster_PBMC.pdf", dir = paths$plot_cellannotation, 
          width = 14, height = 6)


# =========================================================
#           Final consensus assignment
# =========================================================

message("Defining final consensus cell type labels...")
message("Using Azimuth PBMC as primary reference for final assignment...")

final_map <- setNames(
  cluster_annotations$azimuth,
  cluster_annotations$cluster
)

cd_harmony$celltype_final <- unname(final_map[
  as.character(cd_harmony$seurat_clusters)
])

Idents(cd_harmony) <- "celltype_final"


# =========================================================
#             Final consensus visual validation
# =========================================================

message("Generating final annotation plots...")

final_plot <- DimPlot(cd_harmony, reduction = "umap", label = TRUE, repel = TRUE) +
  NoLegend() + 
  ggtitle("Final Cell Type Annotation")

message("Set canonical markers...")

canonical_markers <- c(
  # Pan-T cells
  "CD3D", "CD3E", "TRAC",
  
  # CD4 naïve / memory
  "IL7R", "CCR7", "LEF1", "TCF7", "LTB", "SELL",
  
  # CD8 cytotoxic
  "CD8A", "CD8B",
  "GZMK", "NKG7", "GZMB", "PRF1", "GNLY",
  
  # Treg
  "FOXP3", "IL2RA", "CTLA4", "IKZF2",
  
  # Exhaustion
  "PDCD1", "LAG3", "TIGIT", "HAVCR2", "TOX",
  
  # Tissue-resident memory
  "ITGAE", "CD69", "CXCR6", "ZNF683", "RUNX3",
  
  # γδ T cells
  "TRDC", "TRGC1", "TRGC2",
  
  # Inflammatory axis (validation only)
  "IFNG", "TNF", "RORC"
)

message("Generating canonical marker DotPlot...")

dotplot_final <- DotPlot(cd_harmony, features = canonical_markers, scale = FALSE) + RotatedAxis() + 
  ggtitle("Canonical marker expression across cell types")

# -------------------
#  Save final plots
# -------------------

save_plot(final_plot, filename = "UMAP_Final_Annotation.png", dir = paths$plot_cellannotation,
          width = 12, height = 10)

save_plot(final_plot, filename = "UMAP_Final_Annotation.pdf", dir = paths$plot_cellannotation,
          width = 12, height = 10)

save_plot(dotplot_final, filename = "DotPlot_CanonicalMarkers.pdf", dir = paths$plot_cellannotation,
          width = 14, height = 10)


# =========================================================
#                Save outputs
# =========================================================

message(" Saving consensus annotation table...")

save_csv(cluster_annotations, filename = "cluster_consensus_annotation.csv", dir = paths$tables)

message("Saving annotated Seurat object...")

save_rds(cd_harmony, filename = "srt_consensus_celltype.rds", dir = paths$objects)

message("Consensus celltype annotation pipeline completed successfully")

# =========================================================
#                  Save session info
# =========================================================

message("Saving session information for consensus celltype annotation...")

save_session_info(filename = "sessionInfo_Consensus_CellType.txt", dir = paths$logs, label = "Consensus celltype annotations")

message("Session information saved at: ", paths$logs)


# =========================================================
#                  Final pipeline message
# =========================================================

message("=================================================")
message("PIPELINE STEP COMPLETED: Consensun Celltype Annotations")
message("=================================================")
