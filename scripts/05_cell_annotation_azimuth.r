# =========================================================
# Azimuth Cell Type Annotation
# Input: Harmony-integrated Seurat object with singleR annotation metadata
# Output: Final Annotated Seurat object + diagnostic plots
# =========================================================


# =========================================================
#                  Libraries & Setup
# =========================================================

source("scripts/00_setup/00_paths.R")
source("scripts/00_setup/01_environment.R")
source("scripts/00_setup/02_io_helpers.R")
source("scripts/00_setup/04_seed.R")

set_seed(1234)

message("=== Start Azimuth Annotation ===")


# =========================================================
#      Load data (Seurat object - SingleR Annotated)
# =========================================================

message("Loading Seurat object...")

cd_harmony <- readRDS(
  file.path(paths$objects, "srt_obj_annotated_SingleR.rds")
)


# =========================================================
#          Load Azimuth referenfeces 
# =========================================================

message("Loading references data for Azimuth...")

pbmc_ref <- readRDS(
  file.path(paths$Azimuth_ref, "ref.Rds")
)

# =========================================================
#                  Initial sanity checks
# =========================================================

stopifnot(inherits(cd_harmony, "Seurat"))

stopifnot("seurat_clusters" %in% colnames(cd_harmony@meta.data))

stopifnot(inherits(pbmc_ref, "Seurat"))

message("Seurat object and reference for Azimuth loaded successfully")


# =========================================================
#                   Join layers 
# =========================================================

if ("RNA" %in% names(cd_harmony@assays)) {
  
  message("Joining RNA layers...")
  
  cd_harmony <- JoinLayers(
    object = cd_harmony,
    assay = "RNA",
    layers = "data")
  
} else {
  
  stop("RNA assay not found")
  
}


# =========================================================
#                  Set active identities
# =========================================================

Idents(cd_harmony) <- "seurat_clusters"

message("Active identity class: seurat_clusters")


# =========================================================
#                   Check default assay
# =========================================================

if (DefaultAssay(cd_harmony) != "RNA") {
  
  message(
    "Changing DefaultAssay from ",
    DefaultAssay(cd_harmony),
    " to RNA"
  )
  
  DefaultAssay(cd_harmony) <- "RNA"
  
} else {
  
  message("Default assay already set to RNA")
  
}


# =========================================================
#                 Check of reference
# =========================================================

message("Azimuth reference reductions.")

print(Reductions(pbmc_ref))

# ---------------------------------------------------------
# Distributions of labels
# ---------------------------------------------------------

message("Distributions of labels in cell type level 1: ")

print(table(pbmc_ref$celltype.l1))

message("Distributions of labels in cell type level 2: ")

print(table(pbmc_ref$celltype.l2))

message("Distributions of labels in cell type level 3: ")

print(table(pbmc_ref$celltype.l3))

# =========================================================
#                  Metadata overview
# =========================================================

message("Preview metadata (cd_harmony):")

print(head(cd_harmony@meta.data))

message("Preview metadata (pbmc_ref):")

print(head(pbmc_ref@meta.data))


# =========================================================
#               DimPlots of Azimuth PBMC cells
# =========================================================

level_celltype <- c("celltype.l1", "celltype.l2", "celltype.l3")

for (celltype in level_celltype) {
    p <- DimPlot(pbmc_ref, reduction = "refUMAP", group.by = celltype, label = TRUE, repel = TRUE, raster = TRUE) + 
    NoLegend() + 
    ggtitle(paste0("Cell Type Level ", gsub("\\D", "", celltype)))
  
  ggsave(
    filename = file.path(paths$plot_azimuth, paste0("PBMC_CellType_Level", gsub("\\D", "", celltype), ".png")),
    plot = p,
    width = 14,
    height = 10,
    dpi = 300
  )
  
  message("Saved plot: ", paste0("PBMC_CellType_Level", gsub("\\D", "", celltype), ".png"))
  
  rm(p)
  
}


# =========================================================
#                Prepare query for Azimuth
# =========================================================
# NOTE: Azimuth is based on SCTransform-derived reference space.
# The query object is prepared as a clean RNA-based Seurat object
# to avoid conflicts with precomputed embeddings (e.g., Harmony).

message("Preparing query object for Azimuth annotation...")

colnames_beforeAzimuth <- colnames(cd_harmony@meta.data)

query_az <- cd_harmony

query_az <- DietSeurat(
  query_az,
  assays = "RNA",
  counts = TRUE,
  data = TRUE,
  scale.data = FALSE
)

query_az@reductions <- list()

message("Merging RNA layers for compatibility with Azimuth workflow...")

query_az <- JoinLayers(query_az)

# ---------------------------------------------------------
# Check before Azimuth
# ---------------------------------------------------------

if (all(colnames(cd_harmony) %in% colnames(query_az)) &&
    all(colnames(query_az) %in% colnames(cd_harmony))) {
  
  message("Query data built correctly: cell names are consistent.")
  
} else {
  
  message("ERROR: cell name mismatch between cd_harmony and query_az.")
  stop("Azimuth aborted due to inconsistent cell identifiers.")
}


# =========================================================
#                     Run Azimuth
# =========================================================

message("Running Azimuth cell type annotation...")

query_az <- Azimuth::RunAzimuth(query_az, reference = paths$Azimuth_ref)

colnames_afterAzimuth <- colnames(query_az@meta.data)

new_meta <- setdiff(colnames_afterAzimuth, colnames_beforeAzimuth)

message("New metadata columns generated by Azimuth:")

print(new_meta)

# =========================================================
#      Transfer Azimuth annotations to main object
# =========================================================
# We explicitly transfer only newly created metadata columns
# to maintain separation between analysis (cd_harmony)
# and annotation (query_az)

message("Transferring Azimuth annotations to main Seurat object...")

for (meta in new_meta) {
  
  cd_harmony@meta.data[[meta]] <- query_az@meta.data[[meta]]
  
  message(paste0("Added metadata column: ", meta))
  
}


# =========================================================
#                    Diagnostic summary
# =========================================================
# Inspect annotation output and label distribution quality.

message("Generating diagnostic summary of Azimuth annotations...")

# ---------------------------------------------------------
# Display unique predicted labels
# ---------------------------------------------------------

message("Found labels for first level")

unique(cd_harmony$predicted.celltype.l1)

print(table(cd_harmony$seurat_clusters, cd_harmony$predicted.celltype.l1), useNA = "ifany")

message("Found labels for second level")

unique(cd_harmony$predicted.celltype.l2, useNA = "ifany")

print(table(cd_harmony$seurat_clusters, cd_harmony$predicted.celltype.l2), useNA = "ifany")

message("Found labels for third level")

unique(cd_harmony$predicted.celltype.l3)

print(table(cd_harmony$seurat_clusters, cd_harmony$predicted.celltype.l3))


# =========================================================
#                Azimuth diagnostic plots
# =========================================================

# ---------------------------------------------------------
# Mapping score plot
# ---------------------------------------------------------

message("Showing the confidence distribution")

mapping.score_plot <- VlnPlot(cd_harmony, features = "mapping.score", group.by = "seurat_clusters", pt.size = 0.1) +
  NoLegend() + 
  ggtitle("Mapping Score")

save_plot(mapping.score_plot, filename = "mapping.score_plot.png", dir = paths$plot_azimuth, width = 18, height = 14)

# ---------------------------------------------------------
# Cluster vs label heatmap
# ---------------------------------------------------------

message("Plotting Cluster vs labels heatmap")
tab_celltype_level2 <- table(
  cd_harmony$seurat_clusters,
  cd_harmony$predicted.celltype.l2
)

open_pdf(filename = "Azimuth_PBMC_Cluster_Annotation_Heatmap.pdf", dir = paths$plot_azimuth)

pheatmap(
  log10(tab_celltype_level2 + 1),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("white", "coral2"))(10),
  main = "Azimuth - Cluster vs labels"
         ) 

close_pdf()


# =========================================================
#             Cluster Annotation Summary Tables
# =========================================================

# ---------------------------------------------------------
# Cluster-label contingency table
# --------------------------------------------------------

message("Creating cluster-label contingency table...")

cluster_annotations_contigency <- table(
  cluster = cd_harmony$seurat_clusters,
  azimuth_Level1 = cd_harmony$predicted.celltype.l1,
  azimuth_Level2 = cd_harmony$predicted.celltype.l2,
  azimuth_Level3 = cd_harmony$predicted.celltype.l3
)

# ---------------------------------------------------------
#  Cluster annotation summaries ( CellType Level 2)
# --------------------------------------------------------

message("Creating CellAtlas cluster annotation summary (Choose only CellType Level2)...")

pbmcAzimuth_cluster_summary <- data.frame(
  cluster = cd_harmony$seurat_clusters,
  label = cd_harmony$predicted.celltype.l2
) |> 
  dplyr::count(cluster, label, name = "cells") |> 
  group_by(cluster) |> 
  mutate(
    percentage = cells / sum(cells) * 100
  ) |> 
  arrange(cluster, desc(cells))

# ---------------------------------------------------------
#  Majority-vote cluster annotations ( CellType Level 2)
# --------------------------------------------------------

message("Creating PBMC majority-vote annotations (Choose only CellType Level2)...")

pbmcAzimuth_majority <- pbmcAzimuth_cluster_summary |> 
  group_by(cluster) |> 
  slice_max(order_by = cells, n = 1, with_ties = FALSE) |> 
  ungroup()


# =========================================================
#                Save outputs
# =========================================================

message("Saving cluster annotation tables...")

save_csv(object = as.data.frame(cluster_annotations_contigency), filename = "Azimuth_Contingency_Table.csv", dir = paths$tables_azimuth)

message("Saving cluster annotation tables...")

save_csv(object = as.data.frame(pbmcAzimuth_cluster_summary), filename = "PBMC_Azimuth_ClusterSummary.csv", dir = paths$tables_azimuth)

message("Saving majority tables...")

save_csv(object = as.data.frame(pbmcAzimuth_majority), filename = "PBMC_Azimuth_MajorityVote.csv", dir = paths$tables_azimuth)

message("Saving annotated Seurat object...")

save_rds(cd_harmony, filename = "srt_obj_annotated_SingleR_Azimuth.rds", dir = paths$objects)

message("Azimiuth annotation pipeline completed successfully")

# =========================================================
#                  Save session info
# =========================================================

message("Saving session information (Azimuth annotation...")

save_session_info(filename = "sessionInfo_SingleR_Azimuth.txt", dir = paths$logs, label = "Azimuth annotations")

message("Session information saved at: ", paths$logs)


# =========================================================
#                  Final pipeline message
# =========================================================

message("=================================================")
message("PIPELINE STEP COMPLETED: Azimuth annotation")
message("=================================================")