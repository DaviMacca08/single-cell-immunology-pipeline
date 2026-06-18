# =========================================================
# SingleR Cell Type Annotation
# Input: Harmony-integrated Seurat object
# Output: Annotated Seurat object + diagnostic plots
# =========================================================


# =========================================================
#                  Libraries & Setup
# =========================================================

source("scripts/00_setup/00_paths.R")
source("scripts/00_setup/01_environment.R")
source("scripts/00_setup/02_io_helpers.R")
source("scripts/00_setup/04_seed.R")

set_seed(1234)

message("=== Start SingleR Annotation ===")


# =========================================================
#           Load data (Seurat object - harmony)
# =========================================================

message("Loading Seurat object...")

cd_harmony <- readRDS(
  file.path(paths$objects, "srt_obj_harmony.rds")
)


# =========================================================
#                  Initial sanity checks
# =========================================================

stopifnot(inherits(cd_harmony, "Seurat"))

stopifnot("seurat_clusters" %in% colnames(cd_harmony@meta.data))

message("Seurat object loaded successfully")


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
#                    Metadata overview
# =========================================================

message("Cluster distribution:")

print(table(Idents(cd_harmony)))

# =========================================================
#               SingleR annotation
# =========================================================

message("Preparing normalized matrix...")

counts_norm <- LayerData(cd_harmony, assays = "RNA", layer = "data")

# ---------------------------------------------------------
# Load reference databases
# ---------------------------------------------------------

message("Loading celldex references...")

CellAtlas_db <- celldex::HumanPrimaryCellAtlasData()

MonacoImmune_db <- celldex::MonacoImmuneData()

# ---------------------------------------------------------
# Run SingleR
# ---------------------------------------------------------

message("Running SingleR: CellAtlas")

CellAtlas_pred <- SingleR(
  test = counts_norm,
  ref = CellAtlas_db,
  labels = CellAtlas_db$label.main
)

message("Running SingleR: MonacoImmune")

MonacoImmune_pred <- SingleR(
  test = counts_norm,
  ref = MonacoImmune_db,
  labels = MonacoImmune_db$label.main
)

# ---------------------------------------------------------
# Annotation QC checks
# ---------------------------------------------------------

message("Checking annotation quality...")

message("CellAtlas NA annotations:")

print(table(is.na(CellAtlas_pred$pruned.labels)))

message("MonacoImmune NA annotations:")

print(table(is.na(MonacoImmune_pred$pruned.labels)))

# ---------------------------------------------------------
# Add SingleR labels to metadata
# ---------------------------------------------------------

cd_harmony$Single.R.labels.CellAtlas <- CellAtlas_pred$pruned.labels[match(rownames(cd_harmony@meta.data), rownames(CellAtlas_pred))]
cd_harmony$Single.R.labels.MonacoImmune <- MonacoImmune_pred$pruned.labels[match(rownames(cd_harmony@meta.data), rownames(MonacoImmune_pred))]


# =========================================================
#                SingleR diagnostic plots
# =========================================================

# ---------------------------------------------------------
# CellAtlas heatmap
# ---------------------------------------------------------

open_pdf(filename = "SingleR_CellAtlas_ScoreHeatmap.pdf", dir = paths$plots_singler, width = 12, height = 10)

plotScoreHeatmap(CellAtlas_pred)

close_pdf()

# ---------------------------------------------------------
# MonacoImmune heatmap
# ---------------------------------------------------------

open_pdf(filename = "SingleR_Monaco_ScoreHeatmap.pdf", dir = paths$plots_singler, width = 12, height = 10)

plotScoreHeatmap(MonacoImmune_pred)

close_pdf()

# ---------------------------------------------------------
# Delta distributions
# ---------------------------------------------------------

open_pdf(filename = "SingleR_CellAtlas_DeltaDistribution.pdf", dir = paths$plots_singler)

plotDeltaDistribution(CellAtlas_pred)

close_pdf()

open_pdf(filename = "SingleR_Monaco_DeltaDistribution.pdf", dir = paths$plots_singler)

plotDeltaDistribution(MonacoImmune_pred)

close_pdf()


# =========================================================
#             Cluster-level annotation heatmaps
# =========================================================

tab_CellAtlas <- table(
  Assigned = CellAtlas_pred$pruned.labels,
  Clusters = cd_harmony$seurat_clusters
)

tab_MonacoImmune <- table(
  Assigned = MonacoImmune_pred$pruned.labels,
  Clusters = cd_harmony$seurat_clusters
)

# ---------------------------------------------------------
# CellAtlas pheatmap
# ---------------------------------------------------------

open_pdf(filename = "CellAtlas_Cluster_Annotation_Heatmap.pdf", dir = paths$plots_singler)

pheatmap(
  log10(tab_CellAtlas + 10),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("white", "coral2"))(10),
  main = "Single R Cell Atlas - Cluster vs labels"
)

close_pdf()

# ---------------------------------------------------------
# MonacoImmune pheatmap
# ---------------------------------------------------------

open_pdf(filename = "MonacoImmune_Cluster_Annotation_Heatmap.pdf", dir = paths$plots_singler)

pheatmap(
  log10(tab_MonacoImmune + 10),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("white", "coral2"))(10),
  main = "Single R Monaco Immune - Cluster vs labels"
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
  cellAtlas = cd_harmony$Single.R.labels.CellAtlas,
  monacoImmune = cd_harmony$Single.R.labels.MonacoImmune
)

# ---------------------------------------------------------
#  Cluster annotation summaries
# --------------------------------------------------------

message("Creating CellAtlas cluster annotation summary..")

CellAtlas_cluster_summary <- data.frame(
  cluster = cd_harmony$seurat_clusters,
  label = cd_harmony$Single.R.labels.CellAtlas
) |> 
  dplyr::count(cluster, label, name = "cells") |> 
  group_by(cluster) |> 
  mutate(
    percentage = cells / sum(cells) * 100
  ) |> 
  arrange(cluster, desc(cells))

message("Creating MonacoImmune cluster annotation summary...")

Monaco_cluster_summary <- data.frame(
  cluster = cd_harmony$seurat_clusters,
  label = cd_harmony$Single.R.labels.MonacoImmune
) |> 
  dplyr::count(cluster, label, name = "cells") |> 
  group_by(cluster) |> 
  mutate(
    percentage = cells / sum(cells) * 100
  ) |> 
  arrange(cluster, desc(cells))


# ---------------------------------------------------------
#  Majority-vote cluster annotations
# --------------------------------------------------------

message("Creating CellAtlas majority-vote annotations...")

CellAtlas_majority <- CellAtlas_cluster_summary |> 
  group_by(cluster) |> 
  slice_max(order_by = cells, n = 1, with_ties = FALSE) |> 
  ungroup()

message("Doing majority vote table for MonacoImmune")

Monaco_majority <- Monaco_cluster_summary |> 
  group_by(cluster) |> 
  slice_max(order_by = cells, n = 1, with_ties = FALSE) |>  
  ungroup()

# =========================================================
#                Save outputs
# =========================================================

message("Saving cluster contigency table...")

save_csv(object = as.data.frame(cluster_annotations_contigency), filename = "SingleR_Contingency_Table.csv", dir = paths$tables_singler)

message("Saving cluster annotation tables...")
cat("CellAtlas")

save_csv(object = as.data.frame(CellAtlas_cluster_summary), filename = "CellAtlas_ClusterSummary.csv", dir = paths$tables_singler)

cat("MonacoImmune")

save_csv(object = as.data.frame(Monaco_cluster_summary), filename = "MonacoImmune_ClusterSummary.csv", dir = paths$tables_singler)

message("Saving majority tables...")
cat("CellAtlas")

save_csv(object = as.data.frame(CellAtlas_majority), filename = "CellAtlas_MajorityVote.csv", dir = paths$tables_singler)

cat("MonacoImmune")

save_csv(object = as.data.frame(Monaco_majority), filename = "MonacoImmune_MajorityVote.csv", dir = paths$tables_singler)

message("Saving annotated Seurat object...")

save_rds(cd_harmony, filename = "srt_obj_annotated_SingleR.rds", dir = paths$objects)

message("SingleR annotation pipeline completed successfully")

# =========================================================
#                  Save session info
# =========================================================

message("Saving session information (SingleR annotation...")

save_session_info(filename = "sessionInfo_SingleR.txt", dir = paths$logs, label = "SingleR annotations")

message("Session information saved at: ", paths$logs)


# =========================================================
#                  Final pipeline message
# =========================================================

message("=================================================")
message("PIPELINE STEP COMPLETED: SingleR Annotation")
message("=================================================")