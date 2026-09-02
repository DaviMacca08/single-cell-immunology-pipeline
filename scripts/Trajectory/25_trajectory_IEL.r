# =========================================================
# scRNA-seq trajectory analysis - IEL
# =========================================================
# Input:
#
# - Seurat object of IEL compartment with consensus cell-type annotation
#
# Output:
#
# - SingleCellExperiment object for trajectory analysis
# - Slingshot lineage inference and principal curves
# - Cell-level pseudotime and lineage weights
# - Cell-type composition across inferred lineages
# - Descriptive comparison of pseudotime and cell-state composition
#   between CD and Control conditions
# - Trajectory visualization and diagnostic plots
#
# Methodological note:
#
# - Trajectory inference is performed using Slingshot
# - Trajectory structure is inferred without splitting the dataset
#   by disease condition
# - CD vs Control comparisons are descriptive only due to strong
#   condition-associated cell-type imbalance
# - Condiments is used only for condition-imbalance diagnostics;
#   no differential trajectory tests are performed
# - tradeSeq is not applied
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


message("=== Starting Trajectory analysis for IEL compartment ===")


# =========================================================
#                  Functions
# =========================================================

# Creare SCE
seurat_to_sce <- function(obj, celltype_col, data_layer) {
  
  # check
  stopifnot(
    all(colnames(obj) %in% colnames(counts_full)),
    all(colnames(obj) %in% colnames(LayerData(obj, assay = "RNA", layer = data_layer)))
  )
  
  counts_mat <- counts_full[, colnames(obj), drop = FALSE]
  data_mat   <- LayerData(obj, assay = "RNA", layer = data_layer)
  
  common_cells <- intersect(colnames(counts_mat), colnames(data_mat))
  
  message("[SCE] Retained ", length(common_cells), " cells from ", ncol(obj), " cells.")
  
  counts_mat   <- counts_mat[, common_cells, drop = FALSE]
  data_mat     <- data_mat[, common_cells, drop = FALSE]
  
  meta <- obj@meta.data[common_cells, , drop = FALSE]
  
  sce <- SingleCellExperiment(
    assays  = list(counts = counts_mat, logcounts = data_mat),
    colData = DataFrame(
      celltype = as.character(meta[[celltype_col]]),
      disease  = as.character(meta[["disease"]]),
      donor    = as.character(meta[["donor"]])
    )
  )
  
  reducedDims(sce) <- list(
    PCA  = Embeddings(obj, "pca")[common_cells, , drop = FALSE],
    UMAP = Embeddings(obj, "umap")[common_cells, , drop = FALSE]
  )
  
  sce
}


# =========================================================
#                  Load data 
# =========================================================

message("[LOAD] Loading Seurat objects...")

iel_cells <- readRDS(file.path(paths$objects_subsetting, "srt_IEL_compartment_celltype.rds"))

full_obj <- readRDS(file.path(paths$objects_global, "srt_obj_merge_harmony.rds"))

full_obj <- JoinLayers(full_obj, assay = "RNA")

# All raw counts

counts_full <- LayerData(full_obj, assay = "RNA", layer = "counts")
rm(full_obj)


# =========================================================
#                  Create SingleCellExperiment
# =========================================================

sce_iel <- seurat_to_sce(iel_cells, "celltype_IEL", data_layer = "data")

# Check the SCE created

stopifnot(
  !is.null(assay(sce_iel, "counts")),
  !is.null(assay(sce_iel, "logcounts")),
  identical(dim(assay(sce_iel, "counts")), dim(assay(sce_iel, "logcounts")))
)


# =========================================================
#                  Imbalance Score
# =========================================================

imb_iel <- imbalance_score(
  Object     = reducedDim(sce_iel, "UMAP"),
  conditions = colData(sce_iel)$disease,
  k = 10, smooth = 10  
)

message("[DIAGNOSTIC] Condition imbalance calculated on the IEL UMAP ", "(k = 10, smooth = 10).")

colData(sce_iel)$imbalance_score <- imb_iel$scaled_scores

df <- as.data.frame(reducedDim(sce_iel, "UMAP"))
df$imbalance <- imb_iel$scaled_scores

# plot

plot_imbalance <- ggplot(df, aes(umap_1, umap_2, col = imbalance)) +
  geom_point(size = 0.5) +
  scale_color_viridis_c() +
  theme_classic() +
  ggtitle("IEL — Condition imbalance across the UMAP (CD vs Control)")

save_plot(plot_imbalance, filename = "ImbalanceScore_IEL.png", dir = paths$plots_trajectory_iel,
          width = 15, height = 15)


# =========================================================
#                  Fit trajectory (Slingshot)
# =========================================================

iel_trajectory_cells <- !sce_iel$celltype %in% c(
  "Epithelial cells / enterocytes",
  "CST3+ LYZ+ macrophages"
)

message("[TRAJ] Excluded from trajectory inference: ", "Epithelial cells / enterocytes and CST3+ LYZ+ macrophages.")

message("[TRAJ] Cells retained for trajectory inference: ", sum(iel_trajectory_cells))

sce_iel_traj <- sce_iel[, iel_trajectory_cells]

# Lineages

lineages_iel <- getLineages(
  reducedDim(sce_iel_traj, "PCA"),        
  clusterLabels = colData(sce_iel_traj)$celltype,
  start.clus = "Naïve CD4",
  omega = TRUE
)

open_png(filename = "lineages_IEL.png", dir = paths$plots_trajectory_iel, 
         width = 2400, height = 2400)

par(mar = c(5, 5, 4, 2) + 0.1)
plot(
  reducedDim(sce_iel_traj, "PCA"),
  col = as.factor(sce_iel_traj$celltype),
  pch = 16,
  cex = 0.4,
  main = "IEL — Slingshot lineage structure"
)

lines(
  SlingshotDataSet(lineages_iel),
  type = "lineages",
  lwd = 3
)

close_png()

message("[TRAJ] Lineage topology fitted successfully.")

# Curves

curves_iel <- getCurves(
  lineages_iel,
  approx_points = 150
)

# Project curves onto UMAP for interpretable plotting only

curves_iel_umap <- embedCurves(curves_iel, reducedDim(sce_iel_traj, "UMAP"))

message("[TRAJ] Number of lineages recovered: ", length(slingLineages(curves_iel)))
print(slingLineages(curves_iel))


# =========================================================
#       Pseudotime extraction and core visualization
# =========================================================

pseudo_iel   <- slingPseudotime(curves_iel, na = TRUE)
weights_iel  <- slingCurveWeights(curves_iel)
n_lineages   <- ncol(pseudo_iel)

# Overall pseudotime (weighted average across lineages) 

overall_pt <- rowSums(pseudo_iel * weights_iel, na.rm = TRUE) / rowSums(weights_iel, na.rm = TRUE)
colData(sce_iel_traj)$pseudotime_overall <- overall_pt

df_umap <- as.data.frame(reducedDim(sce_iel_traj, "UMAP"))
df_umap$pseudotime <- overall_pt
df_umap$celltype   <- colData(sce_iel_traj)$celltype
df_umap$disease    <- colData(sce_iel_traj)$disease

# plot

plot_pseudotime_overall <- ggplot(df_umap, aes(umap_1, umap_2, col = pseudotime)) +
  geom_point(size = 0.5) +
  scale_color_viridis_c(option = "plasma") +
  theme_classic() +
  ggtitle("IEL — overall pseudotime (weighted average across lineages)")

save_plot(plot_pseudotime_overall, filename = "Pseudotime_overall_IEL.png", dir = paths$plots_trajectory_iel,
          width = 15, height = 15)

# Per-lineage pseudotime plots (one panel per lineage, cells not assigned to that lineage are left grey)

plot_pseudotime_per_lineage <- lapply(seq_len(n_lineages), function(i) {
  
  d <- df_umap
  d$pt_lineage <- pseudo_iel[, i]
  
  p <- ggplot(d, aes(umap_1, umap_2, col = pt_lineage)) +
    geom_point(size = 0.5) +
    scale_color_viridis_c(option = "plasma", na.value = "grey85") +
    theme_classic() +
    ggtitle(paste0("IEL — pseudotime, lineage ", i))
  
  save_plot(p, filename = paste0("pseudotime_IEL_lineage_", i, ".png"), dir = paths$plots_trajectory_iel,
            width = 15, height = 15) 
  
})

# Celltype composition per lineage (cell assigned to a lineage if curve weight > 0)
# Diagnostic lineage membership defined by positive Slingshot curve weight.
# This is used for descriptive cell-type composition only, not as a hard biological lineage assignment.

lineage_membership <- apply(weights_iel, 2, function(w) w > 0)
for (i in seq_len(n_lineages)) {
  message("[TRAJ] Lineage ", i, " — celltype composition:")
  print(table(colData(sce_iel_traj)$celltype[lineage_membership[, i]]))
}


# =========================================================
#      Diagnostic condition overlay (visual only)
# =========================================================

# Baseline: condition on UMAP

plot_condition_umap <- ggplot(df_umap, aes(umap_1, umap_2, col = disease)) +
  geom_point(size = 0.5, alpha = 0.7) +
  theme_classic() +
  ggtitle("IEL — Condition distribution across the UMAP")

save_plot(plot_condition_umap, filename  = "Condition_umap_IEL.png", dir = paths$plots_trajectory_iel,
          width = 15, height = 15)

# Imbalance score projected on the fitted curves 

df_umap$imbalance <- colData(sce_iel_traj)$imbalance_score <- imb_iel$scaled_scores[iel_trajectory_cells]

plot_imbalance_on_traj <- ggplot(df_umap, aes(umap_1, umap_2, col = imbalance)) +
  geom_point(size = 0.5) +
  scale_color_viridis_c() +
  theme_classic() +
  ggtitle("IEL — imbalance score on trajectory-fitted cells (CD vs Control)")

save_plot(plot_imbalance_on_traj, filename = "Imbalance_Score_fitted_IEL.png", dir = paths$plots_trajectory_iel,
          width = 15, height = 15)

# Pseudotime density by condition, per lineage 
# Descriptive only: no formal statistical comparison of pseudotime distributions is performed.

df_pt_long <- do.call(rbind, lapply(seq_len(n_lineages), function(i) {
  data.frame(
    pseudotime = pseudo_iel[, i],
    lineage    = paste0("Lineage ", i),
    disease    = colData(sce_iel_traj)$disease
  )
}))

df_pt_long <- df_pt_long[!is.na(df_pt_long$pseudotime), ]

plot_pt_density_by_condition <- ggplot(df_pt_long, aes(pseudotime, fill = disease)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~lineage, scales = "free") +
  theme_classic() +
  ggtitle("IEL — pseudotime distribution by condition (descriptive only)")

save_plot(plot_pt_density_by_condition, filename = "Pseudotime_distribution_condition.png",
          dir = paths$plots_trajectory_iel, width = 15, height = 15)

# Compositional profile along pseudotime bins, per lineage: % CD vs Control in successive bins of the path 
# Descriptive compositional profile only.

n_bins <- 10

df_pt_long$bin <- ave(
  df_pt_long$pseudotime, df_pt_long$lineage,
  FUN = function(x) cut(x, breaks = n_bins, labels = FALSE)
)

comp_profile <- aggregate(
  disease ~ lineage + bin, data = df_pt_long,
  FUN = function(x) mean(x == "CD")
)

names(comp_profile)[3] <- "prop_CD"

plot_composition_profile <- ggplot(comp_profile, aes(bin, prop_CD)) +
  geom_line() + geom_point() +
  geom_hline(yintercept = mean(colData(sce_iel_traj)$disease == "CD"),
             linetype = "dashed", col = "grey50") +
  facet_wrap(~lineage) +
  ylim(0, 1) +
  theme_classic() +
  labs(x = "Pseudotime bin (early -> late)", y = "Proportion CD",
       title = "IEL — CD proportion along pseudotime (Dashed = Global cell-level CD proportion)")

save_plot(plot_composition_profile, filename = "CompositionProfile.png", dir = paths$plots_trajectory_iel,
          width = 15, height = 15)


# =========================================================
#   Descriptive marker trends along pseudotime (no model fit)
# =========================================================

marker_panel <- c("GZMA", "GZMB", "GZMK", "NKG7", "PRF1",
                  "ITGAE", "CD69", "PDCD1", "LAG3", "HAVCR2")
marker_panel <- intersect(marker_panel, rownames(sce_iel_traj))
message("[TRAJ] Descriptive marker panel found in data: ", paste(marker_panel, collapse = ", "))

n_bins <- 10
marker_trend_list <- lapply(seq_len(n_lineages), function(i) {
  
  pt_i <- pseudo_iel[, i]
  keep <- !is.na(pt_i)
  if (sum(keep) < 20) return(NULL)   # skip lineages with too few assigned cells to bin meaningfully
  
  bin_i <- cut(pt_i[keep], breaks = n_bins, labels = FALSE)
  expr_i <- as.matrix(logcounts(sce_iel_traj)[marker_panel, keep, drop = FALSE])
  
  df <- as.data.frame(t(expr_i))
  df$bin <- bin_i
  df_long <- reshape2::melt(df, id.vars = "bin", variable.name = "gene", value.name = "expression")
  df_mean <- aggregate(expression ~ bin + gene, data = df_long, FUN = mean)
  df_mean$lineage <- paste0("Lineage ", i)
  df_mean
})

marker_trend_df <- do.call(rbind, marker_trend_list)

plot_marker_trends <- ggplot(marker_trend_df, aes(bin, expression, color = gene)) +
  geom_line(linewidth = 0.8) + geom_point(size = 1) +
  facet_wrap(~lineage, scales = "free_y") +
  theme_classic() +
  labs(x = "Pseudotime bin (early -> late)", y = "Mean logcounts",
       title = "IEL — cytotoxicity, tissue-residency and inhibitory-state markers along pseudotime",
       subtitle = "Descriptive only (binned mean expression; no model fit)")

save_plot(plot_marker_trends, filename = "MarkerTrends_along_pseudotime.png",
          dir = paths$plots_trajectory_iel, width = 20, height = 15)


# =========================================================
#            Final per-lineage summary table
# =========================================================

lineage_summary <- do.call(rbind, lapply(seq_len(n_lineages), function(i) {
  members <- lineage_membership[, i]
  dz <- colData(sce_iel_traj)$disease[members]
  lineage_end_cluster <- tail(slingLineages(curves_iel)[[i]], 1)
  data.frame(
    lineage           = paste0("Lineage ", i),
    n_cells           = sum(members),
    prop_CD           = round(mean(dz == "CD"), 3),
    lineage_end_cluster = lineage_end_cluster
  )
}))

print(lineage_summary)


# =========================================================
#                  Save outputs
# =========================================================

message("[SAVE] Saving IEL trajectory outputs...")

save_rds(curves_iel, filename = "slingshot_curves_IEL.rds", dir = paths$objects_trajectory)
save_rds(sce_iel_traj, filename = "sce_IEL_trajectory.rds", dir = paths$objects_trajectory)
save_rds(pseudo_iel, filename = "pseudotime_IEL.rds", dir = paths$objects_trajectory)
save_rds(weights_iel, filename = "curve_weights_IEL.rds", dir = paths$objects_trajectory)
save_csv(lineage_summary, filename = "IEL_lineage_summary.csv", dir = paths$objects_trajectory)

message("[SAVE] Slingshot curves saved.")
message("[SAVE] SCE trajectory object saved.")
message("[SAVE] Pseudotime matrix saved.")
message("[SAVE] Curve weights saved.")
message("[SAVE] Lineage summary table saved.")


# =========================================================
#                  Save session info
# =========================================================

message("[OUTPUT] Saving session information...")

save_session_info(filename = "sessionInfo_Trajectory_IEL.txt", dir = paths$logs,
                  label = "Trajectory Analysis - IEL compartment")

message("[OUTPUT] Session information saved to: ", paths$logs)


# =========================================================
#                  Final pipeline message
# =========================================================

message("=================================================")
message("[PIPELINE] IEL trajectory analysis completed.")
message("[PIPELINE] Method: Slingshot descriptive trajectory inference.")
message("[PIPELINE] Condition comparisons are descriptive only.")
message("[PIPELINE] Condiments: imbalance diagnostic only; no differential trajectory testing.")
message("[PIPELINE] tradeSeq: not performed.")
message("=================================================")

