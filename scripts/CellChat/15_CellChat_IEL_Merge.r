# =========================================================
# CellChat cell-cell communication analysis - IEL merge
# =========================================================
# Input:
#   - cellchat_IEL_CTRL.rds
#   - cellchat_IEL_CD.rds
#
# Output:
#   - CellChat merged object (CTRL vs CD)
#   - Comparative visualizations
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

message("Starting CellChat merge for IEL compartment (Control vs CD)...")


# =========================================================
#               Load individual CellChat objects
# =========================================================

cellchat_iel_ctrl <- readRDS(file.path(paths$objects_cellchat, "cellchat_IEL_CTRL.rds"))
cellchat_iel_cd <- readRDS(file.path(paths$objects_cellchat, "cellchat_IEL_CD.rds"))


# =========================================================
#        Align cell group levels across conditions
# =========================================================

group_new <- union(levels(cellchat_iel_ctrl@idents), levels(cellchat_iel_cd@idents))

message("Union of cell groups across conditions: ", length(group_new))

missing_ctrl <- setdiff(group_new, levels(cellchat_iel_ctrl@idents))
missing_cd   <- setdiff(group_new, levels(cellchat_iel_cd@idents))

if (length(missing_ctrl) > 0) {
  message("Groups lifted into CTRL (absent originally): ", paste(missing_ctrl, collapse = ", "))
}
if (length(missing_cd) > 0) {
  message("Groups lifted into CD (absent originally): ", paste(missing_cd, collapse = ", "))
}

cellchat_iel_ctrl <- liftCellChat(cellchat_iel_ctrl, group_new)
cellchat_iel_cd   <- liftCellChat(cellchat_iel_cd, group_new)


# =========================================================
#                  Merge CellChat objects
# =========================================================

message("Merging CellChat objects...")

object_list <- list(Control = cellchat_iel_ctrl, CD = cellchat_iel_cd)

cellchat_merged <- mergeCellChat(object_list, add.names = names(object_list))

message("Merge completed. Datasets: ", paste(names(object_list), collapse = ", "))


# =========================================================
# Compare number of interactions between conditions
# =========================================================

message("Generating global comparison plots...")

open_png(filename = "CellChat_IEL_compareInteractions_count.png",dir = paths$plots_cellchat_iel_merge,
         width = 800, height = 1000)

compareInteractions(cellchat_merged, show.legend = FALSE, group = c(1,2), measure = "count", title.name = "Number of inferred interactions (IEL)\nCTRL vs CD")

close_png()

open_png(filename = "CellChat_IEL_compareInteractions_weight.png", dir = paths$plots_cellchat_iel_merge,
         width = 800, height = 1000)

compareInteractions(cellchat_merged, show.legend = FALSE, group = c(1,2), measure = "weight", title.name = "Interaction strength (IEL)\nCTRL vs CD")

close_png()

# =========================================================
#         Differential network (CD vs Control)
# =========================================================

message("Generating differential network visualizations...")

open_png(filename = "CellChat_IEL_diffInteraction_count.png", dir = paths$plots_cellchat_iel_merge,
         width = 1800, height = 2000)

par(mar = c(1,1,8,1), xpd = TRUE)
netVisual_diffInteraction(cellchat_merged, weight.scale = TRUE, measure = "count", 
                          title.name = "Differential number of interactions\nCD vs Control",
                          vertex.label.cex = 0.6)

close_png()

open_png(filename = "CellChat_IEL_diffInteraction_weight.png", dir = paths$plots_cellchat_iel_merge,
         width = 1800, height = 2000)

par(mar = c(1,1,8,1), xpd = TRUE)
netVisual_diffInteraction(cellchat_merged, weight.scale = TRUE, measure = "weight",
                          title.name = "Differential interaction strength\nCD vs Control",
                          vertex.label.cex = 0.6)

close_png()


# =========================================================
#         Differential heatmap (CD vs Control)
# =========================================================

message("Generating differential communication heatmaps...")

open_png(filename = "CellChat_IEL_diffHeatmap_count.png", dir = paths$plots_cellchat_iel_merge,
         width = 1800, height = 1800)

netVisual_heatmap(cellchat_merged, measure = "count",   title.name = "Differential number of interactions\nCD vs Control")

close_png()

open_png(filename = "CellChat_IEL_diffHeatmap_weight.png", dir = paths$plots_cellchat_iel_merge,
         width = 1800, height = 1800)

netVisual_heatmap(cellchat_merged, measure = "weight",   title.name = "Differential interaction strength\nCD vs Control")

close_png()


# =========================================================
#     Global communication networks visualization
# =========================================================

message("Generating global communication networks...")

group_size_ctrl <- as.numeric(table(object_list[[1]]@idents))
group_size_cd <- as.numeric(table(object_list[[2]]@idents))

weight_max <- getMaxWeight(object_list, slot.name = c("idents", "count"), attribute = c("count"))

# -----------------------
# Control condition
# -----------------------

open_png(filename = "CellChat_IEL_globalNetwork_Control.png", dir = paths$plots_cellchat_iel_merge,
  width = 2200, height = 2200)

par(mar = c(1, 1, 8, 1), xpd = TRUE)

netVisual_circle(
  object_list[[1]]@net$count,
  vertex.weight = group_size_ctrl,
  weight.scale = TRUE,
  label.edge = FALSE,
  vertex.label.cex = 0.8,
  edge.weight.max = weight_max[1],
  title.name = "Number of interactions - Control"
)

close_png()


# -----------------------
# Crohn's disease condition
# -----------------------

open_png(filename = "CellChat_IEL_globalNetwork_CD.png", dir = paths$plots_cellchat_iel_merge,
  width = 2200, height = 2200)

par(mar = c(1, 1, 8, 1), xpd = TRUE)

netVisual_circle(
  object_list[[2]]@net$count,
  vertex.weight = group_size_cd,
  weight.scale = TRUE,
  label.edge = FALSE,
  vertex.label.cex = 0.8,
  edge.weight.max = weight_max[1],
  title.name = "Number of interactions - CD"
)

close_png()


# =========================================================
#                Save merged object
# =========================================================

message("Saving CellChat object merged for IEL Compartment")

stopifnot("joint" %in% names(cellchat_merged@idents))

message(
  "Final merged CellChat object (Control vs CD) contains ",
  length(levels(cellchat_merged@idents$joint)),
  " cell groups."
)

save_rds(cellchat_merged, "cellchat_IEL_merged.rds", dir = paths$objects_cellchat)

message("Merged CellChat object saved.")


# =========================================================
#                  Save session info
# =========================================================

message("[OUTPUT] Saving session information...")

save_session_info(filename = "sessionInfo_CellChat_IEL_merge.txt", dir = paths$logs, label = "CellChat merge - IEL compartment (Control vs CD)")

message("[OUTPUT] Session information saved to: ", paths$logs)


# =========================================================
#                  Final pipeline message
# =========================================================

message("=================================================")
message("[PIPELINE] CellChat merge completed for IEL compartment.")
message("=================================================")