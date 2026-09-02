# =========================================================
# CellChat cell-cell communication analysis - LPL merge
# =========================================================
# Input:
#   - cellchat_LPL_CTRL.rds
#   - cellchat_LPL_CD.rds
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

message("Starting CellChat merge for LPL compartment (Control vs CD)...")


# =========================================================
#               Load individual CellChat objects
# =========================================================

cellchat_lpl_ctrl <- readRDS(file.path(paths$objects_cellchat, "cellchat_LPL_CTRL.rds"))
cellchat_lpl_cd <- readRDS(file.path(paths$objects_cellchat, "cellchat_LPL_CD.rds"))


# =========================================================
#        Align cell group levels across conditions
# =========================================================

group_new <- union(levels(cellchat_lpl_ctrl@idents), levels(cellchat_lpl_cd@idents))

message("Union of cell groups across conditions: ", length(group_new))

missing_ctrl <- setdiff(group_new, levels(cellchat_lpl_ctrl@idents))
missing_cd   <- setdiff(group_new, levels(cellchat_lpl_cd@idents))

if (length(missing_ctrl) > 0) {
  message("Groups lifted into CTRL (absent originally): ", paste(missing_ctrl, collapse = ", "))
}
if (length(missing_cd) > 0) {
  message("Groups lifted into CD (absent originally): ", paste(missing_cd, collapse = ", "))
}

cellchat_lpl_ctrl <- liftCellChat(cellchat_lpl_ctrl, group_new)
cellchat_lpl_cd   <- liftCellChat(cellchat_lpl_cd, group_new)


# =========================================================
#                  Merge CellChat objects
# =========================================================

message("Merging CellChat objects...")

object_list <- list(Control = cellchat_lpl_ctrl, CD = cellchat_lpl_cd)

cellchat_merged <- mergeCellChat(object_list, add.names = names(object_list))

message("Merge completed. Datasets: ", paste(names(object_list), collapse = ", "))


# =========================================================
# Compare number of interactions between conditions
# =========================================================

message("Generating global comparison plots...")

open_png(filename = "CellChat_LPL_compareInteractions_count.png",dir = paths$plots_cellchat_lpl_merge,
         width = 800, height = 1000)

compareInteractions(cellchat_merged, show.legend = FALSE, group = c(1,2), measure = "count", title.name = "Number of inferred interactions (LPL)\nCTRL vs CD")

close_png()

open_png(filename = "CellChat_LPL_compareInteractions_weight.png", dir = paths$plots_cellchat_lpl_merge,
         width = 800, height = 1000)

compareInteractions(cellchat_merged, show.legend = FALSE, group = c(1,2), measure = "weight", title.name = "Interaction strength (LPL)\nCTRL vs CD")

close_png()


# =========================================================
#         Differential network (CD vs Control)
# =========================================================

message("Generating differential network visualizations...")

open_png(filename = "CellChat_LPL_diffInteraction_count.png", dir = paths$plots_cellchat_lpl_merge,
         width = 1800, height = 2000)

par(mar = c(1,1,8,1), xpd = TRUE)
netVisual_diffInteraction(cellchat_merged, weight.scale = TRUE, measure = "count", 
                          title.name = "Differential number of interactions\nCD vs Control",
                          vertex.label.cex = 0.6)

close_png()

open_png(filename = "CellChat_LPL_diffInteraction_weight.png", dir = paths$plots_cellchat_lpl_merge,
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

open_png(filename = "CellChat_LPL_diffHeatmap_count.png", dir = paths$plots_cellchat_lpl_merge,
         width = 1800, height = 1800)

netVisual_heatmap(cellchat_merged, measure = "count",   title.name = "Differential number of interactions\nCD vs Control")

close_png()

open_png(filename = "CellChat_LPL_diffHeatmap_weight.png", dir = paths$plots_cellchat_lpl_merge,
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

open_png(filename = "CellChat_LPL_globalNetwork_Control.png", dir = paths$plots_cellchat_lpl_merge,
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

open_png(filename = "CellChat_LPL_globalNetwork_CD.png", dir = paths$plots_cellchat_lpl_merge,
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

message("Saving CellChat object merged for LPL Compartment")

stopifnot("joint" %in% names(cellchat_merged@idents))

message(
  "Final merged CellChat object (Control vs CD) contains ",
  length(levels(cellchat_merged@idents$joint)),
  " cell groups."
)

save_rds(cellchat_merged, "cellchat_LPL_merged.rds", dir = paths$objects_cellchat)

message("Merged CellChat object saved.")


# =========================================================
#                  Save session info
# =========================================================

message("[OUTPUT] Saving session information...")

save_session_info(filename = "sessionInfo_CellChat_LPL_merge.txt", dir = paths$logs, label = "CellChat merge - LPL compartment (Control vs CD)")

message("[OUTPUT] Session information saved to: ", paths$logs)


# =========================================================
#                  Final pipeline message
# =========================================================

message("=================================================")
message("[PIPELINE] CellChat merge completed for LPL compartment.")
message("=================================================")