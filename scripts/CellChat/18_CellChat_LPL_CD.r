# =========================================================
# CellChat cell-cell communication analysis - LPL CD
# =========================================================
# Input:
#   - Seurat object of LPL compartment Crohn's disease with consensus annotation
#
# Output:
#   - CellChat object
#   - Ligand-receptor interaction analysis
#   - Signaling pathway and network visualization plots
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

message("Starting CellChat analysis for LPL compartment Crohn's disease...")


# =========================================================
#                   Load input data
# =========================================================

message("Loading Seurat object...")

lpl_cells <- readRDS(
  file.path(paths$objects_subsetting, "srt_LPL_compartment_celltype.rds")
)

DefaultAssay(lpl_cells) <- "RNA"

lpl_cells_cd <- subset(lpl_cells, disease == "CD")

data_input <- GetAssayData(lpl_cells_cd, assay = "RNA", layer = "data")
labels <- Idents(lpl_cells_cd)
sample_id <- as.factor(lpl_cells_cd$sample_id)
disease <- lpl_cells_cd$disease

# Create metadata fot cellchat object

meta <- data.frame(
  samples = sample_id,
  group = labels,
  disease = disease,
  row.names = names(labels)
)


# =========================================================
#                    Input checks
# =========================================================

message("Performing input checks...")

stopifnot(
  ncol(lpl_cells_cd) ==
    ncol(LayerData(
      lpl_cells_cd,
      assay = "RNA",
      layer = "data"
    ))
)

stopifnot(
  identical(colnames(data_input), rownames(meta)),
  ncol(data_input) == nrow(meta),
  !any(is.na(meta$group)),
  !any(is.na(meta$samples)),
  !any(is.na(meta$disease))
)

message(
  "Input validated: ",
  ncol(data_input), " cells, ",
  nrow(data_input), " genes, ",
  length(unique(meta$group)), " cell groups, ",
  length(unique(meta$samples)), " samples."
)

message("Cell group distribution by disease:")

celltype_disease_table <- table(
  CellType = Idents(lpl_cells),
  Disease = lpl_cells$disease
)

print(celltype_disease_table)

message("Checking cell groups missing in CD subset...")

missing_groups_cd <- setdiff(
  levels(Idents(lpl_cells)),
  levels(Idents(lpl_cells_cd))
)

if(length(missing_groups_cd) > 0){
  message(
    "Cell groups absent in CD: ",
    paste(missing_groups_cd, collapse = ", "),
    ". Tot: ", length(missing_groups_cd)
  )
}


# =========================================================
#               Create CellChat object
# =========================================================

message("Creating CellChat object...")

cellchat_obj <- createCellChat(
  object = data_input, 
  meta = meta, 
  group.by = "group"
)

message("CellChat object created successfully.")


# =========================================================
#      Load the ligand-receptor interactions database
# =========================================================

message("Loading CellChat human interaction database...")

CellChatDB_human <- CellChatDB.human
cellchat_obj@DB <- CellChatDB_human

message(
  "Loaded CellChat database with ",
  nrow(CellChatDB_human$interaction),
  " ligand-receptor interactions."
)

showDatabaseCategory(CellChatDB_human)


# =========================================================
#      Subset and pre-processing the expression data
# =========================================================

message("Subsetting expression data...")

cellchat_obj <- subsetData(cellchat_obj)

message("Pre-processing expression data...")

cellchat_obj <- identifyOverExpressedGenes(cellchat_obj)
cellchat_obj <- identifyOverExpressedInteractions(cellchat_obj)

message("Compute the communication probability...")

cellchat_obj <- computeCommunProb(cellchat_obj)

message(
  "Number of inferred ligand-receptor communications: ",
  nrow(cellchat_obj@net$prob)
)

message("Filter out the cell-cell communication with few number of cells...")

cellchat_obj <- filterCommunication(cellchat_obj, min.cells = 10)

message("Infer the cell-cell communication at a signaling pathway level...")

cellchat_obj <- computeCommunProbPathway(cellchat_obj)

message("Calculated the aggregated cell-cell communication network...")

cellchat_obj <- aggregateNet(cellchat_obj)


# =========================================================
# Visualize aggregated cell-cell communication network
# =========================================================

message("Visualize the aggregated cell-cell communication network...")

group_size <- as.numeric(table(cellchat_obj@idents))

message("Number of cell groups: ",length(group_size))

# -----------------------
# Number of interactions
# -----------------------

open_png(filename = "CellChat_LPL_CD_global_network_count.png", dir = paths$plots_cellchat_lpl,
         width = 1800, height = 2000)

par(mar = c(1,1,8,1), xpd = TRUE)

netVisual_circle(
  cellchat_obj@net$count,
  vertex.weight = group_size,
  weight.scale = TRUE,
  label.edge = FALSE,
  vertex.label.cex = 0.6,
  edge.weight.max = max(cellchat_obj@net$count),
  title.name = "Number of interactions"
)

close_png()

# -----------------------
# Interaction strength
# -----------------------

open_png(filename = "CellChat_LPL_CD_global_network_strength.png", dir = paths$plots_cellchat_lpl,
         width = 1800, height = 2200)

par(mar = c(1,1,8,1), xpd = TRUE)

netVisual_circle(
  cellchat_obj@net$weight,
  vertex.weight = group_size,
  weight.scale = TRUE,
  label.edge = FALSE,
  vertex.label.cex = 0.6,
  edge.weight.max = max(cellchat_obj@net$weight),
  title.name = "Interaction strength"
)

close_png()


# =========================================================
# Visualize outgoing communication per cell group
# =========================================================

message("Generating outgoing communication network panel...")

mat <- cellchat_obj@net$weight

open_png(filename = "CellChat_LPL_CD_outgoing_networks_all_groups.png", dir = paths$plots_cellchat_lpl,
         width = 3000, height = 3000)

par(mfrow = c(5,4), mar = c(1,1,3,1), xpd = TRUE)

for (i in seq_len(nrow(mat))) {
  
  cell_group <- rownames(mat)[i]
  
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  
  mat2[i, ] <- mat[i, ]
  
  netVisual_circle(
    mat2,
    vertex.weight = group_size,
    weight.scale = TRUE,
    edge.weight.max = max(mat),
    label.edge = FALSE,
    vertex.label.cex = 0.5,
    title.name = cell_group
  )
}

close_png()

message("Outgoing communication network panel saved.")


# =========================================================
# Visualize incoming communication per cell group
# =========================================================

message("Generating incoming communication network panel...")

mat_t <- t(cellchat_obj@net$weight)

open_png(filename = "CellChat_LPL_CD_incoming_networks_all_groups.png", dir = paths$plots_cellchat_lpl,
         width = 3000, height = 3000)

par(mfrow = c(5,4), mar = c(1,1,3,1), xpd = TRUE)

for (i in seq_len(nrow(mat_t))) {
  
  cell_group <- rownames(mat_t)[i]
  
  mat2 <- matrix(0, nrow = nrow(mat_t), ncol = ncol(mat_t), dimnames = dimnames(mat_t))
  
  mat2[i, ] <- mat_t[i, ]
  
  netVisual_circle(
    mat2,
    vertex.weight = group_size,
    weight.scale = TRUE,
    edge.weight.max = max(mat_t),
    label.edge = FALSE,
    vertex.label.cex = 0.5,
    title.name = cell_group
  )
}

close_png()

message("Incoming communication network panel saved.")


# =========================================================
# Explore signaling pathway activity
# =========================================================

message("Exploring inferred signaling pathways...")

pathways_all <- cellchat_obj@netP$pathways

message("Number of inferred signaling pathways: ", length(pathways_all))

# ---------------------------------------------------------
# Rank signaling pathways by overall information flow
# ---------------------------------------------------------

open_png(filename = "CellChat_LPL_CD_rankNet.png", dir = paths$plots_cellchat_lpl,
         width = 1200, height = 1600)

rankNet(cellchat_obj, mode = "single", measure = "weight", stacked = FALSE, do.stat = FALSE)

close_png()

# ---------------------------------------------------------
# Visualize each signaling pathway
# ---------------------------------------------------------

message("Generating pathway-specific visualizations...")

for (pw in pathways_all) {
  
  open_png(filename = paste0("CellChat_LPL_CD_", pw, "_circle.png"), dir = paths$plots_cellchat_lpl_path, 
           width = 1400, height = 1400)
  
  par(mar = c(1,1,4,1), xpd = TRUE)
  
  netVisual_aggregate(cellchat_obj, signaling = pw, layout = "circle")
  
  close_png()
  
  p <- netAnalysis_contribution(cellchat_obj, signaling = pw)
  
  save_plot(p, filename = paste0("CellChat_LPL_CD_", pw, "_contribution.png"), dir = paths$plots_cellchat_lpl_path)
  
  rm(p)
  
}

message("Pathway visualizations saved.")


# =========================================================
# Signaling roles 
# =========================================================

# ----------------------------
# Signaling role scatter plot
# ----------------------------

cellchat_obj <- netAnalysis_computeCentrality(cellchat_obj, slot.name = "netP")

open_png(filename = "CellChat_LPL_CD_signalingRole_scatter.png", dir = paths$plots_cellchat_lpl,
         width = 1400, height = 1200)

p <- netAnalysis_signalingRole_scatter(cellchat_obj) + 
  ggplot2::ggtitle("Cell signaling roles: outgoing vs incoming activity") +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      hjust = 0.5,
      face = "bold",
      size = 16
    )
  )

print(p)
rm(p)

close_png()

# ----------------------------
# Signaling role heatmaps
# ----------------------------

open_png(filename = "CellChat_LPL_CD_signalingRole_heatmap_outgoing.png", dir = paths$plots_cellchat_lpl,
         width = 1000, height = 1500)

netAnalysis_signalingRole_heatmap(cellchat_obj, pattern = "outgoing")

grid::grid.text(
  "Outgoing signaling roles",
  y = grid::unit(0.98, "npc"),
  gp = grid::gpar(fontsize = 18, fontface = "bold")
)

close_png()

open_png(filename = "CellChat_LPL_CD_signalingRole_heatmap_incoming.png", dir = paths$plots_cellchat_lpl,
         width = 1000, height = 1500)

netAnalysis_signalingRole_heatmap(cellchat_obj, pattern = "incoming")

grid::grid.text(
  "Incoming signaling roles",
  y = grid::unit(0.98, "npc"),
  gp = grid::gpar(fontsize = 18, fontface = "bold")
)

close_png()


# =========================================================
#                Save full Cellchat object
# =========================================================

message("Saving CellChat object for LPL CD")

message(
  "Final CellChat object for LPL CD contains ",
  length(levels(cellchat_obj@idents)),
  " cell groups and ",
  length(cellchat_obj@netP$pathways),
  " signaling pathways."
)

save_rds(cellchat_obj, "cellchat_LPL_CD.rds", dir = paths$objects_cellchat)


# =========================================================
#                  Save session info
# =========================================================

message("[OUTPUT] Saving session information...")

save_session_info(filename = "sessionInfo_CellChat_LPL_CD.txt", dir = paths$logs, label = "CellChat cell-cell communication analysis - LPL compartment (Crohn's disease)")

message("[OUTPUT] Session information saved to: ", paths$logs)


# =========================================================
#                  Final pipeline message
# =========================================================

message("=================================================")
message("[PIPELINE] CellChat analysis completed successfully for LPL compartment (Crohn's disease condition).")
message("=================================================")


