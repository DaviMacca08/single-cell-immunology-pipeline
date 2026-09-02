# =========================================================
# CellChat cell-cell communication analysis - IEL pathway comparison
# =========================================================
# Input:
#   - cellchat_IEL_merged.rds
#
# Output:
#   - Comparative pathway-level analysis
#   - Differential signaling analyses
#   - Selected pathway visualizations
#   - Ligand-receptor interaction analysis
# =========================================================


# =========================================================
#                  Libraries & Setup
# =========================================================

source("Setup_Environment/00_paths.R")
source("Setup_Environment/01_environment.R")
source("Setup_Environment/02_io_helpers.R")
source("Setup_Environment/03_checks.R")
source("Setup_Environment/04_seed.R")

reticulate::py_require(c("umap-learn", "numpy<2.5"))

set_seed(1234)

message("Starting CellChat pathway comparison for IEL compartment (Control vs CD)...")


# =========================================================
# Load merged and condition-specific CellChat objects
# =========================================================

message("Loading merged and condition-specific CellChat objects...")

cellchat_merged <- readRDS(file.path(paths$objects_cellchat, "cellchat_IEL_merged.rds"))
cellchat_iel_ctrl <- readRDS(file.path(paths$objects_cellchat, "cellchat_IEL_CTRL.rds"))
cellchat_iel_cd   <- readRDS(file.path(paths$objects_cellchat, "cellchat_IEL_CD.rds"))

message("CellChat objects loaded successfully.")

# ---------------------------------------
# Align cell groups across conditions
# ---------------------------------------

message("Aligning cell groups across conditions...")

group_new <- levels(cellchat_merged@idents$joint)

message("Number of cell groups in merged object: ",length(group_new))

cellchat_iel_ctrl <- liftCellChat(cellchat_iel_ctrl, group_new)
cellchat_iel_cd   <- liftCellChat(cellchat_iel_cd, group_new)

stopifnot(identical(levels(cellchat_iel_ctrl@idents), levels(cellchat_iel_cd@idents)))

message("Cell group levels aligned successfully.")

# ---------------------------------------
# Compute signaling network centrality
# ---------------------------------------

message("Computing signaling network centrality for Control and CD...")

cellchat_iel_ctrl <- netAnalysis_computeCentrality(cellchat_iel_ctrl, slot.name = "netP")
cellchat_iel_cd <- netAnalysis_computeCentrality(cellchat_iel_cd, slot.name = "netP")

message("Signaling network centrality computed for both conditions.")

# --------------------------------------------
# Reconstruct condition-specific object list
# --------------------------------------------

object_list <- list(Control = cellchat_iel_ctrl, CD = cellchat_iel_cd)

message("Condition-specific CellChat object list created: ",paste(names(object_list), collapse = ", "))

message("Number of cell groups per condition: ", length(levels(object_list[[1]]@idents)))


# =========================================================
#       Rank signaling pathways by information flow
# =========================================================

message("Ranking signaling pathways by information flow...")

# -----------------------
# Stacked comparison
# -----------------------

open_png(filename = "CellChat_IEL_rankNet_comparison_stacked.png", dir = paths$plots_cellchat_iel_merge,
         width = 1200, height = 1800)

p <- rankNet(cellchat_merged, mode = "comparison", stacked = TRUE, do.stat = TRUE) +
  ggplot2::labs(
    title = "Signaling pathway information flow",
    subtitle = "CTRL vs CD"
  ) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      hjust = 0.5,
      face = "bold",
      size = 18
    ),
    plot.subtitle = ggplot2::element_text(
      hjust = 0.5,
      size = 14
    )
  )

print(p)
rm(p)

close_png()

# -----------------------
# Non-stacked comparison
# -----------------------

open_png(filename = "CellChat_IEL_rankNet_comparison_notstacked.png", dir = paths$plots_cellchat_iel_merge,
         width = 1200, height = 1800)

p <- rankNet(cellchat_merged, mode = "comparison", stacked = FALSE, do.stat = TRUE) +
  ggplot2::labs(
    title = "Signaling pathway information flow",
    subtitle = "CTRL vs CD"
  ) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      hjust = 0.5,
      face = "bold",
      size = 18
    ),
    plot.subtitle = ggplot2::element_text(
      hjust = 0.5,
      size = 14
    )
  )

print(p)
rm(p)

close_png()

message("Signaling pathway ranking plots saved.")


# =========================================================
#        Differential signaling changes
# =========================================================

message("Assessing differential signaling changes between conditions...")

celltypes_of_interest <- c(
  "FOXP3+ IL2RA+ Treg",
  "CST3+ LYZ+ macrophages",
  "TH17-like CD4 T cells",
  "CD39+ tissue-resident CD8 T cells",
  "Epithelial cells / enterocytes",
  "GZMK+ effector memory CD8 T cells",  
  "Effector T cell mixed state"          
)

stopifnot(all(celltypes_of_interest %in% levels(cellchat_merged@idents$joint)))

for (ct in celltypes_of_interest) {
  
  message("Generating signaling changes plot for: ", ct)
  
  fname <- gsub("[^A-Za-z0-9]+", "_", ct)
  
  open_png(filename = paste0("CellChat_IEL_signalingChanges_", fname, ".png"), dir = paths$plots_cellchat_iel_merge_celltype_interes,
           width = 1200, height = 1000)
  
  p <- netAnalysis_signalingChanges_scatter(cellchat_merged, idents.use = ct)
 
  print(p)
  rm(p)
  
  close_png()
}


# =========================================================
#      Functional similarity between conditions
# =========================================================

set_seed(1234)

message("Computing pairwise functional similarity between conditions...")

cellchat_merged <- computeNetSimilarityPairwise(cellchat_merged, type = "functional")

message("Computing functional network embedding...")

cellchat_merged <- netEmbedding(cellchat_merged, type = "functional", umap.method = "uwot")

message("Clustering functional network patterns...")

cellchat_merged <- netClustering(cellchat_merged,type = "functional")

message("Generating functional network similarity visualizations...")


# -----------------------
# Pairwise functional embedding
# -----------------------

open_png(filename = "CellChat_IEL_embeddingPairwise_functional.png", dir = paths$plots_cellchat_iel_merge,
  width = 1200,height = 1000)

netVisual_embeddingPairwise(cellchat_merged, type = "functional", title = "Functional similarity of signaling networks\nCTRL vs CD")

close_png()


# -----------------------
# Functional similarity ranking
# -----------------------

open_png(filename = "CellChat_IEL_rankSimilarity_functional.png",dir = paths$plots_cellchat_iel_merge,
  width = 1200, height = 1000)

rankSimilarity(cellchat_merged, type = "functional", title =   "Functional similarity ranking of signaling pathways\nCTRL vs CD")

close_png()

message("Functional network similarity analysis completed.")


# =========================================================
#           Visualize selected signaling pathways
# =========================================================

message("Generating visualizations for selected signaling pathways...")

pathways_selected <- c("MHC-I", "MHC-II", "LIGHT", "MIF", "GALECTIN", "CDH1")

for (pw in pathways_selected) {
  
  fname <- gsub("[^A-Za-z0-9]+", "_", pw)
  
  message("Processing signaling pathway: ", pw)
  
  # ------------------------
  # Comparative bubble plot 
  # ------------------------
  
  open_png(filename = paste0("CellChat_IEL_bubble_", fname, ".png"),
           dir = paths$plots_cellchat_iel_merge_pathways_bubble, width = 3000, height = 2800)
  
  p <- netVisual_bubble(cellchat_merged, signaling = pw,
                        sources.use = celltypes_of_interest,
                        targets.use = celltypes_of_interest,
                        comparison = c(1, 2), angle.x = 45)
  print(p)
  
  close_png()
  
  # ---------------------------------------------------------
  # Ligand-receptor gene expression
  # ---------------------------------------------------------
  
  open_png(filename = paste0("CellChat_IEL_geneExpr_", fname, ".png"),
           dir = paths$plots_cellchat_iel_merge_pathways_violin, width = 2000, height = 1400)
  
  p <- plotGeneExpression(cellchat_merged, signaling = pw, split.by = "datasets",
                          colors.ggplot = TRUE)
  print(p)
  
  close_png()
  
  # ---------------------------------------------------------
  # Chord diagrams by condition
  # ---------------------------------------------------------
  
  for (i in seq_along(object_list)) {
    
    cond_name <- names(object_list)[i]
    
    if (!(pw %in% object_list[[i]]@netP$pathways)) {
      message("  Pathway ", pw, " absent in ", cond_name, " - skipping chord diagram")
      next
    }
    
    message("Generating chord diagram for ", cond_name, "...")
    
    open_png(filename = paste0("CellChat_IEL_chord_", fname, "_", cond_name, ".png"),
             dir = paths$plots_cellchat_iel_merge_pathways_chord, width = 1600, height = 1600)
    
    par(mar = c(1,1,4,1), xpd = TRUE)
    netVisual_aggregate(object_list[[i]], signaling = pw, layout = "chord",
                        title.name = paste0(pw, " - ", cond_name))
    
    close_png()
  }
}

message("Selected pathway visualizations completed.")


# =========================================================
# Hierarchy plots for selected signaling pathways
# =========================================================

message("Generating hierarchy plots for selected signaling pathways...")

vertex_epithelial   <- which(group_new == "Epithelial cells / enterocytes")
vertex_macrophages  <- which(group_new == "CST3+ LYZ+ macrophages")
vertex_cd8_trm      <- which(group_new == "CD39+ tissue-resident CD8 T cells")
vertex_gd_cytotoxic <- which(group_new == "Cytotoxic γδ T cells")
vertex_th17cd4      <- which(group_new == "TH17-like CD4 T cells")
vertex_treg         <- which(group_new == "FOXP3+ IL2RA+ Treg")

pathway_receiver_map <- list(
  "CDH1"   = c(vertex_cd8_trm, vertex_gd_cytotoxic),  # CD103+ recipient populations; epithelial cells are the main ligand source
  "LIGHT"  = c(vertex_epithelial, vertex_macrophages), # Recipient populations expressing LTBR
  "MHC-II" = c(vertex_th17cd4, vertex_treg),           # CD4+ recipient populations; macrophages act as the sender/APC population
  "MIF"    = c(vertex_macrophages, vertex_epithelial),  # Recipient populations expressing CD74
  "MHC-I"  = c(vertex_cd8_trm, vertex_gd_cytotoxic)    # Recipient populations; receiver structure shifts between conditions 
)

stopifnot(all(lengths(pathway_receiver_map) >= 2))

for (pw in names(pathway_receiver_map)) {
  
  message("Processing hierarchy plot for pathway: ", pw)
  
  fname <- gsub("[^A-Za-z0-9]+", "_", pw)
  vr <- pathway_receiver_map[[pw]]
  
  for (i in seq_along(object_list)) {
    
    cond_name <- names(object_list)[i]
    
    if (!(pw %in% object_list[[i]]@netP$pathways)) {
      message(
        "  Pathway ", pw, " absent in ", cond_name,
        " - skipping hierarchy plot"
      )
      next
    }
    
    open_png(filename = paste0("CellChat_IEL_hierarchy_", fname, "_", cond_name, ".png"),
      dir = paths$plots_cellchat_iel_merge_pathways_hierarchy, width = 3600, height = 1600)
    
    par(mar = c(1,1,4,1), xpd = TRUE)
    
    netVisual_aggregate(object_list[[i]], signaling = pw, layout = "hierarchy", vertex.receiver = vr)

    grid::grid.text(
      paste0(pw, " signaling hierarchy - ", cond_name),
      y = grid::unit(0.98, "npc"),
      gp = grid::gpar(
        fontsize = 18,
        fontface = "bold"
      )
    )
    
    close_png()
  }
}

message("Hierarchy plots completed.")


# =========================================================
# Communication pattern recognition
# =========================================================

message("Starting communication pattern analysis...")

foreach::registerDoSEQ()

library(NMF)

# ---------------------------------------------------------
# Select optimal number of patterns
# ---------------------------------------------------------

set_seed(1234)

message("Selecting the optimal number of communication patterns...")

for (i in seq_along(object_list)){
  
  cond_name <- names(object_list)[i]
  
  message("Evaluating communication patterns for: ", cond_name)
  
  # Outgoing communication patterns
  
  message("  Selecting K for outgoing communication...")
  
  open_png(filename = paste0("CellChat_IEL_", cond_name, "_selectK_outgoing.png"), dir = paths$plots_cellchat_iel_merge_pathways_pattern,
           width = 1200, height = 800)
  
  p <- selectK(object_list[[i]], pattern = "outgoing")
  
  print(p)
  rm(p)
  
  close_png()
  
  message("  Outgoing pattern selection plot saved.")
  
  # Incoming communication patterns
  
  open_png(filename = paste0("CellChat_IEL_", cond_name, "_selectK_incoming.png"), dir = paths$plots_cellchat_iel_merge_pathways_pattern, 
           width = 1200, height = 800)
  
  message("  Selecting K for incoming communication...")
  
  p <- selectK(object_list[[i]], pattern = "incoming")
  
  print(p)
  rm(p)
  
  close_png()
  
  message("  Incoming pattern selection plot saved.")
}


message("Communication pattern selection completed.")

# ---------------------------------------------------------
# Identify outgoing communication patterns
# ---------------------------------------------------------

K_out <- 4

future::plan("sequential")

set_seed(1234)

message("Identifying outgoing communication patterns in Control...")

open_png("CellChat_IEL_Control_pattern_outgoing.png", dir = paths$plots_cellchat_iel_merge_pathways_pattern,
         width = 1300, height = 900)

object_list[[1]] <- identifyCommunicationPatterns(
  object_list[[1]],
  pattern = "outgoing",
  k = K_out
)

close_png()

message("Identifying outgoing communication patterns in CD...")

open_png("CellChat_IEL_CD_pattern_outgoing.png", dir = paths$plots_cellchat_iel_merge_pathways_pattern,
         width = 1300, height = 700)

object_list[[2]] <- identifyCommunicationPatterns(
  object_list[[2]],
  pattern = "outgoing",
  k = K_out
)

close_png()

# ---------------------------------------------------------
# Identify incoming communication patterns
# ---------------------------------------------------------

K_in <- 6

future::plan("sequential")

set_seed(1234)

message("Identifying incoming communication patterns in Control...")

open_png("CellChat_IEL_Control_pattern_incoming.png", dir = paths$plots_cellchat_iel_merge_pathways_pattern,
         width = 1300, height = 700)

object_list[[1]] <- identifyCommunicationPatterns(
  object_list[[1]],
  pattern = "incoming",
  k = K_in
)

close_png()

message("Identifying incoming communication patterns in CD...")

open_png("CellChat_IEL_CD_pattern_incoming.png", dir = paths$plots_cellchat_iel_merge_pathways_pattern,
         width = 1300, height = 700)

object_list[[2]] <- identifyCommunicationPatterns(
  object_list[[2]],
  pattern = "incoming",
  k = K_in
)

close_png()

message("Communication pattern recognition completed.")

# ---------------------------------------------------------
# Visualize outgoing communication patterns
# ---------------------------------------------------------

library("ggalluvial")

for (i in seq_along(object_list)) {
  
  cond_name <- names(object_list)[i]
  
  open_png(filename = paste0("CellChat_IEL_", cond_name, "_communicationPatterns_outgoing_river.png"), dir = paths$plots_cellchat_iel_merge_pathways_pattern,
           width = 1800, height = 1400)
  
  p <- netAnalysis_river(object_list[[i]], pattern = "outgoing")
  
  print(p)
  rm(p)
  
  close_png()
  
  
  open_png(filename = paste0("CellChat_IEL_", cond_name, "_communicationPatterns_outgoing_dot.png"), dir = paths$plots_cellchat_iel_merge_pathways_pattern,
           width = 1800, height = 1400)
  
  p <- netAnalysis_dot(object_list[[i]], pattern = "outgoing")
  
  print(p)
  rm(p)
  
  close_png()
}

# ---------------------------------------------------------
# Visualize incoming communication patterns
# ---------------------------------------------------------

for (i in seq_along(object_list)) {
  
  cond_name <- names(object_list)[i]
  
  open_png(filename = paste0("CellChat_IEL_", cond_name, "_communicationPatterns_incoming_river.png"), dir = paths$plots_cellchat_iel_merge_pathways_pattern,
           width = 1800, height = 1400)
  
  p <- netAnalysis_river(object_list[[i]], pattern = "incoming")
  
  print(p)
  rm(p)
  
  close_png()
  
  
  open_png(filename = paste0("CellChat_IEL_", cond_name, "_communicationPatterns_incoming_dot.png"), dir = paths$plots_cellchat_iel_merge_pathways_pattern,
           width = 1800, height = 1400)
  
  p <- netAnalysis_dot(object_list[[i]], pattern = "incoming")
  
  print(p)
  rm(p)
  
  close_png()
}


# =========================================================
#  Signaling role networks for selected pathways
# =========================================================

message("Generating signaling role networks for selected pathways...")

for (pw in pathways_selected) {
  
  fname <- gsub("[^A-Za-z0-9]+", "_", pw)
  
  message("Processing signaling role network for: ", pw)
  
  for (i in seq_along(object_list)) {
    
    cond_name <- names(object_list)[i]
    
    if (!(pw %in% object_list[[i]]@netP$pathways)) {
      message(
        "  Pathway ", pw, " absent in ", cond_name,
        " - skipping signaling role network"
      )
      next
    }
    
    open_png(filename = paste0("CellChat_IEL_signalingRoleNetwork_", fname, "_", cond_name, ".png"),
      dir = paths$plots_cellchat_iel_merge_pathways_sigrole,
      width = 900, height = 800)
    
    netAnalysis_signalingRole_network(object_list[[i]], signaling = pw)
    
    grid::grid.text(
      paste0(
        pw,
        " signaling roles - ",
        cond_name
      ),
      y = grid::unit(0.98, "npc"),
      gp = grid::gpar(
        fontsize = 18,
        fontface = "bold"
      )
    )
    
    close_png()
  }
}

message("Signaling role network analysis completed.")


# =========================================================
#                Save final object
# =========================================================

message("Finalizing CellChat analysis objects...")

stopifnot(
  length(object_list) == 2,
  identical(names(object_list), c("Control", "CD"))
)

message("Saving final CellChat objects for IEL compartment...")

# ---------------------------------------------------------
# Save merged object
# ---------------------------------------------------------

stopifnot("joint" %in% names(cellchat_merged@idents))

message(
  "Final merged CellChat object (Control vs CD) contains ",
  length(levels(cellchat_merged@idents$joint)),
  " cell groups."
)

save_rds(cellchat_merged, "cellchat_IEL_merged.rds", dir = paths$objects_cellchat)

message("Merged CellChat object saved.")

# ---------------------------------------------------------
# Save condition-specific objects
# ---------------------------------------------------------

stopifnot(identical(names(object_list), c("Control", "CD")))

save_rds(object_list[[1]], "cellchat_IEL_Control_final.rds", dir = paths$objects_cellchat)

save_rds(object_list[[2]], "cellchat_IEL_CD_final.rds", dir = paths$objects_cellchat)

message("Condition-specific CellChat objects saved.")


# =========================================================
# Final analysis summary
# =========================================================

message("CellChat IEL analysis summary:")

message("  Cell groups: ",length(levels(cellchat_merged@idents$joint)))

message("  Conditions: ",paste(names(object_list), collapse = ", "))

message("  Selected signaling pathways: ", paste(pathways_selected, collapse = ", "))

message("  Outgoing communication patterns: K = ",K_out)

message("  Incoming communication patterns: K = ", K_in)


# =========================================================
# Save session information
# =========================================================

message("[OUTPUT] Saving session information...")

save_session_info(filename = "sessionInfo_CellChat_IEL_comparison.txt", dir = paths$logs, 
                  label = "CellChat cell-cell communication analysis - IEL compartment (Control vs CD)")

message("[OUTPUT] Session information saved to: ",paths$logs)


# =========================================================
# Final pipeline message
# =========================================================

message("=================================================")
message("[PIPELINE] CellChat analysis completed successfully for IEL compartment.")
message("[PIPELINE] Conditions compared: Control vs CD.")
message("[PIPELINE] Global, differential, pathway-level, and communication pattern analyses completed.")
message("=================================================")

