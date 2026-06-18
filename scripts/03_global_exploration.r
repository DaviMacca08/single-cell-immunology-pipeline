# =========================================================
# Bulk integration workflow (Merge + PCA + Harmony)
# Input: filtered single-sample Seurat objects
# Output: integrated Seurat object + QC plots
# =========================================================


# =========================================================
#                  Libraries & Setup
# =========================================================

source("scripts/00_setup/00_paths.R")
source("scripts/00_setup/01_environment.R")
source("scripts/00_setup/02_io_helpers.R")
source("scripts/00_setup/03_checks.R")
source("scripts/00_setup/04_seed.R")

set_seed(1234)

message("=== Starting Marker Discovery Analysis ===")


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
#                 Global dimensional plots
# =========================================================

message("Generating Dimension plots...")

meta_vars <- c("sample_id", "donor", "disease", "condition", "compartment")

clusters <- DimPlot(cd_harmony, reduction = "umap", group.by = "seurat_clusters", label = TRUE)

DimPlot<- list()

for (v in meta_vars){
  
  # Check metadata exists
  if (!v %in% colnames(cd_harmony@meta.data)) {
    
    warning(paste(v, "not found in metadata"))
    next
  }
  
  # Check metadata has more than 1 level
  n_groups <- length(unique(cd_harmony@meta.data[[v]]))
  
  if (n_groups < 2) {
    
    warning(paste(v, "has less than 2 groups"))
    next
  }
  
  message("Generating Dimension plots of ", v )
  
  DimPlot[[v]] <- DimPlot(cd_harmony, group.by = v, reduction = "umap") 
  
  p <- DimPlot[[v]] + clusters
  
  # save plot
  save_plot(object = p, filename = paste0(v, "_vs_Clusters",".pdf"), dir = paths$plots_global,
            width = 14, height = 6)
  
  message("Saved plot: ", v, " vs Clusters")

  rm(p)
  
}

rm(DimPlot)

# =========================================================
#                 Find cluster markers
# =========================================================

message("Running FindAllMarkers using MAST...")

table(Idents(cd_harmony))

AllMarkers <- FindAllMarkers(
  cd_harmony, 
  logfc.threshold = 0.25,
  min.pct = 0.1, 
  only.pos = TRUE,
  test.use = "MAST"
)

# ---------------------------------
# Marker sanity checks
# ---------------------------------

stopifnot(length(unique(AllMarkers$cluster)) == length(unique(cd_harmony$seurat_clusters)))

message("Marker detection completed successfully")

# ---------------------------------
# Find Top 5 markers per cluster
# ---------------------------------

top_markers <- AllMarkers |> 
  filter(p_val_adj < 0.05,
         pct.1 > 0.25,
         pct.1 > pct.2,
         !grepl("^RPL|^RPS|^MT-|^HSP|FOS|JUN", gene)) |> 
  mutate(score = avg_log2FC * (pct.1 - pct.2)) |> 
  group_by(cluster) |> 
  slice_max(order_by = score, n = 5) 


# =========================================================
#               Top markers per cluster
# =========================================================

message("Selecting top markers...")

top_markers <- AllMarkers |> 
  filter(p_val_adj < 0.05,
         pct.1 > 0.25,
         pct.1 > pct.2,
         !grepl("^RPL|^RPS|^MT-|^HSP|FOS|JUN", gene)) |> 
  mutate(score = avg_log2FC * (pct.1 - pct.2)) |> 
  group_by(cluster) |> 
  slice_max(order_by = score, n = 5) 

top_genes<- unique(top_markers$gene)

check_features(top_genes, cd_harmony)


# =========================================================
#                 Canonical literature markers
# =========================================================

canonical_markers <- c(
  # T cells
  "CD3D","CD3E","TRAC",
  # CD4 / naïve / memory
  "CD4","IL7R","CCR7","LEF1","TCF7","SELL",
  # CD8 cytotoxic
  "CD8A","CD8B","NKG7","GZMK","GZMB","PRF1","GNLY",
  # Treg
  "FOXP3","IL2RA","CTLA4",
  # exhaustion (CD relevance)
  "PDCD1","LAG3","TIGIT","TOX","HAVCR2",
  # tissue resident / IEL vs LPL
  "ITGAE","CD69","ZNF683","CXCR6","RUNX3",
  # IEL / γδ / innate-like
  "TRDC","TRGC1","TRGC2","KLRD1",
  # inflammation (Crohn axis)
  "IL17A","RORC","IL23R","IFNG","TNF"
)

feature_check <- check_features(canonical_markers, cd_harmony)

final_markers <- intersect(top_genes, canonical_markers)


# =========================================================
#             Biological interpretation plots
# =========================================================

message("Generating biological interpretation plots...")

# ---------------------------------------------------------
# DotPlots
# ---------------------------------------------------------

dot_canonical_compartment <- DotPlot(cd_harmony, features = canonical_markers, group.by = "compartment", scale = FALSE) + RotatedAxis() +
  ggtitle("Canonical Immune Marker Expression Across Tissue Compartments")

dot_canonical_disease <- DotPlot(cd_harmony, features = canonical_markers, group.by = "disease", scale = FALSE) + RotatedAxis() + 
  ggtitle("Canonical Immune Marker Expression Across Disease States")

dot_canonical_clusters <- DotPlot(cd_harmony, features = canonical_markers, scale = TRUE) + RotatedAxis() +
  ggtitle("Canonical Immune Marker Expression Across Seurat Clusters")

save_plot(dot_canonical_compartment, filename = "DotPlot_Canonical_Compartment.pdf", dir = paths$plots_global,  width = 16, height = 8)

save_plot(dot_canonical_disease, filename = "DotPlot_Canonical_Disease.pdf", dir = paths$plots_global, width = 16, height = 8)

save_plot(dot_canonical_clusters, filename = "DotPlot_Canonical_Clusters.pdf", dir = paths$plots_global, width = 16,  height = 8)

# ---------------------------------------------------------
# Violin plots
# ---------------------------------------------------------

vln_compartment <- VlnPlot(cd_harmony, features = c("CD8A", "CCR7", "GNLY") , group.by = "compartment",  ncol = 3) + 
  ggtitle("Representative T-Cell Marker Expression Across Tissue Compartments")

vln_disease <- VlnPlot(cd_harmony, features = c("IL7R","PDCD1", "IL17A") , group.by = "disease",  ncol = 3) + 
  ggtitle("Representative Immune Marker Expression Across Disease States")

vln_clusters <- VlnPlot(cd_harmony, features = c("IL7R", "NKG7", "ITGAE") , group.by = "seurat_clusters",  ncol = 3, pt.size = 0.1) +
  ggtitle("Representative Marker Expression Across Seurat Clusters")

save_plot(vln_compartment, filename = "VlnPlot_Compartment.pdf", dir = paths$plots_global, width = 15, height = 6)

save_plot(vln_disease, filename = "VlnPlot_Disease.pdf", dir = paths$plots_global, width = 15, height = 6)

save_plot(vln_clusters, filename = "VlnPlot_Clusters.pdf", dir = paths$plots_global, width = 15, height = 6)

# ---------------------------------------------------------
# Feature plots
# ---------------------------------------------------------

feature_sets <- list(
  Naive_vs_Cytotoxic = c("IL7R","NKG7","ITGAE","LEF1"), # Naive and Cytotoxic T cells
  Cytotoxic = c( "NKG7","GZMK","GZMB","PRF1"), # Cytotoxic / Effector T cells
  Tissue_Resident = c("ITGAE","CXCR6","ZNF683","RUNX3"), # Tissue Resident / IEL - LPL
  Treg = c("FOXP3","IL2RA","CTLA4","TIGIT") # Treg
)

plot_titles <- c(
  Naive_vs_Cytotoxic = "Naive/Memory vs Cytotoxic T-cell Programs",
  Cytotoxic = "Cytotoxic Effector T-cell Markers",
  Tissue_Resident = "Tissue-Resident Memory (TRM) Cell Markers",
  Treg = "Regulatory T-cell (Treg) Markers"
)

for (plot_name in names(feature_sets)) {
  
  current_features <- feature_sets[[plot_name]]
  
  p <- FeaturePlot(cd_harmony, features = current_features, reduction = "umap", cols = c("lightgray", "red")) +
  plot_annotation(title = plot_titles[[plot_name]])
  
  save_plot(p, filename = paste0("FeaturePlot_", plot_name, ".png"), dir = paths$plots_global, width = 12, height = 10 )
  
  rm(p)
  
}


# =========================================================
#              Data-driven marker visualization
# =========================================================

# ---------------------------------------------------------
# DotPlots
# ---------------------------------------------------------

message("Generating top marker visualizations...")

dot_topgenes_compartment <- DotPlot(cd_harmony, features = top_genes, group.by = "compartment", scale = FALSE) + RotatedAxis() +
  ggtitle("Top Cluster Marker Expression Across Tissue Compartments")

dot_topgenes_disease <- DotPlot(cd_harmony, features = top_genes, group.by = "disease", scale = FALSE) + RotatedAxis() +
  ggtitle("Top Cluster Marker Expression Across Disease Groups")

dot_topgenes_clusters <- DotPlot(cd_harmony, features = top_genes, scale = TRUE) + RotatedAxis() +
  ggtitle("Top Cluster Markers Across Seurat Clusters")

save_plot(dot_topgenes_compartment, filename =  "DotPlot_TopGenes_Compartment.pdf", dir = paths$plots_global,  width = 16, height = 8)

save_plot(dot_topgenes_disease, filename = "DotPlot_TopGenes_Disease.pdf", dir = paths$plots_global, width = 16, height = 8)

save_plot(dot_topgenes_clusters, filename = "DotPlot_TopGenes_Clusters.pdf", dir = paths$plots_global, width = 16,  height = 8)

# ---------------------------------------------------------
# Heatmap
# ---------------------------------------------------------

message("Generating heatmap...")

heatmap_topgenes_obj <- subset(cd_harmony, downsample = 200)

plot_heatmap <- DoHeatmap(heatmap_topgenes_obj, features = top_genes) + 
  theme(axis.text.y = element_text(size = 5)) + 
  NoLegend() +
  ggtitle("Cluster-Specific Marker Gene Expression Heatmap")

save_plot(plot_heatmap, filename = "Heatmap_TopGenes.pdf", dir = paths$plots_global, width = 14, height = 12)

rm(heatmap_topgenes_obj)


# =========================================================
#           Consensus markers (data x literature)
# =========================================================

consensus_dotplot <- DotPlot(cd_harmony, features = final_markers, scale = TRUE) + RotatedAxis() +
  ggtitle("Consensus Marker Expression Across Clusters")

consensus_featureplot <- FeaturePlot(cd_harmony, features = final_markers, reduction = "umap",cols = c("lightgray", "red")) +
  plot_annotation(title = "Spatial Distribution of Consensus Markers on UMAP")

save_plot(consensus_dotplot, filename = "DotPlot_Consensus.pdf", dir = paths$plots_global, width = 14, height = 8)

save_plot(consensus_featureplot, filename = "FeaturePlot_Consensus.pdf", dir = paths$plots_global, width = 16, height = 12)


# =========================================================
#                  Save outputs
# =========================================================

save_csv(AllMarkers, filename = "AllMarkers_clusters.csv", dir = paths$tables_markers)

message("Global exploration pipeline completed successfully")

# =========================================================
#                  Save session info
# =========================================================

message("Saving session information (Global Exploration stage)...")

save_session_info(filename = "sessionInfo_GlobalExploration.txt", dir = paths$logs, label = "Global Exploration stage")

message("Session information saved at: ", paths$logs)


# =========================================================
#                  Final pipeline message
# =========================================================

message("=================================================")
message("PIPELINE STEP COMPLETED: Global Exploration")
message("=================================================")