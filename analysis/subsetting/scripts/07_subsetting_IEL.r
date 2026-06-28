# =========================================================
# Subsetting cells from IEL
# =========================================================
# Consensus Cell Type Annotation for IEL compartment
# Input:
#   - Seurat object with consensus annotation (PBMC Database)
#
# Output:
#   - IEL seurat object with consensus annotation (After reclustering - No Harmony integration)
#   - IEL cluster marker genes table
#   - Diagnostic plots
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

message("Starting subsetting cells from IEL compartment...")


# =========================================================
#                   Load input data
# =========================================================

message("Loading Seurat object...")

cd_harmony <- readRDS(
  file.path(paths$objects_global, "srt_consensus_celltype.rds")
)


# =========================================================
#                  Initial sanity checks
# =========================================================

message(" Running sanity checks...")

stopifnot(inherits(cd_harmony, "Seurat"))
stopifnot("seurat_clusters" %in% colnames(cd_harmony@meta.data))
stopifnot("celltype_final" %in% colnames(cd_harmony@meta.data))


# =========================================================
#          Inspect Cell Type Distribution
# =========================================================

message(" Cell type abundance across tissue compartments")

table(
  cd_harmony$celltype_final,
  cd_harmony$compartment
)

message("Cell type abundance across samples")

table(
  cd_harmony$sample_id,
  cd_harmony$celltype_final
)

message("Disease status distribution across compartments")

table(
  cd_harmony$disease,
  cd_harmony$compartment
)


# =========================================================
#           Calculate Cell Type Composition
# =========================================================

message("Calculating cell type frequencies for each sample...")

composition_table <- cd_harmony@meta.data |> 
  group_by(
    sample_id,
    donor,
    disease,
    compartment,
    celltype_final
  ) |> 
  summarise(n = n(), .groups = "drop") |> 
  group_by(sample_id) |> 
  mutate(
    frequency = n / sum(n),
    percentage = frequency * 100
  )

# -----------------------
#  Stacked Barplot
# -----------------------

stacked_barplot <- ggplot(
  composition_table,
  aes(
    x = sample_id,
    y = percentage,
    fill = celltype_final
  )
) +
  geom_col() +
  theme_bw() +
  labs(
    x = "Sample",
    y = "Cell proportion (%)",
    fill = "Cell type",
    title = "Cell Type Composition Across Samples"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

# -----------------------
#  Celltype Proportion Boxplot
# -----------------------

celltype_proportion_boxplot <- ggplot(
  composition_table,
  aes(
    x = disease,
    y = percentage,
    fill = disease
  )
) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.7
  ) +
  geom_jitter(
    width = 0.15,
    size = 2
  ) +
  scale_fill_manual(
    values = c(
      "Control" = "#4DBBD5",
      "CD"      = "#E64B35"
    )
  ) + 
  facet_wrap(~ celltype_final, scales = "fixed") +
  theme_linedraw() +
  labs(
    x = NULL,
    y = "Cell proportion (%)",
    title = "Cell Type Percentage Across Conditions"
  ) +
  theme(
    legend.position = "none"
  )


# -----------------------
#  Barplot Compartment
# -----------------------

celltype_compartment_boxplot <- ggplot(
  composition_table,
  aes(
    x = compartment,
    y = percentage,
    fill = compartment
  )
) +
  geom_boxplot(
    outlier.shape = NA,
    alpha = 0.7
  ) +
  geom_jitter(
    width = 0.15,
    size = 2
  ) +
  scale_fill_manual(
    values = c(
      "IEL" = "#00A087",
      "LPL" = "#3C5488"
    )
  ) + 
  facet_wrap(~ celltype_final, scales = "fixed") +
  theme_linedraw() +
  labs(
    x = NULL,
    y = "Cell proportion (%)",
    title = "Cell Type Proportions in IEL and LPL Compartments"
  ) +
  theme(
    legend.position = "none"
  )

# -------------------
#  Save final plots
# -------------------

save_plot(stacked_barplot, filename = "cell_type_composition_per_sample.png", dir = paths$plots_IEL_compartment, 
          width = 10, height = 6)
rm(stacked_barplot)

save_plot(celltype_proportion_boxplot, filename = "CellType_Abundance_By_Disease.png", dir = paths$plots_IEL_compartment, 
          width = 10, height = 6)
rm(celltype_proportion_boxplot)

save_plot(celltype_compartment_boxplot, filename = "CellType_Abundance_By_Compartment.png", dir = paths$plots_IEL_compartment, 
          width = 10, height = 6)
rm(celltype_compartment_boxplot)


# =========================================================
#               Subsetting IEL cells  
# =========================================================

message("Subsetting all cells from IEL compartment...")

iel_cells <- subset(
  cd_harmony,
  subset = compartment == "IEL"
)

DefaultAssay(iel_cells) <- "RNA"

# -------------------
paste0("IEL compartment: ", ncol(iel_cells), " cells")
# -------------------

message("Normalizing data...")

iel_cells <- NormalizeData(iel_cells, normalization.method = "LogNormalize", scale.factor = 10000)

message("Finding HVGs...")

iel_cells <- FindVariableFeatures(iel_cells, selection.method = "vst", nfeatures = 3000)

message("Scaling data...")

iel_cells <- ScaleData(
  iel_cells,
  features = VariableFeatures(iel_cells),
  vars.to.regress = "nCount_RNA"
)

message("Running PCA...")

iel_cells <- RunPCA(
  iel_cells,
  features = VariableFeatures(iel_cells)
)

print(iel_cells[["pca"]], dims = 1:15, nfeatures = 10)

message("Elbow Plot...")

print(ElbowPlot(iel_cells, ndims = 50) + ggtitle("Elbow Plot"))

message("Run UMAP...")

iel_cells <- RunUMAP(iel_cells, dims = 1:30, reduction = "pca")


# =========================================================
#                   Visualize data 
# =========================================================

message("Dimplot of subset dataset...")

sampleid_dimplot <- DimPlot(iel_cells, reduction = "umap", group.by = "sample_id", label.size = 5)

disease_dimplot <- DimPlot(iel_cells, reduction = "umap", group.by = "disease", label.size = 5)

condition_dimplot <- DimPlot(iel_cells, reduction = "umap", group.by = "condition", label.size = 5)

# ----------------
# Save plots
#-----------------

save_plot(sampleid_dimplot, filename = "DimPlot_IEL_sampleID.pdf", dir = paths$plots_IEL_compartment,
          width = 15, height = 16)

save_plot(disease_dimplot, filename = "DimPlot_IEL_disease.pdf", dir = paths$plots_IEL_compartment,
          width = 15, height = 16)

save_plot(condition_dimplot, filename = "DimPlot_IEL_condition.pdf", dir = paths$plots_IEL_compartment,
          width = 15, height = 16)


# =========================================================
#                  Find Clusters of IEL dataset
# =========================================================

message("Find neighboors...")

iel_cells <- FindNeighbors(iel_cells, dims = 1:30, reduction = "pca")

message("Find clusters...")

iel_cells <- FindClusters(iel_cells, resolution = 0.4)

paste0("Found total cluster: ", length(unique(iel_cells$seurat_clusters)))

print("--- Number of cells per cluster ---\n")

table(Idents(iel_cells))

Idents(iel_cells) <- "seurat_clusters"


# =========================================================
#                 Find cluster markers
# =========================================================

message("Running FindAllMarkers using MAST...")

iel_cells_AllMarkers <- FindAllMarkers(
  iel_cells, 
  logfc.threshold = 0.25,
  min.pct = 0.1, 
  only.pos = TRUE,
  test.use = "MAST"
)

top_10markers <- iel_cells_AllMarkers |> 
  filter(p_val_adj < 0.05,
         pct.1 > 0.25,
         abs(avg_log2FC) > 0.25,
         !grepl("^RPL|^RPS|^MT-|^HSP|FOS|JUN", gene)) |> 
  mutate(score = avg_log2FC * (pct.1 - pct.2)) |> 
  group_by(cluster) |> 
  slice_max(order_by = score, n = 10) 

check_features(top_10markers$gene, iel_cells)


# =========================================================
#                 Diagnosis of results
# =========================================================

# ---------------------------------------------------------
# DimPlot of identified clusters
# ---------------------------------------------------------

message("Visualize clusters...")

clusters <- DimPlot(iel_cells, reduction = "umap", repel = TRUE, label = TRUE, label.size = 5) 

save_plot(clusters, filename = "DimPlot_IEL_clusters.pdf", dir = paths$plots_IEL_compartment,
          width = 15, height = 11)

# ---------------------------------------------------------
# Heatmap
# ---------------------------------------------------------

message("Generating heatmap of 10 top maerkers per clusters...")

genes_ordered <- top_10markers |>
  arrange(cluster, desc(avg_log2FC)) |>
  pull(gene)

iel_cells_heatmap <- iel_cells 

iel_cells_heatmap <- ScaleData(
  iel_cells_heatmap,
  features = genes_ordered
)
plot_heatmap <- DoHeatmap(iel_cells_heatmap, features = genes_ordered, raster = FALSE) + 
  theme(axis.text.y = element_text(size = 12, face = "bold")) + 
  NoLegend() +
  ggtitle("Cluster-Specific Marker Gene Expression Heatmap")

rm(iel_cells_heatmap)

save_plot(plot_heatmap, filename = "Heatmap_TopGenes_IEL.png", dir = paths$plots_IEL_compartment, width = 20, height = 15)

# ---------------------------------------------------------
# DimPLots
# ---------------------------------------------------------

iel_cluster_DimPlot <- DimPlot(iel_cells, reduction = "umap", label = TRUE, repel = TRUE)

iel_sampleid_DimPlot <- DimPlot(iel_cells, reduction = "umap", group.by = "sample_id")

iel_disease_DimPlot <- DimPlot(iel_cells, reduction = "umap", group.by = "disease")

save_plot(iel_cluster_DimPlot, filename = "DimPlot_IEL_clusters.pdf", dir = paths$plots_IEL_compartment,
          width = 14, height = 6)

save_plot(iel_sampleid_DimPlot, filename = "DimPlot_IEL_sampleID.pdf", dir = paths$plots_IEL_compartment,
          width = 14, height = 6)

save_plot(iel_disease_DimPlot, filename = "DimPlot_IEL_disease.pdf", dir = paths$plots_IEL_compartment,
          width = 14, height = 6)

# ---------------------------------------------------------
# DotPlot
# ---------------------------------------------------------

markers <- c(
  # CD8 T cells (core identity + IEL residency + effector)
  "CD8A","CD8B","ITGAE","CD69","CXCR6",
  "GZMB","GZMK","PRF1","IFNG",
  "KLRG1","FCRL6","KLF2",
  # CD4 T cells (core identity + naive/memory)
  "CD4","CCR7","SELL","TCF7","LEF1",
  # TH17 lineage (paper-specific core)
  "RORC","CCR6","RORA","IL17A","IL23R",
  # Treg lineage (canonical suppressive program)
  "FOXP3","IL2RA","CTLA4","IL10",
  # γδ T / innate-like cytotoxic program 
  "TRDC","TRGC1","TRGC2",
  # TFH lineage (germinal center–like program)
  "CXCL13","BCL6","CXCR5",
  # activation / exhaustion-like state (IEL activated CD8 + CD4)
  "PDCD1","TIGIT","ICOS",
  # proliferation
  "MKI67"
)

p <- DotPlot(iel_cells, features = markers, scale = TRUE) + RotatedAxis() + 
  plot_annotation(title = "IEL cell subset marker expression across annotated cell populations")

save_plot(p, filename = "DotPlot_IEL_Markers.pdf", dir = paths$plots_IEL_compartment,
          width = 18, height = 10)
rm(p)


# =========================================================
#             Feature Plot for cell type discovery
# =========================================================

# ---------------------------------------------------------
# T cell identity markers
# ---------------------------------------------------------

p <- FeaturePlot(iel_cells, features = c("CD3D","TRAC"), reduction = "umap", cols = c("lightgray", "red")) +
  plot_annotation(title = "T cell identity markers")

save_plot(p, filename = "FeaturePlot_Identity.pdf", dir = paths$plots_IEL_compartment,
          width = 14, height = 10)
rm(p)

# ---------------------------------------------------------
# CD4/CD8 lineage markers
# ---------------------------------------------------------

p <- FeaturePlot(iel_cells, features = c("CD4","CD8A","CD8B"), reduction = "umap", cols = c("lightgray", "red")) + 
  plot_annotation(title = "CD4/CD8 lineage markers")

save_plot(p, filename = "FeaturePlot_CD4vsCD8.pdf", dir = paths$plots_IEL_compartment,
          width = 14, height = 10)
rm(p)

# ---------------------------------------------------------
# CD8 IEL residency / effector program
# ---------------------------------------------------------

p <- FeaturePlot(iel_cells,features = c(
  "ITGAE","CD69","CXCR6",
  "KLRG1","GZMB","GZMK",
  "PRF1","IFNG","FCRL6","KLF2"),
  reduction = "umap", cols = c("lightgrey", "red")) +
  plot_annotation(title = "CD8 IEL residency & effector program")

save_plot(p, filename = "FeaturePlot_CD8.png", dir = paths$plots_IEL_compartment,
          width = 18, height = 15)
rm(p)

# ---------------------------------------------------------
# CD4 functional lineages
# ---------------------------------------------------------

p <- FeaturePlot(iel_cells,features = c(
  "CCR7","SELL","TCF7","LEF1",
  "RORC","CCR6","RORA","IL17A","IL23R",
  "FOXP3","IL2RA","CTLA4","IL10",
  "CXCL13","BCL6","CXCR5"),
  reduction = "umap",cols = c("lightgrey", "red")) +
  plot_annotation(title = "CD4 T cell functional programs")

save_plot(p, filename = "FeaturePlot_CD4.png", dir = paths$plots_IEL_compartment,
          width = 20, height = 15)
rm(p)

# ---------------------------------------------------------
# activation / proliferation
# ---------------------------------------------------------

p <- FeaturePlot(iel_cells, features = c(
  "PDCD1","TIGIT","ICOS",
  "TOX","BATF","MAF",
  "MKI67"),
  reduction = "umap", cols = c("lightgrey", "red")) +
  plot_annotation("Activation & proliferation")

save_plot(p, filename = "FeaturePlot_Activ_Prol.png", dir =paths$plots_IEL_compartment,
          width = 18, height = 15)
rm(p)

# -----------------
# Contamination
# -----------------

p <- FeaturePlot(iel_cells,features = c("MS4A1","CD79A","LYZ","MKI67", "EPCAM"),
  reduction = "umap", cols = c("lightgrey", "red")) +
  plot_annotation(title = "Contamination / non-T cells")

save_plot(p, filename = "FeaturePlot_contamination.png", dir = paths$plots_IEL_compartment,
          width = 16, height = 14)
rm(p)

# =========================================================
#                 Manual Annotation
# =========================================================

message("Loading manual cell type annotation table...")

celltype <- read.csv(file = file.path(paths$tables, "IEL_celltype/IEL_celltype_v2.0.csv"))

message("Mapping cluster IDs to manual cell type annotations...")

new_ids <- celltype$final_annotation

message("Renaming Seurat identities using manual annotations...")

names(new_ids) <- celltype$cluster

iel_cells <- RenameIdents(iel_cells, new_ids)

iel_cells$celltype_IEL <- Idents(iel_cells)

message("Generating UMAP colored by manual cell type annotations and seurat clusters...")

dimplot_celltype <- DimPlot(iel_cells, reduction = "umap", label = TRUE, repel = TRUE, label.size = 4) +
  plot_annotation(title = "IEL cell compartment") +
  xlab("UMAP 1") +
  ylab("UMAP 2")

dimplot_clusters <- DimPlot(iel_cells, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE, label.size = 5) + 
  plot_annotation(title = "IEL Clusters") +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  NoLegend() 


# -----------------
# Save Final Plot
# -----------------

save_plot(dimplot_celltype, filename = "DimPlot_IEL_celltype.png", dir = paths$plots_IEL_compartment,
          width = 18, height = 16)

save_plot(dimplot_clusters, filename = "DimPlot_IEL_clusters.pdf", dir = paths$plots_IEL_compartment,
          width = 18, height = 16)


# =========================================================
#                Save outputs
# =========================================================

message("Saving annotated IEL compartment object...")

save_rds(iel_cells, filename = "srt_IEL_compartment_celltype.rds", dir = paths$objects_subsetting)

message("Saving IEL cluster marker genes table...")

save_csv(iel_cells_AllMarkers, filename = "AllMarkers_clusters_Cell_IEL.csv", dir = paths$tables_markers)

message("IEL compartment preprocessing and subsetting pipeline completed successfully.")

# =========================================================
#                  Save session info
# =========================================================

message("Saving session information for IEL consensus cell type annotation...")

save_session_info(filename = "sessionInfo_Consensus_CellType_IEL.txt", dir = paths$logs, label = "IEL consensus cell type annotation")

message("Session information saved at: ", paths$logs)


# =========================================================
#                  Final pipeline message
# =========================================================

message("=================================================")
message("PIPELINE STEP COMPLETED: Consensun Celltype Annotations for IEL compartment")
message("=================================================")

