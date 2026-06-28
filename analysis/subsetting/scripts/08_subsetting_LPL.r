# =========================================================
# Subsetting cells from LPL
# =========================================================
# Consensus Cell Type Annotation for LPL compartment
# Input:
#   - Seurat object with consensus annotation (PBMC Database)
#
# Output:
#   - LPL seurat object with consensus annotation (After reclustering - No Harmony integration)
#   - LPL cluster marker genes table
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

message("Starting subsetting cells from LPL compartment...")


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

save_plot(stacked_barplot, filename = "cell_type_composition_per_sample.png", dir = paths$plots_LPL_compartment, 
          width = 10, height = 6)
rm(stacked_barplot)

save_plot(celltype_proportion_boxplot, filename = "CellType_Abundance_By_Disease.png", dir = paths$plots_LPL_compartment, 
          width = 10, height = 6)
rm(celltype_proportion_boxplot)

save_plot(celltype_compartment_boxplot, filename = "CellType_Abundance_By_Compartment.png", dir = paths$plots_LPL_compartment, 
          width = 10, height = 6)
rm(celltype_compartment_boxplot)


# =========================================================
#               Subsetting LPL cells  
# =========================================================

message("Subsetting all cells from LPL compartment...")

lpl_cells <- subset(
  cd_harmony,
  subset = compartment == "LPL"
)

DefaultAssay(lpl_cells) <- "RNA"

# -------------------
paste0("LPL compartment: ", ncol(lpl_cells), " cells")
# -------------------

message("Normalizing data...")

lpl_cells <- NormalizeData(lpl_cells, normalization.method = "LogNormalize", scale.factor = 10000)

message("Finding HVGs...")

lpl_cells <- FindVariableFeatures(lpl_cells, selection.method = "vst", nfeatures = 3000)

message("Scaling data...")

lpl_cells <- ScaleData(
  lpl_cells,
  features = VariableFeatures(lpl_cells),
  vars.to.regress = "nCount_RNA"
)

message("Running PCA...")

lpl_cells<- RunPCA(
  lpl_cells,
  features = VariableFeatures(lpl_cells)
)

print(lpl_cells[["pca"]], dims = 1:15, nfeatures = 10)

message("Elbow Plot...")

print(ElbowPlot(lpl_cells, ndims = 50) + ggtitle("Elbow Plot"))

message("Run UMAP...")

lpl_cells <- RunUMAP(lpl_cells, dims = 1:30, reduction = "pca")


# =========================================================
#                   Visualize data 
# =========================================================

message("Dimplot of subset dataset...")

sampleid_dimplot <- DimPlot(lpl_cells, reduction = "umap", group.by = "sample_id", label.size = 5)

disease_dimplot <- DimPlot(lpl_cells, reduction = "umap", group.by = "disease", label.size = 5)

condition_dimplot <- DimPlot(lpl_cells, reduction = "umap", group.by = "condition", label.size = 5)

# ----------------
# Save plots
#-----------------

save_plot(sampleid_dimplot, filename = "DimPlot_LPL_sampleID.pdf", dir = paths$plots_LPL_compartment,
          width = 15, height = 16)

save_plot(disease_dimplot, filename = "DimPlot_LPL_disease.pdf", dir = paths$plots_LPL_compartment,
          width = 15, height = 16)

save_plot(condition_dimplot, filename = "DimPlot_LPLL_condition.pdf", dir = paths$plots_LPL_compartment,
          width = 15, height = 16)


# =========================================================
#                  Find Clusters of LPL dataset
# =========================================================

message("Find neighboors...")

lpl_cells <- FindNeighbors(lpl_cells, dims = 1:30, reduction = "pca")

message("Find clusters...")

lpl_cells <- FindClusters(lpl_cells, resolution = 0.5)

paste0("Found total cluster: ", length(unique(lpl_cells$seurat_clusters)))

print("--- Number of cells per cluster ---\n")

table(Idents(lpl_cells))

Idents(lpl_cells) <- "seurat_clusters"


# =========================================================
#                 Find cluster markers
# =========================================================

message("Running FindAllMarkers using MAST...")

lpl_cells_AllMarkers <- FindAllMarkers(
  lpl_cells, 
  logfc.threshold = 0.25,
  min.pct = 0.1, 
  only.pos = TRUE,
  test.use = "MAST"
)

top_10markers <- lpl_cells_AllMarkers |> 
  filter(p_val_adj < 0.05,
         pct.1 > 0.25,
         abs(avg_log2FC) > 0.25,
         !grepl("^RPL|^RPS|^MT-|^HSP|FOS|JUN", gene)) |> 
  mutate(score = avg_log2FC * (pct.1 - pct.2)) |> 
  group_by(cluster) |> 
  slice_max(order_by = score, n = 30) 

check_features(top_10markers$gene, lpl_cells)


# =========================================================
#                 Diagnosis of results
# =========================================================

# ---------------------------------------------------------
# DimPlot of identified clusters
# ---------------------------------------------------------

message("Visualize clusters...")

clusters <- DimPlot(lpl_cells, reduction = "umap", repel = TRUE, label = TRUE, label.size = 5) 

save_plot(clusters, filename = "DimPlot_LPL_clusters.pdf", dir = paths$plots_LPL_compartment,
          width = 15, height = 11)

# ---------------------------------------------------------
# Heatmap
# ---------------------------------------------------------

message("Generating heatmap of 10 top maerkers per clusters...")

genes_ordered <- top_10markers |>
  arrange(cluster, desc(avg_log2FC)) |>
  pull(gene)

lpl_cells_heatmpap <- lpl_cells

lpl_cells_heatmpap <- ScaleData(
  lpl_cells_heatmpap,
  features = genes_ordered
)

plot_heatmap <- DoHeatmap(lpl_cells_heatmpap, features = genes_ordered, raster = FALSE) + 
  theme(axis.text.y = element_text(size = 12, face = "bold")) + 
  NoLegend() +
  ggtitle("Cluster-Specific Marker Gene Expression Heatmap")

rm(lpl_cells_heatmpap)

save_plot(plot_heatmap, filename = "Heatmap_TopGenes_LPL.png", dir = paths$plots_LPL_compartment, width = 20, height = 15)


# ---------------------------------------------------------
# DimPLots
# ---------------------------------------------------------

lpl_cluster_DimPlot <- DimPlot(lpl_cells, reduction = "umap", label = TRUE, repel = TRUE)

lpl_sampleid_DimPlot <- DimPlot(lpl_cells, reduction = "umap", group.by = "sample_id")

lpl_disease_DimPlot <- DimPlot(lpl_cells, reduction = "umap", group.by = "disease")

save_plot(lpl_cluster_DimPlot, filename = "DimPlot_LPL_clusters.pdf", dir = paths$plots_LPL_compartment,
          width = 14, height = 6)

save_plot(lpl_sampleid_DimPlot, filename = "DimPlot_LPL_sampleID.pdf", dir = paths$plots_LPL_compartment,
          width = 14, height = 6)

save_plot(lpl_disease_DimPlot, filename = "DimPlot_LPL_disease.pdf", dir = paths$plots_LPL_compartment,
          width = 14, height = 6)

# ---------------------------------------------------------
# DotPlot
# ---------------------------------------------------------

markers <- c(
  # CD8 T cells (effector + LP CD8 clusters 1/2/8 analogs)
  "CD8A","CD8B",
  "EOMES","TBX21",
  "KLRG1","CCL3","CCL4","CCL5",
  "IFNG","PRF1","GZMB","GZMA","GZMK","GZMH",
  # CD8 tissue-associated / ITGA1+ CD160+ / NK-like CD8 (cluster 2 & 8 LP)
  "CD160","ITGA1",
  "XCL1","XCL2",
  "KLRC1","KLRC2","KLRD1","KIR2DL4",
  "ENTPD1",
  "TRDC",
  # CD4 T cells (core)
  "CD4","CCR7","SELL","LEF1",
  # TH17 (LP quiescent + activated fraction)
  "CCR6","RORA","RORC",
  "KLRB1","LTB","CCL20",
  "IL17A","IL23R","IL26",
  "CXCR6",
  "LGALS3","GPR65",
   # Treg (FOXP3+ suppressive program)
  "FOXP3","IL2RA","CTLA4","IL10",
  "BATF","LAIR2","TNFRSF4",
  "TNFRSF9","GPX1","GLRX",
  # TFH (CXCR5+ TOX2 axis, split CXCL13+ subset)
  "CXCR5","TOX2","MAF",
  "CXCL13","PDCD1","BTLA","CD200",
  "ICOS",
  # Naive / TCM
  "TCF7","KLF2",
  # Stress / heat-shock CD4 cluster (LP-specific)
  "HSPA1A","HSPA1B","DNAJB1","JUN","TNF","IL2",
  # Proliferation 
  "MKI67","STMN1","TUBA1B","TUBB",
   # Activation / exhaustion shared axis
  "TIGIT"
)

p <- DotPlot(lpl_cells, features = markers, scale = TRUE) + RotatedAxis() + 
  plot_annotation(title = "LPL cell subset marker expression across annotated cell populations") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

save_plot(p, filename = "DotPlot_LPL_Markers.png", dir = paths$plots_LPL_compartment,
          width = 20, height = 10)
rm(p)


# =========================================================
#             Feature Plot for cell type discovery
# =========================================================

# ---------------------------------------------------------
# T cell identity markers
# ---------------------------------------------------------

p <- FeaturePlot(lpl_cells, features = c("CD3D","TRAC"), reduction = "umap", cols = c("lightgray", "red")) +
  plot_annotation(title = "T cell identity markers")

save_plot(p, filename = "FeaturePlot_Identity.png", dir = paths$plots_LPL_compartment,
          width = 14, height = 10)

rm(p)

# ---------------------------------------------------------
# CD4/CD8 lineage markers
# ---------------------------------------------------------

p <- FeaturePlot(lpl_cells, features = c("CD4","CD8A","CD8B"), reduction = "umap", cols = c("lightgray", "red")) + 
  plot_annotation(title = "CD4/CD8 lineage markers")

save_plot(p, filename = "FeaturePlot_CD4vsCD8.png", dir = paths$plots_LPL_compartment, 
          width = 14, height = 10)

rm(p)

# ---------------------------------------------------------
# CD8 LP effector + cytotoxic programs
# ---------------------------------------------------------

p <- FeaturePlot(lpl_cells,features = c(
  "KLRG1", "EOMES", "TBX21",
  "CCL3", "CCL4", "CCL5",
  "IFNG", "PRF1", "GZMB",
  "GZMA", "GZMH", "GZMK"),
  reduction = "umap", cols = c("lightgrey", "red")) +
  plot_annotation(title = "CD8 LP effector cytotoxic program")

save_plot(p, filename = "FeaturePlot_CD8_effector.png", dir = paths$plots_LPL_compartment,
          width = 20, height = 17)

rm(p)

# ---------------------------------------------------------
# CD8 CD160 / ITGA1 / NK-like LP program
# ---------------------------------------------------------

p <- FeaturePlot(lpl_cells,features = c(
  "CD160", "ITGA1", "XCL1",
  "XCL2", "KLRC1", "KLRC2",
  "KLRD1", "KIR2DL4", "TRDC", "ENTPD1"),
                 reduction = "umap",
                 cols = c("lightgrey", "red")) +
  plot_annotation(title = "CD8 LP tissue-associated NK-like program")

save_plot(p, filename = "FeaturePlot_CD8_NK_like.png",
          dir = paths$plots_LPL_compartment,
          width = 20, height = 17)

rm(p)

# ---------------------------------------------------------
# CD4 functional lineages
# ---------------------------------------------------------

p <- FeaturePlot(lpl_cells, features = c(
  # Naive
  "CCR7","SELL","TCF7","LEF1","KLF2",
  # TH17 LP
  "CCR6","RORA","RORC",
  "KLRB1","LTB","CCL20", "IL17A","IL23R","IL26", "LGALS3","GPR65",
  # Treg
  "FOXP3","IL2RA",
  "CTLA4","IL10",
  "BATF", "LAIR2",
  "TNFRSF4","TNFRSF9",
  "GPX1","GLRX",
  # TFH
  "CXCR5", "TOX2", "CXCL13",
  "PDCD1", "BTLA", "CD200", "ICOS"),
  reduction = "umap",
  cols = c("lightgrey", "red")) +
  plot_annotation(title = "CD4 T cell functional programs")

save_plot(p, filename = "FeaturePlot_CD4_lineages.png", dir = paths$plots_LPL_compartment,
          width = 22, height = 25)

rm(p)

# ---------------------------------------------------------
# Activation / TFH-Treg shared state / stress
# ---------------------------------------------------------

p <- FeaturePlot(lpl_cells,features = c(
  "PDCD1", "TIGIT", "ICOS",
  "TOX", "TOX2", "MAF",
  # LP stress cluster
  "HSPA1A", "HSPA1B", "DNAJB1",
  "JUN", "TNF", "IL2"),
  reduction = "umap",
  cols = c("lightgrey", "red")) +
  plot_annotation(title = "Activation, TFH-Treg state and stress program")

save_plot(p, filename = "FeaturePlot_Activation_Stress.png",
          dir = paths$plots_LPL_compartment,
          width = 20, height = 16)

rm(p)

# ---------------------------------------------------------
# Proliferation / cell cycle
# ---------------------------------------------------------

p <- FeaturePlot(lpl_cells, features = c(
  "MKI67", "STMN1", "TUBB", "TUBA1B"), 
  reduction = "umap",
  cols = c("lightgrey", "red")) +
  plot_annotation(title = "Proliferating T cells")

save_plot(p, filename = "FeaturePlot_Proliferation.png",
          dir = paths$plots_LPL_compartment,
          width = 16, height = 12)

rm(p)

# ---------------------------------------------------------
# Contamination / non-T cells
# ---------------------------------------------------------

p <- FeaturePlot(lpl_cells,features = c(
  "MS4A1", "CD79A", "LYZ", "EPCAM"),
  reduction = "umap",
  cols = c("lightgrey", "red")) +
  plot_annotation(title = "Contamination / non-T cells")


save_plot(p, filename = "FeaturePlot_contamination.png", dir = paths$plots_LPL_compartment,
          width = 16, height = 14)

rm(p)


# =========================================================
#                 Manual Annotation
# =========================================================

message("Loading manual cell type annotation table...")

celltype <- read.csv(file.path(paths$tables, "LPL_celltype/LPL_celltype_v2.0.csv"))

message("Mapping cluster IDs to manual cell type annotations...")

new_ids <- celltype$final_annotation

message("Renaming Seurat identities using manual annotations...")

names(new_ids) <- celltype$cluster

lpl_cells <- RenameIdents(lpl_cells, new_ids)

lpl_cells$celltype_LPL <- Idents(lpl_cells)

message("Generating UMAP colored by manual cell type annotations and seurat clusters...")

dimplot_celltype <- DimPlot(lpl_cells, reduction = "umap", label = TRUE, repel = TRUE, label.size = 4) +
  plot_annotation(title = "LPL cell compartment") + 
  xlab("UMAP 1") +
  ylab("UMAP 2")

dimplot_clusters <- DimPlot(lpl_cells, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE, label.size = 5) + 
  plot_annotation(title = "LPL Clusters") +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  NoLegend() 


# -----------------
# Save Final Plot
# -----------------

save_plot(dimplot_celltype, filename = "DimPlot_LPL_celltype.png", dir = paths$plots_LPL_compartment,
          width = 18, height = 16)

save_plot(dimplot_clusters, filename = "DimPlot_LPL_clusters.png", dir = paths$plots_LPL_compartment,
          width = 18, height = 16)


# =========================================================
#                Save outputs
# =========================================================

message("Saving annotated LPL compartment object...")

save_rds(lpl_cells, filename = "srt_LPL_compartment_celltype.rds", dir = paths$objects_subsetting)

message("Saving LPL cluster marker genes table...")

save_csv(lpl_cells_AllMarkers, filename = "AllMarkers_clusters_Cell_LPL.csv", dir = paths$tables_markers)

message("LPL compartment preprocessing and subsetting pipeline completed successfully.")

# =========================================================
#                  Save session info
# =========================================================

message("Saving session information for LPL consensus cell type annotation...")

save_session_info(filename = "sessionInfo_Consensus_CellType_LPL.txt", dir = paths$logs, label = "LPL consensus cell type annotation")

message("Session information saved at: ", paths$logs)


# =========================================================
#                  Final pipeline message
# =========================================================

message("=================================================")
message("PIPELINE STEP COMPLETED: Consensun Celltype Annotations for LPL compartment")
message("=================================================")


