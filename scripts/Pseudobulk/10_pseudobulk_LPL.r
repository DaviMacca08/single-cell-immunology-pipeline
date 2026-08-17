# =========================================================
# Pseudobulk differential expression analysis - LPL
# =========================================================
# Input:
#   - Seurat object of LPL compartment with consensus annotation
#
# Output:
#   - Pseudobulk count matrix
#   - DESeq2 objects and differential expression results
#   - MA plots, volcano plots and heatmaps
#   - Session information log
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

message("Starting pseudobulk analysis for LPL compartment...")


# =========================================================
#                  Function
# =========================================================

# ------------------------
# ComplexHeatmap function
# ------------------------

com_heatmap_table <- function(mat, filename, dir, title){
  
  # Basic input check
  stopifnot(!is.null(mat))
  stopifnot(is.matrix(mat) || is.data.frame(mat))
  
  mat <- as.matrix(mat)
  
  if (nrow(mat) == 0 || ncol(mat) == 0) {
    stop("Matrix is empty.")
  }
  
  # NA and Inf values check
  if (any(is.na(mat))) {
    warning("Matrix contains NA values. Replacing with 0.")
    mat[is.na(mat)] <- 0
  }
  
  if (any(is.infinite(mat))) {
    warning("Matrix contains Inf values. Replacing with 0.")
    mat[is.infinite(mat)] <- 0
  }
  
  # Value range check
  if (min(mat) < 0) {
    warning("Matrix contains negative values. Heatmap scale may be inappropriate for proportions.")
  }
  
  if (max(mat) > 1) {
    warning("Values > 1 detected. Check if matrix is properly normalized (expected proportions?).")
  }
  
  
  open_pdf(filename = filename, dir = dir, width = 16, height = 14)
  
  ht <- Heatmap(
    mat,
    name = "prop",
    cluster_columns = FALSE,
    column_title = title,
    border = TRUE,
    rect_gp = grid::gpar(col = "grey6", lwd = 1.5),
    col = colorRamp2(
      c(0, 0.25, 0.5, 1),
      c("#f7fbff", "#6baed6", "#3182bd", "#08306b")
    ),
    
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid::grid.text(
        sprintf("%.2f", mat[i, j]),
        x, y,
        gp = grid::gpar(fontface = "bold", fontsize = 20)
      )
    }
  )
  
  draw(ht)
  
  close_pdf()
  
  message("Saved plot: Heatmap of ", title)
  
}


# =========================================================
#                   Load input data
# =========================================================

message("Loading Seurat object...")

lpl_cells <- readRDS(
  file.path(paths$objects_subsetting, "srt_LPL_compartment_celltype.rds")
)


# =========================================================
#                  Initial sanity checks
# =========================================================

message(" Running sanity checks...")

stopifnot(inherits(lpl_cells, "Seurat"))
stopifnot("seurat_clusters" %in% colnames(lpl_cells@meta.data))
stopifnot("celltype_LPL" %in% colnames(lpl_cells@meta.data))


# =========================================================
#                Ispection of the dataset
# =========================================================

celltype_plot_lpl <- DimPlot(lpl_cells, reduction = "umap", group.by = "celltype_LPL", label = TRUE, repel = TRUE) + 
  xlab("UMAP 1") + 
  ylab("UMAP 2")

message("Combine celltype in major groups base on clusters...")

Idents(lpl_cells) <- "seurat_clusters"

lpl_cells$celltype_pb <- dplyr::case_when(
  
  lpl_cells$seurat_clusters %in% c(0, 5, 7, 9, 15, 16, 18) ~ "LPL-Cytotoxic-TRM-like" , 
  
  lpl_cells$seurat_clusters == 2 ~ "LPL-TH17",
  
  lpl_cells$seurat_clusters == 8  ~ "LPL-TFH",
  
  lpl_cells$seurat_clusters %in% c(6,11) ~ "LPL-Treg",
  
  lpl_cells$seurat_clusters == 17 ~ "LPL-Cycling-T-cells",
  
  lpl_cells$seurat_clusters %in% c(1, 3, 4, 10, 12, 13, 14, 19) ~ "LPL-Other-CD4-T-cells",
  
  lpl_cells$seurat_clusters == 20 ~ "Epithelial-contamination",
  
  TRUE ~ NA_character_  
  
)

stopifnot(!any(is.na(lpl_cells$celltype_pb)))

# ---------------------
# DimPlot Major Groups
# ---------------------

celltype_plot_pb <- DimPlot(lpl_cells, reduction = "umap", group.by = "celltype_pb", label = TRUE, repel = TRUE) +
  xlab("UMAP 1") + 
  ylab("UMAP 2") 

save_plot(celltype_plot_lpl + celltype_plot_pb, filename = "DimPlot_LPL_MajorGroup.png", dir = paths$plots_pb_LPL,
          width = 20, height = 10)


# =========================================================
#                Design Diagnostics Plots
# =========================================================

# ---------------------
# Table celltype_pb x disease
# ---------------------

tab_CellxDisease <- (prop.table(table(lpl_cells$celltype_pb, lpl_cells$disease), margin = 2))

com_heatmap_table(mat = tab_CellxDisease, filename = "HeatMap_LPL_CellxDisease.pdf", 
                  title = "Cell type composition per disease (proportions) - LPL Compartment", dir = paths$plots_pb_LPL)

# ------------------------------
# Table celltype_pb x condition
# ------------------------------

tab_CellxCondition <- (prop.table(table(lpl_cells$celltype_pb, lpl_cells$condition), margin = 2))

com_heatmap_table(mat = tab_CellxCondition, filename = "HeatMap_LPL_CellxCondition.pdf", 
                  title = "Cell type composition per condition (proportions) - LPL Compartment", dir = paths$plots_pb_LPL)

# ------------------------------
# Table celltype_pb x donor
# ------------------------------

tab_CellxDonor <- prop.table(table(lpl_cells$celltype_pb, lpl_cells$donor), margin = 2)

com_heatmap_table(mat = tab_CellxDonor, filename = "HeatMap_LPL_CellxDonor.pdf", 
                  title = "Cell type composition per donor (proportions) - LPL Compartment", dir = paths$plots_pb_LPL)

# ------------------------------
# Table condition x donor
# ------------------------------

tab_ConditionxDonor <- prop.table(table(lpl_cells$condition, lpl_cells$donor), margin = 2)

com_heatmap_table(mat = tab_ConditionxDonor, filename = "HeatMap_LPL_ConditionxDonor.pdf", 
                  title = "Cell type condition per donor (proportions) - LPL Compartment", dir = paths$plots_pb_LPL)

# ------------------------------
# Table condition x disease
# ------------------------------

tab_ConditionxDisease <- (prop.table(table(lpl_cells$condition, lpl_cells$disease), margin = 2))

com_heatmap_table(mat = tab_ConditionxDisease, filename = "HeatMap_LPL_ConditionxDisease.pdf", 
                  title = "Cell type condition per disease (proportions) - LPL Compartment", dir = paths$plots_pb_LPL)


# =============================================================
#   QC: minimum cell/donor requirement per celltype_pb x disease
# =============================================================

min_cells_per_sample     <- 10  # Minimum cells required for a reliable donor x cell type pseudobulk sample
min_donors_per_disease   <- 2   # Minimum biological replicates per disease group for reliable DESeq2 dispersion estimation


cell_counts <- as.data.frame(
  table(
    celltype_pb = lpl_cells$celltype_pb,
    donor       = lpl_cells$donor,
    disease     = lpl_cells$disease
  )
)


cell_counts_present <- cell_counts[cell_counts$Freq > 0, ]


# Identify pseudobulk samples with insufficient cell numbers

low_cell_samples <- cell_counts_present[
  cell_counts_present$Freq < min_cells_per_sample,
]


if (nrow(low_cell_samples) > 0) {
  
  message(
    "The following pseudobulk samples (celltype_pb x donor x disease) contain fewer than ",
    min_cells_per_sample,
    " cells:"
  )
  
  print(low_cell_samples)
  
}


# Count donors passing the minimum cell threshold for each cell type and disease group

donor_summary <- cell_counts_present[
  cell_counts_present$Freq >= min_cells_per_sample,
] |> 
  dplyr::group_by(celltype_pb, disease) |> 
  dplyr::summarise(
    n_donors = dplyr::n_distinct(donor),
    .groups = "drop"
  )


print(donor_summary)



# Cell types eligible for differential expression analysis:
# only cell types with sufficient donor representation in all disease groups

eligible_celltypes <- donor_summary |> 
  dplyr::group_by(celltype_pb) |> 
  dplyr::summarise(
    n_groups = dplyr::n_distinct(disease),
    min_donors = min(n_donors),
    .groups = "drop"
  ) |> 
  dplyr::filter(
    n_groups == dplyr::n_distinct(lpl_cells$disease),
    min_donors >= min_donors_per_disease
  ) |> 
  dplyr::pull(celltype_pb) |> 
  as.character()


message(
  "Eligible cell types for DESeq2 analysis (>= ",
  min_donors_per_disease,
  " donors/group and >= ",
  min_cells_per_sample,
  " cells/sample): ",
  paste(eligible_celltypes, collapse = ", ")
)


# =============================================================
#   Prepare and run DESeq2 for DE across cell types clusters
# =============================================================

message("Mandatory condition: ")
message("Almost 2 donor per celltype_pb")
message("Balance per donor and disease")

celltype_DE <- c("LPL-Cytotoxic-TRM-like", "LPL-Treg")

# check
if (!all(celltype_DE %in% eligible_celltypes)) {
  
  stop(
    "The following cell types in celltype_DE do not meet the minimum cell/donor requirements: ",
    paste(setdiff(celltype_DE, eligible_celltypes), collapse = ", "),
    ". Review celltype_DE selection or pseudobulk aggregation."
  )
  
} else {
  
  message("QC passed: all selected cell types meet the minimum cell/donor requirements for DESeq2 analysis.")
  
}

# subset
lpl_cells_pb <- subset(
  lpl_cells,
  subset =  celltype_pb %in% celltype_DE
)

message("Construct correct id for aggregate counts...")

lpl_cells_pb$pseudobulk_id <- paste(
  lpl_cells_pb$celltype_pb,
  lpl_cells_pb$donor,
  lpl_cells_pb$disease,
  sep = "-"
)

message("Aggregate sample...")

DefaultAssay(lpl_cells_pb) <- "RNA"

cts <- AggregateExpression(
  lpl_cells_pb,
  group.by = "pseudobulk_id",
  assays = "RNA",
  slot = "counts",
  return.seurat = FALSE
)

message("Get counts data...")

cts <- cts$RNA

cts_db <- as.data.frame(cts)

message("Prepare correct colData for DESeq2...")

colData <- data.frame(samples = colnames(cts_db))

meta_sample_id <- strsplit(colData$samples, "-")

colData$celltype <- factor(sapply(meta_sample_id, function(x) paste(x[1:(length(x)-2)], collapse = "-")))
colData$donor <- factor(str_extract(colData$samples, "(?<=-)[0-9]+(?=-)"))
colData$condition <- factor(str_extract(colData$samples, "(Control|CD)$"))

rownames(colData) <- colData$samples

stopifnot(all(colData$samples == colnames(cts_db)))

# ---------------
# Run DESeq2 
# ---------------

message(
  "Cell types selected for differential expression analysis: ",
  paste(celltype_DE, collapse = ", ")
)

dds_list <- list()
deg_list <- list()
dds_obj_list <- list() 
res_obj_list <- list()

table(colData$celltype, colData$condition)
table(colData$celltype, colData$donor)

for (ct in unique(colData$celltype)) {
  
  message("Processing: ", ct)
  
  keep <- colData$celltype == ct
  sub_meta <- colData[keep, ]
  
  print(table(sub_meta$condition))
  print(table(sub_meta$donor))
  
  dds <- DESeqDataSetFromMatrix(
    countData = cts[, keep],
    colData = sub_meta,
    design = ~ condition
  )
  
  keep_genes <- rowSums(counts(dds) >= 10) >= 2
  dds <- dds[keep_genes, ]
  
  print("Design OK")
  
  dds$condition <- relevel(dds$condition, ref = "Control")
  
  dds <- DESeq(dds)
  
  print("DESeq OK")
  
  res <- results(dds, contrast=c("condition","CD","Control"))
  
  res_shrunk <- lfcShrink(dds, coef = "condition_CD_vs_Control", type = "apeglm")
  
  dds_list[[ct]] <- as.data.frame(res)
  
  deg_list[[ct]] <- as.data.frame(res_shrunk)
  
  dds_obj_list[[ct]]  <- dds
  
  res_obj_list[[ct]]  <- res_shrunk 
}


# =========================================================
#                DEG Plots
# =========================================================

message("Generating DEG visualizations (volcano plots, MA plots, and heatmaps)...")

padj_cutoff <- 0.05
lfc_cutoff  <- 1
top_degs_n  <- 30   # Number of genes displayed in heatmaps
top_labels_n <- 15  # Number of genes labelled in volcano plots

condition_palette <- c("Control" = "#0072B2","CD" = "#D55E00")

# Qualitative palette for donors (variable number of samples)
all_donors <- levels(colData$donor)
donor_palette <- setNames(
  RColorBrewer::brewer.pal(max(3, length(all_donors)), "Dark2")[seq_along(all_donors)],
  all_donors
)


for (ct in celltype_DE) {
  
  message("Generating DEG plots for: ", ct)
  
  ct_display <- gsub("-", " ", ct) # Improve cell type labels in plots
  
  stopifnot(ct %in% names(deg_list))
  stopifnot(ct %in% names(dds_obj_list))
  
  res_df <- deg_list[[ct]]
  res_df$gene <- rownames(res_df)
  res_df <- res_df[!is.na(res_df$padj), ]
  
  # ------------------------------------------
  # Volcano Plot
  # ------------------------------------------
  
  res_df$significance <- "Not Significant"
  res_df$significance[res_df$padj < padj_cutoff & res_df$log2FoldChange >  lfc_cutoff] <- "Up in CD"
  res_df$significance[res_df$padj < padj_cutoff & res_df$log2FoldChange < -lfc_cutoff] <- "Down in CD"
  
  n_up   <- sum(res_df$significance == "Up in CD")
  n_down <- sum(res_df$significance == "Down in CD")
  
  
  message(sprintf("[%s] Significant DEGs: %d upregulated in CD, %d downregulated in CD",ct, n_up, n_down))  
  
  genes_to_label <- res_df[res_df$significance != "Not Significant", ]
  genes_to_label <- head(genes_to_label[order(genes_to_label$padj), ], top_labels_n)
  
  volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = significance)) +
    geom_point(size = 1.5, alpha = 0.7) +
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", color = "grey40") +
    scale_color_manual(
      values = c(
        "Up in CD"        = "#CE2915",
        "Down in CD"      = "#0096FF",
        "Not Significant" = "grey75"
      )
    ) +
    ggrepel::geom_text_repel(
      data = genes_to_label,
      aes(label = gene),
      size = 3, max.overlaps = 20, color = "black",
      segment.size = 0.2
    ) +
    labs(
      title = paste0("Volcano plot: CD vs Control (", ct_display, ")"),
      x = expression(log[2] ~ fold~change),
      y = expression(-log[10] ~ adjusted~italic(P)),
      color = NULL
    ) +
    theme_bw()
  
  save_plot(
    volcano_plot,
    filename = paste0("Volcano_", gsub("-", "", ct), "_CDvsControl.png"),
    dir = paths$plots_pb_LPL,
    width = 8, height = 6
  )
  
  # ------------------------------------------
  # MA plot (LFC shrinkage diagnostic)
  # ------------------------------------------
  
  open_png(
    filename = paste0("MAplot_", gsub("-", "", ct), "_CDvsControl.png"),
    dir = paths$plots_pb_LPL, width = 1000, height = 800
  )
  plotMA(res_obj_list[[ct]], ylim = c(-7, 7), alpha = padj_cutoff, cex = 1.5,
         main = paste0("MA plot: CD vs Control (", ct_display, ") — shrunk LFC"))
  abline(h = c(-lfc_cutoff, lfc_cutoff), col = "grey70", lty = 2)
  close_png()
  
  # ------------------------------------------
  # Heatmap of top DEGs (VST Z-score)
  # ------------------------------------------
  
  sig_deg <- subset(res_df, padj < padj_cutoff & abs(log2FoldChange) > lfc_cutoff)
  sig_deg <- sig_deg[order(sig_deg$padj), ]
  
  # Rank genes by adjusted p-value and select top DEGs for visualization
  
  ranked <- res_df[!is.na(res_df$padj), ]
  ranked <- ranked[order(ranked$padj), ]
  
  top_genes <- head(
    ranked$gene,
    min(top_degs_n, nrow(ranked))
  )
  
  heatmap_note <- ifelse(
    sum(top_genes %in% sig_deg$gene) < length(top_genes),
    " — including non-significant top-ranked genes",
    ""
  )
  
  n_show <- length(top_genes)
  
  dds_ct <- dds_obj_list[[ct]]
  vsdata_ct <- tryCatch(
    vst(dds_ct, blind = FALSE),
    error = function(e) varianceStabilizingTransformation(dds_ct, blind = FALSE)
  )
  vsdata_mat <- assay(vsdata_ct)
  
  heatmap_matrix <- vsdata_mat[top_genes, , drop = FALSE]
  heatmap_matrix_scaled <- t(scale(t(heatmap_matrix)))
  
  # ------------------------------------------
  # Row annotation: DEG significance status
  # ------------------------------------------
  
  gene_status <- factor(
    ifelse(top_genes %in% sig_deg$gene, "Significant", "Top ranked"),
    levels = c("Significant", "Top ranked")
  )
  
  status_palette <- c("Significant" = "#D55E00","Top ranked" = "#D9D9D9"
  )  
  row_ann <- rowAnnotation(
    DEG = gene_status,
    col = list(DEG = status_palette),
    annotation_legend_param = list(DEG = list(title = "DEG status"))
  )
  
  sub_meta_ct <- as.data.frame(colData(dds_ct))
  
  ann_col <- HeatmapAnnotation(
    Condition = sub_meta_ct$condition,
    Donor     = sub_meta_ct$donor,
    col = list(
      Condition = condition_palette,
      Donor     = donor_palette
    )
  )
  
  open_png(
    filename = paste0("Heatmap_TopDEG_", gsub("-", "", ct), "_CDvsControl.png"),
    dir = paths$plots_pb_LPL, width = 1400, height = 1600
  )
  
  ht_deg <- ComplexHeatmap::Heatmap(
    heatmap_matrix_scaled,
    name = "Z-score",
    col = circlize::colorRamp2(c(-2, 0, 2),c("#5E3C99", "white", "#E66101")),
    top_annotation = ann_col,
    left_annotation = row_ann,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_row_names = TRUE,
    show_column_names = FALSE,
    row_names_gp = grid::gpar(fontsize = 7),
    column_title = paste0("Top ", n_show, " DEGs — ", ct_display, " (CD vs Control)", heatmap_note),
    heatmap_legend_param = list(title = "Row\nZ-score")
  )
  
  draw(ht_deg)
  close_png()
  
  message("Saved DEG plots for: ", ct)
  
}


# =========================================================
#                Save outputs
# =========================================================

message("Saving pseudobulk analysis outputs...")

save_rds(dds_list, filename = "results_list_LPL.rds", dir = paths$objects_pseudobulk)

save_rds(deg_list, filename = "deg_list_shrunk_LPL.rds", dir = paths$objects_pseudobulk)

save_csv(cts_db, filename = "aggregated_counts_LPL.csv", dir = paths$objects_pseudobulk)

message("Pseudobulk outputs for LPL compartment saved in: ", paths$objects_pseudobulk)


# =========================================================
#                  Save session info
# =========================================================

message("Saving session information (pseudobulk DESeq2 analysis for LPL compartment)...")

save_session_info(filename = "sessionInfo_pseudobulk_DESeq2_LPL.txt", dir = paths$logs, label = "Pseudobulk DESeq2 analysis for LPL compartment")

message("Session information saved at: ", paths$logs)


# =========================================================
#                  Final pipeline message
# =========================================================

message("=================================================")
message("PIPELINE STEP COMPLETED: Pseudobulk differential expression for LPL compartment")
message("=================================================")