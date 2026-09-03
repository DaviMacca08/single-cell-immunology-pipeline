# =========================================================
# Project      : scRNA-Seq project - Bulk RNA-seq Analysis Pipeline
# Dataset      : GSE193677 (GPL16791 Illumina HiSeq 2500)
# Samples      : Ileum biopsies, CD and Control
# Script       : Quality Control, Exploratory Data Analysis and Differential Expression Analysis (Sensitivity analysis)
# Description  : Load raw counts and GEO metadata, validate
#                sample annotation and patient-level replication,
#                define the primary analysis cohort, and perform
#                limma-voom differential expression analysis.
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

message("=== Starting GSE193677 RNA-seq sensitivity analysis ===")

# =========================================================
#                  Load data 
# =========================================================

message("[LOAD] Loading raw count matrix...")

counts_raw <- read.delim(gzfile("scripts/Bulk_RNAseq/Raw_data/GSE193677_MSCCR_Biopsy_counts.txt.gz"), 
                         row.names = 1,
                         check.names = FALSE,
                         header = TRUE,
                         sep = "")


message(
  "[LOAD] Count matrix dimensions: ",
  nrow(counts_raw), " genes x ",
  ncol(counts_raw), " samples."
)


# =========================================================
#                  Load GEO metadata
# =========================================================

gse_id <- "GSE193677"

message("[LOAD] Downloading GEO metadata for ", gse_id, "...")

gse_list <- getGEO(gse_id, GSEMatrix = TRUE, AnnotGPL = FALSE)

if (is.list(gse_list) && !is(gse_list, "ExpressionSet")) {
  message("getGEO() returned ", length(gse_list), " ExpressionSet object(s).")
  gse <- gse_list[[1]]
} else {
  gse <- gse_list
}

meta <- pData(gse)
rownames(meta) <- meta$geo_accession

message(
  "[LOAD] GEO metadata dimensions: ",
  nrow(meta), " samples x ",
  ncol(meta), " metadata variables."
)


# =========================================================
#                  Parse sample metadata
# =========================================================

message("[META] Parsing sample-level metadata...")

# Patient identifier
meta$patient_id <- sub("^MSCCR_reGRID_([0-9]+)_Biopsy_.*", "\\1", meta$title)

# Sample identifier
meta$sample <- sub(",.*$", "", meta$title)

# Anatomical region
meta$region  <- sub("^regionre: ",    "", meta$characteristics_ch1.2)

# Disease status
meta$disease <- sub("^ibd_disease: ", "", meta$characteristics_ch1.4)

# Inflammation status
meta$infl    <- sub("^typere: ",      "", meta$characteristics_ch1.5)


# =========================================================
#                  Metadata validation
# =========================================================

message("[CHECK] Validating parsed metadata...")

# Check for parsing failures
if (any(
  is.na(meta$patient_id) |
  meta$patient_id == "" |
  is.na(meta$region) |
  meta$region == "" |
  is.na(meta$disease) |
  meta$disease == "" |
  is.na(meta$infl) |
  meta$infl == ""
)) {
  stop(
    "[CHECK] Missing or empty values detected in parsed metadata."
  )
}

# Check expected disease categories
message("[CHECK] Disease categories:")

print(table(meta$disease, useNA = "ifany"))

# Check expected anatomical regions
message("[CHECK] Anatomical regions:")

print(table(meta$region, useNA = "ifany"))

# Check inflammation categories
message("[CHECK] Inflammation categories:")

print(table(meta$infl, useNA = "ifany"))


# =========================================================
#                  Count matrix / metadata matching
# =========================================================

message("[CHECK] Checking count matrix and GEO metadata matching...")

count_samples <- colnames(counts_raw)
meta_samples <- meta$sample

# Check duplicated sample identifiers
if (anyDuplicated(count_samples) > 0) {
  stop("[CHECK] Duplicated sample identifiers found in count matrix.")
}

if (anyDuplicated(meta_samples) > 0) {
  stop("[CHECK] Duplicated sample identifiers found in GEO metadata.")
}

# Samples present in counts but missing from metadata
counts_not_in_meta <- setdiff(count_samples, meta_samples)

# Samples present in metadata but missing from counts
meta_not_in_counts <- setdiff(meta_samples, count_samples)

if (length(counts_not_in_meta) > 0) {
  
  stop("[CHECK] ", length(counts_not_in_meta), " count-matrix samples are missing from GEO metadata.")
  
}

if (length(meta_not_in_counts) > 0) {
  
  warning(
    "[CHECK] ",
    length(meta_not_in_counts),
    " GEO metadata samples are not present in the count matrix."
  )
  
}

# Restrict metadata to samples available in the count matrix
meta <- meta[match(count_samples, meta$sample), , drop = FALSE]

# Final ordering check
if (!identical(colnames(counts_raw), meta$sample)) {
  stop(
    "[CHECK] Count matrix columns and metadata rows are not aligned."
  )
}

message(
  "[CHECK] Count matrix and metadata are correctly matched and ordered."
)


# =========================================================
#                  Define Ileum CD/Control cohort
# =========================================================

message("[SUBSET] Selecting ileal CD and Control samples...")

sub_ileum <- meta[meta$region == "Ileum" & meta$disease %in% c("CD", "Control"), , drop = FALSE]

message("[SUBSET] Ileum CD/Control cohort: ", nrow(sub_ileum), " biopsies.")

message("[CHECK] Ileum disease distribution:")

print(table(sub_ileum$disease))

message("[CHECK] Ileum disease x inflammation distribution:")

print(with(sub_ileum, table(disease, infl)))


# =========================================================
#                  Patient-level replication check
# =========================================================

message("[CHECK] Assessing number of ileal biopsies per patient...")

biopsies_per_patient <- table(sub_ileum$patient_id)

print(table(biopsies_per_patient))

dup_patients <- names(biopsies_per_patient[biopsies_per_patient > 1])

message("[CHECK] Patients contributing >1 ileal biopsy: ", length(dup_patients))


# =========================================================
#                  Inspect repeated patients
# =========================================================

if (length(dup_patients) > 0) {
  
  message("[CHECK] Inspecting repeated ileal biopsies from the same patients...")
  
  repeated_ileal <- sub_ileum[
    sub_ileum$patient_id %in% dup_patients,
    c("geo_accession", "patient_id", "region", "disease", "infl", "title"),
    drop = FALSE
    ]
  
  print(repeated_ileal)
  
  message("[CHECK] Disease x inflammation status for repeated patients:")
  
  print(with(repeated_ileal, table(patient_id, disease, infl)))
  
} else {
  
  message("[CHECK] No patients with repeated ileal biopsies.")
  
}


# =========================================================
#            Define sensitivity analysis subset
# =========================================================

message("[SUBSET] Defining sensitivity cohort: ", "all ileal CD/Control biopsies (I + NonI).")

sub_ileum_all <- sub_ileum

message("[SUBSET] Sensitivity cohort: ", nrow(sub_ileum_all)," biopsies.")

message("[CHECK] Sensitivity cohort disease distribution:")

print(table(sub_ileum_all$disease))

message("[CHECK] Sensitivity cohort disease x inflammation distribution:")

print(with(sub_ileum_all, table(disease, infl)))


# =========================================================
#                  Sensitivity cohort replication check
# =========================================================

message("[CHECK] Assessing repeated patients in sensitivity cohort...")

sensitivity_biopsies_per_patient <- table(sub_ileum_all$patient_id)

message(
  "[CHECK] Unique patients in sensitivity cohort: ",
  length(unique(sub_ileum_all$patient_id))
)

message(
  "[CHECK] Patients contributing >1 ileal biopsy: ",
  sum(sensitivity_biopsies_per_patient > 1)
)


# =========================================================
#  Remove pseudoreplicated patients (sensitivity cohort)
# =========================================================

dup_patients_sensitivity <- names(sensitivity_biopsies_per_patient[sensitivity_biopsies_per_patient > 1])

sub_ileum_all_clean <- sub_ileum_all[!sub_ileum_all$patient_id %in% dup_patients_sensitivity, , drop = FALSE]

message("[CHECK] Sensitivity cohort after removing pseudoreplicated patients: ", nrow(sub_ileum_all_clean), " biopsies.")

if (any(table(sub_ileum_all_clean$patient_id) > 1)) {
  stop("[CHECK] Pseudoreplication still present in sensitivity cohort after filtering.")
}

message("[CHECK] Sensitivity cohort disease distribution (clean):")
print(table(sub_ileum_all_clean$disease))


# =========================================================
#                  Sensitivity cohort summary
# =========================================================

message("[SUMMARY] Sensitivity analysis cohort:")

sensitivity_summary <- data.frame(
  Cohort = "Ileum_I_NonI_CD_vs_Control",
  CD = sum(sub_ileum_all_clean$disease == "CD"),
  Control = sum(sub_ileum_all_clean$disease == "Control"),
  Total = nrow(sub_ileum_all_clean),
  Unique_Patients = length(
    unique(sub_ileum_all_clean$patient_id)
  )
)

print(sensitivity_summary)


# =========================================================
#             Define sensitivity count matrices
# =========================================================

message("[COUNTS] Creating sensitivity-specific count matrices...")

counts_sensitivity <- counts_raw[, sub_ileum_all_clean$sample, drop = FALSE]


# =========================================================
#                  Sensitivity cohort alignment
# =========================================================

message("[CHECK] Validating sensitivity count matrix alignment...")

meta_sensitivity <- sub_ileum_all_clean[match(colnames(counts_sensitivity), sub_ileum_all_clean$sample), , drop = FALSE]

if (!identical(colnames(counts_sensitivity), meta_sensitivity$sample)) {
  stop(
    "[CHECK] Sensitivity count matrix and metadata are not aligned."
  )
}

message(
  "[CHECK] Sensitivity counts and metadata are correctly aligned."
)


# =========================================================
#             Count matrix QC (Before limma-voom)
# =========================================================

message("[QC] Checking sensitivity count matrix...")

stopifnot(ncol(counts_sensitivity) == nrow(meta_sensitivity))

stopifnot(nrow(counts_sensitivity) > 1000)

stopifnot(!any(duplicated(rownames(counts_sensitivity))))

stopifnot(!any(duplicated(colnames(counts_sensitivity))))

stopifnot(!anyNA(as.matrix(counts_sensitivity)))

stopifnot(all(is.finite(as.matrix(counts_sensitivity))))

stopifnot(all(counts_sensitivity >= 0))

# Verify sample order

stopifnot(all(colnames(counts_sensitivity) == meta_sensitivity$sample))

message(
  "[QC] Sensitivity count matrix: ",
  nrow(counts_sensitivity),
  " genes x ",
  ncol(counts_sensitivity),
  " samples."
)


# =========================================================
#             limma-voom: Sample metadata
# =========================================================

message("[LIMMA] Building sample metadata...")

colData <- data.frame(
  row.names = meta_sensitivity$sample,
  condition = factor(meta_sensitivity$disease, levels = c("Control", "CD"))
)

message("[META] Sample distribution:")

print(table(colData$condition))

stopifnot(all(colnames(counts_sensitivity) == rownames(colData)))

stopifnot(length(unique(colData$condition)) >= 2)


# =========================================================
#             edgeR DGEList
# =========================================================

message("[LIMMA] Creating edgeR DGEList...")

dge <- edgeR::DGEList(counts = counts_sensitivity)

message(
  "[LIMMA] DGEList created: ",
  nrow(dge),
  " genes x ",
  ncol(dge),
  " samples."
)


# =========================================================
#             Filter lowly expressed genes
# =========================================================

# Lowly expressed genes are removed before normalization and differential expression analysis using edgeR::filterByExpr()

message("[LIMMA] Filtering lowly expressed genes...")

keep_genes <- edgeR::filterByExpr(dge, group = colData$condition)

message(
  "[LIMMA] Genes retained: ",
  sum(keep_genes),
  " / ",
  nrow(dge),
  " (",
  round(100 * mean(keep_genes), 1),
  "%)."
)

dge <- dge[keep_genes, , keep.lib.sizes = FALSE]

stopifnot(nrow(dge) > 1000)


# =========================================================
#             TMM normalization
# =========================================================

# TMM normalization accounts for compositional differences in library sizes across samples

message("[LIMMA] Calculating TMM normalization factors...")

dge <- edgeR::calcNormFactors(dge, method = "TMM")

message("[LIMMA] Normalization factors summary:")

print(summary(dge$samples$norm.factors))


# =========================================================
#             Design matrix
# =========================================================

message("[LIMMA] Building design matrix...")

design <- model.matrix(
  ~ 0 + condition,
  data = colData
)

colnames(design) <- levels(colData$condition)

stopifnot(all(colnames(design) == c("Control", "CD")))


# =========================================================
#             voomLmFit transformation
# =========================================================

message("[LIMMA] Applying voomLmFit with empirical sample quality weights...")

open_png(filename = "voom_mean_variance_trend.png", dir = paths$plots_bulk_rnaseq_sensitivity_qc,
         width = 1000, height = 800)

fit <- edgeR::voomLmFit(
  counts = dge,
  design = design, 
  sample.weights = TRUE,
  plot = TRUE,
  keep.EList = TRUE)

close_png()

message("[LIMMA] voomLmFit completed (sample quality weights estimated).")

# Store voom-transformed log2-CPM values for downstream QC and EDA
logcpm <- as.matrix(fit$EList$E)

# Empirical sample quality weights estimated by voomLmFit
# These are continuous weights used by the model and are not equivalent to binary outlier removal
sample_weights_df <- data.frame(
  sample = rownames(fit$targets),
  sample_weight = fit$targets$sample.weight,
  group = colData$condition
)

message("[QC] Sample quality weights summary:")
print(summary(sample_weights_df$sample_weight))


# =========================================================
#                  Quality Control 
# =========================================================

group_palette <- setNames(RColorBrewer::brewer.pal(3, "Set1")[1:2], levels(colData$condition))
sample_colors <- group_palette[as.character(colData$condition)]

# Distribution of empirical sample quality weights

sample_weights_df$sample <- factor(sample_weights_df$sample,
                                   levels = sample_weights_df$sample[order(sample_weights_df$sample_weight)])

weights_plot <- ggplot(sample_weights_df, aes(x = sample_weight, fill = group)) +
  geom_histogram(bins = 40, alpha = 0.8, position = "identity") +
  scale_fill_manual(values = group_palette) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey30") +
  labs(title = "Distribution of empirical sample quality weights (GSE193677, Sensitivity)",
       x = "Sample weight (voomLmFit)", y = "Count") +
  theme_bw()

save_plot(weights_plot, filename = "sample_quality_weights_distribution.png",
          dir = paths$plots_bulk_rnaseq_sensitivity_qc, width = 8, height = 6)

# Distribution of voom log2-CPM values

open_png(filename = "boxplot_logCPM_distributions.png", dir = paths$plots_bulk_rnaseq_sensitivity_qc,
         width = 1600, height = 1300)

par(mar = c(8, 4, 4, 2))
boxplot(logcpm,
        main = "Distribution of voom log2-CPM per sample (GSE193677, Sensitivity)",
        ylab = "log2-CPM",
        col = sample_colors,
        las = 2, cex.axis = 0.5, outline = FALSE)
legend("topright", legend = levels(colData$condition), fill = group_palette, bty = "n")

close_png()

open_png(filename = "density_plot_logCPM.png", dir = paths$plots_bulk_rnaseq_sensitivity_qc, 
         width = 1200, height = 800)

plot(density(logcpm[, 1]), col = sample_colors[1], lwd = 1,
     main = "Density of voom log2-CPM across samples",
     xlab = "log2-CPM", ylim = c(0, 0.5))
for (i in 2:ncol(logcpm)) {
  lines(density(logcpm[, i]), col = sample_colors[i], lwd = 0.8)
}
legend("topright", legend = levels(colData$condition), col = group_palette, lwd = 2, bty = "n")

close_png()

# Library size per sample (post-filtering, TMM-normalized)

lib_df <- data.frame(
  sample = rownames(dge$samples),
  lib_size = dge$samples$lib.size * dge$samples$norm.factors,
  group = colData$condition
)

libsize_plot <- ggplot(lib_df, aes(x = reorder(sample, lib_size), y = lib_size, fill = group)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = group_palette) +
  labs(title = "Effective library size per sample (GSE193677, Sensitivity)", x = "", y = "TMM-normalized library size") +
  theme_bw() + 
  theme(axis.title.y = element_text(size = 0.5))

save_plot(libsize_plot, filename = "library_size.png", dir = paths$plots_bulk_rnaseq_sensitivity_qc, width = 11, height = 15)


# =========================================================
#               Outlier Detection (QC)
# =========================================================

sample_medians <- matrixStats::colMedians(logcpm, na.rm = TRUE)
med_center <- median(sample_medians, na.rm = TRUE)
med_iqr <- IQR(sample_medians, na.rm = TRUE)
outliers_median <- colnames(logcpm)[abs(sample_medians - med_center) > 3 * med_iqr]

message(sprintf("[QC] Median-logCPM rule flagged %d sample(s).", length(outliers_median)))

# Compute within-group correlations separately for each condition
# to avoid interpreting biological differences between CD and Control
# as technical sample-quality differences.

outliers_corr <- character(0)
corr_diag_df <- data.frame()

for (grp in levels(colData$condition)) {
  
  grp_samples <- rownames(colData)[colData$condition == grp]
  cor_mat_grp <- cor(logcpm[, grp_samples], method = "pearson", use = "pairwise.complete.obs")
  diag(cor_mat_grp) <- NA
  avg_cor_grp <- rowMeans(cor_mat_grp, na.rm = TRUE)
  
  thr_grp <- mean(avg_cor_grp, na.rm = TRUE) - 3 * sd(avg_cor_grp, na.rm = TRUE)
  flagged_grp <- names(avg_cor_grp)[avg_cor_grp < thr_grp]
  
  outliers_corr <- c(outliers_corr, flagged_grp)
  
  corr_diag_df <- rbind(corr_diag_df, data.frame(
    sample = names(avg_cor_grp),
    group = grp,
    avg_within_group_cor = avg_cor_grp
  ))
  
  message(sprintf("[QC] Group %s: correlation rule flagged %d/%d sample(s).",
                  grp, length(flagged_grp), length(grp_samples)))
}

outlier_samples <- union(outliers_median, outliers_corr)

message(sprintf("[QC] Total unique outliers (descriptive only): %d sample(s).", length(outlier_samples)))

qc_outliers_df <- data.frame(
  sample = colnames(logcpm),
  group = colData$condition,
  reason_median_logCPM = colnames(logcpm) %in% outliers_median,
  reason_low_within_group_correlation =
    colnames(logcpm) %in% outliers_corr,
  sample_weight =
    fit$targets$sample.weight,
  stringsAsFactors = FALSE
)

qc_outliers_df <- merge(qc_outliers_df, corr_diag_df[, c("sample", "avg_within_group_cor")], by = "sample")

message(
  "[QC] Outlier assessment completed. Flagged samples were retained; ",
  "sample quality weights were estimated by voomLmFit and incorporated ",
  "into the model."
)

# Median-logCPM rule: visual check

open_png(filename = "outlier_median_logCPM.png", dir = paths$plots_bulk_rnaseq_sensitivity_qc,
         width = 1000, height = 700)

plot(sample_medians, col = ifelse(colnames(logcpm) %in% outliers_median, "red", "grey40"),
     pch = 16, main = "Sample median log2-CPM (outlier rule: median ± 3xIQR)",
     xlab = "Sample index", ylab = "Median log2-CPM")
abline(h = med_center + c(-3, 3) * med_iqr, col = "red", lty = 2)

close_png()

# Within-group correlation rule: visual check

corr_check_plot <- ggplot(corr_diag_df, aes(x = group, y = avg_within_group_cor, color = group)) +
  geom_jitter(width = 0.15, alpha = 0.6) +
  scale_color_manual(values = group_palette) +
  labs(title = "Average within-group sample correlation (GSE193677, Sensitivity)",
       x = "", y = "Average within-group Pearson correlation") +
  theme_bw()

save_plot(corr_check_plot, filename = "outlier_within_group_correlation.png",
          dir = paths$plots_bulk_rnaseq_sensitivity_qc, width = 7, height = 6)


# =========================================================
#           Exploratory Data Analysis (EDA) 
# =========================================================

# PCA

pca_res <- prcomp(t(logcpm), center = TRUE, scale. = FALSE)
pca_var <- pca_res$sdev^2 / sum(pca_res$sdev^2)
pca_var_percent <- round(100 * pca_var, 1)

pca_df <- data.frame(
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2],
  group = colData$condition,
  outlier = colnames(logcpm) %in% outlier_samples
)

pca_plot_group <- ggplot(pca_df, aes(PC1, PC2, color = group, shape = outlier)) +
  geom_point(size = 2.2, alpha = 0.85) +
  xlab(paste0("PC1 (", pca_var_percent[1], "%)")) +
  ylab(paste0("PC2 (", pca_var_percent[2], "%)")) +
  scale_color_manual(values = group_palette) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 4)) +
  ggtitle("PCA of GSE193677 samples (Sensitivity: Ileum I + NonI, CD vs Control)") +
  theme_bw()

save_plot(pca_plot_group, filename = "pca_plot.png", dir = paths$plots_bulk_rnaseq_sensitivity_eda)

# Scree plot

var_df <- data.frame(
  PC = factor(paste0("PC", 1:10), levels = paste0("PC", 1:10)),
  variance = pca_var[1:10] * 100,
  cumvar = cumsum(pca_var[1:10]) * 100
)

open_png(filename = "pca_scree_plot.png", dir = paths$plots_bulk_rnaseq_sensitivity_eda,
         width = 1000, height = 700)

par(mar = c(5, 5, 4, 2))

# Bar plot: percentage of variance explained by each principal component
bar_centers <- barplot(var_df$variance, names.arg = var_df$PC,
                       ylim = c(0, max(var_df$variance) * 1.2),
                       ylab = "% variance explained",
                       main = "PCA Scree Plot (Top 10 PCs)", col = "grey70")

par(new = TRUE)

# Overlay cumulative explained variance on a secondary y-axis
plot(
  bar_centers, var_df$cumvar,
  type = "b",
  axes = FALSE,
  xlab = "", ylab = "",
  xlim = par("usr")[1:2],
  ylim = c(0, 100),
  col = "red", pch = 16
)

# Add secondary axis for cumulative explained variance
axis(side = 4)
mtext("Cumulative variance (%)", side = 4, line = 3)

close_png()

# Hierarchical clustering

sample_dist <- dist(t(logcpm), method = "euclidean")
hc <- hclust(sample_dist, method = "average")

open_png(filename = "hierarchical_clustering.png", dir = paths$plots_bulk_rnaseq_sensitivity_eda, 
         width = 2200, height = 900)
par(mar = c(5, 4, 4, 2))

plot(hc, labels = FALSE,
     main = "Hierarchical clustering of samples (Euclidean, average linkage)",
     xlab = "", sub = "")
ordered_colors <- group_palette[as.character(colData$condition[hc$order])]
for (i in seq_along(hc$order)) {
  axis(1, at = i, labels = FALSE, col.ticks = ordered_colors[i], lwd.ticks = 3)
}
legend("top", legend = names(group_palette), fill = group_palette, bty = "n")

close_png()

# Quantitative cross-check: assess whether unsupervised clustering is broadly consistent with the CD/Control grouping

clusters_k2 <- cutree(hc, k = 2)
cluster_vs_disease <- table(cluster = clusters_k2, disease = colData$condition)

message("[EDA] Cross-tabulation of unsupervised clusters (k = 2) and disease status:")

message("[EDA] Fisher's exact test for cluster-disease association:")

print(fisher.test(cluster_vs_disease))


# Sample-sample correlation heatmap

cor_mat <- cor(logcpm, method = "pearson", use = "pairwise.complete.obs")

annotation_col <- data.frame(Group = colData$condition)
rownames(annotation_col) <- colnames(logcpm)

open_png(filename = "sample_correlation_heatmap.png", dir = paths$plots_bulk_rnaseq_sensitivity_eda,
         width = 1400, height = 1400)

pheatmap::pheatmap(
  cor_mat,
  annotation_col = annotation_col,
  annotation_colors = list(Group = group_palette),
  show_rownames = FALSE, show_colnames = FALSE,
  main = "Sample-sample Pearson correlation (GSE193677, Sensitivity)"
)

close_png()

# Sample-distance heatmap (1 - correlation)

sample_dist_1mcor <- as.dist(1 - cor_mat)

open_png(filename = "sample_distance_heatmap.png", dir = paths$plots_bulk_rnaseq_sensitivity_eda,
         width = 1400, height = 1400)

pheatmap::pheatmap(
  as.matrix(sample_dist_1mcor),
  annotation_col = annotation_col,
  annotation_colors = list(Group = group_palette),
  show_rownames = FALSE, show_colnames = FALSE,
  main = "Sample distance heatmap (1 - Pearson correlation)"
)

close_png()


# =========================================================
#             CD vs Control contrast
# =========================================================

message("[LIMMA] Building CD vs Control contrast...")

contrast_matrix <- limma::makeContrasts(
  CD_vs_Control = CD - Control,
  levels = design
)

print(contrast_matrix)

fit2 <- limma::contrasts.fit(fit, contrast_matrix)

# Positive logFC values indicate higher expression in CD relative to Control.


# =========================================================
#             Empirical Bayes moderation
# =========================================================

# Apply robust empirical Bayes moderation to stabilize gene-wise variance estimates.

message("[LIMMA] Applying empirical Bayes moderation...")

fit2 <- limma::eBayes(fit2, robust = TRUE)


# =========================================================
#             Differential expression results
# =========================================================

# Define significant DEGs using Benjamini-Hochberg FDR < 0.05 and an absolute log2 fold-change threshold of 1.

message("[LIMMA] Extracting differential expression results...")

deg_table <- limma::topTable(
  fit2,
  coef = "CD_vs_Control",
  number = Inf,
  adjust.method = "BH",
  sort.by = "P"
)

deg_table$gene_id <- rownames(deg_table)

message("[LIMMA] Differential expression completed.")

message("[LIMMA] Genes tested: ", nrow(deg_table))


# =========================================================
#             DEG summary
# =========================================================

message("[SUMMARY] Differential expression summary:")

deg_summary <- data.frame(
  Comparison = "CD_vs_Control",
  Genes_tested = nrow(deg_table),
  Significant_FDR_0.05 = sum(deg_table$adj.P.Val < 0.05, na.rm = TRUE),
  Up_FDR_0.05 = sum(deg_table$adj.P.Val < 0.05 & deg_table$logFC > 0, na.rm = TRUE),
  Down_FDR_0.05 = sum(deg_table$adj.P.Val < 0.05 & deg_table$logFC < 0, na.rm = TRUE)
)

print(deg_summary)


# =========================================================
#             P-value distribution diagnostic
# =========================================================

message("[DEG] Generating p-value distribution histogram...")

open_png(filename = "pvalue_histogram.png", dir = paths$plots_bulk_rnaseq_sensitivity_deg,
         width = 1000, height = 700)

hist(
  deg_table$P.Value,
  breaks = 50,
  main = "P-value distribution",
  sub = "CD vs Control | GSE193677 | Sensitivity cohort",
  xlab = "Nominal p-value",
  ylab = "Number of genes"
)

abline(v = 0.05, lty = 2, lwd = 2)

close_png()


# =========================================================
#            Gene annotation
# =========================================================

message("[DEG] Inspecting gene ID format before annotation...")

# Gene identifiers in the DEG table are Ensembl gene IDs.

keytype_used <- "ENSEMBL"  

message("[DEG] Annotating genes using keytype = ", keytype_used, "...")

gene_annot <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = deg_table$gene_id,
  columns = c("SYMBOL", "ENTREZID", "GENENAME"),
  keytype = keytype_used
)

message("[DEG] Removing genes without any annotation...")

gene_annot <- gene_annot[
  !(is.na(gene_annot$SYMBOL) & is.na(gene_annot$ENTREZID) & is.na(gene_annot$GENENAME)),
]

message("[DEG] Removing duplicated key mappings (one annotation per gene)...")

gene_annot <- gene_annot[!duplicated(gene_annot[[keytype_used]]), ]

message(sprintf("[DEG] Retrieved annotations for %d unique genes.", nrow(gene_annot)))

deg_table <- merge(
  deg_table,
  gene_annot,
  by.x = "gene_id",
  by.y = keytype_used,
  all.x = TRUE
)

deg_table <- deg_table[order(deg_table$adj.P.Val), ]

message("[DEG] Gene annotation completed.")


# =========================================================
#             Significant DEG summary (annotated)
# =========================================================

logcf_cutoff <- 1
padj_cutoff  <- 0.05

sig_deg <- subset(
  deg_table,
  adj.P.Val < padj_cutoff & abs(logFC) > logcf_cutoff
)

n_before <- nrow(sig_deg)

sig_deg <- sig_deg[!(is.na(sig_deg$SYMBOL) & is.na(sig_deg$ENTREZID) & is.na(sig_deg$GENENAME)), ]

message(sprintf(
  "[DEG] Removed %d significant gene(s) with no annotation at all (SYMBOL/ENTREZID/GENENAME all NA).",
  n_before - nrow(sig_deg)
))

message(sprintf(
  "[DEG] Significant genes (FDR < %.2f, |log2FC| > %.1f): %d upregulated, %d downregulated (%d total).",
  padj_cutoff,
  logcf_cutoff,
  sum(sig_deg$logFC > 0),
  sum(sig_deg$logFC < 0),
  nrow(sig_deg)
))


# =========================================================
#           Plotting DEGs
# =========================================================

message("[DEG] Generating volcano plot...")

deg_table$significance <- "Not Significant"
deg_table$significance[deg_table$adj.P.Val < padj_cutoff & deg_table$logFC > logcf_cutoff] <- "Up in CD"
deg_table$significance[deg_table$adj.P.Val < padj_cutoff & deg_table$logFC < -logcf_cutoff] <- "Down in CD"

volcano_plot_deg <- ggplot(deg_table, aes(x = logFC, y = -log10(adj.P.Val), color = significance)) +
  geom_point(size = 1.5, alpha = 0.7) +
  geom_vline(xintercept = c(-logcf_cutoff, logcf_cutoff), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", color = "grey40") +
  scale_color_manual(values = c("Up in CD" = "#CE2915", "Down in CD" = "#0096FF", "Not Significant" = "grey75")) +
  labs(
    title = "Volcano plot: CD vs Control (GSE193677, Ileum I + NonI, Sensitivity)",
    x = expression(log[2] ~ fold~change),
    y = expression(-log[10] ~ adjusted~italic(P)),
    color = NULL
  ) +
  theme_bw()

save_plot(volcano_plot_deg, filename = "volcano_plot.png", dir = paths$plots_bulk_rnaseq_sensitivity_deg,
          width = 8, height = 6)

message("[DEG] Generating MA plot...")

open_png(filename = "MA_plot_CD_vs_Control.png", dir = paths$plots_bulk_rnaseq_sensitivity_deg, 
         width = 1000, height = 800)

limma::plotMA(fit2, coef = "CD_vs_Control",
              status = deg_table$significance[match(rownames(fit2), deg_table$gene_id)],
              values = c("Up in CD", "Down in CD"), col = c("firebrick", "steelblue"),
              main = "MA plot: CD vs Control (GSE193677, Sensitivity)")
abline(h = 0, col = "grey40", lty = 2)

close_png()

# Visualize the top significant annotated DEGs using row-scaled voom log2-CPM values.

message("[DEG] Generating heatmap of top DEGs...")

top_degs <- 100

deg_annotated <- sig_deg[!is.na(sig_deg$SYMBOL) & sig_deg$SYMBOL != "", ]
deg_annotated <- deg_annotated[order(deg_annotated$adj.P.Val), ]
deg_annotated <- deg_annotated[!duplicated(deg_annotated$SYMBOL), ]

message(sprintf("[DEG] %d unique annotated genes available after deduplication.", nrow(deg_annotated)))

if (nrow(deg_annotated) < top_degs) {
  warning(sprintf(
    "[DEG] Only %d significant annotated genes available (requested %d) - heatmap will show all %d.",
    nrow(deg_annotated), top_degs, nrow(deg_annotated)
  ))
  top_degs <- nrow(deg_annotated)
}

top_ids <- head(deg_annotated$gene_id, top_degs)
top_symbols <- head(deg_annotated$SYMBOL, top_degs)

heatmap_matrix <- logcpm[top_ids, , drop = FALSE]
heatmap_matrix_scaled <- t(scale(t(heatmap_matrix)))
rownames(heatmap_matrix_scaled) <- make.unique(top_symbols)

ann_col <- ComplexHeatmap::HeatmapAnnotation(
  Group = colData$condition,
  col = list(Group = group_palette)
)

open_png(filename = "Heatmap_TopDegs_CD_vs_Control.png", dir = paths$plots_bulk_rnaseq_sensitivity_deg,
         width = 1600, height = 1800)

ComplexHeatmap::Heatmap(
  heatmap_matrix_scaled,
  name = "Z-score",
  top_annotation = ann_col,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = FALSE,
  row_names_gp = grid::gpar(fontsize = 7),
  column_title = sprintf("Top %d Differentially Expressed Genes (CD vs Control)", top_degs),
  row_title = sprintf("Top %d DEGs", top_degs),
  column_title_gp = grid::gpar(fontsize = 12, fontface = "bold"),
  heatmap_legend_param = list(title = "Row\nZ-score")
)

close_png()


# =========================================================
#                Save outputs
# =========================================================

message("[OUTPUT] Saving voom-limma model objects...")

save_rds(fit, filename = "voomLmFit_object_sensitivity.rds", dir = paths$objects_bulkrnaseq)
save_rds(fit2, filename = "limma_fit2_sensitivity.rds", dir = paths$objects_bulkrnaseq)

message("[OUTPUT] Saving differential expression results...")

save_csv(deg_table, filename = "DEGs_GSE193677_sensitivity.csv", dir = paths$tables_bulkrnaseq)

message("[OUTPUT] Saving background gene list for enrichment analyses...")

# Define the ORA background as all genes retained after expression filtering and included in the differential expression model.

background <- rownames(fit$EList$E)

save_csv(background, filename = "background_GSE193677_sensitivity.csv", dir = paths$tables_bulkrnaseq)

message("[OUTPUT] Saving outliers tables...")

save_csv(qc_outliers_df, filename = "qc_outliers_sensitivity.csv", dir = paths$tables_bulkrnaseq)


# =========================================================
#                  Save session info
# =========================================================

message("[OUTPUT] Saving session information...")

save_session_info(filename = "sessionInfo_DEG_GSE193677_sensitivity.txt", dir = paths$logs,
                  label = "Differential expression analysis - GSE193677 (Sensitivity: Ileum I + NonI, CD vs Control)")

message("[OUTPUT] Session information saved to: ", paths$logs)


# =========================================================
#                  Final pipeline message
# =========================================================

message("=================================================")
message("[PIPELINE] Differential expression analysis completed successfully for GSE193677 (Sensitivity: Ileum I + NonI, CD vs Control).")
message("=================================================")




