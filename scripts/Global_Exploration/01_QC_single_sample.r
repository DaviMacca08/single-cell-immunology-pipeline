# =========================================================
# Single-sample scRNA-seq preprocessing workflow
# Quality control, filtering, and normalization
# This script represents the per-sample processing stage
# prior to multi-sample integration and batch correction
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

message("=== Starting QC pipeline ===")

message("Set the correct directory where the raw data is saved...")


# =========================================================
#                  Load sample
# =========================================================

message("Loading 10X data...")

sample <- Read10X(data.dir = file.path("GRCh38"))

if (is.null(sample)) {
  stop("ERROR: Read10X failed. Check input directory.")
}

message("Creating Seurat object...")

srt_obj <- CreateSeuratObject(counts = sample, min.cells = 3, min.features = 200)

check_seurat(srt_obj)

message("Seurat object created successfully")


# =========================================================
#                  Metadata
# =========================================================

message("Adding metadata...")

srt_obj$sample_id <- "GSM4766849"
srt_obj$compartment <- "IEL"
srt_obj$condition <- "mixed_inflamed_non_inflamed"
srt_obj$disease <- "CD"
srt_obj$donor <- "1818"
srt_obj$origin <- "unknown"

sample_id <- srt_obj$sample_id

message(paste("Processing sample:", unique(srt_obj$sample_id)))


# =========================================================
#      Create sample-specific QC plot directory
# =========================================================

sample_plot_dir <- file.path(paths$plots_qc, unique(srt_obj$sample_id))

if (!dir.exists(sample_plot_dir)) {
  dir.create(sample_plot_dir, recursive = TRUE)
  message("Created directory: ", sample_plot_dir)
} else {
  message("Directory already exists: ", sample_plot_dir)
}


# =========================================================
#                  QC metrics
# =========================================================

message("Computing mitochondrial percentage...")

srt_obj[["percent.mt"]] <- PercentageFeatureSet(srt_obj, pattern = "^MT-")

message("Generating QC plots...")

vln_p <- VlnPlot(srt_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),ncol = 3)

save_plot(vln_p, filename = "01-QC_ViolinPlot.png", dir = sample_plot_dir)

scat1 <- FeatureScatter(srt_obj, "nCount_RNA", "percent.mt", pt.size = 0.3)
scat2 <- FeatureScatter(srt_obj, "nCount_RNA", "nFeature_RNA", pt.size = 0.3)

save_plot(scat1 + scat2, filename = "02-QC_ScatterPlots.png",dir = sample_plot_dir)



# =========================================================
#                  Density QC plots
# =========================================================

message("Generating density plots...")

density_plot <- function(obj, feature, color){
  
  stopifnot(inherits(obj, "Seurat"))
  stopifnot(feature %in% colnames(obj[[]]))
  
  p <- ggplot(obj[[]], aes(x = .data[[feature]])) +
    geom_density(fill = color, alpha = 0.7, linewidth = 0.8) +
    coord_cartesian(xlim = c(0, max(obj[[feature]][,1], na.rm = TRUE))) +
    labs(
      title = paste("Distribution of", feature),
      x = feature,
      y = "Density"
    ) +
    theme_classic(base_size = 14) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"), axis.title = element_text(face = "bold"))
  
  return(p)
  
}

ridge1 <- density_plot(srt_obj, "nFeature_RNA", "coral2")
ridge2 <- density_plot(srt_obj, "nCount_RNA", "lightblue")
ridge3 <- density_plot(srt_obj, "percent.mt", "lightgreen")

save_plot(ridge1 + ridge2 + ridge3, filename = "03-QC_DensityPlot.png",dir = sample_plot_dir)


# =========================================================
#                  Filtering
# =========================================================

message("Filtering low-quality cells...")
message(paste("Cells before filtering:", ncol(srt_obj)))

srt_obj <- subset(srt_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 4500 & nCount_RNA < 33000 & percent.mt < 25)

message(paste("Cells after filtering:", ncol(srt_obj)))

if (ncol(srt_obj) == 0) {
  stop("ERROR: No cells left after filtering. Check thresholds.")
}

vln_filtered <- VlnPlot(srt_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

save_plot(vln_filtered, filename = "04-QC_ViolinPlot_filtered.png", dir = sample_plot_dir)


# =========================================================
#                  Save filtered object
# =========================================================

message("Saving filtered object...")

save_rds(srt_obj, paste0(unique(srt_obj$sample_id), "_Filtered.rds"), dir = paths$filtered)

message(paste("QC and preprocessing completed for sample:", unique(srt_obj$sample_id)))

message("QC pipeline completed successfully!")

# =========================================================
#                  Session information
# =========================================================

message("Saving session information...")

save_session_info(dir = paths$logs, filename = paste0(unique(srt_obj$sample_id), "_sessionInfo_QC.txt"), label = "QC stage")

message("Session info saved at: ", paths$logs)

# =========================================================
#                  Final pipeline message
# =========================================================

message("=== END OF PIPELINE: single-sample QC completed ===")
