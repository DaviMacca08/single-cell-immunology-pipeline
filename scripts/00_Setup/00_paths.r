# =========================================================
# Project paths (single source of truth)
# =========================================================

base_dir <- "/Users/davidemaccarrone/Desktop/Bioinformatics/MyProjects🤞🏻/scRNA-Seq/v3.0-github"

paths <- list(
  
  base = base_dir,
  
  # data
  filtered = file.path(base_dir, "filtered_SeuratObjects"),
  
  # objects
  objects = file.path(base_dir, "objects"),
  
  # plots
  plots = file.path(base_dir, "plots"),
  plots_qc = file.path(base_dir, "plots/QC_by_sample"),
  plots_integration_hvg = file.path(base_dir, "plots/Integration_HVG"),
  plots_global = file.path(base_dir, "plots/Global_Exploration"),
  plot_singler = file.path(base_dir, "plots/SingleR"),
  plot_azimuth = file.path(base_dir, "plots/Azimuth_PBMC"),
  plot_cellannotation = file.path(base_dir, "plots/CellAnnotation"),
  plots_IEL_compartment = file.path(base_dir, "plots/IEL_compartment"),
  plots_LPL_compartment = file.path(base_dir, "plots/LPL_compartment"),
  
  # tables
  tables = file.path(base_dir, "tables"),
  tables_singler = file.path(base_dir, "tables/singler"),
  tables_azimuth = file.path(base_dir, "tables/azimuth"),
  
  # logs
  logs = file.path(base_dir, "logs"),
  
  # Azimuth ref
  Azimuth_ref = file.path(base_dir, "Azimuth_PBMC_ref/")
)

# Ensure directories exist
invisible(lapply(paths, function(x) {
  if (!dir.exists(x)) dir.create(x, recursive = TRUE, showWarnings = FALSE)
}))