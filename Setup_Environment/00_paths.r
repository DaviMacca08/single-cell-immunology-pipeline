# =========================================================
# Project paths (single source of truth)
# =========================================================

base_dir <- "/Users/davidemaccarrone/Desktop/Bioinformatics/MyProjects🤞🏻/scRNA-Seq/v2.0_ongoing"

paths <- list(
  
  base = base_dir,
  
  # data
  raw = file.path(base_dir, "data/raw"),
  filtered = file.path(base_dir, "data/filtered_SeuratObjects"),
  metadata = file.path(base_dir, "data/metadata"),
  
  # Azimuth ref
  Azimuth_ref = file.path(base_dir, "Azimuth_PBMC_ref/"),
  
  # objects
  objects = file.path(base_dir, "objects"),
  objects_global = file.path(base_dir, "objects/global"),
  objects_subsetting = file.path(base_dir, "objects/subsetting"),
  objects_pseudobulk = file.path(base_dir, "objects/pseudobulk"),
  
  # analysis
  analysis = file.path(base_dir, "analysis"),
  analysis_global = file.path(base_dir, "analysis/global_exploration"),
  analysis_subsetting = file.path(base_dir, "analysis/subsetting"),
  analysis_pseudobulk = file.path(base_dir, "analysis/pseudobulk"),
  
  # results
  results = file.path(base_dir, "results"),
  
  # plots inside results
  plots = file.path(base_dir, "results/plots"),
  
    # plots quality control per sample
    plots_qc = file.path(base_dir, "results/plots/QC_by_sample"),
    
    # plots integration Harmony
    plots_integration_hvg = file.path(base_dir, "results/plots/Integration_HVG"),
    
    # plots global exploration
    plots_global = file.path(base_dir, "results/plots/Global_Exploration"),
    plots_singler = file.path(base_dir, "results/plots/Global_Exploration/SingleR"),
    plots_azimuth = file.path(base_dir, "results/plots/Global_Exploration/Azimuth_PBMC"),
    plots_cellannotation = file.path(base_dir, "results/plots/Global_Exploration/CellAnnotation"),
    
    # plots IEL compartment
    plots_IEL_compartment = file.path(base_dir, "results/plots/IEL_compartment"),
    
    # plots LPL compartment
    plots_LPL_compartment = file.path(base_dir, "results/plots/LPL_compartment"),
    
    # plots pseudobulk for IEL compartment
    plots_pb_IEL = file.path(base_dir, "results/plots/PseudoBulk_IEL"),
  
    # plots pseudobulk for LPL compartment
    plots_pb_LPL = file.path(base_dir, "results/plots/PseudoBulk_LPL"),
    
  # tables inside results
  tables = file.path(base_dir, "results/tables"),
  
  # tables for markers
  tables_markers = file.path(base_dir, "results/tables/Markers"),
  
    # tables for singleR annotations
    tables_singler = file.path(base_dir, "results/tables/Singler"),
    
    # tables for azimuth annotations
    tables_azimuth = file.path(base_dir, "results/tables/Azimuth"),
  
  # logs inside results
  logs = file.path(base_dir, "results/logs")
  
)

# Ensure directories exist
invisible(lapply(paths, function(x) {
  if (!dir.exists(x)) dir.create(x, recursive = TRUE, showWarnings = FALSE)
}))
