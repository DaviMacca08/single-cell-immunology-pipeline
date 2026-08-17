# =========================================================
# Project paths (single source of truth)
# =========================================================

base_dir <- "/Users/davidemaccarrone/Desktop/Bioinformatics/MyProjects🤞🏻/scRNA-Seq/UPDATE"

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
  objects_gsea = file.path(base_dir, "objects/gsea"),
  objects_cellchat = file.path(base_dir, "objects/cellchat"),
  objects_liana = file.path(base_dir, "objects/liana"),
  
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
  
    # plots GSEA for IEL compartment
    plots_gsea_IEL = file.path(base_dir, "results/plots/GSEA_IEL"),
  
    # plots GSEA for LPL compartment
    plots_gsea_LPL = file.path(base_dir, "results/plots/GSEA_LPL"),
  
    # plots CellChat for IEL compartment
    plots_cellchat_iel = file.path(base_dir, "results/plots/CellChat_IEL"),
    plots_cellchat_iel_path = file.path(base_dir, "results/plots/CellChat_IEL/pathways"), # folder for plot of pathways
    plots_cellchat_iel_merge = file.path(base_dir, "results/plots/CellChat_IEL/merge"), # folder for plot of merged obj
    plots_cellchat_iel_merge_celltype_interes = file.path(base_dir, "results/plots/CellChat_IEL/merge/ct_interest"), # folder for plot specific celltypes
    plots_cellchat_iel_merge_pathways = file.path(base_dir, "results/plots/CellChat_IEL/merge/pathways"), # folder for plot final pathways
    plots_cellchat_iel_merge_pathways_bubble = file.path(base_dir, "results/plots/CellChat_IEL/merge/pathways/bubble"), # folder for bubble plot final pathways
    plots_cellchat_iel_merge_pathways_violin = file.path(base_dir, "results/plots/CellChat_IEL/merge/pathways/violin"), # folder for violin plot final pathways
    plots_cellchat_iel_merge_pathways_chord = file.path(base_dir, "results/plots/CellChat_IEL/merge/pathways/chord"), # folder for chord plot final pathways
    plots_cellchat_iel_merge_pathways_hierarchy = file.path(base_dir, "results/plots/CellChat_IEL/merge/pathways/hierarchy"), # folder for hierarchy plot final pathways
    plots_cellchat_iel_merge_pathways_pattern = file.path(base_dir, "results/plots/CellChat_IEL/merge/pathways/pattern"), # folder for pattern plot final pathways
    plots_cellchat_iel_merge_pathways_sigrole = file.path(base_dir, "results/plots/CellChat_IEL/merge/pathways/signaling_role"), # folder for signaling role plot final pathways
   
   # plots LIANA for LPL compartment
    plots_liana_iel = file.path(base_dir, "results/plots/LIANA_IEL"),
  
  # plots CellChat for LPL compartment
    plots_cellchat_lpl = file.path(base_dir, "results/plots/CellChat_LPL"),
    plots_cellchat_lpl_path = file.path(base_dir, "results/plots/CellChat_LPL/pathways"), # folder for plot of pathways
    plots_cellchat_lpl_merge = file.path(base_dir, "results/plots/CellChat_LPL/merge"), # folder for plot of merged obj
    plots_cellchat_lpl_merge_celltype_interes = file.path(base_dir, "results/plots/CellChat_LPL/merge/ct_interest"), # folder for plot specific celltypes
    plots_cellchat_lpl_merge_pathways = file.path(base_dir, "results/plots/CellChat_LPL/merge/pathways"), # folder for plot final pathways
    plots_cellchat_lpl_merge_pathways_bubble = file.path(base_dir, "results/plots/CellChat_LPL/merge/pathways/bubble"), # folder for bubble plot final pathways
    plots_cellchat_lpl_merge_pathways_violin = file.path(base_dir, "results/plots/CellChat_LPL/merge/pathways/violin"), # folder for violin plot final pathways
    plots_cellchat_lpl_merge_pathways_chord = file.path(base_dir, "results/plots/CellChat_LPL/merge/pathways/chord"), # folder for chord plot final pathways
    plots_cellchat_lpl_merge_pathways_hierarchy = file.path(base_dir, "results/plots/CellChat_LPL/merge/pathways/hierarchy"), # folder for pattern plot final pathways
    plots_cellchat_lpl_merge_pathways_pattern = file.path(base_dir, "results/plots/CellChat_LPL/merge/pathways/pattern"), # folder for hierarchy plot final pathways
    plots_cellchat_lpl_merge_pathways_sigrole = file.path(base_dir, "results/plots/CellChat_LPL/merge/pathways/signaling_role"), # folder for signaling role plot final pathways
  
    # plots LIANA for LPL compartment
    plots_liana_lpl = file.path(base_dir, "results/plots/LIANA_LPL"),
  
  # tables inside results
  tables = file.path(base_dir, "results/tables"),
  
  # tables for markers
  tables_markers = file.path(base_dir, "results/tables/Markers"),
  
    # tables for singleR annotations
    tables_singler = file.path(base_dir, "results/tables/Singler"),
    
    # tables for azimuth annotations
    tables_azimuth = file.path(base_dir, "results/tables/Azimuth"),
  
    # tables GSEA for IEL compartment
    tables_gsea_IEL = file.path(base_dir, "results/tables/IEL_GSEA"), 
  
    # tables GSEA for LPL compartment
    tables_gsea_LPL = file.path(base_dir, "results/tables/LPL_GSEA"),
  
  # logs inside results
  logs = file.path(base_dir, "results/logs")
  
)

# Ensure directories exist
invisible(lapply(paths, function(x) {
  if (!dir.exists(x)) dir.create(x, recursive = TRUE, showWarnings = FALSE)
}))
