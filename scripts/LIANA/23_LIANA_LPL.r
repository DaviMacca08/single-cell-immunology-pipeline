# =========================================================
# LIANA cell-cell communication analysis - LPL
# =========================================================
# Input:
# - Seurat object of LPL compartment with consensus annotation
# - Global Seurat object containing raw counts
# - Disease metadata: Control and CD
#
# Output:
# - LIANA ligand-receptor interaction analysis
# - Ligand-receptor interaction scores
# - Differential cell-cell communication analysis between CD and Control
# - Visualization of significant ligand-receptor interactions
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

message("Starting LIANA analysis for LPL compartment...")


# =========================================================
#                   Load input data
# =========================================================

message("Loading LPL Seurat object...")

lpl_cells <- readRDS(file.path(paths$objects_subsetting, "srt_LPL_compartment_celltype.rds"))

message("LPL object loaded: ", nrow(lpl_cells), " genes × ", ncol(lpl_cells)," cells")
message("Loading global Seurat object...")

full_obj <- readRDS(file.path(paths$objects_global, "srt_obj_merge_harmony.rds"))

message("Global object loaded: ", nrow(full_obj), " genes × ", ncol(full_obj)," cells")


# =========================================================
# Validate input data
# =========================================================

message("Validating LPL normalized expression layer...")

stopifnot(dim(LayerData(lpl_cells,
                        assay = "RNA",
                        layer = "data")) == dim(lpl_cells))

message("  ✓ LPL data layer dimensions validated")

stopifnot("celltype_LPL" %in% colnames(lpl_cells@meta.data))

message("  ✓ LPL cell-type annotation found")

stopifnot("disease" %in% colnames(lpl_cells@meta.data))

message("  ✓ Disease metadata found")


# =========================================================
# Prepare expression matrices
# =========================================================

message("Preparing expression matrices for LIANA...")

message("Joining raw count layers in global object...")

full_obj <- JoinLayers(full_obj, assay = "RNA")

message("  ✓ Raw count layers joined")

# Extract normalized expression from LPL object

message("Extracting normalized expression from LPL data layer...")

logcounts_mat <- LayerData(lpl_cells, assay = "RNA", layer = "data")

message("  ✓ Normalized expression extracted: ", nrow(logcounts_mat), " genes × ", ncol(logcounts_mat)," cells")

# Extract raw counts from global object

message("Extracting raw counts from global object...")

counts_full <- LayerData(full_obj, assay = "RNA", layer = "counts")

message("  ✓ Raw counts extracted: ", nrow(counts_full), " genes × ", ncol(counts_full), " cells")

# Subset raw counts to LPL cells

message("Subsetting raw counts to LPL cells...")

counts_mat <- counts_full[, colnames(lpl_cells), drop = FALSE]

message("  ✓ LPL raw counts extracted: ", nrow(counts_mat), " genes × ", ncol(counts_mat)," cells")


# =========================================================
# Validate expression matrices
# =========================================================

message("Validating expression matrices...")

stopifnot(identical(colnames(counts_mat), colnames(logcounts_mat)))

message("  ✓ Cell identities and ordering are identical")

stopifnot(identical(rownames(counts_mat), rownames(logcounts_mat)))

message("  ✓ Gene identities and ordering are identical")

stopifnot(dim(counts_mat) == dim(logcounts_mat))

message("  ✓ Counts and logcounts dimensions are identical")

message("Expression matrices validated: ", nrow(logcounts_mat), " genes × ", ncol(logcounts_mat)," cells")


# =========================================================
# Prepare cell metadata
# =========================================================

message("Preparing LPL cell metadata...")

meta_lpl <- lpl_cells@meta.data

message("  ✓ Metadata extracted for ", nrow(meta_lpl), " cells")

message("  ✓ Cell types identified: ", length(unique(lpl_cells$celltype_LPL)))

message("  ✓ Disease groups: ",paste(unique(lpl_cells$disease), collapse = ", "))


# =========================================================
# Build SingleCellExperiment object
# =========================================================

message("Building SingleCellExperiment object...")

lpl_sce <- SingleCellExperiment(assays = list(logcounts = logcounts_mat, counts = counts_mat), colData = meta_lpl)

# Set cell identities

colLabels(lpl_sce) <- lpl_cells$celltype_LPL

message("  ✓ SingleCellExperiment created: ", nrow(lpl_sce), " genes × ", ncol(lpl_sce)," cells")


# =========================================================
# Validate SingleCellExperiment object
# =========================================================

message("Validating SingleCellExperiment object...")

stopifnot(identical(assayNames(lpl_sce),c("logcounts", "counts")))

message("  ✓ counts and logcounts assays detected")

stopifnot(ncol(lpl_sce) == ncol(lpl_cells))

message("  ✓ Cell number validated")

stopifnot(nrow(lpl_sce) == nrow(lpl_cells))

message("  ✓ Gene number validated")

stopifnot(length(colLabels(lpl_sce)) == ncol(lpl_sce))

message("  ✓ Cell labels validated")

message("SingleCellExperiment validation completed successfully")


# =========================================================
# Inspect LIANA resources and methods
# =========================================================

message("Inspecting available LIANA resources...")

show_resources()

message("Inspecting available LIANA methods...")

show_methods()


# =========================================================
# Run LIANA
# =========================================================

message("Running LIANA for LPL compartment...")

liana_test_lpl <- liana_wrap(
  lpl_sce,
  idents_col = "celltype_LPL",
  method = c("natmi", "connectome", "logfc", "sca", "cellphonedb"),
  resource = "Consensus")

message("LIANA analysis completed successfully.")

# Aggregate results 

liana_test_lpl_aggregate <- rank_aggregate(liana_test_lpl)

message("LIANA aggregated results dimensions: ", nrow(liana_test_lpl_aggregate), " rows × ", ncol(liana_test_lpl_aggregate)," columns")


# =========================================================
# Inspect aggregated results
# =========================================================

message("Inspecting aggregated LIANA results...")

print(colnames(liana_test_lpl_aggregate))

message(
  "Aggregated LIANA results: ",
  nrow(liana_test_lpl_aggregate),
  " rows × ",
  ncol(liana_test_lpl_aggregate),
  " columns"
)


# =========================================================
# Global LIANA heatmap and dotplot
# =========================================================

message("Generating global LPL communication frequency heatmap...")

liana_trunc <- liana_test_lpl_aggregate |> filter(specificity_rank <= 0.05)


open_png(filename = "LIANA_global_specificity_Heatmap_LPL.png", dir = paths$plots_liana_lpl,
         width = 3000, height = 1900)

heat_freq(liana_trunc,
          row_title = "Sender (Cell type)",
          column_title = "Receiver (Cell type)",
          name = "Number of LR interactions",
          font_size = 8) 

close_png()

message("Generating specificity-focused LIANA dotplot...")

# 01 - Broad senders → Activated effector Treg (MHC-I focus)
open_png(filename = "LIANA_LPL_activated_effector_treg_dotplot.png", dir = paths$plots_liana_lpl,
         width = 2000, height = 1800)

liana_test_lpl_aggregate |>
  arrange(specificity_rank) |> 
  liana_dotplot(target_groups = "Activated effector Treg",
                ntop = 20,
                magnitude = "magnitude_rank",
                specificity = "specificity_rank",
                invert_magnitude = TRUE,
                invert_specificity = TRUE,
                y.label = "Ligand → Receptor",
                size.label = "Interaction specificity",
                colour.label = "Expression magnitude",
                size_range = c(3, 9)) + 
  labs(
    title = "LPL communication: interactions targeting Activated effector Treg",
    subtitle = "Ranked by specificity | Size: specificity | Colour: magnitude",
    x = NULL,
    y = "Ligand → Receptor"
  ) +
  scale_x_discrete(position = "bottom") +
  theme(
    strip.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 11),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
    axis.text.y = element_text(size = 9)
  )

close_png()

# 02 - Innate-like/TRM CD8 → T follicular helper, CD4 TRM-like, entrambe le Treg
open_png(filename = "LIANA_LPL_innatelike_trm_cd8_to_il2_targets_dotplot.png", dir = paths$plots_liana_lpl,
         width = 2500, height = 1800)

liana_test_lpl_aggregate |>
  arrange(magnitude_rank) |> 
  liana_dotplot(source_groups = "Innate-like / TRM CD8",
                target_groups = c("T follicular helper",
                                  "CD4 TRM-like activated memory T cells",
                                  "Activated effector Treg",
                                  "FOXP3+ activated Treg"),
                ntop = 30,
                magnitude = "magnitude_rank",
                specificity = "specificity_rank",
                invert_magnitude = TRUE,
                invert_specificity = TRUE,
                y.label = "Ligand → Receptor",
                size.label = "Interaction specificity",
                colour.label = "Expression magnitude",
                size_range = c(3, 9)) + 
  labs(
    title = "LPL communication: Innate-like/TRM CD8 → IL2-axis targets",
    subtitle = "Ranked by expression magnitude | Size: specificity | Colour: magnitude",
    x = NULL,
    y = "Ligand → Receptor"
  ) +
  scale_x_discrete(position = "bottom") +
  theme(
    strip.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 11),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10)
  )

close_png()

# 03 - SIRP core: Cycling T cells, CD4 TRM-like, Activated effector Treg (tutte le combinazioni)
open_png(filename = "LIANA_LPL_sirp_core_populations_dotplot.png", dir = paths$plots_liana_lpl,
         width = 2000, height = 1800)

liana_test_lpl_aggregate |>
  arrange(specificity_rank) |> 
  liana_dotplot(source_groups = c("Cycling T cells",
                                  "CD4 TRM-like activated memory T cells",
                                  "Activated effector Treg"),
                target_groups = c("Cycling T cells",
                                  "CD4 TRM-like activated memory T cells",
                                  "Activated effector Treg"),
                ntop = 20,
                magnitude = "magnitude_rank",
                specificity = "specificity_rank",
                invert_magnitude = TRUE,
                invert_specificity = TRUE,
                y.label = "Ligand → Receptor",
                size.label = "Interaction specificity",
                colour.label = "Expression magnitude",
                size_range = c(3, 9)) + 
  labs(
    title = "LPL communication: SIRP-core populations (Cycling, CD4 TRM-like, Activated effector Treg)",
    subtitle = "Ranked by specificity | Size: specificity | Colour: magnitude",
    x = NULL,
    y = "Ligand → Receptor"
  ) +
  scale_x_discrete(position = "bottom") +
  theme(
    strip.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 11),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
    axis.text.y = element_text(size = 9)
  )

close_png()

# 04 - MHC-II axis: T follicular helper (receiver stabile) + Treg (receiver sostituito)
open_png(filename = "LIANA_LPL_mhcii_tfh_treg_targets_dotplot.png", dir = paths$plots_liana_lpl,
         width = 2000, height = 1800)

liana_test_lpl_aggregate |>
  arrange(specificity_rank) |> 
  liana_dotplot(target_groups = c("T follicular helper",
                                  "Activated effector Treg",
                                  "FOXP3+ activated Treg"),
                ntop = 20,
                magnitude = "magnitude_rank",
                specificity = "specificity_rank",
                invert_magnitude = TRUE,
                invert_specificity = TRUE,
                y.label = "Ligand → Receptor",
                size.label = "Interaction specificity",
                colour.label = "Expression magnitude",
                size_range = c(3, 9)) + 
  labs(
    title = "LPL communication: T follicular helper and Treg populations (MHC-II axis)",
    subtitle = "Ranked by specificity | Size: specificity | Colour: magnitude",
    x = NULL,
    y = "Ligand → Receptor"
  ) +
  scale_x_discrete(position = "bottom") +
  theme(
    strip.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 11),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
    axis.text.y = element_text(size = 9)
  )

close_png()


# ============================================================
# Global LPL communication - chord plot
# ============================================================

message("Preparing LIANA global LPL chord plot...")

# Filter interactions based on specificity
liana_trunc <- liana_test_lpl_aggregate |> 
  filter(specificity_rank <= 0.05, magnitude_rank <= 0.25)

# Cell types selected for visualization

celltypes_of_interest_chord1 <- c(
  "Activated effector Treg",              
  "Activated cytotoxic CD8 T cells",       
  "CD4 TRM-like activated memory T cells"
)

celltypes_of_interest_chord2 <- c(
  "Naive CD4 T cells",                     
  "FOXP3+ activated Treg",                 
  "Th17-like activated CD4 T cells"        
)

# ------------------------------------------------------------
# Chord plot 1
# ------------------------------------------------------------

message("Generating chord diagram 1...")

open_png(filename = "LIANA_LPL_global_chord1.png", dir = paths$plots_liana_lpl,
         width = 11000, height = 8000)

chord_freq(
  liana_trunc,
  source_groups = celltypes_of_interest_chord1,
  target_groups = celltypes_of_interest_chord1,
  cex = 2,
  facing = "clockwise",
  adj = c(0, 0.5)
)

close_png()

# ------------------------------------------------------------
# Chord plot 2
# ------------------------------------------------------------

message("Generating chord diagram 2...")

open_png(filename = "LIANA_LPL_global_chord2.png", dir = paths$plots_liana_lpl,
         width = 11000, height = 8000
)

chord_freq(
  liana_trunc,
  source_groups = celltypes_of_interest_chord2,
  target_groups = celltypes_of_interest_chord2,
  cex = 2,
  facing = "clockwise",
  adj = c(0, 0.5)
)

close_png()

message("Individual chord diagrams generated.")


# ============================================================
# CD vs Control comparison - LPL compartment
# ============================================================

message("Starting LIANA comparison between CD and Control in the LPL compartment...")

# Compare CD vs CTRL in LPL compartment
liana_test_lpl_disease <- liana_bysample(
  sce = lpl_sce,
  sample_col = "disease",
  idents_col = "celltype_LPL",
  aggregate_how = "both",
  method = c("natmi", "connectome", "logfc", "sca", "cellphonedb"),
  resource = "Consensus",
  inplace = FALSE,
  verbose = TRUE
)

message("LIANA CD vs Control analysis completed.")

# Check available disease groups
message("Available disease groups: ", paste(names(liana_test_lpl_disease), collapse = ", "))

if (!all(c("Control", "CD") %in% names(liana_test_lpl_disease))) {
  stop("Both Control and CD groups are required for the comparison.")
}

message("✓ Both Control and CD groups are available.")

names(liana_test_lpl_disease)

# ------------------------------------------------------------
# Separate Control and CD results
# ------------------------------------------------------------

message("Extracting Control and CD LIANA results...")

liana_ctrl <- liana_test_lpl_disease[["Control"]]

liana_cd <- liana_test_lpl_disease[["CD"]]

message("Control results extracted.")
message("CD results extracted.")

# Summary of available results
message("Control result dimensions: ", paste(dim(liana_ctrl), collapse = " x "))
message("CD result dimensions: ", paste(dim(liana_cd), collapse = " x "))


# ============================================================
# CD vs Control - adequately represented LPL populations
# ============================================================

message("Generating CD vs Control LIANA dotplots for balance populations...")

tier_lpl <- table(lpl_cells$celltype_LPL, lpl_cells$disease) |> 
  as.data.frame() |> 
  pivot_wider(names_from = Var2, values_from = Freq) |> 
  rename(celltype = Var1) |> 
  mutate(tier = case_when(
    CD == 0 | Control == 0 ~ "condition_exclusive",
    pmin(CD, Control) < 30 ~ "imbalanced",   
    TRUE ~ "fully_comparable"
  ))

fully_comparable_types <- tier_lpl |> filter(tier == "fully_comparable") |>  pull(celltype)

# ------------------------------------------------------------
# CD
# ------------------------------------------------------------

p_mhc2_cd <- liana_cd |>
  arrange(specificity_rank) |>
  liana_dotplot(
    target_groups = fully_comparable_types,
    source_groups = fully_comparable_types,
    ntop = 20,
    magnitude = "magnitude_rank",
    specificity = "specificity_rank",
    invert_magnitude = TRUE,
    invert_specificity = TRUE,
    y.label = "Ligand → Receptor",
    size.label = "Interaction specificity",
    colour.label = "Expression magnitude",
    size_range = c(3, 9)
  ) +
  labs(
    title = "CD",
    subtitle = "Ranked by specificity | Size: specificity | Colour: magnitude",
    x = NULL,
    y = "Ligand → Receptor"
  ) +
  scale_x_discrete(position = "bottom") +
  theme(
    strip.text.x = element_text(
      angle = 90,
      hjust = 0.5,
      vjust = 0.5,
      size = 11
    ),
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1,
      size = 10
    )
  )

# ------------------------------------------------------------
# Control
# ------------------------------------------------------------

p_mhc2_ctrl <- liana_ctrl |>
  arrange(specificity_rank) |>
  liana_dotplot(
    target_groups =  fully_comparable_types,
    source_groups = fully_comparable_types,
    ntop = 20,
    magnitude = "magnitude_rank",
    specificity = "specificity_rank",
    invert_magnitude = TRUE,
    invert_specificity = TRUE,
    y.label = "Ligand → Receptor",
    size.label = "Interaction specificity",
    colour.label = "Expression magnitude",
    size_range = c(3, 9)
  ) +
  labs(
    title = "Control",
    subtitle = "Ranked by specificity | Size: specificity | Colour: magnitude",
    x = NULL,
    y = "Ligand → Receptor"
  ) +
  scale_x_discrete(position = "bottom") +
  theme(
    strip.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 11),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10)
  )

# ------------------------------------------------------------
# Combine and save
# ------------------------------------------------------------

message("Combining CD and Control TH17/Treg dotplots...")

open_png(filename = "LIANA_LPL_CD_vs_Control_balanced_populations_dotplot.png", dir = paths$plots_liana_lpl,
         width = 6000, height = 3000)

p_mhc2_cd + p_mhc2_ctrl

close_png()

message("CD vs Control dotplot saved.")


# =======================================================================
# Targeted inspection - pathway deep-dive LPL (MHC-I, IL2, SIRP, MHC-II)
# =======================================================================

message("Inspecting MHC-I interactions targeting Activated effector Treg...")

bind_rows(
  liana_cd  |> filter(target == "Activated effector Treg") |> mutate(condition = "CD"),
  liana_ctrl |> filter(target == "Activated effector Treg") |> mutate(condition = "Control")
) |> 
  filter(str_detect(ligand.complex, "HLA-A|HLA-B|HLA-C|HLA-E|B2M")) |> 
  select(condition, target, ligand.complex, receptor.complex, magnitude_rank, specificity_rank) |> 
  arrange(condition, specificity_rank) |> 
  print(n = Inf)

message("Inspecting IL2 interactions from Innate-like / TRM CD8...")

bind_rows(
  liana_cd   |> filter(source == "Innate-like / TRM CD8") |> mutate(condition = "CD"),
  liana_ctrl |> filter(source == "Innate-like / TRM CD8") |> mutate(condition = "Control")
) |> 
  filter(str_detect(ligand.complex, "^IL2$")) |> 
  select(condition, source, target, ligand.complex, receptor.complex, magnitude_rank, specificity_rank) |> 
  arrange(condition, target) |> 
  print(n = Inf)

message("Inspecting SIRPG → CD47 among SIRP-core populations...")

sirp_core_lpl <- c("Cycling T cells", "CD4 TRM-like activated memory T cells", "Activated effector Treg")

bind_rows(
  liana_cd |> filter(source %in% sirp_core_lpl, target %in% sirp_core_lpl) |> mutate(condition = "CD"),
  liana_ctrl |> filter(source %in% sirp_core_lpl, target %in% sirp_core_lpl) |> mutate(condition = "Control")
) |> 
  filter(str_detect(ligand.complex, "SIRPG"), str_detect(receptor.complex, "CD47")) |> 
  select(condition, source, target, ligand.complex, receptor.complex, magnitude_rank, specificity_rank) |> 
  arrange(condition, source, target) |> 
  print(n = Inf)

message("Inspecting MHC-II (HLA-D family) targeting T follicular helper and Treg populations...")

bind_rows(
  liana_cd  |> filter(target %in% c("T follicular helper", "Activated effector Treg", "FOXP3+ activated Treg")) |> mutate(condition = "CD"),
  liana_ctrl |> filter(target %in% c("T follicular helper", "Activated effector Treg", "FOXP3+ activated Treg")) |> mutate(condition = "Control")
) |> 
  filter(str_detect(ligand.complex, "HLA-D")) |> 
  select(condition, target, ligand.complex, receptor.complex, magnitude_rank, specificity_rank) |> 
  arrange(target, condition) |> 
  print(n = Inf)


# =========================================================
#                Save LIANA object
# =========================================================

message("Saving LIANA object for LPL compartment...")

save_rds(liana_test_lpl, "LIANA_raw_LPL.rds", dir = paths$objects_liana)

message("Saving LIANA by-sample results for LPL...")

save_rds(liana_test_lpl_disease, "LIANA_LPL_bysample_disease_aggregated.rds", dir = paths$objects_liana)


# =========================================================
#                  Save session info
# =========================================================

message("[OUTPUT] Saving session information...")

save_session_info(filename = "sessionInfo_LIANA_LPL.txt", dir = paths$logs, label = "LIANA cell-cell communication analysis - LPL compartment")

message("[OUTPUT] Session information saved to: ", paths$logs)


# =========================================================
#                  Final pipeline message
# =========================================================

message("=================================================")
message("[PIPELINE] LIANA analysis completed successfully for LPL compartment.")
message("=================================================")



