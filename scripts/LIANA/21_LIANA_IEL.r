# =========================================================
# LIANA cell-cell communication analysis - IEL
# =========================================================
# Input:
# - Seurat object of IEL compartment with consensus annotation
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

message("Starting LIANA analysis for IEL compartment...")


# =========================================================
#                   Load input data
# =========================================================

message("Loading IEL Seurat object...")

iel_cells <- readRDS(file.path(paths$objects_subsetting, "srt_IEL_compartment_celltype.rds"))

message("IEL object loaded: ", nrow(iel_cells), " genes × ", ncol(iel_cells)," cells")
message("Loading global Seurat object...")

full_obj <- readRDS(file.path(paths$objects_global, "srt_obj_merge_harmony.rds"))

message("Global object loaded: ", nrow(full_obj), " genes × ", ncol(full_obj)," cells")


# =========================================================
# Validate input data
# =========================================================

message("Validating IEL normalized expression layer...")

stopifnot(dim(LayerData(iel_cells,
                        assay = "RNA",
                        layer = "data")) == dim(iel_cells))

message("  ✓ IEL data layer dimensions validated")

stopifnot("celltype_IEL" %in% colnames(iel_cells@meta.data))

message("  ✓ IEL cell-type annotation found")

stopifnot("disease" %in% colnames(iel_cells@meta.data))

message("  ✓ Disease metadata found")


# =========================================================
# Prepare expression matrices
# =========================================================

message("Preparing expression matrices for LIANA...")

message("Joining raw count layers in global object...")

full_obj <- JoinLayers(full_obj, assay = "RNA")

message("  ✓ Raw count layers joined")

# Extract normalized expression from IEL object

message("Extracting normalized expression from IEL data layer...")

logcounts_mat <- LayerData(iel_cells, assay = "RNA", layer = "data")

message("  ✓ Normalized expression extracted: ", nrow(logcounts_mat), " genes × ", ncol(logcounts_mat)," cells")

# Extract raw counts from global object

message("Extracting raw counts from global object...")

counts_full <- LayerData(full_obj, assay = "RNA", layer = "counts")

message("  ✓ Raw counts extracted: ", nrow(counts_full), " genes × ", ncol(counts_full), " cells")

# Subset raw counts to IEL cells

message("Subsetting raw counts to IEL cells...")

counts_mat <- counts_full[, colnames(iel_cells), drop = FALSE]

message("  ✓ IEL raw counts extracted: ", nrow(counts_mat), " genes × ", ncol(counts_mat)," cells")


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

message("Preparing IEL cell metadata...")

meta_iel <- iel_cells@meta.data

message("  ✓ Metadata extracted for ", nrow(meta_iel), " cells")

message("  ✓ Cell types identified: ", length(unique(iel_cells$celltype_IEL)))

message("  ✓ Disease groups: ",paste(unique(iel_cells$disease), collapse = ", "))


# =========================================================
# Build SingleCellExperiment object
# =========================================================

message("Building SingleCellExperiment object...")

iel_sce <- SingleCellExperiment(assays = list(logcounts = logcounts_mat, counts = counts_mat), colData = meta_iel)

# Set cell identities

colLabels(iel_sce) <- iel_cells$celltype_IEL

message("  ✓ SingleCellExperiment created: ", nrow(iel_sce), " genes × ", ncol(iel_sce)," cells")


# =========================================================
# Validate SingleCellExperiment object
# =========================================================

message("Validating SingleCellExperiment object...")

stopifnot(identical(assayNames(iel_sce),c("logcounts", "counts")))

message("  ✓ counts and logcounts assays detected")

stopifnot(ncol(iel_sce) == ncol(iel_cells))

message("  ✓ Cell number validated")

stopifnot(nrow(iel_sce) == nrow(iel_cells))

message("  ✓ Gene number validated")

stopifnot(length(colLabels(iel_sce)) == ncol(iel_sce))

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

message("Running LIANA for IEL compartment...")

liana_test_iel <- liana_wrap(
  iel_sce,
  idents_col = "celltype_IEL",
  method = c("natmi", "connectome", "logfc", "sca", "cellphonedb"),
  resource = "Consensus")

message("LIANA analysis completed successfully.")

# Aggregate results 

liana_test_iel_aggregate <- rank_aggregate(liana_test_iel)

message("LIANA aggregated results dimensions: ", nrow(liana_test_iel_aggregate), " rows × ", ncol(liana_test_iel_aggregate)," columns")


# =========================================================
# Inspect aggregated results
# =========================================================

message("Inspecting aggregated LIANA results...")

print(colnames(liana_test_iel_aggregate))

message(
  "Aggregated LIANA results: ",
  nrow(liana_test_iel_aggregate),
  " rows × ",
  ncol(liana_test_iel_aggregate),
  " columns"
)


# =========================================================
# Global LIANA heatmap and dotplot
# =========================================================

message("Generating global IEL communication frequency heatmap...")

liana_trunc <- liana_test_iel_aggregate |> filter(specificity_rank <= 0.05)


open_png(filename = "LIANA_global_specificity_Heatmap_IEL.png", dir = paths$plots_liana_iel,
         width = 2000, height = 1900)

heat_freq(liana_trunc,
          row_title = "Sender (Cell type)",
          column_title = "Receiver (Cell type)",
          name = "Number of LR interactions") 

close_png()

message("Generating specificity-focused LIANA dotplot...")

# 01-Epithelial to CD39 TRM
open_png(filename = "LIANA_IEL_epithelial_to_cd39_trm_dotplot.png", dir = paths$plots_liana_iel,
         width = 1500, height = 1500)

liana_test_iel_aggregate |>
  arrange(specificity_rank) |> 
  liana_dotplot(source_groups = "Epithelial cells / enterocytes",
                target_groups = "CD39+ tissue-resident CD8 T cells",
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
    title = "IEL communication: Epithelial cells → CD39+ TRM",
    subtitle = "Ranked by interaction specificity | Size: specificity | Colour: magnitude",
    x = NULL,
    y = "Ligand → Receptor"
  )

close_png()

# 02-CD39+ TRM -> CST3+ LYZ+ macrophages
open_png(filename = "LIANA_IEL_cd39_trm_to_cst3_lyz_macrophages_dotplot.png", dir = paths$plots_liana_iel,
         width = 1500, height = 1500)

liana_test_iel_aggregate |>
  arrange(magnitude_rank) |> 
  liana_dotplot(source_groups = "CD39+ tissue-resident CD8 T cells",
                target_groups = "CST3+ LYZ+ macrophages",
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
    title = "IEL communication: CD39+ TRM → CST3+ LYZ+ macrophages",
    subtitle = "Ranked by expression magnitude | Size: specificity | Colour: magnitude",
    x = NULL,
    y = "Ligand → Receptor"
  )

close_png()

# 03- T cell targets
open_png(filename = "LIANA_IEL_cytotoxic_tcell_targets_dotplot.png", dir = paths$plots_liana_iel,
         width = 5000, height = 5000)

liana_test_iel_aggregate |> 
  arrange(magnitude_rank) |> 
  liana_dotplot(
    target_groups = c("CD39+ tissue-resident CD8 T cells",
                      "GZMK+ effector memory CD8 T cells",
                      "Terminal effector CD8 T cells",
                      "NK-like cytotoxic T cells",
                      "Cytotoxic γδ T cells"),
    ntop = 30,
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
    title = "IEL communication: interactions targeting cytotoxic T-cell populations",
    subtitle = "Ranked by expression magnitude | Size: specificity | Colour: magnitude",
    x = NULL,
    y = "Ligand → Receptor"
  ) +
  scale_x_discrete(position = "bottom") +
  theme(
    strip.text.x = element_text(
      angle = 90, hjust = 0.5, vjust = 0.5, size = 11),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10)
  )

close_png()      

# 04-TH17 effector" -> "FOXP3+ IL2RA+ Treg
open_png(filename = "LIANA_IEL_th17_treg_targeted_communication_dotplot.png", dir = paths$plots_liana_iel,
         width = 5000, height = 5000)

liana_test_iel_aggregate |> 
  arrange(specificity_rank) |> 
  liana_dotplot(
    target_groups = c("TH17 effector", "FOXP3+ IL2RA+ Treg"),
    ntop = 30, magnitude = "magnitude_rank", 
    specificity = "specificity_rank",
    invert_magnitude = TRUE,
    invert_specificity = TRUE,
    y.label = "Ligand → Receptor",
    size.label = "Interaction specificity",
    colour.label = "Expression magnitude",
    size_range = c(3, 9)
  ) +
  labs(
    title = "IEL communication: interactions targeting TH17 and Treg cells",
    subtitle = "Ranked by interaction specificity | Size: specificity | Colour: magnitude",
    x = NULL,
    y = "Ligand → Receptor"
  ) +
  scale_x_discrete(position = "bottom") +
  theme(
    strip.text.x = element_text(
      angle = 90, hjust = 0.5, vjust = 0.5, size = 11),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10)
  )

close_png()


# ============================================================
# Global IEL communication - chord plot
# ============================================================

message("Preparing LIANA global IEL chord plot...")

# Filter interactions based on specificity
liana_trunc <- liana_test_iel_aggregate |> 
  filter(specificity_rank <= 0.05, magnitude_rank <= 0.25)

# Cell types selected for visualization

celltypes_of_interest_chord1 <- c(
  "TH17-like CD4 T cells",
  "GZMK+ effector memory CD8 T cells",
  "Effector T cell mixed state"
)

celltypes_of_interest_chord2 <- c(
  "FOXP3+ IL2RA+ Treg",
  "CST3+ LYZ+ macrophages",
  "CD39+ tissue-resident CD8 T cells",
  "Epithelial cells / enterocytes"
)

# ------------------------------------------------------------
# Chord plot 1
# ------------------------------------------------------------

message("Generating chord diagram 1...")

open_png(filename = "LIANA_IEL_global_chord1.png", dir = paths$plots_liana_iel,
         width = 11000, height = 8000
)

chord_freq(
  liana_trunc,
  source_groups = celltypes_of_interest_chord1,
  target_groups = celltypes_of_interest_chord1,
  cex = 2.5,
  facing = "clockwise",
  adj = c(0, 0.5)
)

close_png()

# ------------------------------------------------------------
# Chord plot 2
# ------------------------------------------------------------

message("Generating chord diagram 2...")

open_png(filename = "LIANA_IEL_global_chord2.png", dir = paths$plots_liana_iel,
         width = 11000, height = 8000
)

chord_freq(
  liana_trunc,
  source_groups = celltypes_of_interest_chord2,
  target_groups = celltypes_of_interest_chord2,
  cex = 2.5,
  facing = "clockwise",
  adj = c(0, 0.5)
)

close_png()

message("Individual chord diagrams generated.")


# ============================================================
# CD vs Control comparison - IEL compartment
# ============================================================

message("Starting LIANA comparison between CD and Control in the IEL compartment...")

# Compare CD vs CTRL in IEL compartment
liana_test_iel_disease <- liana_bysample(
  sce = iel_sce,
  sample_col = "disease",
  idents_col = "celltype_IEL",
  aggregate_how = "both",
  method = c("natmi", "connectome", "logfc", "sca", "cellphonedb"),
  resource = "Consensus",
  inplace = FALSE,
  verbose = TRUE
)

message("LIANA CD vs Control analysis completed.")

# Check available disease groups
message("Available disease groups: ", paste(names(liana_test_iel_disease), collapse = ", "))

if (!all(c("Control", "CD") %in% names(liana_test_iel_disease))) {
  stop("Both Control and CD groups are required for the comparison.")
}

message("✓ Both Control and CD groups are available.")

names(liana_test_iel_disease)

# ------------------------------------------------------------
# Separate Control and CD results
# ------------------------------------------------------------

message("Extracting Control and CD LIANA results...")

liana_ctrl <- liana_test_iel_disease[["Control"]]

liana_cd <- liana_test_iel_disease[["CD"]]

message("Control results extracted.")
message("CD results extracted.")

# Summary of available results
message("Control result dimensions: ", paste(dim(liana_ctrl), collapse = " x "))
message("CD result dimensions: ", paste(dim(liana_cd), collapse = " x "))


# ============================================================
# CD vs Control - adequately represented IEL populations
# ============================================================

message("Generating CD vs Control LIANA dotplots for balance populations...")

tier_iel <- table(iel_cells$celltype_IEL, iel_cells$disease) |> 
  as.data.frame() |> 
  pivot_wider(names_from = Var2, values_from = Freq) |> 
  rename(celltype = Var1) |> 
  mutate(tier = case_when(
    CD == 0 | Control == 0 ~ "condition_exclusive",
    pmin(CD, Control) < 30 ~ "imbalanced",   
    TRUE ~ "fully_comparable"
  ))

fully_comparable_types <- tier_iel |> filter(tier == "fully_comparable") |>  pull(celltype)


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

open_png(filename = "LIANA_IEL_CD_vs_Control_balanced_populations_dotplot.png", dir = paths$plots_liana_iel,
         width = 6000, height = 3000)

p_mhc2_cd + p_mhc2_ctrl

close_png()

message("CD vs Control dotplot saved.")


# ============================================================
# Targeted inspection of biologically relevant interactions
# ============================================================

message("Inspecting HLA-D interactions targeting TH17 and Treg cells...")

bind_rows(
  liana_cd |> filter(target %in% c("TH17 effector","FOXP3+ IL2RA+ Treg")) |> 
    mutate(condition="CD"),
  liana_ctrl |>  filter(target %in% c("TH17 effector","FOXP3+ IL2RA+ Treg")) |> 
    mutate(condition="Control")
) |> 
  filter(str_detect(ligand.complex, "HLA-D")) |> 
  select(condition, target, ligand.complex, receptor.complex, magnitude_rank, specificity_rank) %>%
  arrange(target, condition) |> 
  print(n = Inf)

message("Inspecting epithelial CDH1 → ITGAE interactions...")

bind_rows(
  liana_cd   |> filter(source=="Epithelial cells / enterocytes", target=="CD39+ tissue-resident CD8 T cells") |> 
    mutate(condition="CD"),
  liana_ctrl |>  filter(source=="Epithelial cells / enterocytes", target=="CD39+ tissue-resident CD8 T cells") |> mutate(condition="Control")
) |> 
  filter(str_detect(ligand.complex,"CDH1"), str_detect(receptor.complex,"ITGAE")) |> 
  select(condition, magnitude_rank, specificity_rank)

message("Inspecting TNFSF14 → LTBR interactions between CD39+ TRM and macrophages...")

liana_cd |> 
  filter(source=="CD39+ tissue-resident CD8 T cells", target=="CST3+ LYZ+ macrophages",
         str_detect(ligand.complex,"TNFSF14"), str_detect(receptor.complex,"LTBR")) |> 
  select(magnitude_rank, specificity_rank)


# =========================================================
#                Save LIANA object
# =========================================================

message("Saving LIANA object for IEL compartment...")

save_rds(liana_test_iel, "LIANA_raw_IEL.rds", dir = paths$objects_liana)

message("Saving LIANA by-sample results for IEL...")

save_rds(liana_test_iel_disease, "LIANA_IEL_bysample_disease_aggregated.rds", dir = paths$objects_liana)


# =========================================================
#                  Save session info
# =========================================================

message("[OUTPUT] Saving session information...")

save_session_info(filename = "sessionInfo_LIANA_IEL.txt", dir = paths$logs, label = "LIANA cell-cell communication analysis - IEL compartment")

message("[OUTPUT] Session information saved to: ", paths$logs)


# =========================================================
#                  Final pipeline message
# =========================================================

message("=================================================")
message("[PIPELINE] LIANA analysis completed successfully for IEL compartment.")
message("=================================================")


