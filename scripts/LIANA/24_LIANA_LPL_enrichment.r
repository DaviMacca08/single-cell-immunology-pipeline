# ============================================================
# LIANA: Functional enrichment of ligand-receptor hits
#
# LPL compartment
# CD vs Control comparison
# Ligand- and receptor-side enrichment
# ============================================================

# =========================================================
#                  Libraries & Setup
# =========================================================

source("Setup_Environment/00_paths.R")
source("Setup_Environment/01_environment.R")
source("Setup_Environment/02_io_helpers.R")
source("Setup_Environment/03_checks.R")
source("Setup_Environment/04_seed.R")

set_seed(1234)

message("Starting LIANA functional enrichment analysis for LPL compartment...")

# =========================================================
#                   Load input data
# =========================================================

message("Loading LIANA LPL disease-level results...")

liana_test_lpl_disease <- readRDS(file.path(paths$objects_liana, "LIANA_LPL_bysample_disease_aggregated.rds"))
lpl_cells <- readRDS(file.path(paths$objects_subsetting, "srt_LPL_compartment_celltype.rds"))

message("Loading LPL cell Seurat Object...")

liana_ctrl <- liana_test_lpl_disease[["Control"]]
liana_cd <- liana_test_lpl_disease[["CD"]]


# =========================================================
#                   Functions
# =========================================================

# 01-ORA

run_liana_ora <- function(liana_res, liana_universe, background_total, celltypes) {
  
  # --------------------------------------------------------------------------
  # 1. Filter LIANA interactions
  # --------------------------------------------------------------------------
  
  # Retain highly specific interactions and restrict the analysis to the selected cell types.
  liana_form <- liana_res |>
    select(source, target, ligand.complex, receptor.complex, specificity_rank) |>
    filter(
      specificity_rank <= 0.05,
      source %in% celltypes,
      target %in% celltypes)
  
  # --------------------------------------------------------------------------
  # 2. Extract unique ligands by source cell type
  # --------------------------------------------------------------------------
  
  # Each ligand is counted once per source cell type.
  liana_ligands <- liana_form |>
    select(source, ligand = ligand.complex) |> 
    distinct() |>
    group_by(source) |> 
    mutate(distinct_hits = n()) |>
    ungroup()
  
  # --------------------------------------------------------------------------
  # 3. Extract unique receptors by target cell type
  # --------------------------------------------------------------------------
  
  # Each receptor is counted once per target cell type.
  liana_receptors <- liana_form |>
    select(target, receptor = receptor.complex) |>
    distinct() |>
    group_by(target) |>
    mutate(distinct_hits = n()) |>
    ungroup()
  
  
  # --------------------------------------------------------------------------
  # 4. Ligand enrichment
  # --------------------------------------------------------------------------
  
  # Map LIANA ligands to the gene-set annotation table 
  ligand_res <- liana_universe |>
    rename(ligand = gene_symbol, geneset = gs_name) |>
    inner_join(liana_ligands, by = "ligand") |>
    group_by(source, geneset) |> 
    mutate(ligands_in_gs = n()) |>
    ungroup() |>
    rowwise() |>
    mutate(pval = phyper(ligands_in_gs - 1, geneset_n, background_total - geneset_n,
                         distinct_hits, lower.tail = FALSE)) |>
    ungroup() |>
    select(source, geneset, pval, ligands_in_gs, distinct_hits) |> distinct() |> # Keep the variables required for downstream interpretation.
    mutate(adj_pval = p.adjust(pval, method = "fdr"),
           GeneRatio = ligands_in_gs / distinct_hits) |>
    arrange(adj_pval)
  
  # --------------------------------------------------------------------------
  # 5. Receptor enrichment
  # --------------------------------------------------------------------------
  
  # Map LIANA receptors to the gene-set annotation table 
  receptor_res <- liana_universe |>
    rename(receptor = gene_symbol, geneset = gs_name) |>
    inner_join(liana_receptors, by = "receptor") |>
    group_by(target, geneset) |> mutate(receptors_in_gs = n()) |> ungroup() |>
    rowwise() |>
    mutate(pval = phyper(receptors_in_gs - 1, geneset_n, background_total - geneset_n,
                         distinct_hits, lower.tail = FALSE)) |>
    ungroup() |>
    select(target, geneset, pval, receptors_in_gs, distinct_hits) |> distinct() |> # Keep the variables required for downstream interpretation.
    mutate(adj_pval = p.adjust(pval, method = "fdr"),
           GeneRatio = receptors_in_gs / distinct_hits) |>
    arrange(adj_pval)
  
  # --------------------------------------------------------------------------
  # 6. Return results
  # --------------------------------------------------------------------------
  list(ligand = ligand_res, receptor = receptor_res)
  
}

# 02-Lollypop plot
liana_ora_lollipop_compare <- function(df, entity_col, db_name = "Hallmark", color = "firebrick",
                                       title = NULL, top_n = 10, fdr_cutoff = 0.05, min_hits = 2) {
  
  # --------------------------------------------------------------------------
  # 1. Filter significant enrichment results
  # --------------------------------------------------------------------------
  
  df_sig <- df |>
    dplyr::filter(!is.na(adj_pval), adj_pval < fdr_cutoff, hits >= min_hits)
  
  # Stop if no significant pathways are available.
  if (nrow(df_sig) == 0) {
    message(sprintf("[%s] No significant pathways in either condition", db_name))
    return(NULL)
  }
  
  # --------------------------------------------------------------------------
  # 2. Select top pathways separately for each cell type and condition
  # --------------------------------------------------------------------------
  
  df_top <- df_sig |>
    dplyr::group_by(.data[[entity_col]], condition) |>
    dplyr::arrange(adj_pval, dplyr::desc(hits)) |>
    dplyr::slice_head(n = top_n) |>
    dplyr::ungroup()
  
  # --------------------------------------------------------------------------
  # 3. Define a common pathway ordering
  # --------------------------------------------------------------------------
  
  geneset_order <- df_top |>
    dplyr::group_by(geneset) |>
    dplyr::summarise(mean_ratio = mean(GeneRatio), .groups = "drop") |>
    dplyr::arrange(mean_ratio) |>
    dplyr::pull(geneset)
  
  
  # --------------------------------------------------------------------------
  # 4. Prepare data for plotting
  # --------------------------------------------------------------------------
  
  df_plot <- df_top |>
    dplyr::mutate(
      FDR = -log10(adj_pval),
      geneset = factor(geneset, levels = geneset_order)
    )
  
  # --------------------------------------------------------------------------
  # 5. Generate comparative lollipop plot
  # --------------------------------------------------------------------------
  
  p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = GeneRatio, y = geneset)) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = GeneRatio, yend = geneset), color = "grey70") +
    ggplot2::geom_point(ggplot2::aes(size = hits, color = FDR)) +
    ggplot2::scale_color_gradient(low = "grey60", high = color, name = "-log10(FDR)") +
    ggplot2::scale_size_continuous(name = "Gene hits") +
    ggplot2::facet_grid(
      rows = ggplot2::vars(condition),
      cols = ggplot2::vars(.data[[entity_col]])
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 8),
      panel.grid.major.y = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", size = 9)
    ) +
    ggplot2::labs(title = title, x = "Gene Ratio", y = "")
  
  # --------------------------------------------------------------------------
  # 6. Return plot and plotting data
  # --------------------------------------------------------------------------
  
  return(list(plot = p, data = df_top))
}

# 03-Clean Hallmark names

clean_hallmark_names <- function(x, prefix = "HALLMARK_") {
  
  if (is.null(x) || nrow(x) == 0)
    return(x)
  
  x$geneset <- gsub(paste0("^", prefix), "", x$geneset)
  x$geneset <- gsub("_", " ", x$geneset)
  x$geneset <- tools::toTitleCase(tolower(x$geneset))
  
  x
}




# ============================================================
# Define adequately represented LPL populations
# ============================================================

message("Assessing LPL population representation across disease conditions...")

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

message("LPL populations retained for CD vs Control comparison: ", length(fully_comparable_types))

message("Excluded populations: ", sum(tier_lpl$tier != "fully_comparable"))


# =========================================================
# Retrieve MSigDB Hallmark gene sets
# =========================================================

message("Retrieving MSigDB Hallmark gene sets...")

hallmark <- msigdbr(species = "Homo sapiens", category = "H") |>
  select(gene_symbol, gs_name)

message("Hallmark gene sets retrieved: ", n_distinct(hallmark$gs_name))


# =========================================================
# Define LIANA Consensus ligand-receptor background universe
# =========================================================

message("Defining LIANA Consensus ligand-receptor background universe...")

sel_res <- select_resource("Consensus")[[1]]

lr_genes_universe <- unique(c(sel_res$source_genesymbol, sel_res$target_genesymbol))

liana_universe <- hallmark |>
  filter(gene_symbol %in% lr_genes_universe) |>
  group_by(gs_name) |> mutate(geneset_n = n()) |> ungroup()

background_total <- length(lr_genes_universe)

message("Unique ligand/receptor genes in LIANA Consensus resource: ", background_total)


# ============================================================
# Run ligand and receptor ORA
# ============================================================

message("Running Hallmark over-representation analysis for LIANA hits...")

ora_cd   <- run_liana_ora(liana_cd,   liana_universe, background_total, fully_comparable_types)
ora_ctrl <- run_liana_ora(liana_ctrl, liana_universe, background_total, fully_comparable_types)

message("ORA completed for CD and Control conditions.")

ora_ligand_compare <- bind_rows(
  ora_cd$ligand   |> mutate(condition = "CD",      hits = ligands_in_gs),
  ora_ctrl$ligand |> mutate(condition = "Control", hits = ligands_in_gs)
)

ora_receptor_compare <- bind_rows(
  ora_cd$receptor   |> mutate(condition = "CD",      hits = receptors_in_gs),
  ora_ctrl$receptor |> mutate(condition = "Control", hits = receptors_in_gs)
)

message("Ligand enrichment results: ", nrow(ora_ligand_compare), " pathway-cell type-condition combinations.")

message("Receptor enrichment results: ", nrow(ora_receptor_compare)," pathway-cell type-condition combinations.")


# ============================================================
# Prepare pathway labels for visualization
# ============================================================

ora_ligand_compare <- clean_hallmark_names(ora_ligand_compare)

ora_receptor_compare <- clean_hallmark_names(ora_receptor_compare)


# ============================================================
# Generate comparative ligand and receptor enrichment plots
# ============================================================

message("Generating comparative ligand enrichment plot...")

res_ligand_compare <- liana_ora_lollipop_compare(
  ora_ligand_compare, entity_col = "source",
  title = "LPL — Ligand-side ORA, CD vs Control (Hallmark)"
)

message("Generating comparative receptor enrichment plot...")

res_receptor_compare <- liana_ora_lollipop_compare(
  ora_receptor_compare, entity_col = "target",color = "blue",
  title = "LPL — Receptor-side ORA, CD vs Control (Hallmark)"
)

message("Saving LIANA enrichment plots...")

open_png(filename = "LIANA_LPL_ORA_ligand_CDvsControl.png", dir = paths$plots_liana_lpl,
         width = 3500, height = 2200)

print(res_ligand_compare$plot)

close_png()

open_png(filename = "LIANA_LPL_ORA_receptor_CDvsControl.png", dir = paths$plots_liana_lpl,
         width = 3500, height = 2200)

print(res_receptor_compare$plot)

close_png()


# =========================================================
#               Save enrichment results
# =========================================================

message("Saving LIANA enrichment results...")

save_rds(ora_cd, "LIANA_LPL_ORA_CD.rds", dir = paths$objects_liana)
save_rds(ora_ctrl, "LIANA_LPL_ORA_Control.rds", dir = paths$objects_liana)

save_rds(ora_ligand_compare, "LIANA_LPL_ORA_ligand_CDvsControl.rds", dir = paths$objects_liana)

save_rds(ora_receptor_compare, "LIANA_LPL_ORA_receptor_CDvsControl.rds", dir = paths$objects_liana)

message("LIANA enrichment results saved to: ", paths$objects_liana)
message("Saving LIANA by-sample results for LPL...")


# =========================================================
#                  Save session info
# =========================================================

message("[OUTPUT] Saving session information...")

save_session_info(filename = "sessionInfo_LIANA_LPL_ORA.txt", dir = paths$logs, label = "LIANA functional enrichment analysis - LPL compartment")

message("[OUTPUT] Session information saved to: ", paths$logs)


# =========================================================
#                  Final pipeline message
# =========================================================


message("=================================================")
message("[PIPELINE] LIANA LPL functional enrichment completed successfully.")
message("[PIPELINE] Ligand-side and receptor-side Hallmark ORA completed.")
message("[PIPELINE] CD vs Control comparison completed.")
message("[OUTPUT] Enrichment results and plots saved.")
message("=================================================")


