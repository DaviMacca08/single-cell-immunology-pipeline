# =========================================================
# Project      : Intestinal Mucosal Immunology - scRNA-seq Pipeline
# Compartment  : IEL (Intraepithelial Lymphocytes)
# Script       : Gene Set Enrichment Analysis (GSEA)
# Description  : Performs pre-ranked Gene Set Enrichment Analysis (GSEA)
#                on pseudobulk DESeq2 differential expression results
#                for each selected cell type using the Wald statistic
#                as the ranking metric. Functional enrichment is
#                performed against GO Biological Process (GO:BP),
#                Reactome pathways, MSigDB Hallmark (H), and
#                MSigDB C7 (Immunologic Signatures).
#
# Input        : results_list_IEL.rds
#                (DESeq2 differential expression results per cell type)
#
# Output       : GSEA result tables (.csv)
#                GSEA result objects (.rds)
#                Lollipop plots (.png)
# =========================================================


# =========================================================
#                     Functions
# =========================================================

# Results and lollipop plot (unchanged from MDS project template)

gsea_lollipop <- function(obj, db_name = "DB", color = "firebrick", title = NULL,
                           top_n = 10, fdr_cutoff = 0.05, min_leading = 5) {

  # 1. Standardize input

  if ("leadingEdge" %in% names(obj)) {

    # fgsea format
    df <- as.data.frame(obj) |>
      dplyr::rename(
        Pathway = pathway,
        p.adjust = padj,
        GeneHit = size
      ) |>
      dplyr::mutate(
        gene = leadingEdge,
        source = "fgsea"
      )

  } else {

    # clusterProfiler format
    df <- obj@result |>
      dplyr::rename(
        Pathway = Description,
        GeneHit = setSize
      ) |>
      dplyr::mutate(
        gene = strsplit(core_enrichment, "/"),
        source = "clusterProfiler"
      )

  }

  # 2. Filter significant

  df <- df |>
    dplyr::filter(!is.na(p.adjust), p.adjust < fdr_cutoff)

  if (nrow(df) == 0) {
    message(sprintf("[%s] No significant pathways", db_name))
    return(NULL)
  }

  # 3. Split NES

  up <- df |> dplyr::filter(NES > 0)
  down <- df |> dplyr::filter(NES < 0)

  # 4. Top selection pathways

  select_top <- function(x) {
    x |>
      dplyr::filter(lengths(gene) >= min_leading) |>
      dplyr::arrange(p.adjust, dplyr::desc(GeneHit)) |>
      dplyr::slice_head(n = top_n)
  }

  up <- select_top(up)
  down <- select_top(down)

  df_top <- dplyr::bind_rows(up, down)

  if (nrow(df_top) == 0) {
    message(sprintf("[%s] No pathways pass the min_leading filter (>= %d genes)", db_name, min_leading))
    return(NULL)
  }

  df_expanded <- df_top |>
    tidyr::unnest(gene)

  # 6. df for plot - one row per pathway

  df_plot <- df_top |>
    dplyr::mutate(
      FDR = -log10(p.adjust),
      Pathway = reorder(Pathway, NES)
    )

  # 7. Lollipop plot

  p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = NES, y = Pathway)) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 0, xend = NES, yend = Pathway),
      color = "grey70"
    ) +
    ggplot2::geom_point(
      ggplot2::aes(size = GeneHit, color = FDR)
    ) +
    ggplot2::scale_color_gradient(
      low = "grey60", high = color, name = "-log10(FDR)"
    ) +
    ggplot2::scale_size_continuous(name = "Leading edge size") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(size = 9),
      panel.grid.major.y = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      title = title,
      x = "Normalized Enrichment Score (NES)",
      y = ""
    )

  return(list(
    plot = p,
    data = df_expanded,
    summary = df_top,
    up = up,
    down = down
  ))
}

# Clean Hallmark names

clean_hallmark_names <- function(x, prefix = "HALLMARK_") {

  if (is.null(x) || nrow(x) == 0)
    return(x)

  x$pathway <- gsub(paste0("^", prefix), "", x$pathway)
  x$pathway <- gsub("_", " ", x$pathway)
  x$pathway <- tools::toTitleCase(tolower(x$pathway))

  x
}

# Map DESeq2 result table (gene symbol rownames) to Entrez IDs and build
# a ranked, deduplicated named vector for pre-ranked GSEA (ranked by `stat`).

build_ranked_list <- function(res_df) {

  res_df$gene <- rownames(res_df)

  mapped <- clusterProfiler::bitr(
    res_df$gene,
    fromType = "SYMBOL",
    toType   = "ENTREZID",
    OrgDb    = org.Hs.eg.db::org.Hs.eg.db
  )

  merged <- dplyr::inner_join(res_df, mapped, by = c("gene" = "SYMBOL"))

  merged <- merged |>
    dplyr::filter(!is.na(stat), is.finite(stat)) |>
    dplyr::group_by(ENTREZID) |>
    dplyr::slice_max(abs(stat), n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::arrange(dplyr::desc(stat))

  ranked <- merged$stat
  names(ranked) <- merged$ENTREZID

  ranked
  
}


# =========================================================
#                  Libraries & Setup
# =========================================================

source("Setup_Environment/00_paths.R")
source("Setup_Environment/01_environment.R")
source("Setup_Environment/02_io_helpers.R")
source("Setup_Environment/03_checks.R")
source("Setup_Environment/04_seed.R")


set_seed(1234)

message("=== Starting IEL compartment GSEA ===")


# =========================================================
#                   Load input data
# =========================================================

message("[LOAD] Reading pseudobulk DEG results (IEL)...")

deg_list <- readRDS(file.path(paths$objects_pseudobulk, "results_list_IEL.rds"))

celltype_DE <- c("IEL-Cytotoxic-TRM-like", "IEL-TH17", "IEL-Cycling-T-cells")

stopifnot(all(celltype_DE %in% names(deg_list)))


# =========================================================
#         Retrieve MSigDB gene sets once (Entrez IDs)
# =========================================================

message("[GSEA] Retrieving MSigDB Hallmark (H) and C7 (ImmunoSigDB) gene sets...")

hallmark_sets_H <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H")
hallmark_sets_C7 <- msigdbr::msigdbr(species = "Homo sapiens", collection = "C7", subcollection = "IMMUNESIGDB")

hm_list_H <- split(hallmark_sets_H$ncbi_gene, hallmark_sets_H$gs_name)
hm_list_C7 <- split(hallmark_sets_C7$ncbi_gene, hallmark_sets_C7$gs_name)


# =========================================================
#          Run GSEA per cell type (celltype_DE loop)
# =========================================================

reactome_list  <- list()
go_list        <- list()
hallmark_list  <- list()
c7_list        <- list()
ranked_list    <- list()

for (ct in celltype_DE) {

  message("=========================================================")
  message("[GSEA] Processing cell type: ", ct)
  message("=========================================================")

  ct_display <- gsub("-", " ", ct)
  ct_tag     <- gsub("-", "", ct)

  res_df <- deg_list[[ct]]
  res_df$gene <- rownames(res_df)

  ranked <- build_ranked_list(res_df)
  ranked_list[[ct]] <- ranked

  message(sprintf("[%s] Ranked gene vector: %d genes (Entrez, deduplicated)", ct, length(ranked)))

  # ---------------------------------------------------
  # Reactome
  # ---------------------------------------------------

  message("[GSEA] Running Reactome for: ", ct)

  reactome <- gsePathway(
    geneList      = ranked,
    organism      = "human",
    pvalueCutoff  = 1,
    pAdjustMethod = "BH",
    minGSSize     = 10,
    maxGSSize     = 500,
    verbose       = TRUE
  )

  reactome_list[[ct]] <- reactome

  res_reactome <- gsea_lollipop(
    reactome, db_name = "Reactome", color = "darkorange3",
    title = paste0("GSEA – Reactome pathways (", ct_display, ")")
  )

  if (!is.null(res_reactome)) {
    save_plot(res_reactome$plot, filename = paste0("GSEA_Reactome_", ct_tag, ".png"),
              dir = paths$plots_gsea_IEL, width = 11, height = 9)
  }

  # ---------------------------------------------------
  # GO Biological Process (+ semantic simplification)
  # ---------------------------------------------------

  message("[GSEA] Running GO:BP for: ", ct)

  go_bp <- gseGO(
    geneList      = ranked,
    OrgDb         = org.Hs.eg.db,
    ont           = "BP",
    keyType       = "ENTREZID",
    pvalueCutoff  = 0.05,
    pAdjustMethod = "BH",
    minGSSize     = 15,
    maxGSSize     = 300,
    eps           = 1e-4,
    verbose       = TRUE
  )

# Reduce GO term redundancy (semantic similarity clustering)
go_bp_simplified <- tryCatch(
  clusterProfiler::simplify(go_bp, cutoff = 0.5, by = "p.adjust", select_fun = min, measure = "Wang"),
  error = function(e) {
    message("[GO:BP] simplify() failed, using unfiltered result: ", conditionMessage(e))
    go_bp
  }
)

  go_list[[ct]] <- go_bp_simplified

  res_go <- gsea_lollipop(
    go_bp_simplified, db_name = "GO:BP", color = "#1B7837",
    title = paste0("GSEA – GO Biological Process (", ct_display, ")")
  )

  if (!is.null(res_go)) {
    save_plot(res_go$plot, filename = paste0("GSEA_GOBP_", ct_tag, ".png"),
              dir = paths$plots_gsea_IEL, width = 12, height = 9)
  }

  # ---------------------------------------------------
  # MSigDB Hallmark (H)
  # ---------------------------------------------------

  message("[GSEA] Running Hallmark (H) for: ", ct)

  hallmarks_H <- fgsea(
    pathways = hm_list_H, 
    stats = ranked, 
    minSize = 15,
    maxSize = 500,
    eps = 0
    )
  
  hallmarks_H_df_clean <- as.data.frame(hallmarks_H)[, !(names(hallmarks_H) %in% "leadingEdge")]
  hallmarks_H <- clean_hallmark_names(hallmarks_H)

  hallmark_list[[ct]] <- hallmarks_H

  res_h <- gsea_lollipop(
    hallmarks_H, db_name = "Hallmark", color = "firebrick",
    title = paste0("GSEA – Hallmark pathways (", ct_display, ")")
  )

  if (!is.null(res_h)) {
    save_plot(res_h$plot, filename = paste0("GSEA_Hallmark_", ct_tag, ".png"),
              dir = paths$plots_gsea_IEL, width = 11, height = 9)
  }

  # ---------------------------------------------------
  # MSigDB C7 (ImmunoSigDB - immunologic signatures)
  # ---------------------------------------------------

  message("[GSEA] Running C7 ImmunoSigDB for: ", ct)

  immunoSigDB_C7 <- fgsea(
    pathways = hm_list_C7,
    stats = ranked,
    minSize = 15,
    maxSize = 500,
    eps = 0
    )
  
  immunoSigDB_C7_df_clean <- as.data.frame(immunoSigDB_C7)[, !(names(immunoSigDB_C7) %in% "leadingEdge")]
  immunoSigDB_C7 <- clean_hallmark_names(immunoSigDB_C7, prefix = "")

  c7_list[[ct]] <- immunoSigDB_C7

  res_c7 <- gsea_lollipop(
    immunoSigDB_C7, db_name = "C7", color = "darkgreen",
    title = paste0("GSEA – ImmunoSigDB (C7) signatures (", ct_display, ")"),
    min_leading = 10
  )

  if (!is.null(res_c7)) {
    save_plot(res_c7$plot, filename = paste0("GSEA_C7_", ct_tag, ".png"),
              dir = paths$plots_gsea_IEL, width = 15, height = 9)
  }

  # ---------------------------------------------------
  # Save per-cell-type result tables (GSEA)
  # ---------------------------------------------------

  save_csv(hallmarks_H_df_clean, filename = paste0("GSEA_hallmark_H_", ct_tag, "_IEL.csv"), dir = paths$tables_gsea_IEL)
  save_csv(immunoSigDB_C7_df_clean, filename = paste0("GSEA_C7_", ct_tag, "_IEL.csv"), dir = paths$tables_gsea_IEL)
  save_csv(as.data.frame(reactome), filename = paste0("GSEA_reactome_", ct_tag, "_IEL.csv"), dir = paths$tables_gsea_IEL)
  save_csv(as.data.frame(go_bp_simplified), filename = paste0("GSEA_GOBP_", ct_tag, "_IEL.csv"), dir = paths$tables_gsea_IEL)

  message("[GSEA] Completed cell type: ", ct)
  
}


# =========================================================
#                Save full result objects
# =========================================================

message("[OUTPUT] Saving GSEA result objects (all cell types, IEL)...")

save_rds(ranked_list,   filename = "GSEA_ranked_lists_IEL.rds",   dir = paths$objects_gsea)
save_rds(reactome_list, filename = "GSEA_reactome_list_IEL.rds",  dir = paths$objects_gsea)
save_rds(go_list,       filename = "GSEA_GOBP_list_IEL.rds",      dir = paths$objects_gsea)
save_rds(hallmark_list, filename = "GSEA_hallmark_list_IEL.rds",  dir = paths$objects_gsea)
save_rds(c7_list,       filename = "GSEA_C7_list_IEL.rds",        dir = paths$objects_gsea)


# =========================================================
#                  Save session info
# =========================================================

message("[OUTPUT] Saving session information...")

save_session_info(filename = "sessionInfo_GSEA_IEL.txt", dir = paths$logs, label = "Gene Set Enrichment Analysis (GSEA) - IEL compartment")

message("[OUTPUT] Session information saved to: ", paths$logs)


# =========================================================
#                  Final pipeline message
# =========================================================

message("=================================================")
message("[PIPELINE] Gene Set Enrichment Analysis completed successfully for IEL compartment (CD vs Control, per cell type).")
message("=================================================")
