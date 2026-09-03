# =========================================================
# Project      : scRNA-Seq project - Bulk RNA-seq Analysis Pipeline
# Dataset      : GSE193677 (GPL16791 Illumina HiSeq 2500)
# Samples      : Ileal biopsy samples from CD patients and healthy controls
# Cohort       : Primary cohort (Ileum, NonI)
# Comparison   : Crohn's disease (CD) vs Control
# Script       : Gene Set Enrichment Analysis (GSEA)
# Description  : Pre-ranked GSEA using gene-level statistics (moderated t)
#                from the limma-voom differential expression analysis.
#                Functional enrichment against GO (BP/MF/CC), KEGG, Reactome,
#                and MSigDB Hallmark (H) and ImmunoSigDB (C7) gene sets.
# =========================================================


# =========================================================
#                     Functions
# =========================================================

# 01-Lollipop function

gsea_lollipop <- function(obj, db_name = "DB", color = "firebrick", title = NULL,
                          top_n = 10, fdr_cutoff = 0.05, min_leading = 10) {
  
  # 1. Standardize input
  
  if ("leadingEdge" %in% names(obj)) {
    
    # fgsea format (raw Bioconductor fgsea output - used for Hallmark H/C7)
    df <- as.data.frame(obj) |>
      dplyr::rename(Pathway = pathway, p.adjust = padj, GeneHit = size) |>
      dplyr::mutate(gene = leadingEdge, source = "fgsea")
    
  } else if (is.data.frame(obj)) {
    
    # Already-combined data.frame (e.g. GO multi-ontology after per-ontology simplify())
    df <- obj |>
      dplyr::rename(Pathway = Description, GeneHit = setSize) |>
      dplyr::mutate(gene = strsplit(core_enrichment, "/"), source = "clusterProfiler_df")
    
  } else {
    
    # clusterProfiler S4 object (single ontology/database - KEGG, Reactome)
    df <- obj@result |>
      dplyr::rename(Pathway = Description, GeneHit = setSize) |>
      dplyr::mutate(gene = strsplit(core_enrichment, "/"), source = "clusterProfiler")
  }
  
  # 2. Filter significant terms
  
  df <- df |> dplyr::filter(!is.na(p.adjust), p.adjust < fdr_cutoff)
  
  if (nrow(df) == 0) {
    message("[LOLLIPOP] No significant pathways.")
    return(NULL)
  }
  
  has_ontology <- "ONTOLOGY" %in% colnames(df)
  
  # 3. Split by direction (NES sign)
  
  up <- df |> dplyr::filter(NES > 0)
  down <- df |> dplyr::filter(NES < 0)
  
  # 4. Top-N selection - per ontology if present, otherwise global
  
  select_top <- function(x) {
    
    x <- x |> dplyr::filter(lengths(gene) >= min_leading)
    
    if (has_ontology) {
      x |>
        dplyr::group_by(ONTOLOGY) |>
        dplyr::arrange(p.adjust, dplyr::desc(GeneHit), .by_group = TRUE) |>
        dplyr::slice_head(n = top_n) |>
        dplyr::ungroup()
    } else {
      x |>
        dplyr::arrange(p.adjust, dplyr::desc(GeneHit)) |>
        dplyr::slice_head(n = top_n)
    }
  }
  
  up <- select_top(up)
  down <- select_top(down)
  
  df_top <- dplyr::bind_rows(up, down)
  
  if (nrow(df_top) == 0) {
    message(sprintf("[LOLLIPOP] No pathways pass the min_leading filter (>= %d genes).", min_leading))
    return(NULL)
  }
  
  df_expanded <- df_top |> tidyr::unnest(gene)
  
  df_plot <- df_top |>
    dplyr::mutate(FDR = -log10(p.adjust), Pathway = reorder(Pathway, NES))
  
  # 5. Lollipop plot
  
  p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = NES, y = Pathway)) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = NES, yend = Pathway), color = "grey70") +
    ggplot2::geom_point(ggplot2::aes(size = GeneHit, color = FDR)) +
    ggplot2::scale_color_gradient(low = "grey60", high = color, name = "-log10(FDR)") +
    ggplot2::scale_size_continuous(name = "Leading edge size") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 9), panel.grid.major.y = ggplot2::element_blank()) +
    ggplot2::labs(title = title, x = "Normalized Enrichment Score (NES)", y = "")
  
  if (has_ontology) {
    p <- p + ggplot2::facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free_y")
  }
  
  return(list(plot = p, data = df_expanded, summary = df_top, up = up, down = down))
}

# 02-Treeplot

add_redundancy_treeplot <- function(enrich_obj, title, out_dir, filename, showCategory = 30, nCluster = 5) {
  
  if (is.null(enrich_obj) || nrow(as.data.frame(enrich_obj)) == 0) {
    message("[REDUNDANCY] Empty enrichment result, skipping treeplot.")
    return(NULL)
  }
  
  n_terms <- nrow(as.data.frame(enrich_obj))
  
  if (n_terms < 3) {
    message(sprintf(
      "[REDUNDANCY] Only %d significant term(s) - treeplot skipped (need >= 3 for meaningful clustering).",
      n_terms
    ))
    return(NULL)
  }
  
  enrich_obj_sim <- enrichplot::pairwise_termsim(enrich_obj)
  
  # nCluster cannot exceed the number of terms actually shown
  n_shown <- min(showCategory, n_terms)
  nCluster_adj <- min(nCluster, n_shown - 1)
  
  message(sprintf(
    "[REDUNDANCY] %d term(s) available, showing %d, clustering into %d group(s).",
    n_terms, n_shown, nCluster_adj
  ))
  
  tp <- enrichplot::treeplot(enrich_obj_sim, showCategory = showCategory, nCluster = nCluster_adj) +
    ggtitle(title)
  
  save_plot(tp, filename = filename, dir = out_dir, width = 12, height = 9)
  
  enrich_obj_sim
}

# 03-Clean MSigDB Hallmark pathway names for plotting 

clean_hallmark_names <- function(x) {
  
  if (is.null(x) || nrow(x) == 0)
    return(x)
  
  x$pathway <- gsub("^HALLMARK_", "", x$pathway)
  x$pathway <- gsub("_", " ", x$pathway)
  x$pathway <- tools::toTitleCase(tolower(x$pathway))
  
  x
}

# 04-Run gseGO separately per ontology (BP/MF/CC) and apply simplify() to each.

run_go_gsea_simplified <- function(gene_ranked, ontologies = c("BP", "MF", "CC"), simplify_cutoff = 0.7) {
  
  results_list <- list()
  
  for (ont in ontologies) {
    
    message(sprintf("[GSEA] Running gseGO for ontology = %s...", ont))
    
    gse_ont <- gseGO(
      geneList = gene_ranked,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = ont,
      minGSSize = 10,
      maxGSSize = 500,
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,   # keep all tested terms, avoid silent internal filtering
      eps = 1e-4,         # Set a lower bound for extremely small P-value estimates
      verbose = TRUE
    )
    
    if (is.null(gse_ont) || nrow(as.data.frame(gse_ont)) == 0) {
      message(sprintf("[GSEA] No enriched terms for ontology = %s.", ont))
      next
    }
    
    message(sprintf("[GSEA] Simplifying redundant %s terms (cutoff = %.1f)...", ont, simplify_cutoff))
    
    gse_ont_simplified <- clusterProfiler::simplify(
      gse_ont, cutoff = simplify_cutoff, by = "p.adjust", select_fun = min, measure = "Wang"
    )
    
    df_ont <- as.data.frame(gse_ont_simplified)
    df_ont$ONTOLOGY <- ont
    
    message(sprintf(
      "[GSEA] %s: %d terms before simplify -> %d after.",
      ont, nrow(as.data.frame(gse_ont)), nrow(df_ont)
    ))
    
    results_list[[ont]] <- df_ont
  }
  
  if (length(results_list) == 0) return(NULL)
  do.call(rbind, results_list)
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

message("=== Starting Gene Set Enrichment Analysis (GSEA) for primary cohort of GSE193677 ===")


# =========================================================
#                  Load data 
# =========================================================

message("[LOAD] Reading differential expression results...")

deg_table <- read.csv(file = file.path(paths$tables_bulkrnaseq, "DEGs_GSE193677_primary.csv"), sep = ",")


# =========================================================
#                  Prepare GSEA ranking
# =========================================================

message("[GSEA] Building pre-ranked gene list (moderated t-statistic)...")

rank <- data.frame(
  score = deg_table$t,
  Entrez = deg_table$ENTREZID
)

rank_df <- rank |>
  filter(!is.na(Entrez))

# One score per Entrez ID: keep the largest-magnitude statistic in case of duplicates
rank_df <- rank_df |>
  group_by(Entrez) |>
  slice_max(abs(score), n = 1, with_ties = FALSE) |>
  ungroup() |>
  arrange(desc(score))

list_allgenes_ranked <- rank_df$score
names(list_allgenes_ranked) <- rank_df$Entrez

list_allgenes_ranked <- list_allgenes_ranked[is.finite(list_allgenes_ranked)]

message(sprintf("[GSEA] Ranked gene list: %d genes.", length(list_allgenes_ranked)))

print(summary(list_allgenes_ranked))


# =========================================================
#                  Run GSEA - GO Database
# =========================================================

message("[GSEA] Running Gene Ontology enrichment analysis (BP/MF/CC)...")

go_df <- run_go_gsea_simplified(list_allgenes_ranked)

res_go <- gsea_lollipop(go_df, db_name = "GO", color = "#1B7837", title = "GSEA enrichment analysis – Primary Cohort - GO terms (BP/MF/CC)")


# =========================================================
#                  Run GSEA - KEGG Database
# =========================================================

message("[GSEA] Running KEGG pathway enrichment analysis...")

kegg <- gseKEGG(
  geneList = list_allgenes_ranked,
  organism = "hsa",
  keyType = "ncbi-geneid",
  minGSSize = 10,
  maxGSSize = 500,
  pAdjustMethod = "BH",
  pvalueCutoff = 1,
  eps = 1e-4,
  verbose = TRUE
)

kegg_df <- as.data.frame(kegg)

res_kegg <- gsea_lollipop(kegg, db_name = "clusterProfiler", color = "#762A83", title = "GSEA enrichment analysis – Primary Cohort - KEGG pathways")

kegg_tree <- add_redundancy_treeplot(
  kegg, title = "KEGG GSEA (redundancy tree) - CD vs Control",
  out_dir = paths$plots_bulk_rnaseq_primary_gsea, filename = "GSEA_kegg_treeplot.png"
)


# =========================================================
#                  Run GSEA - Reactome Database
# =========================================================

message("[GSEA] Running Reactome pathway enrichment analysis...")

reactome <- gsePathway(
  geneList = list_allgenes_ranked,
  organism = "human",
  minGSSize = 10,
  maxGSSize = 500,
  pAdjustMethod = "BH",
  pvalueCutoff = 1,
  verbose = TRUE
)

reactome_df <- as.data.frame(reactome)

res_reactome <- gsea_lollipop(reactome, db_name = "clusterProfiler", color = "darkorange3", title = "GSEA enrichment analysis – Primary Cohort - Reactome pathways")

reactome_tree <- add_redundancy_treeplot(
  reactome, title = "Reactome GSEA (redundancy tree) - CD vs Control",
  out_dir = paths$plots_bulk_rnaseq_primary_gsea, filename = "GSEA_reactome_treeplot.png"
)


# =========================================================
#                  Run GSEA - MSigDB Hallmark 
# =========================================================

message("[GSEA] Retrieving MSigDB Hallmark (H) and ImmunoSigDB (C7) gene sets...")

hallmark_sets_H  <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H")
hallmark_sets_C7 <- msigdbr::msigdbr(species = "Homo sapiens", collection = "C7")

hm_list_H <- split(
  x = hallmark_sets_H$ncbi_gene,
  f = hallmark_sets_H$gs_name
)

hm_list_C7 <- split(
  x = hallmark_sets_C7$ncbi_gene,
  f = hallmark_sets_C7$gs_name
)

# --- Collection H ---

message("[GSEA] Running fgsea on Hallmark (H) gene sets...")

hallmarks_H <- fgsea(
  pathways = hm_list_H,
  stats = list_allgenes_ranked,
  eps = 1e-4
)

hallmarks_H_df <- as.data.frame(hallmarks_H)
hallmarks_H_df_clean <- hallmarks_H_df[, !(names(hallmarks_H_df) %in% "leadingEdge")]

hallmarks_H <- clean_hallmark_names(hallmarks_H)

res_h <- gsea_lollipop(obj = hallmarks_H, db_name = "Hallmark", color = "firebrick", title = "GSEA enrichment analysis – Primary Cohort - Hallmark H pathways")

# --- Collection C7 (Immunologic signatures) ---

message("[GSEA] Running fgsea on ImmunoSigDB (C7) gene sets...")

hallmarks_C7 <- fgsea(
  pathways = hm_list_C7,
  stats = list_allgenes_ranked,
  eps = 1e-4
)

hallmarks_C7_df <- as.data.frame(hallmarks_C7)
hallmarks_C7_df_clean <- hallmarks_C7_df[, !(names(hallmarks_C7_df) %in% "leadingEdge")]

hallmarks_C7 <- clean_hallmark_names(hallmarks_C7)

res_c7 <- gsea_lollipop(obj = hallmarks_C7, db_name = "Hallmark", color = "darkgreen", title = "GSEA enrichment analysis – Primary Cohort - Immunologic Signatures C7 (ImmunoSigDB) pathways")


# =========================================================
#             Saving result tables
# =========================================================

message("[OUTPUT] Saving GSEA result tables...")

save_csv(go_df,               filename = "results_GSE193677_primary_GSEA_go_df.csv",         dir = paths$tables_bulkrnaseq)
save_csv(kegg_df,             filename = "results_GSE193677_primary_GSEA_kegg_df.csv",       dir = paths$tables_bulkrnaseq)
save_csv(reactome_df,         filename = "results_GSE193677_primary_GSEA_reactome_df.csv",   dir = paths$tables_bulkrnaseq)
save_csv(hallmarks_H_df_clean,  filename = "results_GSE193677_primary_GSEA_hallmarks_h_df.csv",  dir = paths$tables_bulkrnaseq)
save_csv(hallmarks_C7_df_clean, filename = "results_GSE193677_primary_GSEA_ImmunoSigDB_c7_df.csv", dir = paths$tables_bulkrnaseq)


# =========================================================
#             Saving final plots 
# =========================================================

message("[OUTPUT] Saving GSEA lollipop plots...")

save_plot(res_go[[1]],       filename = "GSEA_go.png",           dir = paths$plots_bulk_rnaseq_primary_gsea, width = 11, height = 9)
save_plot(res_kegg[[1]],     filename = "GSEA_kegg.png",         dir = paths$plots_bulk_rnaseq_primary_gsea, width = 11, height = 9)
save_plot(res_reactome[[1]], filename = "GSEA_reactome.png",     dir = paths$plots_bulk_rnaseq_primary_gsea, width = 11, height = 9)
save_plot(res_h[[1]],        filename = "GSEA_hallmarks_H.png",  dir = paths$plots_bulk_rnaseq_primary_gsea, width = 11, height = 9)
save_plot(res_c7[[1]],       filename = "GSEA_ImmunoSigDB_C7.png", dir = paths$plots_bulk_rnaseq_primary_gsea, width = 15, height = 9)


# =========================================================
#                  Save session info
# =========================================================

message("[OUTPUT] Saving session information...")

save_session_info(filename = "sessionInfo_GSEA_Primary_GSE193677.txt", dir = paths$logs, label = "Gene Set Enrichment Analysis (GSEA) - Primary cohort of GSE193677")

message("[OUTPUT] Session information saved to: ", paths$logs)


# =========================================================
#                  Final pipeline message
# =========================================================

message("=================================================")
message("[PIPELINE] Gene Set Enrichment Analysis completed successfully for GSE193677 (Ileum NonI, CD vs Control).")
message("=================================================")


