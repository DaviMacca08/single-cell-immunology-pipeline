# =========================================================
# Project      : scRNA-Seq project - Bulk RNA-seq Analysis Pipeline
# Dataset      : GSE193677 (GPL16791 Illumina HiSeq 2500)
# Samples      : Ileal biopsy samples from CD patients and healthy controls
# Cohort       : Sensitivity cohort (Ileum I + NonI)
# Comparison   : Crohn's disease (CD) vs Control
# Script       : Over-Representation Analysis (ORA)
# Description  : Functional over-representation analysis of differentially
#                expressed genes using GO, KEGG, Reactome and MSigDB Hallmark
# =========================================================

# =========================================================
#                     Functions
# =========================================================

# 01-Create enrichment dot plot

plot_enrichment <- function(df, title, color, top_n = 10, wrap_width = 35, fdr_cutoff = 0.05) {
  
  if (is.null(df) || nrow(df) == 0) {
    message("[PLOT] Empty or NULL enrichment result. No plot generated.")
    return(NULL)
  }
  
  df <- as.data.frame(df)
  df <- df[!is.na(df$p.adjust), ]
  
  n_before <- nrow(df)
  df <- df[!is.na(df$p.adjust) & df$p.adjust < fdr_cutoff, ]
  
  message(sprintf(
    "[PLOT] %d / %d terms retained after FDR filter (p.adjust < %.2f).",
    nrow(df), n_before, fdr_cutoff
  ))
  
  if (nrow(df) == 0) {
    message("[PLOT] No enriched terms passed the FDR threshold.")
    return(NULL)
  }
  
  if (!"FoldEnrichment" %in% colnames(df)) {
    
    message("[PLOT] Calculating Fold Enrichment...")
    
    gene_ratio_num <- as.numeric(sub("/.*", "", df$GeneRatio))
    gene_ratio_den <- as.numeric(sub(".*/", "", df$GeneRatio))
    
    bg_ratio_num <- as.numeric(sub("/.*", "", df$BgRatio))
    bg_ratio_den <- as.numeric(sub(".*/", "", df$BgRatio))
    
    df$FoldEnrichment <- (gene_ratio_num / gene_ratio_den) / (bg_ratio_num / bg_ratio_den)
  }
  
  has_set <- "Set" %in% colnames(df)
  has_ontology <- "ONTOLOGY" %in% colnames(df)
  
  if (has_set) {
    df$Set <- factor(df$Set, levels = c("Module", "DEG", "Intersection"))
  } else {
    df$GeneRatioNum <- as.numeric(sub("/.*", "", df$GeneRatio)) /
      as.numeric(sub(".*/", "", df$GeneRatio))
  }
  
  df$Description <- stringr::str_wrap(df$Description, width = wrap_width)
  
  group_vars <- intersect(c("Set", "ONTOLOGY"), colnames(df))
  
  if (length(group_vars) > 0) {
    Top <- df |>
      dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) |>
      dplyr::arrange(p.adjust, .by_group = TRUE) |>
      dplyr::slice_head(n = top_n) |>
      dplyr::ungroup()
  } else {
    Top <- df |>
      dplyr::arrange(p.adjust) |>
      dplyr::slice_head(n = top_n)
  }
  
  Top$Description <- factor(
    Top$Description,
    levels = Top$Description[order(Top$p.adjust)]
  )
  
  x_var <- if (has_set) "Set" else "GeneRatioNum"
  
  stopifnot(x_var %in% colnames(Top))
  
  p <- ggplot(Top, aes(x = .data[[x_var]], y = Description)) +
    geom_point(aes(size = FoldEnrichment, color = p.adjust)) +
    scale_color_gradient(low = color, high = "#2C7BB6", labels = scales::scientific) +
    scale_size_continuous(range = c(2, 6), breaks = pretty(Top$FoldEnrichment, n = 3)) +
    theme_bw() +
    theme(
      axis.text.y = element_text(size = 9),
      axis.text.x = element_text(size = 11, angle = 0, hjust = 0.5),
      strip.text.y = element_text(face = "bold", size = 10)
    ) +
    labs(
      title = title,
      x = if (has_set) "" else "Gene Ratio",
      y = "",
      size = "Fold Enrichment",
      color = "FDR"
    )
  
  if (has_ontology) {
    p <- p + facet_grid(ONTOLOGY ~ ., scales = "free_y", space = "free_y")
  }
  
  message("[PLOT] Top terms plotted: ", nrow(Top))
  
  p
  
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
  
  # nCluster non può eccedere il numero di termini effettivamente mostrati/disponibili
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

# 03-GO Simplified

run_go_ora_simplified <- function(gene_list, universe, ontologies = c("BP", "MF", "CC"),
                                  simplify_cutoff = 0.7) {
  
  results_list <- list()
  
  for (ont in ontologies) {
    
    message(sprintf("[ORA] Running enrichGO for ontology = %s...", ont))
    
    ego_ont <- enrichGO(
      gene = gene_list,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = ont,
      universe = universe,
      pAdjustMethod = "BH",
      minGSSize = 10,
      maxGSSize = 500
    )
    
    if (is.null(ego_ont) || nrow(as.data.frame(ego_ont)) == 0) {
      message(sprintf("[ORA] No enriched terms for ontology = %s.", ont))
      next
    }
    
    message(sprintf("[ORA] Simplifying redundant %s terms (cutoff = %.1f)...", ont, simplify_cutoff))
    
    ego_ont_simplified <- clusterProfiler::simplify(
      ego_ont,
      cutoff = simplify_cutoff,
      by = "p.adjust",
      select_fun = min,
      measure = "Wang"
    )
    
    df_ont <- as.data.frame(ego_ont_simplified)
    df_ont$ONTOLOGY <- ont
    
    message(sprintf(
      "[ORA] %s: %d terms before simplify -> %d after.",
      ont, nrow(as.data.frame(ego_ont)), nrow(df_ont)
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

message("=== Starting over-representation sensitivity analysis for GSE193677 ===")


# =========================================================
#                  Load data 
# =========================================================

message("[LOAD] Reading differential expression results...")

deg_table <- read.csv(file = file.path(paths$tables_bulkrnaseq, "DEGs_GSE193677_sensitivity.csv"), sep = "," )

message("[LOAD] Reading background gene list...")

background <- read.csv(file = file.path(paths$tables_bulkrnaseq, "background_GSE193677_sensitivity.csv"), sep = "," )


# =========================================================
#                  Prepare ORA gene universe
# =========================================================

message("[ORA] Converting gene universe from ENSEMBL to ENTREZID...")

universe_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = background$x,
  keytype = "ENSEMBL",
  columns = c("ENSEMBL", "ENTREZID")
)

list_universe_entrez <- universe_map$ENTREZID |>
  na.omit() |>
  unique()

message(
  "[ORA] Universe: ",
  length(background$x),
  " ENSEMBL IDs -> ",
  length(list_universe_entrez),
  " unique ENTREZID."
)


# =========================================================
#  Prepare DEG and background gene lists for enrichment analysis
# =========================================================

logcf_cutoff <- 1
padj_cutoff  <- 0.05

sig_deg <- subset(
  deg_table,
  adj.P.Val < padj_cutoff & abs(logFC) > logcf_cutoff
)

n_before <- nrow(sig_deg)

sig_deg <- sig_deg[!(is.na(sig_deg$SYMBOL) & is.na(sig_deg$ENTREZID) & is.na(sig_deg$GENENAME)), ]

# Up-regultaed genes

sig_deg_up <- sig_deg |> 
  filter(logFC > 0)

# Down-regultaed genes

sig_deg_down <- sig_deg |> 
  filter(logFC < 0)

message(sprintf(
  "[DEG] Removed %d significant gene(s) lacking all annotation fields (SYMBOL, ENTREZID and GENENAME).",
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

# Create Entrez ID gene lists

list_deg <- unique(sig_deg$ENTREZID[!is.na(sig_deg$ENTREZID)])
list_deg_up <- unique(sig_deg_up$ENTREZID[!is.na(sig_deg_up$ENTREZID)])
list_deg_down <- unique(sig_deg_down$ENTREZID[!is.na(sig_deg_down$ENTREZID)])


# =========================================================
#                  Run ORA - GO Database
# =========================================================

# GO Enrichment - all / up / down DEGs

go_all_df  <- run_go_ora_simplified(list_deg, list_universe_entrez)
go_up_df   <- run_go_ora_simplified(list_deg_up, list_universe_entrez)
go_down_df <- run_go_ora_simplified(list_deg_down, list_universe_entrez)

# Generate enrichment plots

go_all_plot <- plot_enrichment(go_all_df, "GO Enrichment - Sensitivity Cohort: CD vs Control (All genes)", color = "#C44E52")
go_up_plot <- plot_enrichment(go_up_df, "GO Enrichment - Sensitivity Cohort: CD vs Control (Up-regulated genes)", color = "#009E73")
go_down_plot <- plot_enrichment(go_down_df, "GO Enrichment - Sensitivity Cohort: CD vs Control (Down-regulated genes)", color = "#CC79A7")

# Save enrichment plots

save_plot(go_all_plot, filename = "GO_all_genes.png", dir = paths$plots_bulk_rnaseq_sensitivity_ora, width = 10, height = 9)
save_plot(go_up_plot, filename = "GO_Up_genes.png", dir = paths$plots_bulk_rnaseq_sensitivity_ora, width = 10, height = 9)
save_plot(go_down_plot, filename = "GO_Down_genes.png", dir = paths$plots_bulk_rnaseq_sensitivity_ora, width = 10, height = 9)


# =========================================================
#                  Run ORA - KEGG Database
# =========================================================

# KEGG enrichment - all DEGs

kegg_all <- enrichKEGG(
  gene = list_deg,
  organism = "hsa",
  keyType = "ncbi-geneid",
  universe = list_universe_entrez,
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 500
)

# KEGG enrichment - upregulated DEGs

kegg_up <- enrichKEGG(
  gene = list_deg_up,
  organism = "hsa",
  keyType = "ncbi-geneid",
  universe = list_universe_entrez,
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 500
)

# KEGG enrichment - downregulated DEGs

kegg_down <- enrichKEGG(
  gene = list_deg_down,
  organism = "hsa",
  keyType = "ncbi-geneid",
  universe = list_universe_entrez,
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 500
)

# Generate enrichment plots

kegg_all_plot <- plot_enrichment(kegg_all, title = "KEGG Enrichment - Sensitivity Cohort: CD vs Control (All genes)", color = "#C44E52")
kegg_up_plot <- plot_enrichment(kegg_up, title = "KEGG Enrichment - Sensitivity Cohort: CD vs Control (Up-regulated genes)", color = "#009E73")
kegg_down_plot <- plot_enrichment(kegg_down, title = "KEGG Enrichment - Sensitivity Cohort: CD vs Control (Down-regulated genes)", color = "#CC79A7")

# Generate treeplots

kegg_all_tree  <- add_redundancy_treeplot(kegg_all,  "KEGG (redundancy tree) - Sensitivity Cohort: CD vs Control (All genes)",  paths$plots_bulk_rnaseq_sensitivity_ora, "KEGG_all_genes_treeplot.png")
kegg_up_tree   <- add_redundancy_treeplot(kegg_up,   "KEGG (redundancy tree) - Sensitivity Cohort: CD vs Control (Up-regulated genes)", paths$plots_bulk_rnaseq_sensitivity_ora, "KEGG_up_genes_treeplot.png")
kegg_down_tree <- add_redundancy_treeplot(kegg_down, "KEGG (redundancy tree) - Sensitivity Cohort: CD vs Control (Down-regulated genes)", paths$plots_bulk_rnaseq_sensitivity_ora, "KEGG_down_genes_treeplot.png")

# Save enrichment plots

save_plot(kegg_all_plot, filename = "KEGG_all_genes.png", dir = paths$plots_bulk_rnaseq_sensitivity_ora, width = 10, height = 9)
save_plot(kegg_up_plot, filename = "KEGG_up_genes.png", dir = paths$plots_bulk_rnaseq_sensitivity_ora, width = 10, height = 9)


# =========================================================
#                  Run ORA - Reactome Database
# =========================================================

# Reactome enrichment - all DEGs

react_all <- enrichPathway(
  gene = list_deg,
  organism = "human",
  universe = list_universe_entrez,
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 500,
  readable = TRUE
)

# Reactome enrichment - upregulated DEGs

react_up <- enrichPathway(
  gene = list_deg_up,
  organism = "human",
  universe = list_universe_entrez,
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 500,
  readable = TRUE
)

# Reactome enrichment - downregulated DEGs

react_down <- enrichPathway(
  gene = list_deg_down,
  organism = "human",
  universe = list_universe_entrez,
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 500,
  readable = TRUE
)

# Generate enrichment plots

react_all_plot <- plot_enrichment(react_all, title = "Reactome Enrichment - Sensitivity Cohort: CD vs Control (All genes)", color = "#C44E52")
react_up_plot <- plot_enrichment(react_up, title = "Reactome Enrichment - Sensitivity Cohort: CD vs Control (Up-regulated genes)", color = "#009E73")
react_down_plot <- plot_enrichment(react_down, title = "Reactome Enrichment - Sensitivity Cohort: CD vs Control (Down-regulated genes)", color = "#CC79A7")

# Generate treeplots

react_all_tree  <- add_redundancy_treeplot(react_all,  "Reactome (redundancy tree) - Sensitivity Cohort: CD vs Control (All genes)",  paths$plots_bulk_rnaseq_sensitivity_ora, "Reactome_all_genes_treeplot.png")
react_up_tree   <- add_redundancy_treeplot(react_up,   "Reactome (redundancy tree) - Sensitivity Cohort: CD vs Control (Up-regulated genes)", paths$plots_bulk_rnaseq_sensitivity_ora, "Reactome_up_genes_treeplot.png")
react_down_tree <- add_redundancy_treeplot(react_down, "Reactome (redundancy tree) - Sensitivity Cohort: CD vs Control (Down-regulated genes)", paths$plots_bulk_rnaseq_sensitivity_ora, "Reactome_down_genes_treeplot.png")

# Save enrichment plots

save_plot(react_all_plot, filename = "Reactome_all_genes.png", dir = paths$plots_bulk_rnaseq_sensitivity_ora, width = 10, height = 9)
save_plot(react_up_plot, filename = "Reactome_up_genes.png", dir = paths$plots_bulk_rnaseq_sensitivity_ora, width = 10, height = 9)


# =========================================================
#                  Run ORA - MSigDB Hallmark
# =========================================================

# ------------------
# Collection H
# ------------------

message("[ORA] Retrieving MSigDB Hallmark (H) gene sets...")

if ("collection" %in% names(formals(msigdbr::msigdbr))) {
  hallmark_sets <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H")
  hallmark_c7 <- msigdbr::msigdbr(species = "Homo sapiens", collection = "C7")
  message("Collection found.")
} else {
  hallmark_sets <- msigdbr::msigdbr(species = "Homo sapiens", category = "H")
  hallmark_c7 <- msigdbr::msigdbr(species = "Homo sapiens", category = "C7")
  
  message("Category found.")
}

entrez_col <- if ("ncbi_gene" %in% colnames(hallmark_sets)) "ncbi_gene" else "entrez_gene"

message(sprintf("[ORA] Using '%s' as the Entrez ID column.", entrez_col))

# TERM2GENE - H
term2gene_h <- unique(hallmark_sets[, c("gs_name", entrez_col)])
colnames(term2gene_h) <- c("gs_name", "entrez_gene")
term2gene_h$entrez_gene <- as.character(term2gene_h$entrez_gene)

# MSigDB Hallmark enrichment - all DEGs 

hallmark_all_h <- enricher(
  gene = as.character(list_deg),
  TERM2GENE = term2gene_h,
  universe = as.character(list_universe_entrez),
  pAdjustMethod = "BH",
  minGSSize = 1,
  maxGSSize = 5000
)

# Hallmark - upregulated DEGs

hallmark_up_h <- enricher(
  gene = as.character(list_deg_up),
  TERM2GENE = term2gene_h,
  universe = as.character(list_universe_entrez),
  pAdjustMethod = "BH",
  minGSSize = 1,
  maxGSSize = 5000
)

# Hallmark - downregulated DEGs

hallmark_down_h <- enricher(
  gene = as.character(list_deg_down),
  TERM2GENE = term2gene_h,
  universe = as.character(list_universe_entrez),
  pAdjustMethod = "BH",
  minGSSize = 1,
  maxGSSize = 5000
)

clean_hallmark_names <- function(x) {
  if (is.null(x) || nrow(as.data.frame(x)) == 0) return(x)
  x@result$Description <- tolower(gsub("_", " ", sub("^HALLMARK_", "", x@result$ID)))
  x
}

hallmark_all_h  <- clean_hallmark_names(hallmark_all_h)
hallmark_up_h   <- clean_hallmark_names(hallmark_up_h)
hallmark_down_h <- clean_hallmark_names(hallmark_down_h)

# Generate enrichment plots

hallmark_all_plot_h  <- plot_enrichment(hallmark_all_h,  title = "Hallmark H Enrichment - Sensitivity Cohort: CD vs Control (All genes)", color = "#C44E52")
hallmark_up_plot_h   <- plot_enrichment(hallmark_up_h,   title = "Hallmark H Enrichment - Sensitivity Cohort: CD vs Control (Up-regulated genes)", color = "#009E73")
hallmark_down_plot_h <- plot_enrichment(hallmark_down_h, title = "Hallmark H Enrichment - Sensitivity Cohort: CD vs Control (Down-regulated genes)", color = "#CC79A7")

# Save enrichment plots

save_plot(hallmark_all_plot_h,   filename = "Hallmark_H_all_genes.png",   dir = paths$plots_bulk_rnaseq_sensitivity_ora, width = 10, height = 9)
save_plot(hallmark_up_plot_h,   filename = "Hallmark_H_up_genes.png",   dir = paths$plots_bulk_rnaseq_sensitivity_ora, width = 10, height = 9)
save_plot(hallmark_down_plot_h,   filename = "Hallmark_H_down_genes.png",   dir = paths$plots_bulk_rnaseq_sensitivity_ora, width = 10, height = 9)

# ------------------
# Collection C7
# ------------------

# TERM2GENE - C7
term2gene_c7 <- unique(hallmark_c7[, c("gs_name", entrez_col)])
colnames(term2gene_c7) <- c("gs_name", "entrez_gene")
term2gene_c7$entrez_gene <- as.character(term2gene_c7$entrez_gene)

# MSigDB Hallmark enrichment - all DEGs 

hallmark_all_c7 <- enricher(
  gene = as.character(list_deg),
  TERM2GENE = term2gene_c7,
  universe = as.character(list_universe_entrez),
  pAdjustMethod = "BH",
  minGSSize = 1,
  maxGSSize = 5000
)

# Hallmark - upregulated DEGs

hallmark_up_c7 <- enricher(
  gene = as.character(list_deg_up),
  TERM2GENE = term2gene_c7,
  universe = as.character(list_universe_entrez),
  pAdjustMethod = "BH",
  minGSSize = 1,
  maxGSSize = 5000
)

# Hallmark - downregulated DEGs

hallmark_down_c7 <- enricher(
  gene = as.character(list_deg_down),
  TERM2GENE = term2gene_c7,
  universe = as.character(list_universe_entrez),
  pAdjustMethod = "BH",
  minGSSize = 1,
  maxGSSize = 5000
)

hallmark_all_c7  <- clean_hallmark_names(hallmark_all_c7)
hallmark_up_c7  <- clean_hallmark_names(hallmark_up_c7)
hallmark_down_c7 <- clean_hallmark_names(hallmark_down_c7)

# Generate enrichment plots

hallmark_all_plot_c7  <- plot_enrichment(hallmark_all_c7,  title = "Immunologic Signatures C7 Enrichment - Sensitivity Cohort: CD vs Control (All genes)", color = "#C44E52")
hallmark_up_plot_c7   <- plot_enrichment(hallmark_up_c7,   title = "Immunologic Signatures C7 Enrichment - Sensitivity Cohort: CD vs Control (Up-regulated genes)", color = "#009E73")
hallmark_down_plot_c7 <- plot_enrichment(hallmark_down_c7, title = "Immunologic Signatures C7 Enrichment - Sensitivity Cohort: CD vs Control (Down-regulated genes)", color = "#CC79A7")

# Generate treeplots

hallmark_all_c7_tree  <- add_redundancy_treeplot(hallmark_all_c7,  title = "Immunologic Signatures C7 (redundancy tree) - Sensitivity Cohort: CD vs Control (All genes)",  paths$plots_bulk_rnaseq_sensitivity_ora, "ImmunoSigDB_C7_all_genes_treeplot.png")
hallmark_up_c7_tree   <- add_redundancy_treeplot(hallmark_up_c7,   title = "Immunologic Signatures C7 (redundancy tree) - Sensitivity Cohort: CD vs Control (Up-regulated genes)", paths$plots_bulk_rnaseq_sensitivity_ora, "ImmunoSigDB_C7_up_genes_treeplot.png")
hallmark_down_c7_tree <- add_redundancy_treeplot(hallmark_down_c7, title = "Immunologic Signatures C7 (redundancy tree) - Sensitivity Cohort: CD vs Control (Down-regulated genes)", paths$plots_bulk_rnaseq_sensitivity_ora, "ImmunoSigDB_C7_down_genes_treeplot.png")

# Save enrichment plots

save_plot(hallmark_all_plot_c7,  filename = "ImmunoSigDB_C7_all_genes.png",  dir = paths$plots_bulk_rnaseq_sensitivity_ora, width = 10, height = 9)
save_plot(hallmark_up_plot_c7, filename = "ImmunoSigDB_C7_up_genes.png", dir = paths$plots_bulk_rnaseq_sensitivity_ora, width = 10, height = 9)
save_plot(hallmark_down_plot_c7, filename = "ImmunoSigDB_C7_down_genes.png", dir = paths$plots_bulk_rnaseq_sensitivity_ora, width = 10, height = 9)


# =========================================================
#                  Save session info
# =========================================================

message("[OUTPUT] Saving session information...")

save_session_info(filename = "sessionInfo_ORA_sensitivity_GSE193677.txt", dir = paths$logs, label = "Over-Representation Analysis (ORA) - Sensitivity cohort GSE193677")

message("[OUTPUT] Session information saved to:" , paths$logs)


# =========================================================
#                  Final pipeline message
# =========================================================

message("=================================================")
message("[PIPELINE] Over-Representation Analysis completed successfully for sensitivity analysis of GSE193677 (CD vs  Control).")
message("=================================================")



