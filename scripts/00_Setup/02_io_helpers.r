# =========================================================
# IO helpers (plot + file saving utilities)
# =========================================================

# ---------------------------------------------------------
# Save ggplot / patchwork objects
# ---------------------------------------------------------
save_plot <- function(..., filename, dir, width = 12, height = 8, dpi = 300) {
  
  plots <- list(...)
  
  # if (length(plots) == 1) {
  #   final_plot <- plots[[1]]
  # } else {
  #   final_plot <- patchwork::wrap_plots(plots)
  # }
  # 
  ggsave(
    filename = file.path(dir, filename),
    plot = plots,
    width = width,
    height = height,
    dpi = dpi
  )
  
  message("Saved plot: ", filename)
}

# ---------------------------------------------------------
# Save base R graphics (e.g. pheatmap)
# ---------------------------------------------------------
open_pdf <- function(filename, dir, width = 10, height = 8) {
  
  pdf(file.path(dir, filename), width = width, height = height)
  invisible(NULL)
}

close_pdf <- function() {
  dev.off()
}

# ---------------------------------------------------------
# Save R object safely
# ---------------------------------------------------------
save_rds <- function(object, filename, dir) {
  
  saveRDS(object, file.path(dir, filename))
  message("Saved object: ", filename)
}

# ---------------------------------------------------------
# Save csv object safely
# ---------------------------------------------------------

save_csv <- function(object, filename, dir) {
  
  write.csv(x = object, file = file.path(dir, filename), row.names = FALSE)
  
  message("Saved table: ", filename)
}

# ---------------------------------------------------------
# Save session info safely
# ---------------------------------------------------------
save_session_info <- function(filename = "sessionInfo.txt", dir, label = NULL) {
  
  if (!dir.exists(dir)) {
    stop("Directory does not exist: ", dir)
  }
  
  session_file <- file.path(dir, filename)
  
  header <- if (!is.null(label)) {
    paste0("Session information (", label, ")\n")
  } else {
    "Session information\n"
  }
  
  writeLines(
    c(
      header,
      as.character(Sys.time()),
      "\n",
      capture.output(sessionInfo())
    ),
    session_file
  )
  
  message("Session information saved at: ", session_file)
  
  invisible(session_file)
}
