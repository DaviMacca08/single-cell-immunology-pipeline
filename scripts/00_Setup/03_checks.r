# =========================================================
# Validation utilities
# =========================================================

# ---------------------------------------------------------
# Check Seurat object validity
# ---------------------------------------------------------
check_seurat <- function(obj) {
  
  stopifnot(inherits(obj, "Seurat"))
  
  if (!"seurat_clusters" %in% colnames(obj@meta.data)) {
    warning("seurat_clusters not found in metadata")
  }
  else{
    message("OK: seurat_clusters present in metadata")
  }
  
  message("Seurat object validated")
}

# ---------------------------------------------------------
# Check feature availability
# ---------------------------------------------------------
check_features <- function(features, object) {
  
  present <- features[features %in% rownames(object)]
  missing <- setdiff(features, present)
  
  if (length(missing) > 0) {
    warning(
      length(missing), "/", length(features),
      " features missing"
    ) }
    else print("Everything is ok!")
  
  return(invisible(NULL))
}

# ---------------------------------------------------------
# Safe remove object
# ---------------------------------------------------------
safe_rm <- function(object_name) {
  
  if (exists(object_name, envir = .GlobalEnv)) {
    rm(list = object_name, envir = .GlobalEnv)
    gc()
    message("Removed: ", object_name)
  }
}