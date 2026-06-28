# =========================================================
# Reproducibility
# =========================================================

set_seed <- function(seed = 1234) {
  
  set.seed(seed)
  
  message("Seed set to: ", seed)
}