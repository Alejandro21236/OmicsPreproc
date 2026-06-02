#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

# -----------------------------
# Paths
# -----------------------------

input_dir <- "/fs/scratch/PAS2942/Alejandro/RNA"
output_dir <- "/fs/scratch/PAS2942/Alejandro/slides/RNA"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Find TCGA expression matrices
# -----------------------------

matrix_files <- list.files(
  input_dir,
  pattern = "^TCGA-.*_expression_matrix\\.tsv$",
  full.names = TRUE
)

if (length(matrix_files) == 0) {
  stop("No TCGA expression_matrix.tsv files found in: ", input_dir)
}

cat("[INFO] Found", length(matrix_files), "expression matrices\n")

# -----------------------------
# Log-normalization function
# -----------------------------
# Assumes:
#   - first column is named gene
#   - remaining columns are case/sample IDs
#   - values are raw or pre-normalized expression-like counts
#
# Operation:
#   log2(x + 1)
#
# This avoids exploding on zeros because biology already causes enough damage.

log_normalize_matrix <- function(dt) {
  if (!("gene" %in% colnames(dt))) {
    stop("Input matrix does not contain a 'gene' column.")
  }

  gene_col <- dt$gene
  expr_cols <- setdiff(colnames(dt), "gene")

  expr <- as.matrix(dt[, ..expr_cols])
  mode(expr) <- "numeric"

  expr[is.na(expr)] <- 0
  expr[expr < 0] <- 0

  expr_log <- log2(expr + 1)

  out <- data.table(gene = gene_col)
  out <- cbind(out, as.data.table(expr_log))

  setnames(out, c("gene", expr_cols))

  return(out)
}

# -----------------------------
# Process each cohort
# -----------------------------

for (file_path in matrix_files) {
  file_name <- basename(file_path)

  cohort <- sub("_expression_matrix\\.tsv$", "", file_name)

  cat("[INFO] Processing:", cohort, "\n")

  dt <- fread(file_path)

  if (nrow(dt) == 0) {
    warning("[SKIP] Empty matrix: ", file_name)
    next
  }

  if (!("gene" %in% colnames(dt))) {
    warning("[SKIP] Missing gene column: ", file_name)
    next
  }

  out_dt <- log_normalize_matrix(dt)

  out_file <- file.path(
    output_dir,
    paste0(cohort, "_expression_matrix_log2.tsv")
  )

  fwrite(
    out_dt,
    file = out_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  cat("[OK] Wrote:", out_file, "\n")
}

cat("[DONE] RNA log-normalization complete.\n")
