#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(GSVA)
  library(msigdbr)
})

# -----------------------------
# Paths
# -----------------------------

input_dir <- "/fs/scratch/PAS2942/Alejandro/slides/RNA"
output_dir <- "/fs/scratch/PAS2942/Alejandro/slides/Hallmark"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Find log-normalized matrices
# -----------------------------

matrix_files <- list.files(
  input_dir,
  pattern = "^TCGA-.*_expression_matrix_log2\\.tsv$",
  full.names = TRUE
)

if (length(matrix_files) == 0) {
  stop("No log-normalized expression matrices found in: ", input_dir)
}

cat("[INFO] Found", length(matrix_files), "log-normalized matrices\n")

# -----------------------------
# Load Hallmark gene sets
# -----------------------------

cat("[INFO] Loading MSigDB Hallmark gene sets\n")

hallmark_df <- msigdbr(
  species = "Homo sapiens",
  category = "H"
)

hallmark_sets <- split(
  hallmark_df$gene_symbol,
  hallmark_df$gs_name
)

cat("[INFO] Loaded", length(hallmark_sets), "Hallmark gene sets\n")

# -----------------------------
# ssGSEA function
# -----------------------------

run_ssgsea <- function(expr_mat, hallmark_sets) {
  # Compatible with newer GSVA versions
  if ("ssgseaParam" %in% getNamespaceExports("GSVA")) {
    param <- ssgseaParam(
      exprData = expr_mat,
      geneSets = hallmark_sets,
      normalize = TRUE
    )

    scores <- gsva(param)

  } else {
    # Compatible with older GSVA versions
    scores <- gsva(
      expr = expr_mat,
      gset.idx.list = hallmark_sets,
      method = "ssgsea",
      kcdf = "Gaussian",
      abs.ranking = FALSE,
      ssgsea.norm = TRUE,
      verbose = FALSE
    )
  }

  return(scores)
}

# -----------------------------
# Process each cohort
# -----------------------------

for (file_path in matrix_files) {
  file_name <- basename(file_path)

  cohort <- sub("_expression_matrix_log2\\.tsv$", "", file_name)

  cat("[INFO] Processing:", cohort, "\n")

  dt <- fread(file_path)

  if (!("gene" %in% colnames(dt))) {
    warning("[SKIP] Missing gene column: ", file_name)
    next
  }

  if (nrow(dt) == 0) {
    warning("[SKIP] Empty matrix: ", file_name)
    next
  }

  genes <- dt$gene
  expr_cols <- setdiff(colnames(dt), "gene")

  expr_mat <- as.matrix(dt[, ..expr_cols])
  mode(expr_mat) <- "numeric"

  rownames(expr_mat) <- genes

  # Remove duplicated gene symbols by averaging them.
  # Because naturally TCGA files cannot simply behave.
  if (any(duplicated(rownames(expr_mat)))) {
    cat("[INFO] Collapsing duplicated genes by mean:", cohort, "\n")

    expr_dt <- as.data.table(expr_mat)
    expr_dt[, gene := rownames(expr_mat)]

    expr_dt <- expr_dt[
      ,
      lapply(.SD, mean, na.rm = TRUE),
      by = gene
    ]

    genes <- expr_dt$gene
    expr_mat <- as.matrix(expr_dt[, -"gene"])
    mode(expr_mat) <- "numeric"
    rownames(expr_mat) <- genes
  }

  expr_mat[is.na(expr_mat)] <- 0

  common_genes <- intersect(rownames(expr_mat), unique(unlist(hallmark_sets)))

  cat("[INFO]", cohort, "genes overlapping Hallmark sets:", length(common_genes), "\n")

  if (length(common_genes) < 100) {
    warning("[SKIP] Too few Hallmark-overlapping genes for: ", cohort)
    next
  }

  scores <- run_ssgsea(expr_mat, hallmark_sets)

  scores_dt <- as.data.table(scores, keep.rownames = "pathway")

  out_file <- file.path(
    output_dir,
    paste0(cohort, "_Hallmark_ssGSEA_scores.tsv")
  )

  fwrite(
    scores_dt,
    file = out_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  cat("[OK] Wrote:", out_file, "\n")
}

cat("[DONE] Hallmark ssGSEA scoring complete.\n")
