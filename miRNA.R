#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

mirna_dir <- "/fs/scratch/PAS2942/Alejandro/slides/miRNA"
sample_sheet_path <- file.path(mirna_dir, "mirnasam.tsv")
output_dir <- "/fs/scratch/PAS2942/Alejandro/slides"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(sample_sheet_path)) {
  stop("Sample sheet not found: ", sample_sheet_path)
}

sample_sheet <- fread(sample_sheet_path)

required_sample_cols <- c("File ID", "Case ID", "Project ID")
missing_sample_cols <- setdiff(required_sample_cols, colnames(sample_sheet))

if (length(missing_sample_cols) > 0) {
  stop(
    "mirnasam.tsv is missing required columns: ",
    paste(missing_sample_cols, collapse = ", ")
  )
}

sample_sheet <- unique(sample_sheet[, .(`File ID`, `Case ID`, `Project ID`)])

cat("[INFO] Loaded sample sheet with", nrow(sample_sheet), "file mappings\n")

find_mirna_file <- function(file_id) {
  file_dir <- file.path(mirna_dir, file_id)

  if (!dir.exists(file_dir)) {
    warning("[SKIP] Directory not found for File ID: ", file_id)
    return(NA_character_)
  }

  mirna_files <- list.files(
    file_dir,
    pattern = "\\.(txt|tsv|csv|gz)$",
    full.names = TRUE
  )

  mirna_files <- mirna_files[
    !grepl("mirnasam\\.tsv$", basename(mirna_files), ignore.case = TRUE)
  ]

  if (length(mirna_files) == 0) {
    warning("[SKIP] No miRNA file found in directory: ", file_dir)
    return(NA_character_)
  }

  if (length(mirna_files) > 1) {
    warning(
      "[INFO] Multiple miRNA files found for File ID ",
      file_id,
      ". Using first file: ",
      basename(mirna_files[1])
    )
  }

  return(mirna_files[1])
}

read_mirna_for_file_id <- function(file_id, case_id, project_id) {
  mirna_path <- find_mirna_file(file_id)

  if (is.na(mirna_path)) {
    return(NULL)
  }

  cat("[INFO] Reading:", mirna_path, "\n")

  dt <- tryCatch(
    fread(
      mirna_path,
      header = TRUE,
      data.table = TRUE,
      fill = TRUE,
      quote = ""
    ),
    error = function(e) {
      warning("[SKIP] Failed to read miRNA file for File ID ", file_id, ": ", e$message)
      return(NULL)
    }
  )

  if (is.null(dt)) {
    return(NULL)
  }

  required_mirna_cols <- c("miRNA_ID", "reads_per_million_miRNA_mapped")
  missing_mirna_cols <- setdiff(required_mirna_cols, colnames(dt))

  if (length(missing_mirna_cols) > 0) {
    warning(
      "[SKIP] miRNA file missing required columns for File ID ",
      file_id,
      ": ",
      paste(missing_mirna_cols, collapse = ", ")
    )
    return(NULL)
  }

  dt <- dt[
    !is.na(miRNA_ID) &
      miRNA_ID != "" &
      !is.na(reads_per_million_miRNA_mapped)
  ]

  if (nrow(dt) == 0) {
    warning("[SKIP] No usable miRNA rows for File ID: ", file_id)
    return(NULL)
  }

  dt[, reads_per_million_miRNA_mapped := as.numeric(reads_per_million_miRNA_mapped)]

  dt <- dt[!is.na(reads_per_million_miRNA_mapped)]

  if (nrow(dt) == 0) {
    warning("[SKIP] No numeric RPM values for File ID: ", file_id)
    return(NULL)
  }

  dt[reads_per_million_miRNA_mapped < 0, reads_per_million_miRNA_mapped := 0]

  dt[
    ,
    log2_RPM := log2(reads_per_million_miRNA_mapped + 1)
  ]

  dt <- dt[
    ,
    .(
      log2_RPM = mean(log2_RPM, na.rm = TRUE)
    ),
    by = miRNA_ID
  ]

  dt[, `Project ID` := project_id]
  dt[, `Case ID` := case_id]
  dt[, `File ID` := file_id]
  dt[, miRNA_File := basename(mirna_path)]

  dt <- dt[
    ,
    .(
      `Project ID`,
      `Case ID`,
      `File ID`,
      miRNA_File,
      miRNA_ID,
      log2_RPM
    )
  ]

  return(dt)
}

all_long_rows <- list()

project_ids <- sort(unique(sample_sheet$`Project ID`))

cat("[INFO] Found", length(project_ids), "unique projects\n")

for (project_id in project_ids) {
  cat("\n[INFO] Processing project:", project_id, "\n")

  project_sheet <- sample_sheet[`Project ID` == project_id]

  project_long_rows <- list()

  for (i in seq_len(nrow(project_sheet))) {
    file_id <- project_sheet$`File ID`[i]
    case_id <- project_sheet$`Case ID`[i]

    case_dt <- read_mirna_for_file_id(
      file_id = file_id,
      case_id = case_id,
      project_id = project_id
    )

    if (!is.null(case_dt)) {
      project_long_rows[[length(project_long_rows) + 1]] <- case_dt
      all_long_rows[[length(all_long_rows) + 1]] <- case_dt
    }
  }

  if (length(project_long_rows) == 0) {
    warning("[SKIP] No readable miRNA rows for project: ", project_id)
    next
  }

  project_long <- rbindlist(
    project_long_rows,
    use.names = TRUE,
    fill = TRUE
  )

  project_long <- project_long[
    ,
    .(
      log2_RPM = mean(log2_RPM, na.rm = TRUE),
      File_IDs = paste(sort(unique(`File ID`)), collapse = ";")
    ),
    by = .(`Project ID`, `Case ID`, miRNA_ID)
  ]

  setorder(project_long, `Project ID`, `Case ID`, miRNA_ID)

  project_long_file <- file.path(
    output_dir,
    paste0(project_id, "_miRNA_log2RPM_long.tsv")
  )

  fwrite(
    project_long,
    file = project_long_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  project_matrix <- dcast(
    project_long,
    miRNA_ID ~ `Case ID`,
    value.var = "log2_RPM"
  )

  setorder(project_matrix, miRNA_ID)

  project_matrix_file <- file.path(
    output_dir,
    paste0(project_id, "_miRNA_log2RPM_matrix.tsv")
  )

  fwrite(
    project_matrix,
    file = project_matrix_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  cat("[OK] Wrote long file:", project_long_file, "\n")
  cat("[OK] Wrote matrix:", project_matrix_file, "\n")
  cat("[OK]", project_id, "miRNAs:", nrow(project_matrix), "\n")
  cat("[OK]", project_id, "cases:", ncol(project_matrix) - 1, "\n")
}

if (length(all_long_rows) == 0) {
  stop("No readable miRNA data found.")
}

combined_long <- rbindlist(
  all_long_rows,
  use.names = TRUE,
  fill = TRUE
)

combined_long <- combined_long[
  ,
  .(
    log2_RPM = mean(log2_RPM, na.rm = TRUE),
    File_IDs = paste(sort(unique(`File ID`)), collapse = ";")
  ),
  by = .(`Project ID`, `Case ID`, miRNA_ID)
]

setorder(combined_long, `Project ID`, `Case ID`, miRNA_ID)

combined_long_file <- file.path(
  output_dir,
  "ALL_TCGA_miRNA_log2RPM_long.tsv"
)

fwrite(
  combined_long,
  file = combined_long_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

combined_matrix <- dcast(
  combined_long,
  miRNA_ID ~ `Case ID`,
  value.var = "log2_RPM"
)

setorder(combined_matrix, miRNA_ID)

combined_matrix_file <- file.path(
  output_dir,
  "ALL_TCGA_miRNA_log2RPM_matrix.tsv"
)

fwrite(
  combined_matrix,
  file = combined_matrix_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\n[OK] Wrote combined long file:", combined_long_file, "\n")
cat("[OK] Wrote combined matrix:", combined_matrix_file, "\n")
cat("[OK] Combined miRNAs:", nrow(combined_matrix), "\n")
cat("[OK] Combined cases:", ncol(combined_matrix) - 1, "\n")
cat("[DONE] miRNA log2 RPM preprocessing complete.\n")
