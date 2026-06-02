#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

wxs_dir <- "/fs/scratch/PAS2942/Alejandro/slides/WXS"
sample_sheet_path <- file.path(wxs_dir, "WXSam.tsv")
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
    "WXSam.tsv is missing required columns: ",
    paste(missing_sample_cols, collapse = ", ")
  )
}

sample_sheet <- unique(sample_sheet[, .(`File ID`, `Case ID`, `Project ID`)])

cat("[INFO] Loaded sample sheet with", nrow(sample_sheet), "file mappings\n")

# For binary mutation calling, use nonsilent coding/splice mutations.
# This avoids counting Silent, Intron, UTR, RNA, etc. as "mutated" for the model.
mutated_variant_classes <- c(
  "Missense_Mutation",
  "Nonsense_Mutation",
  "Nonstop_Mutation",
  "Splice_Site",
  "Translation_Start_Site",
  "Frame_Shift_Del",
  "Frame_Shift_Ins",
  "In_Frame_Del",
  "In_Frame_Ins"
)

find_maf_file <- function(file_id) {
  file_dir <- file.path(wxs_dir, file_id)

  if (!dir.exists(file_dir)) {
    warning("[SKIP] Directory not found for File ID: ", file_id)
    return(NA_character_)
  }

  maf_files <- list.files(
    file_dir,
    pattern = "\\.maf(\\.gz)?$",
    full.names = TRUE
  )

  if (length(maf_files) == 0) {
    warning("[SKIP] No MAF file found in directory: ", file_dir)
    return(NA_character_)
  }

  if (length(maf_files) > 1) {
    warning(
      "[INFO] Multiple MAF files found for File ID ",
      file_id,
      ". Using first file: ",
      basename(maf_files[1])
    )
  }

  return(maf_files[1])
}

read_wxs_for_file_id <- function(file_id, case_id, project_id) {
  maf_path <- find_maf_file(file_id)

  if (is.na(maf_path)) {
    return(NULL)
  }

  cat("[INFO] Reading:", maf_path, "\n")

  maf_dt <- tryCatch(
    fread(
      maf_path,
      sep = "\t",
      header = TRUE,
      data.table = TRUE,
      fill = TRUE,
      quote = "",
      comment.char = "#"
    ),
    error = function(e) {
      warning("[SKIP] Failed to read MAF for File ID ", file_id, ": ", e$message)
      return(NULL)
    }
  )

  if (is.null(maf_dt)) {
    return(NULL)
  }

  if (nrow(maf_dt) == 0) {
    warning("[SKIP] Empty MAF file for File ID: ", file_id)
    return(NULL)
  }

  required_maf_cols <- c("Hugo_Symbol", "Variant_Classification")
  missing_maf_cols <- setdiff(required_maf_cols, colnames(maf_dt))

  if (length(missing_maf_cols) > 0) {
    warning(
      "[SKIP] MAF file missing required columns for File ID ",
      file_id,
      ": ",
      paste(missing_maf_cols, collapse = ", ")
    )
    return(NULL)
  }

  maf_dt <- maf_dt[
    !is.na(Hugo_Symbol) &
      Hugo_Symbol != "" &
      !is.na(Variant_Classification) &
      Variant_Classification %in% mutated_variant_classes
  ]

  if ("Mutation_Status" %in% colnames(maf_dt)) {
    maf_dt <- maf_dt[
      is.na(Mutation_Status) |
        Mutation_Status == "" |
        Mutation_Status == "Somatic"
    ]
  }

  if (nrow(maf_dt) == 0) {
    return(data.table(
      `Project ID` = character(),
      `Case ID` = character(),
      Gene = character(),
      Mutated = integer()
    ))
  }

  out <- unique(maf_dt[, .(Gene = Hugo_Symbol)])
  out[, `Project ID` := project_id]
  out[, `Case ID` := case_id]
  out[, Mutated := 1L]

  out <- out[, .(`Project ID`, `Case ID`, Gene, Mutated)]

  return(out)
}

all_mutation_rows <- list()

project_ids <- sort(unique(sample_sheet$`Project ID`))

cat("[INFO] Found", length(project_ids), "unique projects\n")

for (project_id in project_ids) {
  cat("\n[INFO] Processing project:", project_id, "\n")

  project_sheet <- sample_sheet[`Project ID` == project_id]
  project_cases <- sort(unique(project_sheet$`Case ID`))

  project_rows <- list()

  for (i in seq_len(nrow(project_sheet))) {
    file_id <- project_sheet$`File ID`[i]
    case_id <- project_sheet$`Case ID`[i]

    case_dt <- read_wxs_for_file_id(
      file_id = file_id,
      case_id = case_id,
      project_id = project_id
    )

    if (!is.null(case_dt) && nrow(case_dt) > 0) {
      project_rows[[length(project_rows) + 1]] <- case_dt
      all_mutation_rows[[length(all_mutation_rows) + 1]] <- case_dt
    }
  }

  if (length(project_rows) == 0) {
    warning("[SKIP] No mutated genes found for project: ", project_id)
    next
  }

  project_long <- rbindlist(
    project_rows,
    use.names = TRUE,
    fill = TRUE
  )

  project_long <- unique(project_long[, .(`Project ID`, `Case ID`, Gene, Mutated)])

  project_matrix <- dcast(
    project_long,
    Gene ~ `Case ID`,
    value.var = "Mutated",
    fill = 0
  )

  missing_cases <- setdiff(project_cases, colnames(project_matrix))

  for (case_id in missing_cases) {
    project_matrix[, (case_id) := 0L]
  }

  case_cols <- sort(setdiff(colnames(project_matrix), "Gene"))

  setcolorder(project_matrix, c("Gene", case_cols))
  setorder(project_matrix, Gene)

  out_file <- file.path(
    output_dir,
    paste0(project_id, "_WXS_binary_mutation_matrix.tsv")
  )

  fwrite(
    project_matrix,
    file = out_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  cat("[OK] Wrote:", out_file, "\n")
  cat("[OK]", project_id, "genes:", nrow(project_matrix), "\n")
  cat("[OK]", project_id, "cases:", length(case_cols), "\n")
}

if (length(all_mutation_rows) == 0) {
  stop("No readable WXS mutation data found.")
}

combined_long <- rbindlist(
  all_mutation_rows,
  use.names = TRUE,
  fill = TRUE
)

combined_long <- unique(combined_long[, .(`Project ID`, `Case ID`, Gene, Mutated)])

all_cases <- sort(unique(sample_sheet$`Case ID`))

combined_matrix <- dcast(
  combined_long,
  Gene ~ `Case ID`,
  value.var = "Mutated",
  fill = 0
)

missing_cases <- setdiff(all_cases, colnames(combined_matrix))

for (case_id in missing_cases) {
  combined_matrix[, (case_id) := 0L]
}

case_cols <- sort(setdiff(colnames(combined_matrix), "Gene"))

setcolorder(combined_matrix, c("Gene", case_cols))
setorder(combined_matrix, Gene)

combined_file <- file.path(
  output_dir,
  "ALL_TCGA_WXS_binary_mutation_matrix.tsv"
)

fwrite(
  combined_matrix,
  file = combined_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\n[OK] Wrote combined matrix:", combined_file, "\n")
cat("[OK] Combined genes:", nrow(combined_matrix), "\n")
cat("[OK] Combined cases:", length(case_cols), "\n")
cat("[DONE] WXS binary mutation matrix preprocessing complete.\n")
