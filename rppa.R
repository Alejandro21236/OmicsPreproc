#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

rppa_dir <- "/fs/scratch/PAS2942/Alejandro/slides/rppa"
sample_sheet_path <- file.path(rppa_dir, "proteinsam.tsv")
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
    "proteinsam.tsv is missing required columns: ",
    paste(missing_sample_cols, collapse = ", ")
  )
}

sample_sheet <- unique(sample_sheet[, .(`File ID`, `Case ID`, `Project ID`)])

cat("[INFO] Loaded sample sheet with", nrow(sample_sheet), "file mappings\n")

read_rppa_for_file_id <- function(file_id, case_id, project_id) {
  file_dir <- file.path(rppa_dir, file_id)

  if (!dir.exists(file_dir)) {
    warning("[SKIP] Directory not found for File ID: ", file_id)
    return(NULL)
  }

  rppa_files <- list.files(
    file_dir,
    pattern = "\\.(txt|tsv|csv|gz)$",
    full.names = TRUE
  )

  rppa_files <- rppa_files[!grepl("proteinsam\\.tsv$", basename(rppa_files))]

  if (length(rppa_files) == 0) {
    warning("[SKIP] No RPPA data file found in directory: ", file_dir)
    return(NULL)
  }

  if (length(rppa_files) > 1) {
    warning(
      "[INFO] Multiple RPPA files found for File ID ",
      file_id,
      ". Using first file: ",
      basename(rppa_files[1])
    )
  }

  rppa_path <- rppa_files[1]

  cat("[INFO] Reading:", rppa_path, "\n")

  dt <- tryCatch(
    fread(
      rppa_path,
      header = TRUE,
      data.table = TRUE,
      fill = TRUE,
      quote = ""
    ),
    error = function(e) {
      warning("[SKIP] Failed to read RPPA file for File ID ", file_id, ": ", e$message)
      return(NULL)
    }
  )

  if (is.null(dt)) {
    return(NULL)
  }

  required_rppa_cols <- c("peptide_target", "protein_expression")
  missing_rppa_cols <- setdiff(required_rppa_cols, colnames(dt))

  if (length(missing_rppa_cols) > 0) {
    warning(
      "[SKIP] RPPA file missing required columns for File ID ",
      file_id,
      ": ",
      paste(missing_rppa_cols, collapse = ", ")
    )
    return(NULL)
  }

  dt <- dt[
    !is.na(peptide_target) &
      peptide_target != "" &
      !is.na(protein_expression)
  ]

  if (nrow(dt) == 0) {
    warning("[SKIP] No usable peptide_target/protein_expression rows for File ID: ", file_id)
    return(NULL)
  }

  dt[, protein_expression := as.numeric(protein_expression)]

  dt <- dt[
    !is.na(protein_expression),
    .(
      protein_expression = mean(protein_expression, na.rm = TRUE)
    ),
    by = peptide_target
  ]

  wide <- dcast(
    dt,
    . ~ peptide_target,
    value.var = "protein_expression"
  )

  wide[, "." := NULL]

  wide[, `Project ID` := project_id]
  wide[, `Case ID` := case_id]
  wide[, `File ID` := file_id]
  wide[, RPPA_File := basename(rppa_path)]

  setcolorder(
    wide,
    c(
      "Project ID",
      "Case ID",
      "File ID",
      "RPPA_File",
      setdiff(colnames(wide), c("Project ID", "Case ID", "File ID", "RPPA_File"))
    )
  )

  return(wide)
}

all_case_rows <- list()

project_ids <- sort(unique(sample_sheet$`Project ID`))

cat("[INFO] Found", length(project_ids), "unique projects\n")

for (project_id in project_ids) {
  cat("\n[INFO] Processing project:", project_id, "\n")

  project_sheet <- sample_sheet[`Project ID` == project_id]

  project_rows <- list()

  for (i in seq_len(nrow(project_sheet))) {
    file_id <- project_sheet$`File ID`[i]
    case_id <- project_sheet$`Case ID`[i]

    case_dt <- read_rppa_for_file_id(
      file_id = file_id,
      case_id = case_id,
      project_id = project_id
    )

    if (!is.null(case_dt)) {
      project_rows[[length(project_rows) + 1]] <- case_dt
      all_case_rows[[length(all_case_rows) + 1]] <- case_dt
    }
  }

  if (length(project_rows) == 0) {
    warning("[SKIP] No readable RPPA rows for project: ", project_id)
    next
  }

  project_out <- rbindlist(
    project_rows,
    use.names = TRUE,
    fill = TRUE
  )

  project_out <- project_out[
    ,
    lapply(.SD, function(x) {
      if (is.numeric(x)) {
        mean(x, na.rm = TRUE)
      } else {
        paste(sort(unique(na.omit(x))), collapse = ";")
      }
    }),
    by = .(`Project ID`, `Case ID`)
  ]

  setorder(project_out, `Project ID`, `Case ID`)

  out_file <- file.path(
    output_dir,
    paste0(project_id, "_RPPA_protein_expression.tsv")
  )

  fwrite(
    project_out,
    file = out_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  cat("[OK] Wrote:", out_file, "\n")
  cat("[OK]", project_id, "cases:", uniqueN(project_out$`Case ID`), "\n")
}

if (length(all_case_rows) == 0) {
  stop("No readable RPPA data found.")
}

combined_out <- rbindlist(
  all_case_rows,
  use.names = TRUE,
  fill = TRUE
)

combined_out <- combined_out[
  ,
  lapply(.SD, function(x) {
    if (is.numeric(x)) {
      mean(x, na.rm = TRUE)
    } else {
      paste(sort(unique(na.omit(x))), collapse = ";")
    }
  }),
  by = .(`Project ID`, `Case ID`)
]

setorder(combined_out, `Project ID`, `Case ID`)

combined_file <- file.path(
  output_dir,
  "ALL_TCGA_RPPA_protein_expression.tsv"
)

fwrite(
  combined_out,
  file = combined_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("\n[OK] Wrote combined file:", combined_file, "\n")
cat("[OK] Total projects:", uniqueN(combined_out$`Project ID`), "\n")
cat("[OK] Total cases:", uniqueN(combined_out$`Case ID`), "\n")
cat("[DONE] RPPA preprocessing complete.\n")
