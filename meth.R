#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

methyl_dir <- "/fs/scratch/PAS2942/Alejandro/slides/methyl"
sample_sheet_path <- file.path(methyl_dir, "methylsam.tsv")
annotation_path <- "/fs/scratch/PAS2942/Alejandro/slides/meth450k.csv"
output_dir <- "/fs/scratch/PAS2942/Alejandro/slides"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(sample_sheet_path)) {
  stop("Sample sheet not found: ", sample_sheet_path)
}

if (!file.exists(annotation_path)) {
  stop("450K annotation file not found: ", annotation_path)
}

sample_sheet <- fread(sample_sheet_path)

required_sample_cols <- c("File ID", "Case ID", "Project ID")
missing_sample_cols <- setdiff(required_sample_cols, colnames(sample_sheet))

if (length(missing_sample_cols) > 0) {
  stop(
    "methylsam.tsv is missing required columns: ",
    paste(missing_sample_cols, collapse = ", ")
  )
}

sample_sheet <- unique(sample_sheet[, .(`File ID`, `Case ID`, `Project ID`)])

cat("[INFO] Loaded sample sheet with", nrow(sample_sheet), "file mappings\n")

cat("[INFO] Reading 450K annotation:", annotation_path, "\n")
annotation_lines <- readLines(annotation_path, warn = FALSE)

header_line <- grep(
  "^IlmnID,",
  annotation_lines
)[1]

if (is.na(header_line)) {
  stop("Could not find Illumina annotation header line starting with 'IlmnID,' in: ", annotation_path)
}

cat("[INFO] Illumina annotation header found at line:", header_line, "\n")

annotation_dt <- fread(
  annotation_path,
  skip = header_line - 1,
  header = TRUE,
  data.table = TRUE,
  fill = TRUE,
  quote = ""
)
required_annotation_cols <- c(
  "IlmnID",
  "UCSC_RefGene_Name",
  "UCSC_RefGene_Group"
)

missing_annotation_cols <- setdiff(required_annotation_cols, colnames(annotation_dt))

if (length(missing_annotation_cols) > 0) {
  stop(
    "meth450k.csv is missing required annotation columns: ",
    paste(missing_annotation_cols, collapse = ", ")
  )
}

promoter_groups <- c("TSS1500", "TSS200", "5'UTR", "1stExon")

annotation_dt <- annotation_dt[
  !is.na(IlmnID) &
    IlmnID != "" &
    !is.na(UCSC_RefGene_Name) &
    UCSC_RefGene_Name != "" &
    !is.na(UCSC_RefGene_Group) &
    UCSC_RefGene_Group != ""
]

probe_gene_rows <- list()

for (i in seq_len(nrow(annotation_dt))) {
  probe_id <- annotation_dt$IlmnID[i]

  gene_names <- unlist(strsplit(annotation_dt$UCSC_RefGene_Name[i], ";", fixed = TRUE))
  gene_groups <- unlist(strsplit(annotation_dt$UCSC_RefGene_Group[i], ";", fixed = TRUE))

  gene_names <- trimws(gene_names)
  gene_groups <- trimws(gene_groups)

  n_pair <- min(length(gene_names), length(gene_groups))

  if (n_pair == 0) {
    next
  }

  tmp <- data.table(
    Probe_ID = probe_id,
    Gene = gene_names[seq_len(n_pair)],
    Gene_Group = gene_groups[seq_len(n_pair)]
  )

  tmp <- tmp[
    Gene != "" &
      Gene_Group %in% promoter_groups
  ]

  if (nrow(tmp) > 0) {
    probe_gene_rows[[length(probe_gene_rows) + 1]] <- tmp
  }
}

probe_gene_map <- rbindlist(
  probe_gene_rows,
  use.names = TRUE,
  fill = TRUE
)

probe_gene_map <- unique(probe_gene_map[, .(Probe_ID, Gene)])

cat("[INFO] Promoter probe-gene mappings:", nrow(probe_gene_map), "\n")
cat("[INFO] Unique promoter probes:", uniqueN(probe_gene_map$Probe_ID), "\n")
cat("[INFO] Unique promoter genes:", uniqueN(probe_gene_map$Gene), "\n")

if (nrow(probe_gene_map) == 0) {
  stop("No promoter probe-gene mappings found.")
}

find_methyl_file <- function(file_id) {
  file_dir <- file.path(methyl_dir, file_id)

  if (!dir.exists(file_dir)) {
    warning("[SKIP] Directory not found for File ID: ", file_id)
    return(NA_character_)
  }

  methyl_files <- list.files(
    file_dir,
    pattern = "\\.(txt|tsv|csv|gz)$",
    full.names = TRUE
  )

  methyl_files <- methyl_files[
    !grepl("methylsam\\.tsv$", basename(methyl_files), ignore.case = TRUE)
  ]

  if (length(methyl_files) == 0) {
    warning("[SKIP] No methylation beta file found in directory: ", file_dir)
    return(NA_character_)
  }

  if (length(methyl_files) > 1) {
    warning(
      "[INFO] Multiple methylation files found for File ID ",
      file_id,
      ". Using first file: ",
      basename(methyl_files[1])
    )
  }

  return(methyl_files[1])
}

read_methyl_for_file_id <- function(file_id, case_id, project_id) {
  methyl_path <- find_methyl_file(file_id)

  if (is.na(methyl_path)) {
    return(NULL)
  }

  cat("[INFO] Reading:", methyl_path, "\n")

  dt <- tryCatch(
    fread(
      methyl_path,
      header = TRUE,
      data.table = TRUE,
      fill = TRUE,
      quote = ""
    ),
    error = function(e) {
      warning("[SKIP] Failed to read methylation file for File ID ", file_id, ": ", e$message)
      return(NULL)
    }
  )

  if (is.null(dt)) {
    return(NULL)
  }

  if (ncol(dt) < 2) {
    warning("[SKIP] Methylation file has fewer than 2 columns for File ID: ", file_id)
    return(NULL)
  }

  probe_col <- colnames(dt)[1]
  beta_col <- colnames(dt)[2]

  setnames(dt, old = c(probe_col, beta_col), new = c("Probe_ID", "Beta_Value"))

  dt <- dt[
    !is.na(Probe_ID) &
      Probe_ID != "" &
      !is.na(Beta_Value)
  ]

  if (nrow(dt) == 0) {
    warning("[SKIP] No usable beta rows for File ID: ", file_id)
    return(NULL)
  }

  dt[, Beta_Value := as.numeric(Beta_Value)]

  dt <- dt[
    !is.na(Beta_Value) &
      Beta_Value >= 0 &
      Beta_Value <= 1
  ]

  if (nrow(dt) == 0) {
    warning("[SKIP] No numeric beta values in [0,1] for File ID: ", file_id)
    return(NULL)
  }

  dt <- merge(
    dt,
    probe_gene_map,
    by = "Probe_ID",
    all = FALSE,
    allow.cartesian = TRUE
  )

  if (nrow(dt) == 0) {
    warning("[SKIP] No promoter probe overlaps for File ID: ", file_id)
    return(NULL)
  }

  gene_level <- dt[
    ,
    .(
      Promoter_Beta_Mean = mean(Beta_Value, na.rm = TRUE),
      Promoter_CpG_Count = uniqueN(Probe_ID)
    ),
    by = Gene
  ]

  gene_level[, `Project ID` := project_id]
  gene_level[, `Case ID` := case_id]
  gene_level[, `File ID` := file_id]
  gene_level[, Methylation_File := basename(methyl_path)]

  gene_level <- gene_level[
    ,
    .(
      `Project ID`,
      `Case ID`,
      `File ID`,
      Methylation_File,
      Gene,
      Promoter_Beta_Mean,
      Promoter_CpG_Count
    )
  ]

  return(gene_level)
}

all_long_rows <- list()

project_ids <- sort(unique(sample_sheet$`Project ID`))

cat("[INFO] Found", length(project_ids), "unique projects\n")

for (project_id in project_ids) {
  cat("\n[INFO] Processing project:", project_id, "\n")

  project_sheet <- sample_sheet[`Project ID` == project_id]

  project_rows <- list()

  for (i in seq_len(nrow(project_sheet))) {
    file_id <- project_sheet$`File ID`[i]
    case_id <- project_sheet$`Case ID`[i]

    case_dt <- read_methyl_for_file_id(
      file_id = file_id,
      case_id = case_id,
      project_id = project_id
    )

    if (!is.null(case_dt)) {
      project_rows[[length(project_rows) + 1]] <- case_dt
      all_long_rows[[length(all_long_rows) + 1]] <- case_dt
    }
  }

  if (length(project_rows) == 0) {
    warning("[SKIP] No readable methylation rows for project: ", project_id)
    next
  }

  project_long <- rbindlist(
    project_rows,
    use.names = TRUE,
    fill = TRUE
  )

  project_long <- project_long[
    ,
    .(
      Promoter_Beta_Mean = mean(Promoter_Beta_Mean, na.rm = TRUE),
      Promoter_CpG_Count = sum(Promoter_CpG_Count, na.rm = TRUE),
      File_IDs = paste(sort(unique(`File ID`)), collapse = ";")
    ),
    by = .(`Project ID`, `Case ID`, Gene)
  ]

  setorder(project_long, `Project ID`, `Case ID`, Gene)

  project_long_file <- file.path(
    output_dir,
    paste0(project_id, "_promoter_methylation_beta_long.tsv")
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
    Gene ~ `Case ID`,
    value.var = "Promoter_Beta_Mean"
  )

  setorder(project_matrix, Gene)

  project_matrix_file <- file.path(
    output_dir,
    paste0(project_id, "_promoter_methylation_beta_matrix.tsv")
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
  cat("[OK]", project_id, "genes:", nrow(project_matrix), "\n")
  cat("[OK]", project_id, "cases:", ncol(project_matrix) - 1, "\n")
}

if (length(all_long_rows) == 0) {
  stop("No readable methylation data found.")
}

combined_long <- rbindlist(
  all_long_rows,
  use.names = TRUE,
  fill = TRUE
)

combined_long <- combined_long[
  ,
  .(
    Promoter_Beta_Mean = mean(Promoter_Beta_Mean, na.rm = TRUE),
    Promoter_CpG_Count = sum(Promoter_CpG_Count, na.rm = TRUE),
    File_IDs = paste(sort(unique(`File ID`)), collapse = ";")
  ),
  by = .(`Project ID`, `Case ID`, Gene)
]

setorder(combined_long, `Project ID`, `Case ID`, Gene)

combined_long_file <- file.path(
  output_dir,
  "ALL_TCGA_promoter_methylation_beta_long.tsv"
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
  Gene ~ `Case ID`,
  value.var = "Promoter_Beta_Mean"
)

setorder(combined_matrix, Gene)

combined_matrix_file <- file.path(
  output_dir,
  "ALL_TCGA_promoter_methylation_beta_matrix.tsv"
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
cat("[OK] Combined genes:", nrow(combined_matrix), "\n")
cat("[OK] Combined cases:", ncol(combined_matrix) - 1, "\n")
cat("[DONE] Promoter methylation beta preprocessing complete.\n")
