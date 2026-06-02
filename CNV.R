#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(GenomicFeatures)
  library(AnnotationDbi)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
})

wgs_dir <- "/fs/scratch/PAS2942/Alejandro/slides/WGS"
sample_sheet_path <- file.path(wgs_dir, "WGSam.tsv")
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
    "WGSam.tsv is missing required columns: ",
    paste(missing_sample_cols, collapse = ", ")
  )
}

sample_sheet <- unique(sample_sheet[, .(`File ID`, `Case ID`, `Project ID`)])

cat("[INFO] Loaded sample sheet with", nrow(sample_sheet), "file mappings\n")

cat("[INFO] Loading hg38 gene annotation\n")

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
genes_gr <- genes(txdb)

gene_ids <- as.character(mcols(genes_gr)$gene_id)

gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = gene_ids,
  column = "SYMBOL",
  keytype = "ENTREZID",
  multiVals = "first"
)

genes_gr$Gene <- as.character(gene_symbols)

genes_gr <- genes_gr[
  !is.na(genes_gr$Gene) &
    genes_gr$Gene != "" &
    as.character(seqnames(genes_gr)) %in% paste0("chr", c(1:22, "X", "Y"))
]

genes_gr <- sort(genes_gr)

cat("[INFO] Annotated genes retained:", length(genes_gr), "\n")

standardize_chromosome <- function(x) {
  x <- as.character(x)
  x <- gsub("^chr", "", x, ignore.case = TRUE)
  x <- gsub("^23$", "X", x)
  x <- gsub("^24$", "Y", x)
  x <- paste0("chr", x)
  return(x)
}

find_segment_file <- function(file_id) {
  file_dir <- file.path(wgs_dir, file_id)

  if (!dir.exists(file_dir)) {
    warning("[SKIP] Directory not found for File ID: ", file_id)
    return(NA_character_)
  }

  seg_files <- list.files(
    file_dir,
    pattern = "\\.(txt|tsv|seg|csv|gz)$",
    full.names = TRUE
  )

  seg_files <- seg_files[
    !grepl("WGSam\\.tsv$", basename(seg_files), ignore.case = TRUE)
  ]

  if (length(seg_files) == 0) {
    warning("[SKIP] No segment file found in directory: ", file_dir)
    return(NA_character_)
  }

  if (length(seg_files) > 1) {
    warning(
      "[INFO] Multiple segment files found for File ID ",
      file_id,
      ". Using first file: ",
      basename(seg_files[1])
    )
  }

  return(seg_files[1])
}

process_one_case <- function(file_id, case_id, project_id) {
  seg_path <- find_segment_file(file_id)

  if (is.na(seg_path)) {
    return(NULL)
  }

  cat("[INFO] Reading:", seg_path, "\n")

  seg_dt <- tryCatch(
    fread(
      seg_path,
      header = TRUE,
      data.table = TRUE,
      fill = TRUE,
      quote = ""
    ),
    error = function(e) {
      warning("[SKIP] Failed to read segment file for File ID ", file_id, ": ", e$message)
      return(NULL)
    }
  )

  if (is.null(seg_dt)) {
    return(NULL)
  }

  required_seg_cols <- c("Chromosome", "Start", "End", "Segment_Mean")
  missing_seg_cols <- setdiff(required_seg_cols, colnames(seg_dt))

  if (length(missing_seg_cols) > 0) {
    warning(
      "[SKIP] Segment file missing required columns for File ID ",
      file_id,
      ": ",
      paste(missing_seg_cols, collapse = ", ")
    )
    return(NULL)
  }

  seg_dt <- seg_dt[
    !is.na(Chromosome) &
      !is.na(Start) &
      !is.na(End) &
      !is.na(Segment_Mean)
  ]

  if (nrow(seg_dt) == 0) {
    warning("[SKIP] No usable segments for File ID: ", file_id)
    return(NULL)
  }

  seg_dt[, Chromosome := standardize_chromosome(Chromosome)]
  seg_dt[, Start := as.integer(Start)]
  seg_dt[, End := as.integer(End)]
  seg_dt[, Segment_Mean := as.numeric(Segment_Mean)]

  seg_dt <- seg_dt[
    Chromosome %in% paste0("chr", c(1:22, "X", "Y")) &
      !is.na(Start) &
      !is.na(End) &
      !is.na(Segment_Mean) &
      End >= Start
  ]

  if (nrow(seg_dt) == 0) {
    warning("[SKIP] No valid autosome/sex-chromosome segments for File ID: ", file_id)
    return(NULL)
  }

  seg_gr <- makeGRangesFromDataFrame(
    seg_dt,
    seqnames.field = "Chromosome",
    start.field = "Start",
    end.field = "End",
    keep.extra.columns = TRUE
  )

  seqlevelsStyle(seg_gr) <- "UCSC"

  overlap_hits <- findOverlaps(
    genes_gr,
    seg_gr,
    ignore.strand = TRUE
  )

  if (length(overlap_hits) == 0) {
    warning("[SKIP] No gene-segment overlaps for File ID: ", file_id)
    return(NULL)
  }

  gene_hit <- genes_gr[queryHits(overlap_hits)]
  seg_hit <- seg_gr[subjectHits(overlap_hits)]

  overlap_width <- width(pintersect(gene_hit, seg_hit))

  gene_cnv_dt <- data.table(
    `Project ID` = project_id,
    `Case ID` = case_id,
    `File ID` = file_id,
    Gene = gene_hit$Gene,
    Segment_Mean = seg_hit$Segment_Mean,
    Overlap_Width = overlap_width
  )

  gene_cnv_dt <- gene_cnv_dt[
    !is.na(Gene) &
      Gene != "" &
      !is.na(Segment_Mean) &
      !is.na(Overlap_Width) &
      Overlap_Width > 0
  ]

  if (nrow(gene_cnv_dt) == 0) {
    warning("[SKIP] No usable gene CNV rows for File ID: ", file_id)
    return(NULL)
  }

  gene_level <- gene_cnv_dt[
    ,
    .(
      Segment_Mean_Gene_Weighted =
        sum(Segment_Mean * Overlap_Width, na.rm = TRUE) /
        sum(Overlap_Width, na.rm = TRUE),
      Gene_Overlap_Bases = sum(Overlap_Width, na.rm = TRUE)
    ),
    by = .(`Project ID`, `Case ID`, `File ID`, Gene)
  ]

  gene_level[
    ,
    Estimated_Absolute_Copy_Number := 2 * (2 ^ Segment_Mean_Gene_Weighted)
  ]

  gene_level <- gene_level[
    ,
    .(
      `Project ID`,
      `Case ID`,
      `File ID`,
      Gene,
      Segment_Mean_Gene_Weighted,
      Estimated_Absolute_Copy_Number,
      Gene_Overlap_Bases
    )
  ]

  return(gene_level)
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

    case_dt <- process_one_case(
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
    warning("[SKIP] No readable WGS CNV rows for project: ", project_id)
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
      Segment_Mean_Gene_Weighted = mean(Segment_Mean_Gene_Weighted, na.rm = TRUE),
      Estimated_Absolute_Copy_Number = mean(Estimated_Absolute_Copy_Number, na.rm = TRUE),
      Gene_Overlap_Bases = sum(Gene_Overlap_Bases, na.rm = TRUE),
      File_IDs = paste(sort(unique(`File ID`)), collapse = ";")
    ),
    by = .(`Project ID`, `Case ID`, Gene)
  ]

  setorder(project_long, `Project ID`, `Case ID`, Gene)

  project_long_file <- file.path(
    output_dir,
    paste0(project_id, "_gene_absolute_CNV_long.tsv")
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
    value.var = "Estimated_Absolute_Copy_Number"
  )

  setorder(project_matrix, Gene)

  project_matrix_file <- file.path(
    output_dir,
    paste0(project_id, "_gene_absolute_CNV_matrix.tsv")
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

if (length(all_case_rows) == 0) {
  stop("No readable WGS CNV data found.")
}

combined_long <- rbindlist(
  all_case_rows,
  use.names = TRUE,
  fill = TRUE
)

combined_long <- combined_long[
  ,
  .(
    Segment_Mean_Gene_Weighted = mean(Segment_Mean_Gene_Weighted, na.rm = TRUE),
    Estimated_Absolute_Copy_Number = mean(Estimated_Absolute_Copy_Number, na.rm = TRUE),
    Gene_Overlap_Bases = sum(Gene_Overlap_Bases, na.rm = TRUE),
    File_IDs = paste(sort(unique(`File ID`)), collapse = ";")
  ),
  by = .(`Project ID`, `Case ID`, Gene)
]

setorder(combined_long, `Project ID`, `Case ID`, Gene)

combined_long_file <- file.path(
  output_dir,
  "ALL_TCGA_gene_absolute_CNV_long.tsv"
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
  value.var = "Estimated_Absolute_Copy_Number"
)

setorder(combined_matrix, Gene)

combined_matrix_file <- file.path(
  output_dir,
  "ALL_TCGA_gene_absolute_CNV_matrix.tsv"
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
cat("[DONE] WGS gene-level estimated absolute CNV preprocessing complete.\n")
