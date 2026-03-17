#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(tools)
})

is_absolute_path <- function(path) {
  grepl("^(?:[A-Za-z]:[/\\\\]|/|\\\\\\\\)", path)
}

script_path <- function() {
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg) > 0) {
    script_file <- sub("^--file=", "", file_arg[[1]])
    return(normalizePath(dirname(script_file), winslash = "/", mustWork = FALSE))
  }

  args_full <- commandArgs(trailingOnly = FALSE)
  flag_index <- which(args_full == "-f")
  if (length(flag_index) > 0) {
    script_file <- args_full[flag_index[[1]] + 1]
    return(normalizePath(dirname(script_file), winslash = "/", mustWork = FALSE))
  }

  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

kegni_root <- script_path()
project_root <- dirname(kegni_root)
repo_root <- dirname(project_root)

resolve_path <- function(path) {
  if (is_absolute_path(path)) {
    return(normalizePath(path, winslash = "/", mustWork = FALSE))
  }
  normalizePath(file.path(kegni_root, path), winslash = "/", mustWork = FALSE)
}

# ============================================================
# extract_kegni_tf_targets.R
#
# Reconstruct KEGNI's native edge scores from saved scg embeddings
# using Trainer.recon(self, z):
#
#   z <- tanh(z)
#   if norm == -1:
#       pred <- z %*% t(z)
#   else:
#       pred <- normalized similarity with p = norm
#
# Then:
#   - remove self edges
#   - sort descending by score
#   - filter to TF -> target edges
#
# Inputs:
#   --inputs_root   Root folder containing per-cell-type TF_list.csv files
#                   default: kegni_inputs
#   --outputs_root  Folder containing KEGNI output embedding files
#                   default: outputs
#   --out_root      Folder for extracted TF-target edge lists
#                   default: tf_target_edges
#   --norm          KEGNI recon norm parameter
#                   default: -1
#
# Expected KEGNI output files per cell type:
#   <cell_type>-scg_embedding.csv
#   or any file matching *<cell_type>*scg_embedding.csv
#
# Output files:
#   tf_target_edges/<cell_type>_tf_target_ranked.csv
#   tf_target_edges/all_celltypes_tf_target_ranked.csv
#   tf_target_edges/extraction_summary.csv
# ============================================================



default_meta <- file.path(project_root, "resources", "CHOOSE-multiome-wt-metadata.csv")
default_inputs_root  <- file.path(kegni_root, "kegni_inputs")
default_outputs_root <- file.path(kegni_root, "outputs")
default_out_root     <- file.path(kegni_root, "tf_target_edges")
default_norm         <- -1

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) == 0) return(default)
  if (idx[length(idx)] == length(args)) {
    stop("Missing value for argument: ", flag)
  }
  args[idx[length(idx)] + 1]
}

meta_file    <- resolve_path(get_arg("--meta",         default_meta))
inputs_root  <- resolve_path(get_arg("--inputs_root",  default_inputs_root))
outputs_root <- resolve_path(get_arg("--outputs_root", default_outputs_root))
out_root     <- resolve_path(get_arg("--out_root",     default_out_root))
norm_value   <- as.numeric(get_arg("--norm", default_norm))

required_files <- c(meta_file)
missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    "The following required input file(s) do not exist:\n",
    paste("  -", missing_files, collapse = "\n")
  )
}

if (!dir.exists(inputs_root)) {
  stop("The inputs root directory does not exist: ", inputs_root)
}

if (!dir.exists(outputs_root)) {
  stop(
    "The outputs root directory does not exist: ", outputs_root,
    "\nPass the KEGNI outputs folder via --outputs_root (for example, KEGNI/outputs)."
  )
}

available_scg <- list.files(
  outputs_root,
  recursive = TRUE,
  full.names = TRUE,
  pattern = "scg_embedding\\.csv$",
  ignore.case = TRUE
)

if (length(available_scg) == 0) {
  stop(
    "No scg embedding files were found under outputs_root: ", outputs_root,
    "\nExpected files matching *scg_embedding.csv."
  )
}

meta_df <- read.csv(meta_file, row.names=1, check.names = FALSE, stringsAsFactors = FALSE)
if (!("celltype_jf" %in% colnames(meta_df))) {
  stop("Metadata file must contain a column named 'celltype_jf'.")
}

cell_types_arg <- get_arg("--cell_types", NULL)
cell_types <- if (is.null(cell_types_arg)) {
  unique(as.character(meta_df$celltype_jf))
} else {
  trimws(unlist(strsplit(cell_types_arg, ",", fixed = TRUE)))
}
cell_types <- cell_types[nzchar(cell_types) & !is.na(cell_types)]

if (length(cell_types) == 0) {
  stop("No cell types were found to process.")
}

dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

message2 <- function(...) cat(..., "\n", sep = "")

read_tf_list <- function(path) {
  if (!file.exists(path)) {
    stop("TF list not found: ", path)
  }

  tf_df <- read_csv(path, show_col_types = FALSE)
  tf_col <- if ("TF" %in% names(tf_df)) "TF" else names(tf_df)[1]

  tf_df %>%
    transmute(TF = toupper(as.character(.data[[tf_col]]))) %>%
    filter(!is.na(TF), TF != "") %>%
    distinct()
}

find_scg_embedding_file <- function(cell_type, outputs_root) {
  files <- list.files(outputs_root, recursive = TRUE, full.names = TRUE)

  candidates <- files[
    grepl("scg_embedding\\.csv$", basename(files), ignore.case = TRUE)
  ]

  if (length(candidates) == 0) return(NULL)

  # Prefer files that mention the cell type in the path or file name
  hits <- candidates[grepl(cell_type, candidates, ignore.case = TRUE)]
  if (length(hits) > 0) return(hits[1])

  # Fallback: if there is exactly one scg embedding, use it
  if (length(candidates) == 1) return(candidates[1])

  NULL
}

read_scg_embedding <- function(path) {
  df <- read_csv(path, show_col_types = FALSE)

  first_col <- names(df)[1]
  gene_names <- toupper(as.character(df[[first_col]]))

  mat <- as.matrix(df[, -1, drop = FALSE])
  storage.mode(mat) <- "numeric"

  rownames(mat) <- gene_names
  mat
}

kegni_recon_scores <- function(z, norm_value = -1) {
  z <- tanh(z)

  if (norm_value == -1) {
    pred <- z %*% t(z)
  } else {
    row_norms <- apply(z, 1, function(x) sum(abs(x)^norm_value)^(1 / norm_value))
    row_norms[row_norms == 0] <- 1
    z_norm <- z / row_norms
    pred <- z_norm %*% t(z_norm)
  }

  pred
}

matrix_to_ranked_edges <- function(score_mat, tf_set, cell_type) {
  df <- as.data.frame(score_mat, check.names = FALSE) %>%
    rownames_to_column("Gene1") %>%
    pivot_longer(cols = -Gene1, names_to = "Gene2", values_to = "EdgeWeight") %>%
    mutate(
      Gene1 = toupper(as.character(Gene1)),
      Gene2 = toupper(as.character(Gene2)),
      EdgeWeight = as.numeric(EdgeWeight)
    ) %>%
    filter(Gene1 != Gene2) %>%
    filter(Gene1 %in% tf_set) %>%
    arrange(desc(EdgeWeight), Gene1, Gene2) %>%
    mutate(
      cell_type = cell_type,
      rank = row_number()
    ) %>%
    select(cell_type, rank, Gene1, Gene2, EdgeWeight)

  df
}

all_edges <- list()
summary_rows <- list()

for (cell_type in cell_types) {
  message2("Processing ", cell_type, "...")

  tf_path <- file.path(inputs_root, cell_type, "TF_list.csv")
  if (!file.exists(tf_path)) {
    warning("Missing TF list for ", cell_type, ": ", tf_path)
    next
  }

  emb_path <- find_scg_embedding_file(cell_type, outputs_root)
  if (is.null(emb_path)) {
    warning("No scg_embedding file found for ", cell_type, " under ", outputs_root)
    summary_rows[[cell_type]] <- tibble(
      cell_type = cell_type,
      scg_embedding_file = NA_character_,
      norm = norm_value,
      n_genes = NA_integer_,
      n_tfs = NA_integer_,
      n_ranked_tf_target_edges = 0L,
      status = "missing_scg_embedding"
    )
    next
  }

  tf_df <- read_tf_list(tf_path)
  tf_set <- tf_df$TF

  z <- read_scg_embedding(emb_path)
  score_mat <- kegni_recon_scores(z, norm_value = norm_value)

  edges <- matrix_to_ranked_edges(score_mat, tf_set = tf_set, cell_type = cell_type)

  out_file <- file.path(out_root, paste0(cell_type, "_tf_target_ranked.csv"))
  write_csv(edges, out_file)

  message2("  embedding: ", emb_path)
  message2("  wrote:     ", out_file, " (", nrow(edges), " edges)")

  all_edges[[cell_type]] <- edges
  summary_rows[[cell_type]] <- tibble(
    cell_type = cell_type,
    scg_embedding_file = emb_path,
    norm = norm_value,
    n_genes = nrow(z),
    n_tfs = length(tf_set),
    n_ranked_tf_target_edges = nrow(edges),
    status = "ok"
  )
}

summary_df <- bind_rows(summary_rows)
write_csv(summary_df, file.path(out_root, "extraction_summary.csv"))

if (length(all_edges) > 0) {
  combined <- bind_rows(all_edges) %>%
    arrange(cell_type, rank)
  write_csv(combined, file.path(out_root, "all_celltypes_tf_target_ranked.csv"))
}

message2("")
message2("Done.")
message2("Summary: ", file.path(out_root, "extraction_summary.csv"))
if (length(all_edges) > 0) {
  message2("Combined: ", file.path(out_root, "all_celltypes_tf_target_ranked.csv"))
}