#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tibble)
  library(stringr)
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
# make_kegni_inputs.R
#
# Purpose:
#   Build KEGNI-ready input folders for the five cell types:
#     ge_in, ctx_ex, ge_npc, ctx_npc, ctx_ip
#
# Inputs:
#   1) Global expression matrix (genes x cells; first column = gene)
#   2) Global metadata table   (first column = cell barcode)
#   3) KEGG KG TSV            (gene1, relation, gene2)
#   4) TRRUST KG TSV          (gene1, relation, gene2)
#   5) TF list CSV            (first column = TF gene symbol, or a column named TF)
#
# Outputs:
#   kegni_inputs/<cell_type>/expression.csv.gz
#   kegni_inputs/<cell_type>/metadata.csv.gz
#   kegni_inputs/<cell_type>/merged_kg.tsv
#   kegni_inputs/<cell_type>/TF_list.csv
#   kegni_inputs/summary.csv
#
# Usage:
#   Rscript make_kegni_inputs.R \
#     --expr ../resources/CHOOSE-multiome-wt-log-norm.csv.gz \
#     --meta ../resources/CHOOSE-multiome-wt-metadata.csv \
#     --kegg ./inputs/KG/KEGG_hs.tsv \
#     --trrust ./inputs/KG/TRRUST_hs.tsv \
#     --tf ../resources/TF_list.csv \
#     --out kegni_inputs
# ============================================================

# -----------------------------
# Defaults
# -----------------------------
default_expr   <- file.path(project_root, "resources", "CHOOSE-multiome-wt-log-norm.csv.gz")
default_meta   <- file.path(project_root, "resources", "CHOOSE-multiome-wt-metadata.csv")
default_kegg   <- file.path(kegni_root, "inputs", "KG", "KEGG_hs.tsv")
default_trrust <- file.path(kegni_root, "inputs", "KG", "TRRUST_hs.tsv")
default_tf     <- file.path(project_root, "resources", "TF_list.csv")
default_out    <- file.path(kegni_root, "kegni_inputs")

# -----------------------------
# Parse command-line args
# -----------------------------
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) == 0) return(default)
  if (idx[length(idx)] == length(args)) {
    stop("Missing value for argument: ", flag)
  }
  args[idx[length(idx)] + 1]
}

expr_file   <- resolve_path(get_arg("--expr",   default_expr))
meta_file   <- resolve_path(get_arg("--meta",   default_meta))
kegg_file   <- resolve_path(get_arg("--kegg",   default_kegg))
trrust_file <- resolve_path(get_arg("--trrust", default_trrust))
tf_file     <- resolve_path(get_arg("--tf",     default_tf))
out_root    <- resolve_path(get_arg("--out",    default_out))

cell_types_arg <- get_arg("--cell_types", NULL)

# -----------------------------
# Validate files
# -----------------------------
required_files <- c(expr_file, meta_file, kegg_file, trrust_file, tf_file)
missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    "The following required input file(s) do not exist:\n",
    paste("  -", missing_files, collapse = "\n")
  )
}

dir.create(out_root, showWarnings = FALSE, recursive = TRUE)

message("Reading expression matrix: ", expr_file)
expr_mat <- read.csv(expr_file, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

# Normalize gene symbols in expression
rownames(expr_mat) <- toupper(rownames(expr_mat))

# Collapse duplicated gene rows if they exist after uppercasing
if (anyDuplicated(rownames(expr_mat)) > 0) {
  message("Detected duplicated gene symbols after uppercasing; collapsing by sum.")
  expr_mat <- expr_mat %>%
    rownames_to_column("gene") %>%
    group_by(gene) %>%
    summarise(across(everything(), ~ sum(.x, na.rm = TRUE)), .groups = "drop") %>%
    column_to_rownames("gene") %>%
    as.data.frame(check.names = FALSE)
}

message("Reading metadata: ", meta_file)
meta_df <- read.csv(meta_file, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
# Organize metadata rows so they're in the same order as expression columns (after intersecting shared cells)
# meta_df <- meta_df[colnames(expr_mat), , drop = FALSE]


# print(meta_df[1:5, 106:110])

if (!("celltype_jf" %in% colnames(meta_df))) {
  stop("metadata file must contain a column named 'celltype_jf'.")
}

# Normalize once so comparisons are robust across R versions and input quirks.
meta_df$celltype_jf <- trimws(as.character(meta_df$celltype_jf))

if (is.null(cell_types_arg)) {
  cell_types <- unique(meta_df$celltype_jf)
} else {
  cell_types <- trimws(unlist(strsplit(cell_types_arg, ",", fixed = TRUE)))
}
cell_types <- unique(cell_types[nzchar(cell_types) & !is.na(cell_types)])

if (length(cell_types) == 0) {
  stop("No cell types were found to process.")
}

# Keep only shared cells and align metadata to expression columns
shared_cells <- intersect(colnames(expr_mat), rownames(meta_df))
if (length(shared_cells) == 0) {
  stop("No overlapping cell barcodes were found between expression columns and metadata row names.")
}

expr_mat <- expr_mat[, shared_cells, drop = FALSE]
meta_df  <- meta_df[shared_cells, , drop = FALSE]

message("Shared cells retained: ", length(shared_cells))

# -----------------------------
# Read TF list
# -----------------------------
message("Reading TF list: ", tf_file)
tf_df <- read_csv(tf_file, show_col_types = FALSE)

if ("TF" %in% colnames(tf_df)) {
  tf_col <- "TF"
} else {
  tf_col <- colnames(tf_df)[1]
}

tf_df <- tf_df %>%
  transmute(TF = toupper(as.character(.data[[tf_col]]))) %>%
  filter(!is.na(TF), TF != "") %>%
  distinct()

tf_set <- unique(tf_df$TF)

# -----------------------------
# Read KG files
# -----------------------------
normalize_relation <- function(x) {
  x_low <- str_to_lower(as.character(x))
  case_when(
    x_low == "activation" ~ "activation",
    x_low == "inhibition" ~ "inhibition",
    x_low == "repression" ~ "repression",
    x_low == "effect"     ~ "effect",
    TRUE                  ~ as.character(x)
  )
}

read_kg <- function(path) {
  read_tsv(
    path,
    col_names = c("gene1", "relation", "gene2"),
    show_col_types = FALSE
  ) %>%
    mutate(
      gene1 = toupper(as.character(gene1)),
      gene2 = toupper(as.character(gene2)),
      relation = normalize_relation(relation)
    ) %>%
    filter(
      !is.na(gene1), gene1 != "",
      !is.na(gene2), gene2 != "",
      !is.na(relation), relation != ""
    ) %>%
    distinct()
}

message("Reading KG files...")
kegg   <- read_kg(kegg_file)
trrust <- read_kg(trrust_file)
merged <- bind_rows(kegg, trrust) %>% distinct()

summary_list <- list()

# -----------------------------
# Build cell-type-specific inputs
# -----------------------------
for (celltype in cell_types) {
  message("Processing cell type: ", celltype)

  celltype_dir <- file.path(out_root, celltype)
  dir.create(celltype_dir, showWarnings = FALSE, recursive = TRUE)

  celltype_cells <- rownames(meta_df)[meta_df$celltype_jf == celltype]

  if (length(celltype_cells) == 0) {
    warning("No cells found for cell type: ", celltype)
    summary_list[[celltype]] <- tibble(
      celltype = celltype,
      n_cells = 0L,
      n_genes = 0L,
      n_kegg_edges = 0L,
      n_trrust_edges = 0L,
      n_merged_edges = 0L,
      n_tfs_total = length(tf_set),
      n_tfs_in_expr = 0L,
      n_tfs_in_kg_gene1 = 0L
    )
    next
  }

  celltype_expr <- expr_mat[, celltype_cells, drop = FALSE]
  celltype_meta <- meta_df[celltype_cells, , drop = FALSE]

  # Keep genes with nonzero total expression in this cell type
  keep_genes <- rowSums(celltype_expr, na.rm = TRUE) > 0
  celltype_expr <- celltype_expr[keep_genes, , drop = FALSE]

  gene_set <- rownames(celltype_expr)

  celltype_kegg <- kegg %>%
    filter(gene1 %in% gene_set, gene2 %in% gene_set)

  celltype_trrust <- trrust %>%
    filter(gene1 %in% gene_set, gene2 %in% gene_set)

  celltype_merged <- merged %>%
    filter(gene1 %in% gene_set, gene2 %in% gene_set)

  # Write expression
  celltype_expr_out <- celltype_expr %>%
    rownames_to_column("gene")
  write_csv(celltype_expr_out, file.path(celltype_dir, "expression.csv.gz"))

  # Write metadata
  celltype_meta_out <- celltype_meta %>%
    rownames_to_column("cell")
  write_csv(celltype_meta_out, file.path(celltype_dir, "metadata.csv.gz"))

  # Write merged KG (headerless TSV, as expected by KEGNI)
  write_tsv(
    celltype_merged,
    file.path(celltype_dir, "merged_kg.tsv"),
    col_names = FALSE
  )

  # Write TF list into each cell-type folder for downstream TF->target extraction
  write_csv(tf_df, file.path(celltype_dir, "TF_list.csv"))

  summary_list[[celltype]] <- tibble(
    celltype = celltype,
    n_cells = length(celltype_cells),
    n_genes = nrow(celltype_expr),
    n_kegg_edges = nrow(celltype_kegg),
    n_trrust_edges = nrow(celltype_trrust),
    n_merged_edges = nrow(celltype_merged),
    n_tfs_total = length(tf_set),
    n_tfs_in_expr = sum(tf_set %in% gene_set),
    n_tfs_in_kg_gene1 = sum(tf_set %in% celltype_merged$gene1)
  )
}

summary_df <- bind_rows(summary_list)
write_csv(summary_df, file.path(out_root, "summary.csv"))

message("\nSummary:")
message(capture.output(print(summary_df)), sep = "\n")

message("\nDone.")
message("KEGNI-ready input sets written to: ", out_root)
