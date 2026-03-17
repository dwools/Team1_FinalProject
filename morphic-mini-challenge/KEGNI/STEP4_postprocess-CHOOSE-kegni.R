#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages(library(pacman))
suppressPackageStartupMessages(p_load(data.table))
suppressPackageStartupMessages(p_load(plyr))
suppressPackageStartupMessages(p_load(dplyr))

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

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) == 0) {
    return(default)
  }

  if (idx[length(idx)] == length(args)) {
    stop("Missing value for argument: ", flag, call. = FALSE)
  }

  args[idx[length(idx)] + 1]
}

default_input <- file.path(kegni_root, "tf_target_edges", "all_celltypes_tf_target_ranked.csv")
default_output <- file.path(project_root, "resources", "postprocessed-kegni-all-tfs-celltypes.csv.gz")

input.file <- resolve_path(get_arg("--input", default_input))
output.file <- resolve_path(get_arg("--out", default_output))

if (!file.exists(input.file)) {
  stop("Input file not found: ", input.file)
}

all.res <- as.data.frame(fread(input.file))

req.cols <- c("cell_type", "rank", "Gene1", "Gene2", "EdgeWeight")
miss.cols <- setdiff(req.cols, colnames(all.res))
if (length(miss.cols) > 0) {
  stop("Missing expected column(s) in ", basename(input.file), ": ",
       paste(miss.cols, collapse = ", "))
}

# Standardize to the GRNBoost2 postprocessing naming convention.
setnames(all.res,
         old = c("cell_type", "Gene1", "Gene2", "EdgeWeight"),
         new = c("cell.type", "TF", "target", "importance"))

# Keep score definition/normalization parallel to GRNBoost2 postprocessing.
# No correlation term is available from KEGNI ranked-edge inputs.
all.res[, "score"] <- all.res[, "importance"]
all.res[, "score"] <- all.res[, "score"] / max(abs(all.res[, "score"]))

# Add per-(cell type, TF) p-values
pval_cutoff <- pnorm(4, lower.tail = FALSE)

all.res <-
  ddply(all.res, .variables = c("cell.type", "TF"),
        .fun = function(df) {
          cols <- colnames(df)
          cols <- cols[!(cols %in% c("cell.type", "TF"))]
          s <- sd(df$importance)
          if (is.na(s) || s == 0) {
            s <- .Machine$double.eps
          }
          df[, "pval"] <- pnorm(df$importance,
                                mean = 0,
                                sd   = s,
                                lower.tail = FALSE)
          df[, "padj"]        <- p.adjust(df[, "pval"])
          df[, "significant"] <- FALSE
          df[df[, "pval"] < pval_cutoff, "significant"] <- TRUE
          df[, c(cols, "pval", "padj", "significant")]
        })

odir <- dirname(output.file)
dir.create(odir, recursive = TRUE, showWarnings = FALSE)
print(output.file)
fwrite(as.data.frame(all.res),
  file      = output.file,
       row.names = FALSE,
       col.names = TRUE,
       quote     = FALSE,
       sep       = ",")