#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

is_absolute_path <- function(path) {
  grepl("^(?:[A-Za-z]:[/\\\\]|/|\\\\\\\\)", path)
}

# This project assumes the directory structure Team1_FinalProject/morphic-mini-challenge/KEGNI. repo_root should be "Team1_FinalProject", project_root should be "Team1_FinalProject/morphic-mini-challenge", and kegni_root should be "Team1_FinalProject/morphic-mini-challenge/KEGNI". This script is assumed to be in the directory "Team1_FinalProject/morphic-mini-challenge/KEGNI". The script_path function should return the full path to this script, regardless of the current working directory. The resolve_path function should take a relative path and resolve it relative to the kegni_root directory, or return an absolute path if given one.

# Find the path to this script, even if the working directory is not the same as the script's directory. This allows the script to be run from any location and still correctly find its inputs and outputs relative to its own location.
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


# Find the path to the directory containing this script (Team1_FinalProject/morphic-mini-challenge/KEGNI).
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
  index <- which(args == flag)
  if (length(index) == 0) {
    return(default)
  }

  last_index <- index[[length(index)]]
  if (last_index >= length(args)) {
    stop("Missing value for argument: ", flag, call. = FALSE)
  }

  args[[last_index + 1]]
}

parse_flag <- function(value, default = FALSE) {
  if (is.null(value)) {
    return(default)
  }

  normalized <- tolower(trimws(as.character(value)))
  if (normalized %in% c("1", "true", "t", "yes", "y")) {
    return(TRUE)
  }
  if (normalized %in% c("0", "false", "f", "no", "n")) {
    return(FALSE)
  }

  stop("Invalid logical value: ", value, call. = FALSE)
}

find_rscript <- function() {
  candidates <- c(
    Sys.which("Rscript"),
    file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  )

  candidates <- unique(candidates[nzchar(candidates)])
  existing <- candidates[file.exists(candidates)]
  if (length(existing) == 0) {
    stop("Rscript is not available.", call. = FALSE)
  }

  existing[[1]]
}

find_python <- function(user_python = NULL) {
  if (!is.null(user_python) && nzchar(trimws(user_python))) {
    return(list(command = trimws(user_python), prefix_args = character()))
  }

  python_bin <- Sys.which("python")
  if (nzchar(python_bin)) {
    return(list(command = python_bin, prefix_args = character()))
  }

  python3_bin <- Sys.which("python3")
  if (nzchar(python3_bin)) {
    return(list(command = python3_bin, prefix_args = character()))
  }

  py_launcher <- Sys.which("py")
  if (nzchar(py_launcher)) {
    return(list(command = py_launcher, prefix_args = c("-3")))
  }

  stop("Python was not found. Provide --python or add Python to PATH.", call. = FALSE)
}

run_command <- function(command, command_args, error_message) {
  output <- system2(
    command,
    args = shQuote(command_args),
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(output, "status")
  if (is.null(status)) {
    status <- 0L
  }
  if (!identical(as.integer(status), 0L)) {
    tail_output <- if (length(output) > 0) {
      paste(tail(output, 20), collapse = "\n")
    } else {
      "<no output captured>"
    }
    stop(
      error_message,
      " Exit status: ", status,
      "\nCommand: ", command,
      "\nArgs: ", paste(command_args, collapse = " "),
      "\nOutput tail:\n", tail_output,
      call. = FALSE
    )
  }
}

expr_file <- resolve_path(get_arg("--expr", "../resources/CHOOSE-multiome-wt-log-norm.csv.gz"))
meta_file <- resolve_path(get_arg("--meta", "../resources/CHOOSE-multiome-wt-metadata.csv"))
kegg_file <- resolve_path(get_arg("--kegg", "./inputs/KG/KEGG_hs.tsv"))
trrust_file <- resolve_path(get_arg("--trrust", "./inputs/KG/TRRUST_hs.tsv"))
tf_file <- resolve_path(get_arg("--tf", "../resources/TF_list.csv"))
out_root <- resolve_path(get_arg("--out", "kegni_inputs"))

n_neighbors <- as.integer(get_arg("--n_neighbors", "30"))
max_steps <- as.integer(get_arg("--max_steps", "1"))
num_hidden <- as.integer(get_arg("--num_hidden", "256"))
num_heads <- as.integer(get_arg("--num_heads", "4"))
num_layers <- as.integer(get_arg("--num_layers", "2"))
lambda_kge <- as.numeric(get_arg("--lambda_kge", "1"))
device <- as.integer(get_arg("--device", "-1"))
run_eval <- parse_flag(get_arg("--run_eval", "0"), default = FALSE)
ground_truth_dir <- resolve_path(get_arg("--ground_truth_dir", "./data/ground_truth"))
python_config <- find_python(get_arg("--python", ""))

old_wd <- getwd()
setwd(kegni_root)
on.exit(setwd(old_wd), add = TRUE)

train_script <- file.path(kegni_root, "train.py")
step1_script <- file.path(kegni_root, "STEP1_make_kegni_inputs.R")
step3_script <- file.path(kegni_root, "STEP3_extract_kegni_tf_targets.R")
step4_script <- file.path(kegni_root, "STEP4_postprocess-CHOOSE-kegni.R")
step5_script <- file.path(kegni_root, "STEP5_apply-viper-evaluation_kegni.R")
tf_target_edges_dir <- file.path(kegni_root, "tf_target_edges")
postprocessed_output <- file.path(project_root, "resources", "postprocessed-kegni-all-tfs-celltypes.csv.gz")

required_files <- c(
  train_script,
  step1_script,
  step3_script,
  step4_script,
  step5_script,
  expr_file,
  meta_file,
  kegg_file,
  trrust_file,
  tf_file
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop(
    "The following required file(s) do not exist:\n",
    paste("  -", missing_files, collapse = "\n"),
    call. = FALSE
  )
}

meta_df <- read.csv(meta_file, row.names=1, check.names = FALSE, stringsAsFactors = FALSE)
if (!("celltype_jf" %in% colnames(meta_df))) {
  stop("Metadata file must contain a 'celltype_jf' column.", call. = FALSE)
}

cell_types_arg <- get_arg("--cell_types", NULL)
cell_types <- if (is.null(cell_types_arg)) {
  unique(as.character(meta_df$celltype_jf))
} else {
  trimws(unlist(strsplit(cell_types_arg, ",", fixed = TRUE)))
}
cell_types <- cell_types[nzchar(cell_types) & !is.na(cell_types)]

if (length(cell_types) == 0) {
  stop("No cell types were found to process.", call. = FALSE)
}

cat("==================================================\n")
cat("Step 1: Building KEGNI input folders\n")
cat("==================================================\n")

run_command(
  command = find_rscript(),
  command_args = {
    step1_args <- c(
    step1_script,
    "--expr", expr_file,
    "--meta", meta_file,
    "--kegg", kegg_file,
    "--trrust", trrust_file,
    "--tf", tf_file,
    "--out", out_root
    )

    if (!is.null(cell_types_arg) && nzchar(trimws(cell_types_arg))) {
      step1_args <- c(step1_args, "--cell_types", cell_types_arg)
    }

    step1_args
  },
  error_message = "Failed while building KEGNI input folders."
)

cat("\n")
cat("==================================================\n")
cat("Step 2: Running KEGNI for each cell type\n")
cat("==================================================\n")

for (cell_type in cell_types) {
  input_expr <- file.path(out_root, cell_type, "expression.csv.gz")
  input_kg <- file.path(out_root, cell_type, "merged_kg.tsv")

  if (!file.exists(input_expr)) {
    cat(sprintf("Skipping %s: missing %s\n", cell_type, input_expr))
    next
  }

  if (!file.exists(input_kg)) {
    cat(sprintf("Skipping %s: missing %s\n", cell_type, input_kg))
    next
  }

  cat("\n")
  cat("----------------------------------------------\n")
  cat(sprintf("Running KEGNI for cell type: %s\n", cell_type))
  cat("----------------------------------------------\n")

  command_args <- c(
    python_config$prefix_args,
    train_script,
    "--input", input_expr,
    "--data_path", input_kg,
    "--n_neighbors", as.character(n_neighbors),
    "--max_steps", as.character(max_steps),
    "--num_hidden", as.character(num_hidden),
    "--num_heads", as.character(num_heads),
    "--num_layers", as.character(num_layers),
    "--lambda_kge", as.character(lambda_kge),
    "--device", as.character(device)
  )

  if (run_eval) {
    if (dir.exists(ground_truth_dir)) {
      cat(sprintf("RUN_EVAL=1 and ground truth directory exists; adding --eval for %s\n", cell_type))
      command_args <- c(command_args, "--eval")
    } else {
      cat(sprintf("RUN_EVAL=1 but GROUND_TRUTH_DIR does not exist; running without --eval for %s\n", cell_type))
    }
  }

  # Snapshot outputs/ before this training run so we can identify new files afterwards
  outputs_dir <- file.path(kegni_root, "outputs")
  embedding_pattern <- "-(?:kgg|relation|scg)_embedding\\.csv$"

  run_command(
    command = python_config$command,
    command_args = command_args,
    error_message = sprintf("KEGNI training failed for cell type: %s.", cell_type)
  )

  # Rename the newly created embedding CSVs to include the cell type prefix
  # change working directory to outputs/ to find the new files, then change back after renaming
  setwd(outputs_dir)
  files_to_rename <- c(list.files(pattern = embedding_pattern)) # replace with actual patternfiles_to_rename

  for (file in files_to_rename) {
      if (!grepl("_expression.csv", file)) {
          new_file_name <- paste0(cell_type, "_", file)
          file.rename(file, new_file_name)
          print(new_file_name)
      }
    }
  setwd(kegni_root)
  }

cat("\n")
cat("==================================================\n")
cat("Step 3: Extracting KEGNI TF-target edges\n")
cat("==================================================\n")

run_command(
  command = find_rscript(),
  command_args = {
    step3_args <- c(
    step3_script,
    "--meta", meta_file,
    "--inputs_root", out_root,
    "--outputs_root", file.path(kegni_root, "outputs"),
    "--out_root", tf_target_edges_dir
    )

    if (!is.null(cell_types_arg) && nzchar(trimws(cell_types_arg))) {
      step3_args <- c(step3_args, "--cell_types", cell_types_arg)
    }

    step3_args
  },
  error_message = "Failed while extracting KEGNI TF-target edges."
)

cat("\n")
cat("==================================================\n")
cat("Step 4: Postprocessing KEGNI outputs\n")
cat("==================================================\n")

run_command(
  command = find_rscript(),
  command_args = c(
    step4_script,
    "--input", file.path(tf_target_edges_dir, "all_celltypes_tf_target_ranked.csv"),
    "--out", postprocessed_output
  ),
  error_message = "Failed while postprocessing KEGNI outputs."
)

cat("\n")
cat("==================================================\n")
cat("Step 5: Running VIPER evaluation\n")
cat("==================================================\n")

run_command(
  command = find_rscript(),
  command_args = c(step5_script),
  error_message = "Failed while running VIPER evaluation."
)

cat("\nAll requested KEGNI runs have completed.\n")