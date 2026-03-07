suppressPackageStartupMessages(library(pacman))
suppressPackageStartupMessages(p_load(foreach))
suppressPackageStartupMessages(p_load(parallel))
suppressPackageStartupMessages(p_load(data.table))
suppressPackageStartupMessages(p_load(plyr))
suppressPackageStartupMessages(p_load(dplyr))
suppressPackageStartupMessages(p_load(docstring))
suppressPackageStartupMessages(p_load(decoupleR))
suppressPackageStartupMessages(p_load(this.path))
suppressPackageStartupMessages(p_load(ggpubr))

# Initialization parallelism
num.cores <- detectCores()
if(!is.na(num.cores) && (num.cores > 1)) {
  suppressPackageStartupMessages(p_load("doMC"))
  cat(paste("Registering ", num.cores-1, " cores.\n", sep=""))
  registerDoMC(cores=(num.cores-1))
}

# Initialize variables
# Copy scRNA-seq metadata and data files from Drive to resources folder
root <- dirname(dirname(this.path::this.path()))
resources.dir <- paste0(root, "/resources/")
sc.metadata.file <- paste0(resources.dir, "CHOOSE-sc-wt-and-tf-metadata.csv")
sc.data.file <- paste0(resources.dir, "CHOOSE-sc-wt-and-tf-log-norm.csv.gz")
grnboost2.res.file <- paste0(resources.dir, "postprocessed-grnboost2-all-tfs-celltypes-global-var-genes.csv.gz")
tfs.to.analyze.file <- paste0(resources.dir, "CHOOSE-tf-to-analyze-metadata.csv")

plot.dir <- paste0(root, "/plots/")
dir.create(plot.dir, showWarnings = FALSE)

# Load in predictions to be evaluated. Here from GRNBoost2
grnboost2.res <- as.data.frame(fread(grnboost2.res.file))

# Read in CHOOSE scRNA-seq data and metadata
fread_with_rownames <- function(file) {
  tbl <- as.data.frame(fread(file))
  rownames(tbl) <- tbl$V1
  tbl <- tbl[, !(colnames(tbl) == "V1")]
  tbl
}

sc.data <- fread_with_rownames(sc.data.file)
sc.metadata <- fread_with_rownames(sc.metadata.file)

# Read in the TF / cell type combinations to analyze
tfs.to.analyze <- read.csv(tfs.to.analyze.file, header=TRUE)

score.viper <- function(predictions.tbl, sc.metadata, sc.data, 
                        sc.cell.type.col = "celltype_jf", sc.genotype.col = "gRNA_perturb", sc.wt.status = "CTRL") 
{
  #' Apply VIPER to compute single-cell transcription factor (TF) activities from predicted TF targets.
  #' 
  #' @description Compute VIPER TF activity scores using predicted TF targets and wild type (WT) and TF knockout (KO) scRNA-seq
  #' data.
  #' 
  #' @param predictions.tbl data.frame. A data.frame of cell type-specific TF targets, where each row corresponds to a 
  #' cell type-specific TF target, which is described by columns TF, cell.type, target, and score. TF and target are
  #' gene symbols, appearing as row names in the scRNA-seq data in sc.data. score is the strength of the interaction
  #' and should be in [-1, +1], with positive values corresponding to a positively regulated (induced) target and
  #' negative values corresponding to a negatively regulated (repressed) target.
  #' @param sc.metadata data.frame. Metadata describing the single-cell data, where each row is a cell. Row names
  #' correspond to the cell id. The cell's type is provided in column sc.cell.type.col and the KO'ed TF is in
  #' column sc.genotype.col. WT cells have gRNA_perturb set to CTRL.
  #' @param sc.data data.frame. scRNA-seq data, with row names gene symbols and column names cell ids. The data
  #' are log normalized.
  #' @param sc.cell.type.col. character string providing the name of the column in sc.metadata holding the cell type.
  #' @param sc.genotype.col. character string providing the name of the column in sc.metadata holding WT vs TF KO status.
  #' @param sc.wt.status. character string providing the status of WT cells in the sc.genotype.col
  #' @return A data.frame with VIPER TF activity scores for each cell. Each row is a cell, with VIPER TF activity score in column score.
  #' The TF is in column TF. The cell's type is in column cell.type and its WT or TF KO status
  #' is in column condition (as WT or KO, respectively). 
  
  ddply(predictions.tbl,
        .variables = c("TF", "cell.type"),
        .parallel = FALSE,
        .fun = function(df) {
          tf <- df[1, "TF"]
          ct <- df[1, "cell.type"]
          
          ko.flag <- (sc.metadata[, sc.genotype.col] == tf) & (sc.metadata[, sc.cell.type.col] == ct)
          wt.flag <- (sc.metadata[, sc.genotype.col] == sc.wt.status) & (sc.metadata[, sc.cell.type.col] == ct)
          ko.mat <- sc.data[, ko.flag]
          wt.mat <- sc.data[, wt.flag]
          
          network <- data.frame("source" = as.character(tf), "target" = as.character(df$target), "mor" = df$score)
          # VIPER seems to crash if we don't have at least two TFs, so create a dummy.
          dummy.network <- network
          dummy.network$source <- "dummy"
          network <- rbind(network, dummy.network)
          
          combined.mat <- cbind(ko.mat, wt.mat)
          vp <- run_viper(combined.mat, network, .source = "source", .target = "target", .mor = "mor", pleiotropy = FALSE, method = "none")
          
          vp <- subset(vp, source != "dummy")
          vp <- merge(vp, sc.metadata[ko.flag | wt.flag, c(sc.genotype.col), drop=FALSE], by.x = c("condition"), by.y = c("row.names"))
          vp <- vp[, !(colnames(vp) %in% c("condition", "statistic", "source", "p_value"))]
          colnames(vp)[colnames(vp) == sc.genotype.col] <- "perturbation"
          vp$condition <- vp$perturbation
          vp$perturbation[vp$perturbation == sc.wt.status] <- "none"
          vp$condition[vp$condition != sc.wt.status] <- "KO"
          vp$condition[vp$condition == sc.wt.status] <- "WT"
          vp <- cbind(vp, TF = tf, cell.type = ct)
          vp
        })
}

compare.viper.distributions <- function(viper.df) {
  #' Compare VIPER TF activity scores between WT and TK KO cells.
  #' 
  #' @description Apply a one-sided Wilcoxon test of the hypothesis that TF activity scores are higher in WT than in TF KO cells.
  #' 
  #' @param viper.df data.frame. A data.frame with VIPER TF activity scores for each cell, e.g., as returned by score.viper.
  #' Each row is a cell, with VIPER TF activity score in column score. The TF is in column TF.
  #' The cell's type is in column cell.type and its WT or TF KO status is in column condition (as WT or KO, respectively).
  #' @return A data.frame with the Wilcoxon significance (column pval) that the VIPER activity score is higher for WT
  #' cells than for TF KO cells for the corresponding cell.type and TF.
  stats <- 
    ddply(viper.df, .variables = c("cell.type", "TF"), 
          .fun = function(df) {
            ret <- wilcox.test(x = subset(df, condition == "WT")$score, 
                               y = subset(df, condition == "KO")$score, alternative="greater")
            data.frame(pval = ret$p.value)
          })
  stats <- stats[order(stats$pval, decreasing=FALSE),]
  stats
}

tmp <- tfs.to.analyze[, c("TF", "celltype_jf")]
colnames(tmp) <- c("TF", "cell.type")
tmp <- merge(grnboost2.res, tmp)
gb2.preds <- subset(tmp, significant == TRUE)

# Compute VIPER TF activity scores for predictions using WT and TF KO scRNA-seq data
gb2.viper <- score.viper(gb2.preds, sc.metadata, sc.data)

# Compare the distribution (over single cells) of TF activity scores, expecting that 
# TF activity will be higher in WT cells than in KO cells where the TF is disrupted.
gb2.viper.pvals <- compare.viper.distributions(gb2.viper)
print(gb2.viper.pvals)

# Create a plot of the distributions of TF activity scores, stratified by WT or KO condition.
# g <- ggplot(data = gb2.viper, aes(x = condition, y = score)) + geom_violin()
# Use ggpubr to include stats on plot
g <- ggviolin(data = gb2.viper, "condition", "score") + stat_compare_means(method.args = list(alternative = "greater"), ref.group = "KO")
g <- g + facet_wrap(cell.type ~ TF, scales = "free_y")

ofile <- paste0(plot.dir, "CHOOSE-TF-viper-distributions.png")
png(ofile)
print(g)
d <- dev.off()
print(ofile)