# Gene Regulatory Network Inference & VIPER Evaluation

This repository evaluates the performance of different gene regulatory network (GRN) inference algorithms and control baselines using single-cell RNA sequencing data. 

Each algorithm and baseline control has its own standalone R Jupyter Notebook located in its respective directory. These notebooks handle the end-to-end process: loading the data, running the specific inference/control method, and applying VIPER for protein activity evaluation.

## :open_file_folder: Project Structure

Ensure your repository is structured as follows. All standalone notebooks will read from the shared `resources` directory.

```text
your_project_folder/
├── GENIE3/
│   └── apply-viper-evaluation_genie3.ipynb                 # GENIE3 inference & VIPER
├── TIGRESS/
│   └── apply-viper-evaluation_tigress.ipynb                # TIGRESS inference & VIPER
├── KEGNI/
│   └── apply-viper-evaluation_kegni.ipynb                  # KEGNI inference & VIPER
├── MOSTEXPRESSED/
│   └── apply-viper-evaluation_control_mostexpressed.ipynb  # Control: Most expressed TFs
├── MOSTVARIABLE/
│   └── apply-viper-evaluation_control_mostvariable.ipynb   # Control: Most variable TFs
└── resources/
    ├── postprocessed-mean-baseline.csv.gz          # Baseline mean metrics
    ├── postprocessed-variance-baseline.csv.gz      # Baseline variance metrics
    ├── CHOOSE-sc-wt-and-tf-log-norm.csv.gz         # Compressed log-normalized expression data
    ├── CHOOSE-sc-wt-and-tf-metadata.csv            # Cell metadata 
    └── CHOOSE-tf-to-analyze-metadata.csv           # Metadata for target transcription factors
```

# IMPORTANT DATA INFORMATION

The data and scripts in this repository were obtained from the following sources:

- Brian White's "morphic-mini-challenge" GitHub Repo: https://github.com/bswhite/morphic-mini-challenge/tree/main

- Ka Yee Yeung's "Data" Google Drive: https://drive.google.com/drive/folders/12tt1rCTQtdDrb04PAOcDoRwOHgb1kd28
  - Data: (dir)
    - *KaYee-notebook-incomplete.ipynb*
    - *README - Description of training data.docx*
      
    - CHOOSE-multiome-expr: (dir)
      - *CHOOSE-multiome-wt-log-norm.csv.gz*
      - *CHOOSE-multiome-wt-metadata.csv*
        
    - CHOOSE-scRNA-seq (dir)
      - *CHOOSE-sc-wt-and-tf-log-norm.csv.gz* (**NOTE**: This file is **NOT** included in this remote repository, as it exceeds 100 MB. It will need to be obtained from the Google Drive above and copied to your local directory of this repo.)
      - *CHOOSE-sc-wt-and-tf-metadata.csv* (**NOTE**: After downloading this CSV file from the Google Drive, its entire header row needed to be shifted 1 column to the right, such that Column A (which should be the data's row names) has no header value.)  
     
(**NOTE**: All of these data files have been relocated to the "morphic-mini-challenge/resources" directory, as this was specified by the reference VIPER script (apply-viper-evaluation.R) that Brian White provided, on which our VIPER evaluation notebooks are based)

## GRNBoost2
### Reminder of what GRNBoost2 does
GRNBoost2 applies gradient boosting machine (GBM) regression to learn a relationship between a (potential target) gene and all input TFs. More specifically, the expression of input TFs is used to predict the expression of each target gene. 

**Remember**: a transcription factor prompts/stimulates the transcription of a specific set of genes (this occurs by the transcription factor *protein* binding with DNA in such a manner that those genes are transcribed into mRNA transcripts).

So, if a transcription factor has been knocked out (i.e., the protein isn't floating around in the cell anymore), the expression of those genes (i.e., the amount of mRNA molecules assembled in concordance with that gene's DNA sequence, and the amount of protein molecules assembled from those mRNA molecules) will drop.

Thus, we assume that those genes will be expressed more in wild-type cells than in cells in their transcription factor is knocked out.

We label putative TF-target interactions as significant if they are above four standard deviations from the mean.

### Input
Among our default data, the *inputs* for GRNBoost2 are our *training data*:
- `CHOOSE-multiome-wt-metadata.csv` (metadata)
- `CHOOSE-multiome-wt-log-norm.csv.gz` (expression data matrix)

In the *metadata*, each row corresponds to a cell, with its cell id the row name.

These cell ids are the column names of the *data matrix*. The cell type is provided in the `celltype_jf` column.

The data are log-normalized (i.e., we applied `NormalizeData(normalization.method = “LogNormalize”, scale.factor=10000)`, which divides the gene (or “feature”) count of each cell by the total counts for that cell, multiplies it by 10,000, and finally adds one and natural log transforms (`log1p`).

We suggest using `fread` from the `data.table` R package for reading large expression matrices.

### Output:
Among our default data, the *outputs* for GRNBoost2 (that are input into VIPER) are:
- `postprocessed-grnboost2-all-tfs-celltypes-global-var-genes.csv.gz`

In `postprocessed-grnboost2-all-tfs-celltypes-global-var-genes.csv.gz`, each row lists a putative TF-target interaction. 

Among other columns, the cell type is in `cell.type`, the transcription factor gene symbol is in `TF`, the target gene symbol is in `target`, the predicted mode of regulation (MoR) is in `score`, and whether that score passes a (GRNBoost2-specific) significance threshold is in `significant`.

MoR (a termed used by our evaluation metric) is between -1 and 1, with the absolute value indicating the strength of the TF target interaction (higher absolute values imply stronger interactions) and the sign indicating the direction of regulation (with positive values indicating TF induction/positive regulation of the target and negative values indicating TF repression/negative regulation of the target).

## VIPER evaluation
### Reminder of what VIPER does
VIPER quantifies the activity of a transcription factor protein via the mRNA expression of the genes that specific transcraption factor regulates.

(Reiterated from above):

_Remember: a transcription factor prompts/stimulates the transcription of a specific set of genes (this occurs by the transcription factor *protein* binding with DNA in such a manner that those genes are transcribed into mRNA transcripts)._

_So, if a transcription factor has been knocked out (i.e., the protein isn't floating around in the cell anymore), the expression of those genes (i.e., the amount of mRNA molecules assembled in concordance with that gene's DNA sequence, and the amount of protein molecules assembled from those mRNA molecules) will drop._

_Thus, we assume that those genes will be expressed more in wild-type cells than in cells in their transcription factor is knocked out._

(End of reiteration)

These expression discrepancies—between wild-type cells and cells with a transcription factor knocked out—are what VIPER computes and plots.

But wait, how is this different from GRNBoost2??

VIPER is where the GRNBoost2 classifier model is evaluated. VIPER's job isn't to identify correlations between transcription factors and gene expressions. VIPER's job is to quantify the accuracy of the results generated by GRNBoost2's machine learning model.

### Input
Among our default data, the inputs for VIPER evaluation are:
- `CHOOSE-sc-wt-and-tf-log-norm.csv.gz`
  - (These are single-cell RNA sequencing data; they represent how much each gene is expressed in each cell type. These are the _evaluation_ data on which to apply our GRNBoost2 classifier that was trained in the GRNBoost2 step.)
- `CHOOSE-sc-wt-and-tf-metadata.csv`
  - (These are single-cell RNA sequencing metadata, providing information on each of the cell types included in the scRNA-seq evaluation data above.)
- `CHOOSE-tf-to-analyze-metadata.csv`
- `postprocessed-grnboost2-all-tfs-celltypes-global-var-genes.csv.gz` (GRNBoost2 output)
  - (These two are combined to represent the relevant GRNBoost2 information. They represent the relevant the pairs (transcription factors, and the genes they stimulate) identified by GRNBoost2. `postprocepostprocessed-grnboost2-all-tfs-celltypes-global-var-genes.csv.gz` is filtered for rows (which represent putative TF-target interactions) where `significant` = `TRUE`, and for the transcription factors listed in `CHOOSE-tf-to-analyze-metadata.csv`.

### Output
VIPER generates a plot showing WT and TF KO distributions for six TF/cell type pairs, which includes the Wilcoxon p-value. See the plot generated from the default data at `mini-morphic-challenge` >> `plots` >> `CHOOSE-TF-viper-distributions.png`).

**NOTE**: only *one* TF/cell-type pair (TF BCL11A in cell type ctx_npc) is significant. This is the most reliable pair, and we recommend you focus on it during your testing (in the image file indicated above, it's the top-middle plot.)
