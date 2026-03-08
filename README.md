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
GRNBoost2 results are inferred from the **CHOOSE multiome** data, specifically the file `postprocessed-grnboost2-all-tfs-celltypes-global-var-genes.csv.gz`.

In `postprocessed-grnboost2-all-tfs-celltypes-global-var-genes.csv.gz`, each row lists a putative TF target. 

Among other columns, the cell type is in `cell.type`, the transcription factor gene symbol is in `TF`, the target gene symbol is in `target`, the predicted mode of regulation (MoR) is in `score`, and whether that score passes a (GRNBoost2-specific) significance threshold is in `significant`.

MoR (a termed used by our evaluation metric) is between -1 and 1, with the absolute value indicating the strength of the TF target interaction (higher absolute values imply stronger interactions) and the sign indicating the direction of regulation (with positive values indicating TF induction/positive regulation of the target and negative values indicating TF repression/negative regulation of the target).

Our default GRNBoost2 outputs are...


