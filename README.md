# dropSplit
dropSplit is designed to identify true cells from droplet-based scRNAseq data. 
It consists of three parts: 
* quality control
* model construction and droplet classification
* summarizing features

dropSplit provides some special QC metrics such as CellEntropy or CellGini which can help identification. 

In general, user can use the predefined parameters in the XGBoost and get the important features that help in cell identification. It also provides a automatic XGBoost hyperparameters-tuning function to optimize the model.

# Installation
```
# install.packages("remotes")
remotes::install_github("zh542370159/dropSplit")
```

# Example
```
library(dropSplit)
counts <- DropletUtils:::simCounts()
colnames(counts) <- paste0("Cell-", 1:ncol(counts))
result <- dropSplit(counts)
head(result$meta_info)
```
