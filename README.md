# dropSplit
dropSplit is designed to identify true cells from droplet-based scRNAseq data. 
It consists of three parts: 
* quality control
* model construction and droplet classification
* summarizing features

dropSplit provides some special droplet QC metrics such as CellEntropy or CellGini which can help identification. 

In general, user can use the predefined parameters in the XGBoost and get the important features that help in cell identification. It also provides a automatic XGBoost hyperparameters-tuning function to optimize the model.

# Installation
```
# install.packages("remotes")
remotes::install_github("zh542370159/dropSplit")
```

# Example
```
library(dropSplit)
# Simulate a counts matrix including 20000 empty droplets, 2000 large cells and 200 small cells.
counts <- simCounts(nempty = 20000, nlarge = 2000, nsmall = 200)
counts_label <- gsub(pattern = "-.*", replacement = "", x = colnames(counts), perl = TRUE)
result <- dropSplit(counts)
head(result$meta_info)

dropSplitClass <- result$meta_info$dropSplitClass
# True positive
sum(counts_label %in% c("LargeCell", "SmallCell") & dropSplitClass == "Cell")
# False negative
sum(counts_label %in% c("LargeCell", "SmallCell") & dropSplitClass != "Cell")
# True negative
sum(counts_label == "Empty" & dropSplitClass != "Cell")
# False positive
sum(counts_label == "Empty" & dropSplitClass == "Cell")

# Various QC metrics plot
qc <- QCPlot(result$meta_info)
qc$RankMSE$Merge
qc$CellEntropy$Merge

# Feature importance plot
pl <- ImportancePlot(result$meta_info, result$train, result$importance_matrix, top_n = 20)
pl$Importance
pl$preDefinedClassExp

```
