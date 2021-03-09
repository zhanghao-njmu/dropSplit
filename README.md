# dropSplit

dropSplit is designed to accurately identify 'Cell' droplets for the droplet-based scRNAseq data. 

It consists of four main steps:
* Pre-define droplets as 'Cell', 'Uncertain', 'Empty' and 'Discarded' droplets according to the RankMSE curve.
* Simulate 'Cell' and 'Uncertain' droplets under a depth of 'Empty' used for model construction and prediction.
* Iteratively buid the model, classify 'Uncertain' droplets, and update the training set use newly predicted 'Empty'.
* Make classification with the optimal model.

dropSplit provides some special droplet QC metrics such as CellEntropy or CellGini and plot function, which can help quickly check quality of the prediction.
    
In general, user can use the predefined parameters in the XGBoost. It also provides a automatic XGBoost hyperparameters-tuning function to optimize the model.

# Installation
```
# install.packages("remotes")
remotes::install_github("zh542370159/dropSplit")
```

# Example
```
# Library -----------------------------------------------------------------
library(dropSplit)

# Simple Simulation ---------------------------------------------------------------
simple_counts <- simSimpleCounts()
true <- strsplit(colnames(simple_counts), "-")
true <- as.data.frame(Reduce(function(x, y) rbind(x, y), true))
colnames(true) <- c("label", "Type", "Cluster", "Cell")
rownames(true) <- colnames(simple_counts)
true_label <- true$label
table(true_label)

## DropSplit ---------------------------------------------------------------
result1 <- CellCalling(simple_counts, method = "dropSplit")

# check the result with true classification
table(true_label,result1$meta_info$dropSplitClass)
table(true_label,result1$meta_info$preDefinedClass)

# check the classification by various QC metrics
qc_class <- QCplot(result1$meta_info,colorBy = "dropSplitClass")
qc_class$CellEntropy$Merge
qc_score <- QCplot(result1$meta_info,colorBy = "dropSplitScore")
qc_score$CellEntropy$Merge

# compare with the true labels
result1$meta_info$true_label <- true_label
qc_true <- QCplot(result1$meta_info,colorBy = "true_label")
qc_true$CellEntropy$Merge

## CellRangerV2 -----------------------------------------------------------
result2 <- CellCalling(simple_counts, method = "CellRangerV2")
table(true_label,result2$meta_info$CellRangerV2Class)

## CellRangerV3 -----------------------------------------------------------
result3 <- CellCalling(simple_counts, method = "CellRangerV3")
table(true_label,result3$meta_info$CellRangerV3Class)

## EmptyDrops -----------------------------------------------------------
result4 <- CellCalling(simple_counts, method = "EmptyDrops")
table(true_label,result4$meta_info$EmptyDropsClass)

## zUMIs -----------------------------------------------------------
result5 <- CellCalling(simple_counts, method = "zUMIs")
table(true_label,result5$meta_info$zUMIsClass)


```
