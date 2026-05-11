# Affymetrix Microarray Analysis using R

## Project Purpose

The purpose of this project is to analyze Affymetrix microarray gene expression data and identify genes that are differentially expressed between biological conditions.

This workflow helps in:

- Understanding gene expression patterns
- Identifying potential biomarker genes
- Comparing disease vs control samples
- Visualizing transcriptomic differences

---

## What is Microarray Analysis?

Microarray analysis is a technique used to measure the expression levels of thousands of genes simultaneously.

Affymetrix microarrays generate raw `.CEL` files which contain probe intensity data for each sample.

These raw files need preprocessing before biological interpretation.

---

## Dataset Used

- Dataset: GSE1004
- Platform: Affymetrix HG_U95Av2
- Data Type: Affymetrix Microarray CEL files

---

## Workflow Used

```text
Download GEO Dataset (from NCBI)
        ↓
Extract CEL Files
        ↓
Read CEL Files in R
        ↓
RMA Normalization
        ↓
Generate Expression Matrix
        ↓
Hierarchical Clustering
        ↓
PCA Analysis
        ↓
Differential Expression Analysis
        ↓
Heatmap Visualization
        ↓
Biomarker Discovery using geNetClassifier
```

---


###  RMA Normalization

RMA (Robust Multi-array Average) normalization was performed.

This step includes:
- Background correction
- Quantile normalization
- Expression summarization

Purpose:
- Remove technical variation
- Make samples comparable

Output:
- Normalized expression matrix

---

###  Hierarchical Clustering

Hierarchical clustering was performed using:

```r
sample_dist <- dist(t(expr))
hc <- hclust(sample_dist)
plot(hc)
```

Purpose:
- Group samples based on similarity
- Identify clustering patterns
- Detect possible outliers

Observation:
- Samples formed distinct clusters indicating biological differences between conditions.

---

### PCA Analysis

Principal Component Analysis (PCA) was performed using:

```r
pca <- prcomp(t(expr))

plot(
  pca$x[,1],
  pca$x[,2],
  col = as.numeric(labels),
  pch = 19
)
```

Purpose:
- Reduce high-dimensional gene expression data
- Visualize sample separation
- Check class discrimination

Observation:
- PCA showed clear separation between sample groups.

---

###  Differential Expression Analysis

Differential expression analysis was performed using the `limma` package.

Purpose:
- Identify significantly expressed genes
- Compare disease vs control samples

Output Columns:
- logFC → log fold change
- P.Value → statistical significance
- adj.P.Val → adjusted p-value
- t → moderated t-statistic

Observation:
- Multiple genes showed strong differential expression.

---

###  Heatmap Visualization

Heatmap visualization was performed using `pheatmap`.

Purpose:
- Visualize expression patterns of top genes
- Observe clustering of samples and genes

Observation:
- Top differentially expressed genes showed distinct expression patterns across samples.

---

###  Biomarker Discovery using geNetClassifier

`geNetClassifier` was used for biomarker discovery and classification.

Purpose:
- Rank genes based on classification power
- Identify biomarker candidate genes
- Perform machine learning-based classification

The package internally uses:
- Support Vector Machine (SVM)

Biomarker candidate genes were identified from:
- Top ranked genes
- Differentially expressed genes
- Classification-associated genes

---

##  geNetClassifier

`geNetClassifier` is a Bioconductor package used for:

- Biomarker discovery
- Gene ranking
- Sample classification
- Machine learning-based analysis

It helps identify genes that best distinguish biological conditions.

---


## Outputs Generated

### Plots

- Hierarchical clustering plot
- PCA plot
- Heatmap
- Volcano plot

### Result Files

- Normalized expression matrix
- Differential expression results
- Biomarker candidate genes
- Saved classifier model

---

## Tools and Packages Used

- R
- Bioconductor
- affy
- limma
- pheatmap
- geNetClassifier

---

## Conclusion

This project demonstrates a complete Affymetrix microarray analysis workflow using R.

The workflow includes:
- Data preprocessing
- Normalization
- Clustering
- PCA analysis
- Differential expression analysis
- Heatmap visualization
- Biomarker discovery using machine learning

The analysis successfully identified genes with significant expression differences and potential biomarker candidates between biological conditions.

