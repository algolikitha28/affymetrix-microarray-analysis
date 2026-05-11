#set the working directory
setwd("C:/Users/HP/Downloads/CEL_files")

#check the samples
list.files()

#install packages
install.packages("BiocManager")

BiocManager::install(c(
  "affy",
  "limma",
  "geNetClassifier"
))

install.packages("pheatmap")

#load packages
library(affy)
library(limma)
library(geNetClassifier)
library(pheatmap)

#read CEL files
data <- ReadAffy()
#check samples
sampleNames(data)

#normalise using rma
eset <- rma(data)

#extract expression matrix
expr <- exprs(eset)
#check dimension
dim(expr)

#view expression matrix
head(expr)

#create sample labels
labels <- factor(c(
  rep("Control",11),
  rep("Disease",11)
))
#check
length(labels)

#hirerachical clustering
sample_dist <- dist(t(expr))
hc <- hclust(sample_dist)
plot(hc)

#heat map clustering
pheatmap(expr)

#pca analysis
pca <- prcomp(t(expr))

plot(
  pca$x[,1],
  pca$x[,2],
  col = as.numeric(labels),
  pch = 19
)

#differential expression analysis
design <- model.matrix(~labels)
fit <- lmFit(expr, design)
fit <- eBayes(fit)

#top differentially expressed genes
results <- topTable(
  fit,
  coef = 2,
  number = 20
)
results

#save results
write.csv(results, "Top_DE_Genes.csv")

#volcano plot
volcanoplot(fit, coef=2)

#train gNetclassifier
model <- geNetClassifier(
  eset = expr,
  sampleLabels = labels
)

#check models
class(model)

#explore models
slotNames(model)

#view classification gens
model@classificationGenes

#view gene rankings
model@genesRanking

#save model
save(
  model,
  file = "geNetClassifier_model.RData"
)