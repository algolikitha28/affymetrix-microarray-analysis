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

#network construction

#install packages
install.packages("igraph")
install.packages("ggraph")
install.packages("tidygraph")

#load 
library(igraph)
library(ggraph)
library(tidygraph)

#select top genes
topgenes <- rownames(results)[1:20]
network_data <- expr[topgenes, ]

#calculate corelation matrix
cor_matrix <- cor(
  t(network_data),
  method = "pearson"
)

#remove weak corelations
cor_matrix[abs(cor_matrix) < 0.8] <- 0

#create network
library(igraph)
network <- graph_from_adjacency_matrix(
  cor_matrix,
  mode = "undirected",
  weighted = TRUE,
  diag = FALSE
)

#plot network
plot(
  network,
  vertex.size = 20,
  vertex.label.cex = 0.8,
  main = "Gene Co-expression Network"
)

#instal annotation package
BiocManager::install("hgu95av2.db")

#load
library(hgu95av2.db)
library(AnnotationDbi)

#Get Probe IDs from Network
probe_ids <- V(network)$name

#map to gene symbols
gene_symbols <- mapIds(
  hgu95av2.db,
  keys = probe_ids,
  column = "SYMBOL",
  keytype = "PROBEID",
  multiVals = "first"
)

#add labels to network
V(network)$label <- gene_symbols

#plot network
plot(
  network,
  vertex.size = 20,
  vertex.label.cex = 0.7,
  main = "Gene Co-expression Network"
)