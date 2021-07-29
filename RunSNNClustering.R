library(dplyr)
library(Seurat)
library(patchwork)
library(mclust)

source("clustering_evaluation.R")

### Step 1: Choose data set from ProcessedDasets:

## ElongatedGaussiansWithBridge3
#mydata <- read.csv(file = 'ElongatedGaussiansWithBridge3_Denoised_2d.csv',header = FALSE)
#mydata <- read.csv(file = 'ElongatedGaussiansWithBridge3_Denoised.csv',header = FALSE)
#true_labels <- read.csv(file = 'ElongatedGaussiansWithBridge3_Denoised_Labels.csv',header = FALSE)

## Balls
#mydata <- read.csv(file = 'Balls_Denoised.csv',header = FALSE)
#true_labels <- read.csv(file = 'Balls_DenoisedLabels.csv',header = FALSE)


## SwissRoll
#mydata <- read.csv(file = 'SwissRoll_Denoised.csv',header = FALSE)
#true_labels <- read.csv(file = 'SwissRoll_DenoisedLabels.csv',header = FALSE)


## GL_Manifold
#mydata <- read.csv(file = 'GL_Manifold_Denoised.csv',header = FALSE)
#true_labels <- read.csv(file = 'GL_Manifold_DenoisedLabels.csv',header = FALSE)


## ElongatedWithBridge
#mydata <- read.csv(file = 'ElongatedWithBridge_Denoised.csv',header = FALSE)
#true_labels <- read.csv(file = 'ElongatedWithBridge_DenoisedLabels.csv',header = FALSE)


## RNAMix1Basic
#mydata <- read.csv(file = 'RNAMix1Basic.csv',header = FALSE)
#true_labels <- read.csv(file = 'RNAMix1Labels.csv',header = FALSE)

##RNAMix2Basic
#mydata <- read.csv(file = 'RNAMix2Basic.csv',header = FALSE)
#true_labels <- read.csv(file = 'RNAMix2Labels.csv',header = FALSE)

## TMLungBasic
#mydata <- read.csv(file = 'TMLungBasic.csv',header = FALSE)
#true_labels <- read.csv(file = 'TMLungLabels.csv',header = FALSE)

## Beta2Basic 
#mydata <- read.csv(file = 'Beta2Basic.csv',header = FALSE)
#true_labels <- read.csv(file = 'Beta2Labels.csv',header = FALSE)

## TMPancBasic
#mydata <- read.csv(file = 'TMPancBasic.csv',header = FALSE)
#true_labels <- read.csv(file = 'TMPancLabels.csv',header = FALSE)

## BaronPancSCT
#mydata <- read.csv(file = 'BaronPancSCT.csv',header = FALSE)
#true_labels <- read.csv(file = 'BaronPancLabels.csv',header = FALSE)

## PBMC3kSCT
#mydata <- read.csv(file = 'PBMC3kSCT.csv',header = FALSE)
#true_labels <- read.csv(file = 'PBMC3klabels.csv',header = FALSE)

## PBMC4k BASIC
mydata <- read.csv(file = 'PBMC4kMAINBasic.csv',header = FALSE)
true_labels <- read.csv(file = 'PBMC4kmainlabels_merged_Tcells.csv',header = FALSE)


## CellMixSngSCT
#mydata <- read.csv(file = 'CellMixSngSCT.csv',header = FALSE)
#true_labels <- read.csv('CellMixSngLabels.csv',header = FALSE)

### Step 3: Create Seurat Object 
rownames(mydata) = 1:nrow(mydata)
Seurat_mydata <- CreateSeuratObject(counts = mydata, project = "mydata", min.cells = 0, min.features = 0)
Seurat_mydata <- ScaleData(Seurat_mydata, features = rownames(mydata))
Seurat_mydata@assays$RNA@scale.data <- as.matrix(mydata)
rownames(mydata) <- paste0("PC_", 1:min(nrow(mydata),40))
#Seurat_mydata <- RunPCA(Seurat_mydata, features = rownames(Seurat_mydata),return.only.var.genes = FALSE, npcs = 2)
Seurat_mydata[["pca"]] <- CreateDimReducObject(embeddings = as.matrix(t(mydata)), key = "PC_",assay='RNA')
Seurat_mydata <- FindNeighbors(Seurat_mydata, dims = 1:min(nrow(mydata),40))
colnames(true_labels)=NULL
rownames(true_labels) <- colnames(x =Seurat_mydata)
Seurat_mydata<-AddMetaData(Seurat_mydata, true_labels, col.name = 'true_labels')

### Step 3: Clustering

## Elongated Gaussians:
#Seurat_mydata <- FindClusters(Seurat_mydata, resolution = 0.03)
#Seurat_labels <- Seurat_mydata$RNA_snn_res.0.03

## Balls:
#Seurat_mydata <- FindClusters(Seurat_mydata, resolution = 0.05)
#Seurat_labels <- Seurat_mydata$RNA_snn_res.0.05

## SwissRoll:
#Seurat_mydata <- FindClusters(Seurat_mydata, resolution = 0.05)
#Seurat_labels <- Seurat_mydata$RNA_snn_res.0.05

## GL_Manifold
#Seurat_mydata <- FindClusters(Seurat_mydata, resolution = 0.05)
#Seurat_labels <- Seurat_mydata$RNA_snn_res.0.05

## ElongatedWithBridge
#Seurat_mydata <- FindClusters(Seurat_mydata, resolution = 0.03)
#Seurat_labels <- Seurat_mydata$RNA_snn_res.0.03


## RNAMix1:
#Seurat_mydata <- FindClusters(Seurat_mydata, resolution = 0.11)
#Seurat_labels <- Seurat_mydata$RNA_snn_res.0.11

## RNAMix2:
# Seurat_mydata <- FindClusters(Seurat_mydata, resolution = 0.15)
# Seurat_labels <- Seurat_mydata$RNA_snn_res.0.15

## Beta:
# Seurat_mydata <- FindClusters(Seurat_mydata, resolution = 0.4)
# Seurat_labels <- Seurat_mydata$RNA_snn_res.0.4

## Baron Panc:
# Seurat_mydata <- FindClusters(Seurat_mydata, resolution = 0.06)
# Seurat_labels <- Seurat_mydata$RNA_snn_res.0.06

## TM Panc:
# Seurat_mydata <- FindClusters(Seurat_mydata, resolution = 0.02)
# Seurat_labels <- Seurat_mydata$RNA_snn_res.0.02

## TM Lung:
#Seurat_mydata <- FindClusters(Seurat_mydata, resolution = 0.2)
#Seurat_labels <- Seurat_mydata$RNA_snn_res.0.2

## PBMC3KSCT:
#Seurat_mydata <- FindClusters(Seurat_mydata, resolution = 0.05)
#Seurat_labels <- Seurat_mydata$RNA_snn_res.0.05

## PBMC4KBASIC:
#Seurat_mydata <- FindClusters(Seurat_mydata, resolution = 0.003)
#Seurat_labels <- Seurat_mydata$RNA_snn_res.0.003


## CellmixSNGSCT:
#Seurat_mydata <- FindClusters(Seurat_mydata, resolution = 0)
#Seurat_labels <- Seurat_mydata$RNA_snn_res.0

### Step 4: Clustering Evaluation
adjustedRandIndex(Seurat_labels, t(true_labels))
head(Idents(Seurat_mydata), 5)
clustering_evaluation(Seurat_labels,t(true_labels))

### Step 5:PCA, UMAP and TSNE plots
Seurat_mydata<-AddMetaData(Seurat_mydata, Seurat_labels, col.name = 'Seurat_labels')
x11()
p1<-DimPlot(Seurat_mydata, reduction = "pca",group.by='Seurat_labels',label=TRUE)
p2<-DimPlot(Seurat_mydata, reduction = "pca",group.by='true_labels')
p1+p2
x11()
Seurat_mydata <- RunUMAP(Seurat_mydata, dims = 1:min(nrow(mydata),40))
u1<-DimPlot(Seurat_mydata, reduction = "umap",group.by='Seurat_labels',label=TRUE)
u2<-DimPlot(Seurat_mydata, reduction = "umap",group.by='true_labels')
u1+u2
x11()
Seurat_mydata <- RunTSNE(Seurat_mydata, dims = 1:min(nrow(mydata),40),check_duplicates = FALSE)
t1<-DimPlot(Seurat_mydata, reduction = "tsne",group.by='Seurat_labels',label=FALSE)
t2<-DimPlot(Seurat_mydata, reduction = "tsne",group.by='true_labels')
t1+t2

