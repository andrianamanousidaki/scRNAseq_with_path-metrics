library(dplyr)
library(Seurat)
library(patchwork)
library(mclust)



source("clustering_evaluation.R")

### Step 1: Choose data set from Original Datasets:

#choose dataset

## ElongatedGaussiansWithBridge3
#mydata <- read.csv(file = 'ElongatedGaussiansWithBridge3_Denoised_2d.csv',header = FALSE)
#mydata <- read.csv(file = 'ElongatedGaussiansWithBridge3_Denoised.csv',header = FALSE)
#true_labels <- read.csv(file = 'ElongatedGaussiansWithBridge3_Denoised_Labels.csv',header = FALSE)

## RNAMix1
#dataset='RNAMix1'
#mydata <- read.csv(file = 'RNAmix1_original.csv')
#true_labels <- read.csv(file = 'rnamix1_original_labels.csv')

##RNAMix2
#dataset='RNAMix2'
#mydata <- read.csv(file = 'RNAmix2_original.csv')
#true_labels <- read.csv(file = 'rnamix2_original_labels.csv')

## TMLung
#dataset='TMLung'
#mydata <- read.csv(file = 'lung_tabula_muris_saver.csv')
#true_labels <- read.csv(file = 'lung_tabula_muris_labels.csv')

## Beta2
#dataset='Beta2'
#mydata <- read.csv(file = 'beta_3_4_10_filtered_saver.csv')
#true_labels <- read.csv(file = 'beta_cell_groups_3_4_10_after_filtering.csv')
#true_labels<-true_labels[,-1]

## TMPanc
#dataset='TMPanc'
#mydata <- read.csv(file = 'pancreas_tabula_muris_saver.csv')
#true_labels <- read.csv(file = 'pancreatic_tabula_muris_labels.csv')

## BaronPancSCT
#dataset='BaronPancSCT'
#mydata <- read.csv(file = 'BaronPancSCT_filtered.csv')
#true_labels <- read.csv(file = 'pancreatic_num_labels_filtered.csv')


## PBMC3kSCT
#dataset='PBMC3kSCT'
#mydata <- read.csv(file = 'pbmc3k_SCT_filtered.csv')
#true_labels <- read.csv(file = 'pbmc3k_labels_filtered.csv')

## PBMC3kSCT
#dataset='PBMC3kSCT_singleR'
#mydata <- read.csv(file = 'pbmc3k_SCT_singleR_filtered.csv')
#true_labels <- read.csv(file = 'pbmc3k_singleRlabels_filtered.csv')

## PBMC4kBASIC
dataset='PBMC4k'
mydata <- read.csv(file = 'subset_pbmc4k.csv')
true_labels <- read.csv(file = 'subset_pbmc4k_labels.csv')


## CellMixSngSCT
#dataset='CellMixSngSCT'
#mydata <- read.csv(file = 'Cellmix_sng_SCT.csv')
#true_labels <- read.csv('cellmix_sng_labels.csv')

### Step 3: Create Seurat Object 

rownames(mydata) = mydata[,1]
mydata<-mydata[,-1]
true_labels<-true_labels[,-1]

Seurat_mydata <- CreateSeuratObject(counts = mydata, project = dataset, min.cells = 0, min.features = 0)
n_hvg=2000
dim=50

if(dataset=='TMPanc' | dataset=='TMLung'){ 
  Seurat_mydata[["percent.mt"]] <- PercentageFeatureSet(Seurat_mydata, pattern = "^mt-")
}else 
 {Seurat_mydata[["percent.mt"]] <-PercentageFeatureSet(Seurat_mydata, pattern = "^MT-")}

if(dataset!='BaronPancSCT' && dataset!='PBMC3kSCT' && dataset!='CellMixSngSCT'&& dataset!='PBMC3kSCT_singleR'){ 
  Seurat_mydata <- NormalizeData(Seurat_mydata)
  Seurat_mydata <- FindVariableFeatures(Seurat_mydata, selection.method = "vst", nfeatures = n_hvg)
  Seurat_mydata <- ScaleData(Seurat_mydata, features = rownames(Seurat_mydata))
  Seurat_mydata<- RunPCA(Seurat_mydata, features = VariableFeatures(object = Seurat_mydata))
  Seurat_mydata <- JackStraw(Seurat_mydata, num.replicate = 100, dims = dim)
  Seurat_mydata <- ScoreJackStraw(Seurat_mydata, dims = 1:dim)
}else 
{Seurat_mydata <- ScaleData(Seurat_mydata, features = rownames(mydata))
Seurat_mydata@assays$RNA@scale.data <- as.matrix(mydata)
Seurat_mydata<- RunPCA(Seurat_mydata, features = rownames(Seurat_mydata))
Seurat_mydata <- JackStraw(Seurat_mydata, num.replicate = 100, dims = dim)
Seurat_mydata <- ScoreJackStraw(Seurat_mydata, dims = 1:dim)
}

#find important pcs
pcid<-which(Seurat_mydata@reductions$pca@jackstraw@overall.p.values[,2]<0.001)
#find neighbors
Seurat_mydata <- FindNeighbors(Seurat_mydata, dims =pcid )
Seurat_mydata<-AddMetaData(Seurat_mydata, true_labels, col.name = 'true_labels')

### Step 3: Clustering

## Elongated Gaussians:
#Seurat_mydata <- FindClusters(Seurat_mydata, resolution = 0.03)
#Seurat_labels <- Seurat_mydata$RNA_snn_res.0.03

## RNAMix1:
#Seurat_mydata <- FindClusters(Seurat_mydata, resolution = 0.5)
#Seurat_labels <- Seurat_mydata$RNA_snn_res.0.5

## RNAMix2:
# Seurat_mydata <- FindClusters(Seurat_mydata, resolution = 0.5)
# Seurat_labels <- Seurat_mydata$RNA_snn_res.0.5

## TM Lung:
#Seurat_mydata <- FindClusters(Seurat_mydata, resolution = 0.5)
#Seurat_labels <- Seurat_mydata$RNA_snn_res.0.5


## Beta:
# Seurat_mydata <- FindClusters(Seurat_mydata, resolution = 0.3)
# Seurat_labels <- Seurat_mydata$RNA_snn_res.0.3

## TM Panc:
# Seurat_mydata <- FindClusters(Seurat_mydata, resolution = 0.3)
# Seurat_labels <- Seurat_mydata$RNA_snn_res.0.3

## Baron Panc:
# Seurat_mydata <- FindClusters(Seurat_mydata, resolution = 0.15)
# Seurat_labels <- Seurat_mydata$RNA_snn_res.0.15

## PBMC3KSCT:
#Seurat_mydata <- FindClusters(Seurat_mydata, resolution = 0.12)
#Seurat_labels <- Seurat_mydata$RNA_snn_res.0.12

## PBMC3KSCT_singleR:
#Seurat_mydata <- FindClusters(Seurat_mydata, resolution = 0.12)
#Seurat_labels <- Seurat_mydata$RNA_snn_res.0.12

## PBMC4K :
Seurat_mydata <- FindClusters(Seurat_mydata, resolution = 0.007)
Seurat_labels <- Seurat_mydata$RNA_snn_res.0.007


## CellmixSNGSCT:
#Seurat_mydata <- FindClusters(Seurat_mydata, resolution = 0.01)
#Seurat_labels <- Seurat_mydata$RNA_snn_res.0.01

### Step 4: Clustering Evaluation
clustering_evaluation(Seurat_labels,t(true_labels))

### Step 5:PCA, UMAP and TSNE plots
Seurat_mydata<-AddMetaData(Seurat_mydata, Seurat_labels, col.name = 'Seurat_labels')

p1<-DimPlot(Seurat_mydata, reduction = "pca",group.by='Seurat_labels',label=TRUE)
p2<-DimPlot(Seurat_mydata, reduction = "pca",group.by='true_labels')
p1+p2

Seurat_mydata <- RunUMAP(Seurat_mydata, dims =pcid)
u1<-DimPlot(Seurat_mydata, reduction = "umap",group.by='Seurat_labels',label=TRUE)
u2<-DimPlot(Seurat_mydata, reduction = "umap",group.by='true_labels')
u1+u2

Seurat_mydata <- RunTSNE(Seurat_mydata, dims = pcid,check_duplicates = FALSE)
t1<-DimPlot(Seurat_mydata, reduction = "tsne",group.by='Seurat_labels',label=FALSE)
t2<-DimPlot(Seurat_mydata, reduction = "tsne",group.by='true_labels')
t1+t2

