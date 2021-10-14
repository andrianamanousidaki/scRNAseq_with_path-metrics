#SCTranform version of data

library(sctransform)
library(Seurat)
library(patchwork)

#Choose dataset ater SAVER imputation

data<-read.csv('pancreatic_saver_no_normalization.csv')
#data<-read.csv('pancreas_tabula_muris_saver.csv')
#data<-read.csv('lung_tabula_muris_saver.csv')
#data<-read.csv('beta_3_4_10_filtered_saver.csv')
#data<-read.csv('cellmix_sng.csv')
#data<-read.csv('pbmc4k_qc_saver.csv')
#data<-read.csv("RNAmix1_original.csv")
#data<-read.csv("RNAmix2_original.csv")

rownames(data)<-data[,1]
data<-data[,-1]
dim(data)

Data <- CreateSeuratObject(counts = data,min.cell=0,min.feat=0)

#For tabula muris data sets: pattern="^mt-"
Data <- PercentageFeatureSet(Data, pattern = "^MT-", col.name = "percent.mt")
Data <- SCTransform(Data, vars.to.regress = "percent.mt")
tosave<-GetAssayData(Data,assay="SCT")
genes<-Data@assays$SCT@var.features

#Choose how to save output data set

write.csv(tosave[genes,],'Baron_Pancreatic_SCT.csv')
#write.csv(tosave[genes,],'TM_Pancreatic_SCT.csv')
#write.csv(tosave[genes,],'TM_Lung_SCT.csv')
#write.csv(tosave[genes,],'Beta_filtered_SCT.csv')
#write.csv(tosave[genes,],'Cellmix_sng_SCT.csv')
#write.csv(tosave[genes,],'pbmc4k_SCT.csv')
#write.csv(tosave[genes,],'rnamix2_SCT.csv')
#write.csv(tosave[genes,],'rnamix1_SCT.csv')






