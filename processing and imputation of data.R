# Preparation for SAVER Imputation

library(SAVER)
library(dplyr)
library(Seurat)
library(patchwork)

##################################################################################
#Lung Tabula Muris - no preparation needed
################################################################################
# Download data file "https://ndownloader.figshare.com/files/10700143"

tm_lung <-read.csv('Lung-counts.csv')
rownames(tm_lung)<-tm_lung[,1]
tm_lung<-tm_lung[,-1]
tm_lung_saver<-saver(as.matrix(tm_lung),ncores=20,size.factor=1)
tm_lung_saver<-tm_lung_saver$estimate
saveRDS(tm_lung_saver, 'lung_saver.rds')

annotation<-read.csv('annotations_facs.csv')
counts<-readRDS("lung_saver.rds")
cell_id<-which( colnames(counts) %in% annotation$cell)
counts<-counts[,cell_id]
names<-colnames(counts)
labels<-cbind(annotation$cell,annotation$free_annotation)
rownames(labels)<-annotation$cell
labels<-labels[names,]
labels<-labels[which(labels[,2] != ""),]
counts<-counts[,rownames(labels)]

write.csv(counts,'lung_tabula_muris_saver.csv')


labels[,2]<-revalue(labels[,2],c("alveolar epithelial type 1 cells, alveolar epithelial type 2 cells, club cells, and basal cells"=
                                 "alveolar ep 1&2,club & basal",
                              "dendritic cells, alveolar macrophages, and interstital macrophages"=
                               "dendritic,alveolar macro,interstital macro"))
num_labels<-revalue(labels[,2],c("lung neuroendocrine cells and unknown cells"="1",
                                 "alveolar ep 1&2,club & basal"  ="2"             ,
                                 "mast cells and unknown immune cells"="3"        ,
                                 "circulating monocytes" = "4"                      ,
                                 "invading monocytes"  = "5"                       ,
                                 "multiciliated cells"  = "6"                      ,
                                 "dendritic,alveolar macro,interstital macro" ="7"))

write.csv(num_labels,"lung_tabula_muris_labels.csv")

cores<-cbind(1:7,c("lung neuroendocrine cells and unknown cells",
                   "alveolar ep 1&2,club & basal"          ,
                   "mast cells and unknown immune cells"      ,
                   "circulating monocytes"                       ,
                   "invading monocytes"                     ,
                   "multiciliated cells"                       ,
                   "dendritic,alveolar macro,interstital macro" ))
write.csv(cores,'num_to_word_tabula_muris_lung.csv')


################################################################################
#Pancreas Tabula Muris - no preparation needed
################################################################################

# Download data file "https://ndownloader.figshare.com/files/10700143"

tm_panc <-read.csv('Pancreas-counts.csv')
rownames(tm_panc)<-tm_panc[,1]
tm_panc<-tm_panc[,-1]
tm_panc_saver<-saver(tm_panc,ncores=20,size.factor=1)
tm_panc_saver<-tm_panc_saver$estimate
saveRDS(tm_panc_saver, 'pancreas_saver.rds')

##labels
annotation<-read.csv('annotations_facs.csv')

counts<-readRDS('pancreas_saver.rds')
cell_id<-which( colnames(counts) %in% annotation$cell)
counts<-counts[,cell_id]
names<-colnames(counts)
labels<-cbind(annotation$cell,annotation$free_annotation)
rownames(labels)<-annotation$cell
labels<-labels[names,]
labels<-labels[which(labels[,2] != ""),]
counts<-counts[,rownames(labels)]
write.csv(counts,'pancreas_tabula_muris_saver.csv')

labels[,2]<-revalue(labels[,2],c("acinar cell"="1",
                                 "beta cell"="2",
                                 "stellate cell"="3",
                                 "pancreatic A cell"="4",
                                 "ductal cell"="5",
                                 "pancreatic D cell"="6",
                                 "pancreatic PP cell"="7"))


write.csv(num_labels,"pancreatic_tabula_muris_labels.csv")




################################################################################
#Beta Simulated Dataset 
################################################################################

beta<-readRDS("beta_sim_34_10_for_saver.rds")
beta_saver<- saver(beta,ncores=20)
saveRDS(beta_saver$estimate,'beta_sim_filtered_saver.rds')
# for labels see Simulationscript.R 

################################################################################
# PBMC3k(not inckuded in the paper)
################################################################################

# pbmc.data <- readRDS('pbmc.data.rds')
# pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
# pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
# pbmc <- subset(pbmc, subset = percent.mt < 20)
# saveRDS(GetAssayData(pbmc),'pbmc3k_for_saver.rds')
# pbmc_data<-readRDS('pbmc3k_qc_saver.rds')
# pbmc_saver<- saver(pbmc_data,ncores=20)
# saveRDS(pbmc_saver$estimate,'pbmc3k_qc_saver.rds')


################################################################################
# RNAmix1 and RNAmix2
################################################################################
#Download data from 'https://github.com/LuyiTian/sc_mixology/tree/master/data'

load('mRNAmix_qc.RData')
RNAmix1<-sce8_qc
l1<-paste(sce8_qc@colData$H2228_prop,sce8_qc@colData$H1975_prop,sce8_qc@colData$HCC827_prop,sep="_")
c1<-sce8_qc@colData$cell_name
l1<-cbind(c1,l1)

RNAmix2<-sce2_qc
l2<-paste(sce2_qc@colData$H2228_prop,sce2_qc@colData$H1975_prop,sce2_qc@colData$HCC827_prop,sep="_")
c2<-sce2_qc@colData$cell_name
l2<-cbind(c2,l2)

rnamix1<-as.Seurat(RNAmix1,data=NULL)
saveRDS(GetAssayData(rnamix1),"rnamix1_filtered.rds")  

rnamix2<-as.Seurat(RNAmix2,data=NULL)
saveRDS(GetAssayData(rnamix2),"rnamix2_filtered.rds")  

rnamix1_saver<-saver(GetAssayData(rnamix1),ncores=20,size.factor=1)
saveRDS(rnamix1_saver$estimate,"RNAmix1_original_saver.rds")

rnamix2_saver<-saver(GetAssayData(rnamix2),ncores=20,size.factor=1)
saveRDS(rnamix2_saver$estimate,"RNAmix2_original_saver.rds")


################################################################################
#Baron's Pancreatic data
################################################################################
#Download data from 'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2230757'

x <- read.csv("GSM2230757_human1_umifm_counts.csv.gz", header = TRUE, check.names = FALSE)
cell_clusters<-x[,2:3]
rownames(cell_clusters)<-make.names(x[,2], unique=TRUE)
cell_clusters[,1]<-rownames(cell_clusters)

x.dat <- t(as.matrix(x[, 4:ncol(x)]))
colnames(x.dat) <- make.names(x[,2], unique=TRUE)
x <- x.dat



filtered_data<-seurat_analysis(counts=x,processing =TRUE,project="Baron's Pancreatic data",
                               normalization = FALSE,clustering =FALSE,
                               signature_genes = FALSE)

data_for_saver<-GetAssayData(filtered_data)
saveRDS(data_for_saver, 'pancreatic_for_saver.rds')

panc_saver<-saver(data_for_saver,ncores=20,size.factor=1)
saveRDS(panc_saver,'pancreatic_saver_no_normalization.rds')


################################################################################
# PBMC4k
################################################################################
#Download data from 'https://support.10xgenomics.com/single-cell-gene-expression/datasets/2.1.0/pbmc4k'
#and SingleR annotation from 'https://github.com/dviraran/SingleR/blob/master/manuscript_figures/FiguresData/SingleR.PBMC.4K.RData'

## Find cell names annotated by SingleR
load('SingleR.PBMC.4K.RData')
singler2 = singler$singler[[2]]
sr_cells<-singler2$SingleR.single$cell.names

#Prepare data set for SAVER
pbmc.data <- Read10X(data.dir = "filtered_gene_bc_matrices/GRCh38/")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc4k", min.cells = 3, min.features = 200)
cell <- colnames(pbmc)
cell<-substr(cell,1,nchar(cell)-2)
id<-which(cell%in%sr_cells)
pbmc <- subset(pbmc, cells=colnames(pbmc)[id])
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = percent.mt < 20)
saveRDS(GetAssayData(pbmc),'pbmc4k_for_saver.rds')

#Imputation
c<-readRDS("pbmc4k_for_saver.rds")
set.seed(1)
saverc<- saver(c,ncores=20)
saveRDS(saverc$estimate,'pbmc4k_qc_saver.rds')

data<-readRDS('pbmc4k_qc_saver.rds')
write.csv(data,'pbmc4k_qc_saver.csv')

#Labels

pbmc_cells<-colnames(data)
pbmc_cells<-substr(pbmc_cells,1,nchar(pbmc_cells)-2)

main_labels<-data.frame(cells=singler2$SingleR.single.main$cell.names,labels=singler2$SingleR.single.main$labels)
id_cells<-c()

for( i in 1:length(pbmc_cells)){
  id_cells[i]<-which(main_labels[,1]== pbmc_cells[i])
}

main_labels<-main_labels[id_cells,]
library(plyr)

num_main_labels<-revalue(main_labels[,2],c("Monocytes"="1","CD8+ T-cells"="2","CD4+ T-cells"="3",
                                           "B-cells"="4","HSC"="5","NK cells"="6","DC"="7"))

main_labels$num_labels=num_main_labels

#Check that colnames in data set are in correct order
identical(as.character(main_labels$cells),as.character(pbmc_cells))
identical(as.character(main_labels$cells),as.character(rownames(main_labels)))
main_labels<-main_labels[pbmc_cells,]

#Word to numerical labels
d<-data.frame(levels=unique(main_labels[,3]),num=1:7)
write.csv(d,row.names = FALSE,'pbmc4k_num_to_main_labels.csv')
write.csv(num_main_labels,row.names = FALSE,'pbmc4k_main_labels.csv')

#Remove cluster 5,7 
l<-read.csv('pbmc4k_main_labels.csv')
identical(as.character(l$x),as.character(main_labels$num_labels))
rownames(l)<-colnames(data)

ce<-which(l$x!=5 & l$x!=7)
data<-data[,rownames(l)[ce]]
colnames(data)<-substr(colnames(data),1,nchar(colnames(data))-2)
write.csv(data,'subset_pbmc4k.csv')

#Merge the two T cell cell types
l$x<-revalue(as.character(l$x),c("3"="2"))

write.csv(l[ce,1],'subset_pbmc4k_labels.csv')




################################################################################
# CellmixSNG
################################################################################
##No imputation was applied to CELLMIXSNG data set 

#Download data from  'https://github.com/LuyiTian/sc_mixology/tree/master/data'

load("sincell_with_class_5cl.RData")

## Find singletons of CELLMIX

data<-sce_sc_10x_5cl_qc
metadata<-read.csv('sc_10x_5cl.metadata.csv.gz')
counts<-data@assays$data$counts
sng_index<-data.frame(rownames(metadata),metadata$demuxlet_cls)
rownames(sng_index)<-sng_index[,1]

counts<-counts[,which(sng_index[,2]=='SNG')]
dim(counts)
write.csv(counts,"cellmix_sng.csv")


## Filtered labels

info<-sce_sc_10x_5cl_qc@colData@listData
cell_labels<-info$cell_line_demuxlet
c_labels<-cbind(colnames(sce_sc_10x_5cl_qc), cell_labels)
rownames(c_labels)<-c_labels[,1]
filtered_labels<-c_labels[which(sng_index[,2]=='SNG'),]
library(plyr)
revallab<-revalue(filtered_labels[,2], c("HCC827"=1, "H1975"=2 , "H838"=3 ,  "H2228"=4 , "A549"=5))
write.csv(revallab,"cellmix_sng_labels.csv")






