##                                 Information                                        ##
##The simulations are part of a project of:Andriana Manousidaki,Anna Little, Yuying Xie
##Author of R script: Andriana Manousidaki
##email:manousid@msu.edu
##Department of Statistics and Probability, MSU
#
#########################################################################################

library(SAVER)
library(dplyr)
library(Seurat)
library(ggplot2)
library(cowplot)


## Create the reference dataset  based on https://github.com/mohuangx/SAVER-paper/blob/master/1-filter_data.R  ##

#Find "GSM2230757_human1_umifm_counts.csv.gz" in https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2230757

x <- read.csv("GSM2230757_human1_umifm_counts.csv.gz", header = TRUE, check.names = FALSE)
cell_clusters<-x[,2:3]
rownames(cell_clusters)<-make.names(x[,2], unique=TRUE)
cell_clusters[,1]<-rownames(cell_clusters)

x.dat <- t(as.matrix(x[, 4:ncol(x)]))
colnames(x.dat) <- make.names(x[,2], unique=TRUE)
x <- x.dat


ercc <- which(grepl("ERCC", rownames(x), ignore.case = TRUE))
rgenes <- which("r_" == substring(rownames(x), 1, 2))
mt <- which("mt-" == substring(rownames(x), 1, 3))
x1 <- x[-c(ercc, rgenes, mt), ]

X11()
plot(density(log10(colSums(x1))))
summary(colSums(x1))
dim(x1)

# no need to filter by cells

# look at mean gene expression
X11()
plot(density(log10(rowMeans(x1))))
abline(v = log10(0.001))
# Remove genes with mean expression less than 0.001

x2 <- x1[rowMeans(x1) >= 0.001, ]
dim(x2)

# Look at nonzero cells
plot(density(log10(rowSums(x2 != 0))))
abline(v = log10(3))
# Remove genes with less than 3 non-zero cells

x3 <- x2[rowSums(x2 != 0) >= 3, ]
dim(x3)


# build reference
lib.size <- colSums(x3)
non.zero.prop <- apply(x3, 1, function(x) sum(x != 0)/length(x))

cells.filt <- which(lib.size > 5000)
genes.filt <- which(non.zero.prop > 0.25)

data.filt <- x3[genes.filt, cells.filt]
cluster_ids <- cell_clusters[cells.filt,]
write.csv(data.filt,"reference_barons.csv")

# Focus on beta cells

set.seed(5)
beta.cells<-cluster_ids[cluster_ids[,2]=="beta",]
filtered<-data.filt[,beta.cells[,1]]

saveRDS(filtered,"subset_reference.rds")
write.csv(filtered,"subset_reference.csv")

###Simulations###

##Based on the simulation method of https://www.nature.com/articles/s41467-018-03405-7#Sec8


###      Scaling of the subset of the reference ##############
baron<-readRDS("subset_reference.rds")


c<-1:dim(baron)[2]
set.seed(5)
mix<-sample(c,dim(baron)[2],replace=FALSE)
sim<-baron

# Data frame sim has shuffled columns ->cells. 
for(i in c){
  sim[,i]= baron[, mix[i]]
  colnames(sim)[i]= colnames(baron)[mix[i]]
}

# Randomly choosing 210 marker genes, that is roughly 10% of the genes
set.seed(1)
m.gene.id<-sample(1:dim(baron)[1], 210, replace = FALSE)

markers<-rownames(sim)[m.gene.id]

#Rearranging the rows so that the first 210 rows correspond to the marker genes
s<-sim
sim<-rbind(sim[m.gene.id,],sim[-m.gene.id,])
s<-rbind(s[m.gene.id,],s[-m.gene.id,])

#Scaling of mean expression
#change values of first cluster
k<-round(dim(sim)[2]/3)

for(i in 1:70){
  set.seed(i)
  s[i,1:k]<-s[i,1:k]+runif(1,3,4)*mean(s[i,1:k])
}

#change values of second cluster

for(i in 71:140){
  set.seed(i)
  s[i,k+1:2*k]<-s[i,k+1:2*k]+runif(1,3,4)*mean(s[i,k+1:2*k])}


#change values of third cluster

for(i in 141:210){
  set.seed(i)
  s[i,(2*k+1):dim(sim)[2]]<-s[i,(2*k+1):dim(sim)[2]]+runif(1,3,4)*mean(s[i,(2*k+1):dim(sim)[2]])
}

#Marker genes and cell groups

mtable<-data.frame(cbind(rownames(sim[1:210,]),rep(c(1,2,3),each=70)))
colnames(mtable)<-c("gene","cluster")
ctable<-data.frame(cbind(colnames(sim),rep(c(1,2,3),times=c(k,k,k+1))))
colnames(ctable)<-c("cell","cluster")

write.csv(mtable,"beta_marker_genes_3_4_10.csv")
write.csv(ctable,"beta_cell_groups_3_4_10.csv")




#Down-sampling

alpha <- rgamma(ncol(s), 10, 100)
data.samp <- t(apply(sweep(s, 2, alpha, "*"), 1, function(x)
  rpois(length(x), x)))

colnames(data.samp) <- colnames(s)
down.zero<-length(which(data.samp==0))/(dim(data.samp)[1]*dim(data.samp)[2])
mean.exp.original<-mean(x3)
mean.exp.down<-mean(data.samp)

#SCALED DATA

saveRDS(s, "scaled_data.rds")
write.csv(s, "scaled_data.csv", quote = FALSE)

## Final simulated data set

write.csv(data.samp, "beta_simulated_data_3_4Gamma_10_100.csv", quote = FALSE)
saveRDS(data.samp,"beta_simulated_3_4Gamma_10_100.rds")



