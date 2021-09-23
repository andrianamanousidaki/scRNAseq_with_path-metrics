library(Linnorm)

#choose dataset
data<-read.csv('pancreatic_saver_no_normalization.csv')
#data<-read.csv('pancreas_tabula_muris_saver.csv')
#data<-read.csv('lung_tabula_muris_saver.csv')
#data<-read.csv('beta_3_4_10_filtered_saver.csv')
#data<-read.csv('sc_10x_5cl_no_preprocessing.csv')
#data<-read.csv('pbmc3k_qc_saver.csv')
#data<-read.csv('cellmix_sng.csv')
#data<-read.csv('pbmc4k_qc_saver.csv')

rownames(data)<-data[,1]
data<-data[,-1]

nor_data<-Linnorm(data)

#Choose output
write.csv(nor_data,'pancreatic_saver_Linnorm.csv')
#write.csv(nor_data,'pancreas_tabula_muris_saver_Linnorm.csv')
#write.csv(nor_data,'lung_tabula_muris_saver_Linnorm.csv')
#write.csv(nor_data,'beta_3_4_10_filtered_saver_Linnorm.csv')
#write.csv(nor_data,'sc_10x_5cl_Linnorm.csv')
#write.csv(nor_data,'pbmc3k_Linnorm.csv')
#write.csv(nor_data,'cellmix_sng_Linnorm.csv')
#write.csv(nor_data,'pbmc4k_Linnorm.csv')
