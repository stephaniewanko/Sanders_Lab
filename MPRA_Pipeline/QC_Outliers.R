#This is the QC metrics and outlier script for the Sanders' lab MPRA pipeline
#last edited: 02/05/2019

##THIS STILL NEED TO BE CLEANED UP/COMMENTED/EDITED! 


library(tidyr)
library(gtools)
setwd('/Users/stephaniewankowicz/Downloads/')
#take in the Raw RNA and DNA counts
DNA<-read.csv('DNA_Rep3_combined_pickle.csv')
RNA<-read.csv('RNA_Rep3_combined_pickle.csv')
DNA_RNA_sum<-read.csv('Rep3_SummaryTable.csv')

print('DNA barcode summary:')
print(paste('Number of DNA barcodes:', nrow(DNA)))
print(paste('Number of Variants captured by DNA barcodes:',length(unique(DNA$Varaint))))
print(paste('Summary of the DNA barcode counts:', summary(DNA$count)))
print(paste('95th quantile of your DNA barcode counts (this is what the histogram will cut off at):', quantile(DNA$count, 0.95)))
#creating histogram
hist(DNA[DNA$count < quantile(DNA$count, 0.95), ]$count)
#DNA_few_barcodes2=setNames(data.frame(matrix(ncol = 2, nrow = 0)), c("Variant", "Barcode_Count"))
#FIX THIS SO IT IS A DF
DNA_few_barcodes=list()
for(i in unique(DNA$Varaint)){ #fix this (speed it up)
  rows=nrow(DNA[DNA$Varaint==i,])
  if (rows<=3){
    #print('bad!')
  DNA_few_barcodes[[i]]<-i
  }
}
print(paste('The number of DNA variants with less than 3 barcodes is:', length(DNA_few_barcodes)))
write.csv(DNA_few_barcodes, 'Rep3_DNA_variants_lessthan3barcodes.csv')

print('RNA barcode summary:')
print(paste('Number of RNA barcodes:', nrow(RNA)))
print(paste('Number of Variants captured by RNA barcodes:',length(unique(RNA$Varaint))))
print(paste('Summary of the RNA barcode counts:', summary(RNA$count)))
print(paste('95th quantile of your RNA barcode counts (this is what the histogram will cut off at):', quantile(RNA$count, 0.95)))
#creating histogram
hist(RNA[RNA$count < quantile(RNA$count, 0.95), ]$count)
RNA_few_barcodes=list()
for(i in unique(RNA$Varaint)){
  rows=nrow(RNA[RNA$Varaint==i,])
  if (rows<=3){
    #print('bad!')
    RNA_few_barcodes[[i]]<-i
  }
}
print(paste('The number of RNA variants with less than 3 barcodes is:', length(RNA_few_barcodes)))
write.csv(RNA_few_barcodes, 'Rep3_RNA_variants_lessthan3barcodes.csv')

#overlap of RNA and DNA variants:
count=0
for(i in RNA_few_barcodes){
  if(i%in% DNA_few_barcodes){
    count=count+1}}
print(paste('Number of ', count))


#summary
Variants<-as.data.frame(DNA_RNA_sum$Variant_x)
colnames(Variants)<-'Variant'
DNA_RNA_sum<-separate(data = DNA_RNA_sum, col = Variant_x, into = c('Type', 'Number','Location', 'description', 'Variant_WT', 'Details', 'unknown'), sep = "_")
DNA_RNA_sum$unknown<-NULL
DNA_RNA_sum<-cbind(DNA_RNA_sum,Variants)
print(paste('Mean RNA/DNA:', mean(DNA_RNA_sum$RNA.DNA)))
RNA.DNA_mean=mean(DNA_RNA_sum$RNA.DNA, na.rm = TRUE)
RNA.DNA_sd=sd(DNA_RNA_sum$RNA.DNA,na.rm = TRUE)#*sqrt((length(DNA_RNA_sum$RNA.DNA)-1)/(length(DNA_RNA_sum$RNA.DNA)))

z_score<-function(i){
  z<-(i-RNA.DNA_mean)/RNA.DNA_sd
  return(z)
}

head(DNA_RNA_sum)
DNA_RNA_sum$z_score<-do.call(rbind, lapply(DNA_RNA_sum$RNA.DNA, z_score))
hist(DNA_RNA_sum$z_score, xlim=c(-3,20), breaks = 100)
DNA_RNA_sum$RNA_low_count<-ifelse(DNA_RNA_sum$Variant%in%RNA_few_barcodes, 1, 0)
DNA_RNA_sum$DNA_low_count<-ifelse(DNA_RNA_sum$Variant%in%DNA_few_barcodes, 1, 0)

low_outliers_DNA_RNA<-DNA_RNA_sum[DNA_RNA_sum$z_score<=-1.5, ]
high_outliers_DNA_RNA<-DNA_RNA_sum[DNA_RNA_sum$z_score>=1.5, ]
high_outliers_DNA_RNA<-high_outliers_DNA_RNA[complete.cases(high_outliers_DNA_RNA), ]
low_outliers_DNA_RNA<-low_outliers_DNA_RNA[complete.cases(low_outliers_DNA_RNA), ]

print('The z-score cut off used: 1.5')

nrow(low_outliers_DNA_RNA)
nrow(high_outliers_DNA_RNA)
#positive and negative
DNA_RNA_pos<-DNA_RNA_sum[grep('Positive*', DNA_RNA_sum$Variant), ]
hist(DNA_RNA_pos$z_score)
summary(DNA_RNA_pos$z_score)
print(paste('The median z-score of the positive controls:', median(DNA_RNA_pos$z_score)))

DNA_RNA_neg<-DNA_RNA_sum[grep('Negative*', DNA_RNA_sum$Variant), ]
hist(DNA_RNA_neg$z_score)
summary(DNA_RNA_neg$z_score)
print(paste('The median z-score of the negative controls:', median(DNA_RNA_neg$z_score)))



count=0
high_outliers_DNA_RNA_wide=data.frame()
high_outliers_DNA_RNA_nodup=data.frame()

for(i in unique(high_outliers_DNA_RNA$description)){
  print(i)
  subset<-high_outliers_DNA_RNA[high_outliers_DNA_RNA$description==i, ]
  if(nrow(subset)==2){
    #head(subset)
    subset_var<-subset[subset$Variant_WT=='Variant', ]
    subset_WT<-subset[subset$Variant_WT=='Wildtype', ]
    #head(subset_WT)
    #head(subset_var)
    colnames(subset_var)<-c('Type_var', 'Numb_var', 'Location', 'Description', 'Variant', 'Details_var', 'DNA_count_var', 'RNA_count_var', 'RNA.DNA_var', 'Variant_var', 'z_score_var', 'RNA_low_count_var', 'DNA_low_count_var')
    colnames(subset_WT)<-c('Type_WT', 'Numb_WT', 'Location', 'Description', 'Wildtype', 'Details_WT', 'DNA_count_WT', 'RNA_count_WT', 'RNA.DNA_WT', 'Variant_WT', 'z_score_WT','RNA_low_count_WT', 'DNA_low_count_WT')
    head(subset_var)
    #tmp<-cbind(subset_var, subset_WT)
    #print(tmp)
    high_outliers_DNA_RNA_wide<-smartbind(high_outliers_DNA_RNA_wide, tmp)
    rm(subset)
    count=count+1
  }
  else{
    high_outliers_DNA_RNA_nodup<-rbind(high_outliers_DNA_RNA_nodup,subset)
  }
}

write.csv(high_outliers_DNA_RNA_wide, 'Rep3_high_outliers_DNA_RNA_wide.csv')
write.csv(high_outliers_DNA_RNA_nodup, 'Rep3_high_outliers_DNA_RNA_nodup.csv')

count=0
low_outliers_DNA_RNA_wide=data.frame()
low_outliers_DNA_RNA_nodup=data.frame()
for(i in unique(low_outliers_DNA_RNA$description)){
  subset<-low_outliers_DNA_RNA[low_outliers_DNA_RNA$description==i, ]
  if(nrow(subset)==2){
    subset_var<-subset[subset$Variant_WT=='Variant', ]
    subset_WT<-subset[subset$Variant_WT=='Wildtype', ]
    colnames(subset_var)<-c('Type_var', 'Numb_var', 'Location', 'Description', 'Variant', 'Details_var', 'DNA_count_var', 'RNA_count_var', 'RNA.DNA_var', 'Variant_var', 'z_score_var', 'RNA_low_count_var', 'DNA_low_count_var')
    colnames(subset_WT)<-c('Type_WT', 'Numb_WT', 'Location', 'Description', 'Wildtype', 'Details_WT', 'DNA_count_WT', 'RNA_count_WT', 'RNA.DNA_WT', 'Variant_WT', 'z_score_WT', 'RNA_low_count_WT', 'DNA_low_count_WT')
    tmp<-cbind(subset_var, subset_WT)
    print(head(tmp))
    low_outliers_DNA_RNA_wide<-cbind(low_outliers_DNA_RNA_wide, tmp)
    rm(subset)
    count=count+1
  }
  else{
    low_outliers_DNA_RNA_nodup<-rbind(low_outliers_DNA_RNA_nodup,subset)
  }
}


###Looking at WT versus mutant scores



write.csv(low_outliers_DNA_RNA_wide, 'Rep3_low_outliers_DNA_RNA_wide.csv')
write.csv(low_outliers_DNA_RNA_nodup, 'Rep3_low_outliers_DNA_RNA_nodup.csv')
