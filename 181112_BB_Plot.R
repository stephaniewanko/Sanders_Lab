#Stephanie Wankowicz
#2018_11_01
'''
This script will take in VEP annotated variants, which include cases and controls, and output a bbplot. 
The script is currently made to look at variants with 5kb of the TSS
Currently, many thing in this script is hardcoded, and will need to be adjusted for other applications.
'''

#packages needed
library(plotly)
require(epitools)

#functions
fixScores <- function(df, scores){
  df[,scores] <- ifelse(df[,scores]=='', NA, df[,scores])
  df[,scores] <- as.numeric(sapply(strsplit(df[,scores], "&"),max))
  new_scores = paste('Var',scores,sep='_')
  df[,new_scores] = df[,scores]
  return (df)
}

preprocess_df<-function(df){
  df = df[grepl("upstream|downstream", df$Consequence),]
  #case/control
  df$isCase <- ifelse( grepl('p1', df$SampleID), 1, 0)
  df[is.na(df)] <- '' # Fix all na values to blank
  print(nrow(df))
  df$ID<-as.character(df$ID)
  print('Annotating Indels')
  df$isIndel <- ifelse(nchar(as.character(do.call(rbind.data.frame, strsplit(df$ID, ':')))[[3]]) == 1 & nchar(as.character(do.call(rbind.data.frame, strsplit(df$ID, ':')))[[4]]) == 1, 0, 1)
  df$variant_loc <- as.numeric(do.call(rbind.data.frame, strsplit(df$ID, ':'))[[2]])
  # Add conservation scores (from Joon's script)
  print('Fixing Conservation Scores')
  df <- fixScores(df, 'phastCons46wayVt')
  df <- fixScores(df, 'phyloP46wayVt')
  df$isPhyloP <- ifelse(!is.na(df$phyloP46wayVt) & df$phyloP46wayVt>=2, 1, 0)
  df$isPhast <- ifelse(!is.na(df$phastCons46wayVt) & df$phastCons46wayVt>=0.2, 1, 0)
  df$isCons <- ifelse(df$isPhast==1 | df$isPhyloP==1, 1, 0)
  print('Fixing ID labels')
  df$ID2 <- paste(df$ID, df$SampleID, sep=':')
  df$DISTANCE<-as.numeric(df$DISTANCE)
  df$DISTANCE<-ifelse(df$Consequence=='downstream_gene_variant',(df$DISTANCE*(-1)),(df$DISTANCE))
  df<-df[df$isCons=='1',] #subset down to only conserved regions
}

bb_plot_preprocessing<-function(df, window_size, distance){
  output = as.data.frame(matrix(nrow=0, ncol=8))
  colnames(output) <- c('Start_Dist', 'End_Dist', 'Mean_Distance','num_P', 'num_S', 'RR', 'Total', 'direction')
  n=1
  num_cases=length(unique(df_con[df_con$isCase==1,]$SampleID))
  num_controls=length(unique(df_con[df_con$isCase==0,]$SampleID))
  neg_dist=distance*(-1)
  end=distance-window_size
  print('Calculaing the number of alterations in each window.')
  for (i in neg_dist:end){
    output[n,]$Start_Dist<-i
    output[n,]$End_Dist<-i+window_size
    subset=df_con[df_con$DISTANCE>i & df_con$DISTANCE<=i+window_size,]
    num_P<-nrow(subset[subset$isCase==1,])
    num_S<-nrow(subset[subset$isCase==0,])
    output[n,]$num_P<-nrow(subset[subset$isCase==1,])
    output[n,]$num_S<-nrow(subset[subset$isCase==0,])
    output[n,]$Total<-nrow(subset)
    rm(subset)
    n=n+1
  }
  output$RR=(output$num_P/num_controls)/(output$num_S/num_cases)
  output$Mean_Distance=(output$Start_Dist+output$End_Dist)/2
  output$Total_P = num_controls
  output$Total_S = num_cases
  output$CI_low_95<-exp(log(output$RR)-(1.96*(sqrt((((output$Total_P-output$num_P)/output$num_P)/output$Total_P)+((output$Total_S-output$num_S)/output$num_S)/output$Total_S))))
  output$CI_high_95<-exp(log(output$RR)+(1.96*(sqrt((((output$Total_P-output$num_P)/output$num_P)/output$Total_P)+((output$Total_S-output$num_S)/output$num_S)/output$Total_S))))
  output$CI_high_80<-exp(log(output$RR)+(1.285*(sqrt((((output$Total_P-output$num_P)/output$num_P)/output$Total_P)+((output$Total_S-output$num_S)/output$num_S)/output$Total_S))))
  output$CI_low_80<-exp(log(output$RR)-(1.285*(sqrt((((output$Total_P-output$num_P)/output$num_P)/output$Total_P)+((output$Total_S-output$num_S)/output$num_S)/output$Total_S))))
  return(output)
}

#create candlesick plot
plot_candlesick<-function(title,df){
  f <- list(family = "Arial", size = 18, color = "#7f7f7f")
  title=list(title=title, titlefont=f)
  x=list(title="Average Distance from TSS", titlefont=f, list(autorange = "reversed"))
  y=list(title="Relative Risk", titlefont=f)
  p <- df %>%
    plot_ly(x = ~rev(Mean_Distance), type="candlestick",
            open = ~rev(CI_high_95), close = ~rev(CI_low_95),
            high = ~rev(CI_high_95), low = ~rev(CI_low_95)) %>%
    add_lines(x = ~rev(Mean_Distance), y = ~rev(RR), line = list(color = 'black', width = 1.25, title='Relative Risk'), inherit = F) %>%
    layout(title = title, xaxis=list(autorange = "reversed"), yaxis=y)
}



#read in data
setwd('/Users/stephaniewankowicz/Dropbox/Stephanie_TSS/')
df0 = read.delim('list_gene.hqDNV_sscwgs.20171115.P231_WGS519_256.hg38.vep_gene_5kb_simpleAnnot_20180924.txt')

#preprocess the dataframe
df_con<-preprocess_df(df0)
#calculate the CI and windows
df_bbplot<-bb_plot_preprocessing(df_con, 200, 5000)
p<-plot_candlesick('5000kb Up and Downstream from TSS',df_bbplot)
htmlwidgets::saveWidget(as_widget(p), "5000kb_TSS.html")
write.csv(df_bbplot, '5000kb_TSS_RR_Table.csv')
