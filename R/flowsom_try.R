library(plyr)
library(tidyr)
setwd("/Users/s1249052/PhD/flow\ data/FlowSOM_Baseline_v_DOD_v_T+6_results/results/cluster_medians")

### list all csv files in working directory, import as list

temp = list.files(pattern="*.csv")
myfiles = lapply(temp, read.csv)

### import random file to extract colnames
thirsty<-read.csv("C-1_301_cluster_medians.csv")

### add column filename
col_names<-colnames(thirsty)
col_names[18]<-"filename"

### populate filename with filename info from wd
data<-setNames(do.call(rbind,Map(`cbind`, 
                           lapply(temp, read.csv), V4=temp)), col_names)

### get rid of laser filter info in column name
colnames(data)<-gsub("\\..*", "", colnames(data))

### get range of values of each column (need to remove string column filename, otherwise matrix will be non-numeric)
data2<-as.matrix(data[,1:17])
ranges<-apply(data2, 2, range)

### values seem to be arcsinh transformed

          





