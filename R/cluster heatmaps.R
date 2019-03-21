library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(ggthemes)
library(plyr)

setwd("C:/users/Florian/Downloads")

data<-read.csv("results.csv", header=T)

#split data frame into dataframes for different time points


data1<-data
data8<-data
data9<-data
data10<-data

data1$Timepoints<-ifelse(data1$Timepoints=="C-1", "Baseline", NA)
data1<-na.omit(data1)

data8$Timepoints<-ifelse(data8$Timepoints=="C+8", "C+8", NA)
data8<-na.omit(data8)

data9$Timepoints<-ifelse(data9$Timepoints=="C+9", "C+9", NA)
data9<-na.omit(data9)

data10$Timepoints<-ifelse(data10$Timepoints=="C+10", "C+10", NA)
data10<-na.omit(data10)


######## add additional row of median values

data1<-rbind(data1, apply(data1[2:12,], 2 ,median))
data8<-rbind(data8, apply(data8[2:12,], 2 ,median))
data9<-rbind(data9, apply(data9[2:12,], 2 ,median))

# delete gunk from column names

colnames(data1)<-substring(colnames(data), 24, 90)
colnames(data8)<-substring(colnames(data8), 24, 90)
colnames(data9)<-substring(colnames(data9), 24, 90)

# make heatmap of median values

cd25<-data8[13,seq(from=1, to=560, by=7)]
cd25<-as.matrix(cd25)

hladr<-t(data8[13,seq(from=2, to=560, by=7)])
cd4<-t(data8[13,seq(from=3, to=560, by=7)])
cd127<-t(data8[13,seq(from=4, to=560, by=7)])
cd45ra<-t(data8[13,seq(from=5, to=560, by=7)])
cd38<-t(data8[13,seq(from=6, to=560, by=7)])
ccr7<-t(data8[13,seq(from=7, to=560, by=7)])



hm8<-data.frame(t(cd25), hladr, cd4, cd127, cd45ra, cd38, ccr7)
colnames(hm8)<-c("CD25", "HLADR", "CD4", "CD127", "CD45RA", "CD38", "CCR7")
rownames(hm8)<-substring(rownames(hm8), 17, 90)

View(hm8) ### hm8 is a data frame


#convert character numbers into real numbers

hm8$CD25<-as.numeric(as.character(hm8$CD25))
hm8$HLADR<-as.numeric(as.character(hm8$HLADR))
hm8$CD4<-as.numeric(as.character(hm8$CD4))
hm8$CD127<-as.numeric(as.character(hm8$CD127))
hm8$CD45RA<-as.numeric(as.character(hm8$CD45RA))
hm8$CD38<-as.numeric(as.character(hm8$CD38))
hm8$CCR7<-as.numeric(as.character(hm8$CCR7))


##############################
#cutoff cd4: 1111
#cutoff hladr: 1000, 10000
#cd38:3000
#cd45ra: 2200
#ccr7: 600
#cd127: 200
#cd25 500, 3000


####### categorise marker expression aspos neg & high

hm8$CD25<-ifelse(hm8$CD25 < 500, "negative", ifelse(hm8$CD25 > 3000, "high", "positive"))
hm8$HLADR<-ifelse(hm8$HLADR < 1000, "negative", ifelse(hm8$HLADR > 10000, "high", "positive"))
hm8$CD4<-ifelse(hm8$CD4 < 3000, "negative", "positive")
hm8$CD127<-ifelse(hm8$CD127 < 200, "negative", ifelse(hm8$CD127 > 10000, "high", "positive"))
hm8$CD45RA<-ifelse(hm8$CD45RA < 3200, "negative", "positive")
hm8$CD38<-ifelse(hm8$CD38 < 3000, "negative", "positive")
hm8$CCR7<-ifelse(hm8$CCR7 < 650, "negative", "positive")



hm8$ClusterID<-rownames(hm8)
map<-melt(hm8, id="ClusterID")
colnames(map)<-c("ClusterID", "Marker", "Median Expression")


yes<-map[order(map$ClusterID ),]

##### reorder markers for figure legibility
yes$Marker<-factor(yes$Marker, levels=c("CD4", "CD45RA", "CCR7", "CD25", "CD127", "CD38", "HLADR" ))

####### good heatmap

windows();ggplot(data=yes, aes(x=yes$Marker, y=yes$ClusterID))+
  geom_tile(aes(fill=yes$`Median Expression`), color="white")+
  scale_fill_manual(values= c("red", "blue", "yellow"))

######### heatmap only of the top 20 clusters   #######

yes$top<-ifelse(yes$ClusterID %in% top20$variable, "top20", NA)

maybs<-na.omit(yes)


windows();ggplot(data=maybs, aes(x=maybs$Marker, y=maybs$ClusterID=="Cluster_55"))+
  geom_tile(aes(fill=maybs$`Median Expression`), color="white")+
  scale_fill_manual(values= c("red", "blue", "yellow"))





top20["CLUSTER_55",]





top20$variable=="CLUSTER55"





c55<-subset(top20, top20$variable=="CLUSTER_55")






















################     garbage code nobody needs       ###########



#determine maxima for each colum

apply(hm8, 2, max)
#######    CD25     HLADR       CD4     CD127    CD45RA      CD38      CCR7 
#       8442.840  5146.382  6570.381  4136.627 20578.840  5196.646  1984.875 




#make new data frame with all fluorescent intensities divided by the maximum

relative<-data.frame(hm8$CD25/8442.84, hm8$HLADR/5146, hm8$CD4/6570, hm8$CD127/4136, hm8$CD45RA/20578, hm8$CD38/5196, hm8$CCR7/1984, rownames(hm8))
colnames(relative)<-c("CD25", "HLADR", "CD4", "CD127", "CD45RA", "CD38", "CCR7", "ClusterID")
relative$ClusterID<-substring(relative$ClusterID, 17, 90)
apply(relative, 2, max)
#remove negative/infinetisimally small values

relative$CD25<-ifelse(relative$CD25<0.00001, 0.00001, relative$CD25)
relative$HLADR<-ifelse(relative$HLADR<0.00001, 0.00001, relative$HLADR)
relative$CD4<-ifelse(relative$CD4<0.00001, 0.00001, relative$CD4)
relative$CD127<-ifelse(relative$CD127<0.00001, 0.00001, relative$CD127)
relative$CD45RA<-ifelse(relative$CD45RA<0.00001, 0.00001, relative$CD45RA)
relative$CD38<-ifelse(relative$CD38<0.00001, 0.00001, relative$CD38)
relative$CCR7<-ifelse(relative$CCR7<0.00001, 0.00001, relative$CCR7)

goodmap<-melt(relative, id="ClusterID")

goodmap<-goodmap[order(goodmap$Marker),]

colnames(goodmap)<-c("ClusterID", "Marker", "Relative Expression")

goodmap$`Relative Expression`<-log(goodmap$`Relative Expression`)

goodmap$`Relative Expression`<-as.numeric(as.character(goodmap$`Relative Expression`))


windows();ggplot(data=goodmap, aes(x=goodmap$Marker, y=goodmap$ClusterID, fill=goodmap$`Relative Expression`))+
  geom_tile()+
  scale_fill_gradient(low="beige", high="red", limits=c(-2, 0), na.value="blue")

apply(relative, 2, min)






