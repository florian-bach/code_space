library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(ggthemes)
library(cowplot)
library(plyr)

setwd("C:/users/Florian/Documents")

cd4<-read.csv("cd4.csv", header=T)
cd4neg<-read.csv("cd4neg.csv", header=T)

cd4$Timepoints<-factor(cd4$Timepoints, levels=c("C-1", "C+8", "C+9", "C+10"))
cd4neg$Timepoints<-factor(cd4neg$Timepoints, levels=c("C-1", "C+8", "C+9", "C+10"))


########## relative  ############

cd4plot<-ggplot(cd4, aes(x=Timepoints, y=georgio..CD4, group=factor(Individuals)))+
  geom_point(aes(color=factor(Individuals)))+
  geom_line(aes(color=factor(Individuals)))+
  theme_bw()+
  theme(axis.line = element_line(colour = "black"))+theme(axis.title.x=element_blank())+
  ylab("CD4+ CD38+ CCR7- Effectors")+
  theme(legend.position="none")+
  scale_y_continuous(limits=c(0, 15), labels = function(x) paste0(x, "%"))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))+
  scale_color_cividis(alpha = 1, begin = 0, end = 1, direction = 1,
                      discrete = TRUE, option = "V")

cd4negplot<-ggplot(cd4neg, aes(x=Timepoints, y=cd4neg., group=factor(Individuals)))+
  geom_point(aes(color=factor(Individuals)))+
  geom_line(aes(color=factor(Individuals)))+
  theme_bw()+
  theme(axis.line = element_line(colour = "black"))+theme(axis.title.x=element_blank())+
  ylab("CD4- CD38+ CCR7- Effectors")+
  theme(legend.position="none")+
  scale_y_continuous(limits=c(0, 12), labels = function(x) paste0(x, "%"))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

cd4hladrplot<-ggplot(cd4, aes(x=Timepoints, y=Medians.of.HLADR, group=factor(Individuals)))+
  geom_point(aes(color=factor(Individuals)))+
  geom_line(aes(color=factor(Individuals)))+
  theme_bw()+
  theme(axis.line = element_line(colour = "black"))+
  theme(axis.title.x=element_blank())+
  ylab("HLADR median on CD4+ CD38+ CCR7- Effectors")+
  scale_y_continuous(limits=c(0, 300))+
    theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))+
  stat_summary(aes(y = Medians.of.HLADR ,group=1), fun.y=median, colour="black", geom="line",group=1)


cd4neghladrplot<-ggplot(cd4neg, aes(x=Timepoints, y=Medians.of.HLADR, group=factor(Individuals)))+
  geom_point(aes(color=factor(Individuals)))+
  geom_line(aes(color=factor(Individuals)))+
  theme_bw()+
  theme(axis.line = element_line(colour = "black"))+
  theme(axis.title.x=element_blank())+
  ylab("HLADR median on CD4- CD38+ CCR7- Effectors")+
  scale_y_continuous(limits=c(0, 300))+
    theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))+
  stat_summary(aes(y = Medians.of.HLADR ,group=1), fun.y=median, colour="black", geom="line",group=1)


########### group ##########

cd4groupplot<-ggplot(cd4, aes(x=Timepoints, y=georgio..CD4, group=factor(Individuals)))+
  geom_point(aes(color=factor(Conditions)))+
  geom_line(aes(color=factor(Conditions)))+
  theme_bw()+
  theme(axis.line = element_line(colour = "black"))+theme(axis.title.x=element_blank())+
  ylab("CD4+ CD38+ CCR7- Effectors")+
  theme(legend.position="none")+
  scale_y_continuous(limits=c(0, 15), labels = function(x) paste0(x, "%"))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

cd4groupnegplot<-ggplot(cd4neg, aes(x=Timepoints, y=cd4neg., group=factor(Individuals)))+
  geom_point(aes(color=factor(Infection)))+
  geom_line(aes(color=factor(Infection)))+
  theme_bw()+
  theme(axis.line = element_line(colour = "black"))+theme(axis.title.x=element_blank())+
  ylab("CD4- CD38+ CCR7- Effectors")+
  theme(legend.position="none")+
  scale_y_continuous(limits=c(0, 12), labels = function(x) paste0(x, "%"))+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))


########### cray.dat ###########


setwd("C:/users/Florian/PhD/flowdata/re_vac063_T_cell_sorts")

data<- read.csv("better cluster percentages.csv", header=T )

primary<-c(806, 808, 818, 820, 825)
secondary<-c(802, 812, 815, 819, 821, 822, 823, 824)

data$group<- ifelse(data$Individuals %in% primary, "Primary", ifelse(data$Individuals %in% secondary, "Secondary", NA))

data<-data[order(data$Individuals, data$Timepoints),]

data<-data.frame(data)

data$FCS.Filename<-NULL

by_person<-split(data, data$Individuals)
 


data802<-data.frame(by_person[1])
data808<-data.frame(by_person[2])
data812<-data.frame(by_person[3])
data815<-data.frame(by_person[4])
data818<-data.frame(by_person[5])
data819<-data.frame(by_person[6])
data820<-data.frame(by_person[7])
data821<-data.frame(by_person[8])#4
data822<-data.frame(by_person[9])
data823<-data.frame(by_person[10])#4
data824<-data.frame(by_person[11])#4
data825<-data.frame(by_person[12])








change802<- data802[3,3:72] / data802[1,3:72] 
change802<-melt(change802)
change802<-na.omit(change802)
change802<-change802[order(change802$value, decreasing=T),]
change802$infection<-"secondary"

change808<- data808[3,3:72] / data808[1,3:72] 
change808<-melt(change808)
change808<-na.omit(change808)
change808<-change808[order(change808$value, decreasing=T),]
change808$infection<-"primary"

change812<- data812[3,3:72] / data812[1,3:72] 
change812<-melt(change812)
change812<-na.omit(change812)
change812<-change812[order(change812$value, decreasing=T),]
change812$infection<-"secondary"


change815<- data815[3,3:72] / data815[1,3:72] 
change815<-melt(change815)
change815<-na.omit(change815)
change815<-change815[order(change815$value, decreasing=T),]
change815$infection<-"secondary"

change818<- data818[3,3:72] / data818[1,3:72] 
change818<-melt(change818)
change818<-na.omit(change818)
change818<-change818[order(change818$value, decreasing=T),]
change818$infection<-"primary"

change819<- data819[3,3:72] / data819[1,3:72] 
change819<-melt(change819)
change819<-na.omit(change819)
change819<-change819[order(change819$value, decreasing=T),]
change819<-change819[order(change819$value, decreasing=T),]
change819$infection<-"secondary"


change820<- data820[3,3:72] / data820[1,3:72] 
change820<-melt(change820)
change820<-na.omit(change820)
change820<-change820[order(change820$value, decreasing=T),]
change820<-change820[order(change820$value, decreasing=T),]
change820$infection<-"primary"

change821<- data821[2,3:72] / data821[1,3:72] 
change821<-melt(change821)
change821<-na.omit(change821)
change821<-change821[order(change821$value, decreasing=T),]
change821<-change821[order(change821$value, decreasing=T),]
change821$infection<-"secondary"

change822<- data822[3,3:72] / data822[1,3:72] 
change822<-melt(change822)
change822<-na.omit(change822)
change822<-change822[order(change822$value, decreasing=T),]
change822<-change822[order(change822$value, decreasing=T),]
change822$infection<-"secondary"

change823<- data823[2,3:72] / data823[1,3:72] 
change823<-melt(change823)
change823<-na.omit(change823)
change823<-change823[order(change823$value, decreasing=T),]
change823<-change823[order(change823$value, decreasing=T),]
change823$infection<-"secondary"

change824<- data824[2,3:72] / data824[1,3:72] 
change824<-melt(change824)
change824<-na.omit(change824)
change824<-change824[order(change824$value, decreasing=T),]
change824<-change824[order(change824$value, decreasing=T),]
change824$infection<-"secondary"

change825<- data825[3,3:72] / data825[1,3:72] 
change825<-melt(change825)
change825<-na.omit(change825)
change825<-change825[order(change825$value, decreasing=T),]
change825<-change825[order(change825$value, decreasing=T),]
change825$infection<-"primary"



top20<-rbind(change802[1:20,1:3], change808[1:20,1:3], change812[1:20,1:3],change815[1:20,1:3],change818[1:20,1:3],change819[1:20,1:3],change820[1:20,1:3],change821[1:20,1:3],change822[1:20,1:3],change823[1:20,1:3],change824[1:20,1:3],change825[1:20,1:3])








bot20<-rbind(tail(change802, n=20), tail(change808, n=20), tail(change812, n=20), tail(change815, n=20), tail(change818, n=20), tail(change819, n=20), tail(change820, n=20), tail(change821, n=20), tail(change822, n=20), tail(change823, n=20), tail(change824, n=20), tail(change825, n=20))

topprim<-rbind(change808[1:20,1:3],change818[1:20,1:3],change820[1:20,1:3],change825[1:20,1:3])

topsec<-rbind(change802[1:20,1:3], change812[1:20,1:3], change815[1:20,1:3], change819[1:20,1:3], change821[1:20,1:3],change822[1:20,1:3],change823[1:20,1:3],change824[1:20,1:3])


top20$ratio<-ifelse(topprim$variable %in% topsec$variable, table(topprim$variable)[,]/table(topsec$variable)[,], "unique")



top20$individual<-substring(top20$variable, 2,4)



top20$variable<-substring(top20$variable, 12, 21)
topprim$variable<-substring(topprim$variable, 12, 90)
topsec$variable<-substring(topsec$variable, 12, 90)
bot20$variable<-substring(bot20$variable, 12, 90)

top20$value<-ifelse(top20$value==Inf, NA, top20$value)
top20<-na.omit(top20)

topprim$value<-ifelse(topprim$value==Inf, NA, topprim$value)
topprim<-na.omit(topprim)

topsec$value<-ifelse(topsec$value==Inf, NA, topsec$value)
topsec<-na.omit(topsec)


clusters<-data.frame(top20$variable, top20$infection, top20$individual)


View(top20)

nice<-table(top20$variable)
nice<-nice[order()]


top20factor<-(as.factor(top20$variable))

colors <- colorRampPalette(brewer.pal(8,"Dark2"))


ggplot(top20)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  geom_bar(aes(x=variable, color=infection, fill=infection), position=position_dodge())

prim<-subset(top20, top20$infection=="primary")
sec<-subset(top20, top20$infection=="secondary")


selecsec<-as.matrix(table(sec$variable))
selecsec<-subset(selecsec, selecsec[,1]>=5)

selecprim<-as.matrix(table(prim$variable))
selecprim[rownames(selecsec),]

ggplot(bot20)+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+geom_bar(aes(x=variable, color=infection, fill=infection), position=position_dodge())

head(top20)

newtop<-data.frame(table(clusters$top20.variable))

newtop$Freq<-ifelse(newtop$Freq==1, NA, newtop$Freq)
newtop<-na.omit(newtop)



prim<-as.data.frame(table(topprim$variable))
prim$infection<-"primary"

zwei<-as.data.frame(table(topsec$variable))
zwei$infection<-"secondary"

drei<-rbind(prim, zwei)
drei<-drei[order(drei$Var1),]

alone<-tail(drei, n=23)


##########################      HEATMAP CODE        #########################







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

windows();ggplot(data=yes, aes(x=yes$Marker, y=yes$ClusterID=="CLUSTER_3"))+
  geom_tile(aes(fill=yes$`Median Expression`), color="white")+
  scale_fill_manual(values= c("red", "blue", "yellow"))


rownames(alone$Var1)<-substring(rownames(alone$var1), 1,10)

alone$var2<-substring(alone[,1], 0, 10) 

alone$var2<-lapply(alone$var2, gsub, pattern = ".", replacement = "")

hot<-subset(yes, yes$ClusterID %in% alone$var2)
rbind(hot, yes$ClusterID=="CLUSTER_1")

windows();ggplot(data=hot, aes(x=hot$Marker, y=hot$ClusterID))+
  geom_tile(aes(fill=hot$`Median Expression`), color="white")+
  scale_fill_manual(values= c("red", "green"))

######### heatmap only of the top 20 clusters   #######

yes$top<-ifelse(yes$ClusterID %in% top20$variable, "top20", NA)

maybs<-na.omit(yes)


windows();ggplot(data=maybs, aes(x=maybs$Marker, y=maybs$ClusterID))+
  geom_tile(aes(fill=maybs$`Median Expression`), color="white")+
  scale_fill_manual(values= c("red", "blue", "yellow"))




