library(tidyr)
library(dplyr)
library(ggplot2)
library(viridis)

#read in data
data <- read.delim("/Users/s1249052/Downloads/perproteinGroups.txt")

#get rid of mouse (if you want to work on mouse proteins just remove the exclamation mark)
plasmodium_only <- subset(data, !grepl("Mus musculus (Mouse)", data$Fasta.headers, fixed=T))

#get rid of first row
plasmodium_only <- plasmodium_only[2:nrow(plasmodium_only),]

#extract plasmo_db codes
plasmodium_only$plasmo_db <- stringr::str_match(plasmodium_only$Fasta.headers, "PBANKA_[0-9]*")[, 1]

#some NAs are produced, not sure why... get rid of them:
plasmodium_only_clean <- plasmodium_only[!is.na(plasmodium_only$plasmo_db),]

#get rid of useless columns
slim_plasmodium <- plasmodium_only_clean[,c(1:12, 33)]

#select first three "morning" samples
morning <- select(slim_plasmodium, c(colnames(slim_plasmodium)[1:3], plasmo_db))

#convert factors to numbers
morning$LFQ.intensity._1 <- as.numeric(as.character(morning$LFQ.intensity._1))
morning$LFQ.intensity._2 <- as.numeric(as.character(morning$LFQ.intensity._2))
morning$LFQ.intensity._3 <- as.numeric(as.character(morning$LFQ.intensity._3))


#remove proteins that are absent
clean_morning <- subset(morning, ifelse(morning$LFQ.intensity._1+morning$LFQ.intensity._2+morning$LFQ.intensity._3==0, FALSE, TRUE))

#make fold change columns
clean_morning$early <- clean_morning$LFQ.intensity._2/clean_morning$LFQ.intensity._1
clean_morning$late <- clean_morning$LFQ.intensity._3/clean_morning$LFQ.intensity._1


#convert to long format needed for ggplotting, get rid of NAs in fold change ()
long_morning <- gather(clean_morning, timepoint, fold_change, c(early, late))

long_morning <- long_morning[!is.na(long_morning$fold_change),]


heatmap <- ggplot(long_morning, aes(x=timepoint, y=plasmo_db, label=plasmo_db, fc=fold_change))+
  geom_tile(aes(fill=fold_change))+
  #scale_fill_viridis(option="A")+
  #scale_fill_brewer(palette = "RdYlGn")
  theme_minimal()+
  labs(fill="Fold Change\nrelative to baseline")+
  theme(legend.title = element_text(hjust=0.5))


ggplotly(heatmap, tooltip = c("label", "fc"))


