#download RNAseqPower, mark both lines of code and run simultaneously
source("https://bioconductor.org/biocLite.R")
biocLite("RNASeqPower")


library("RNASeqPower")
library("tidyr")
libraray("ggplot2")




# depth = coverage
# cv = coefficient of variance ~ reproducibility between replicates
# effect = effect size i.e. fold change you reliably want to detect (1.5 = 50% increase, 2= doubling etc.)
# alpha = p value (two sided test for difference)
# power = power

results<-rnapower(depth=70, cv=.2, effect=c(1.5, 1.75, 2, 4, 6),
                  alpha= .05, power=c(0.7,0.8, 0.9, 0.95, 0.99))

results<-data.frame(results)
results$effect<-rownames(results)

# the 1:5 term needs to be changed, according to how many versions of power you're testing (right now it's 5)..
# if you're testing two it should be 1:2 and so on

results2<-gather(results, power, sample_size, 1:5)

ggplot(results2, aes(x=power, y=sample_size, group=factor(effect)))+
  geom_point(aes(color=factor(effect)))+
  #ggtitle("70x coverage, variance 0.2")+
  geom_line(aes(color=factor(effect)))+
  scale_y_continuous(breaks=seq(0,20, by=2))+
  geom_hline(aes(yintercept=5), colour = "red", linetype = "dashed", size=1)+
  geom_hline(aes(yintercept=4), colour = "green", linetype = "dashed", size=1)+
 get

  