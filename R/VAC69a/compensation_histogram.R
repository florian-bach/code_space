library(ggplot2)

comp_mat <- read.csv("/Users/s1249052/PhD//flow_data/vac69a/FlorianBACHspillMat.csv", header = T, row.names=1)
dim(comp_mat)
comp_mat <- dplyr::select(comp_mat, -Pd102Di, -Pd104Di, -Pd105Di, -Pd106Di, -Pd108Di, -Pd110Di, -Cd111Di, -Cd112Di, -Cd113Di, -Cd116Di)
dim(comp_mat)

lin_comp_mat <- unlist(comp_mat)
names(lin_comp_mat) <- NULL

clean_comp_mat <- subset(lin_comp_mat, lin_comp_mat !=1)
#clean_comp_mat <- subset(clean_comp_mat, clean_comp_mat !=0)

ggplot()+
  geom_histogram(aes(x=clean_comp_mat*100), fill="maroon", binwidth = 0.25)+
  theme_minimal()+
  scale_y_log10()+
  stat_bin(binwidth= 0.25, geom="text", aes(x=clean_comp_mat*100, label=..count..), vjust=-1)+
  xlab("Compensation in %")+
  ylab("Count")+
  scale_x_continuous(breaks = seq(1,7))
