library(tidyr)
library(dplyr)
library(ggplot2)

omiq_data <- read.csv("~/postdoc/stanford/cytometry/spectral/MUSICAL/big_experiment/ungated/ungated MUSICAL-ungated workflow-Export Statistics copy-250304_2323.csv")


omiq_data%>%
  filter(control==TRUE)%>%
  pivot_longer(cols=starts_with("Q"), names_to="gate", values_to = "freq")%>%
  mutate("normed"=grepl("*Norm", gate))%>%
  mutate(gate_kind=case_when(grepl("*B515", gate)~"CD14_CD7",
                             grepl("*CD7", gate)~"CD14_CD7",
                             grepl("*506", gate)~"CD20_IgD",
                             grepl("*CD20", gate)~"CD20_IgD",
                             grepl("*V450", gate)~"CD16_56",
                             grepl("*CD16", gate)~"CD16_56",
                             grepl("*BV750", gate)~"CCR7_CD45RA",
                             grepl("*CCR7", gate)~"CCR7_CD45RA")
         )%>%
  mutate(quadrant=case_when(grepl("*Q1", gate)~"Q1",
                            grepl("*Q2", gate)~"Q2",
                            grepl("*Q3", gate)~"Q3",
                            grepl("*Q4", gate)~"Q4"))%>%
  ggplot(., aes(x=normed, y=freq, fill=experiment))+
    scale_y_continuous(labels = scales::label_percent())+
    geom_boxplot(outliers=F)+
    facet_wrap(~gate_kind+quadrant, scales="free")+
    theme_minimal()
  
  