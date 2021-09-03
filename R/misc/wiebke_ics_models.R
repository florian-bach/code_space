library(tidyr)
library(dplyr)
library(ggplot2)

cytokines <- read.csv("~/Documents/wiebke_cytokines.csv")

cytokines$timepoint <- factor(cytokines$timepoint, levels=c("baseline", "diagnosis", "C28", "C90"))
cytokines$volunteer <- as.character(cytokines$volunteer)

model <- glm(IFN~timepoint+volunteer, data=cytokines, family = "binomial", weights = )

summary(model)
# 
# dod_contrast <- t(matrix(c(0,1,0,0,0,0,0)))
# c28_contrast <- t(matrix(c(0,0,1,0,0,0,0)))
# c90_contrast <- t(matrix(c(0,0,0,1,0,0,0)))
# 
# linear_test <- multcomp::glht(model, dod_contrast)
# summary(linear_test)
# 
# linear_test <- multcomp::glht(model, c28_contrast)
# summary(linear_test)
# 
# linear_test <- multcomp::glht(model, c90_contrast)
# summary(linear_test)


ggplot(cytokines, aes(x=timepoint, y=IFN, group=volunteer))+
  geom_point(aes(colour=volunteer))+
  geom_line(aes(colour=volunteer))+
  ylab("% IFN producing Cells")+
  theme_minimal()
