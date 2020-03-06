library(ggplot2)
library(tidyr)
library(dplyr)
library(gtools)
library(cowplot)

`%!in%` = Negate(`%in%`)

setwd("/Users/s1249052/PhD/clinical_data/vac69b")
my_paired_palette <- c("#FB9A99","#E31A1C","#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C")

##########   biochem findings   #################

data <- read.csv("biochem.csv", header=T)

data_no_ae <- select(data, -c(colnames(data)[grep("_ae",colnames(data) ,fixed=T)]))

long_data <- gather(data_no_ae, biochem, value, colnames(data_no_ae)[6:13])

long_data <- filter(long_data, long_data$flo_timepoint!="extra")

alt <- subset(long_data, long_data$flo_timepoint=="T6"&long_data$biochem=="alt")


biochem_plot <- ggplot(long_data, aes(x=factor(flo_timepoint, levels=c("Baseline", "C7", "C14", "DoD", "T6", "C96")), y=value, group=trial_number))+
  geom_point(aes(color=factor(trial_number)))+
  geom_line(aes(color=factor(trial_number)))+
  facet_wrap(~biochem, scales="free")+
  theme_minimal()+
  theme(axis.title = element_blank(),
        legend.title = element_blank())

ggsave("biochem.png", biochem_plot, height = 6, width=10)


data <- read.csv("haem.csv", header=T)


data_no_ae <- select(data, -c(colnames(data)[grep("_ae",colnames(data) ,fixed=T)]))

long_data <- gather(data_no_ae, haem, value, colnames(data_no_ae)[7:13])

long_data <- filter(long_data, long_data$flo_timepoint!="extra")

haem_plot <- ggplot(long_data, aes(x=factor(flo_timepoint, levels=c("Baseline", "C7", "C14", "DoD", "T6", "C96")), y=value, group=trial_number))+
  geom_point(aes(color=factor(trial_number)))+
  geom_line(aes(color=factor(trial_number)))+
  facet_wrap(~haem, scales="free")+
  theme_minimal()+
  theme(axis.title = element_blank(),
        legend.title = element_blank())

ggsave("haem.png", haem_plot, height = 6, width=10)




data <- read.csv("/Users/s1249052/PhD/clinical_data/alt_vivax.csv")
first_data <- subset(data, data$N_infection=="1")

fit_Null <- lme4::lmer(alt ~ 0 + (1|Volunteer), data=data) #124.3561
fit_T <- lme4::lmer(alt ~ T_cell_activation+(1|Volunteer), data=data)#103.7344
fit_T_V <- lme4::lmer(alt ~ T_cell_activation+Volunteer+(1|Volunteer), data=data) #45.26805
fit_T_N <- lme4::lmer(alt ~ T_cell_activation+N_infection+(1|Volunteer), data=data) #96.52253
fit_T_V_N <- lme4::lmer(alt ~ T_cell_activation+Volunteer+N_infection+(1|Volunteer), data=data) #34.44452


modelz <- c(fit_Null, fit_T, fit_T_V, fit_T_N)
aic <- lapply(modelz, function(x){AIC(x)})


ggplot(data, aes(x=T_cell_activation, y=alt))+
  geom_point(aes(color=factor(Volunteer), shape=factor(N_infection)))+
  geom_smooth(method = "lm")+
  theme_minimal()+
  #theme(axis.title = element_blank(),
   #     legend.title = element_blank())+
  labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                     "Intercept =",signif(fit$coef[[1]],5 ),
                     " Slope =",signif(fit$coef[[2]], 5),
                     " P =",signif(summary(fit)$coef[2,4], 5)))
