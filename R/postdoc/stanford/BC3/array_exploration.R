library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)


igg_data <- read.csv("~/postdoc/stanford/clinical_data/BC3/igg)NormalizedData.csv")
# igg3_data <- read.csv("~/postdoc/stanford/clinical_data/BC3/igg3_NormalizedData.csv")

demographic_data <- readxl::read_excel("~/postdoc/stanford/clinical_data/BC3/adi_demographics.xlsx")
#add an s so that the id values match the antibody dataset
demographic_data$id <- paste("s", demographic_data$id, sep="")

long_data <- igg_data %>%
  rename("antigen"=ID)%>%
  pivot_longer(cols = colnames(igg_data)[5:104], names_to = "id", values_to = "concentration")

data_summary <- long_data %>%
  group_by(antigen)%>%
  summarise("min"=min(concentration),
            "max"=max(concentration),
            "mean"=mean(concentration),
            "median"=median(concentration),
            "sd"=sd(concentration))


# scale everything to unit variance
long_data$scaled_concentration <- (long_data$concentration-data_summary$mean[match(long_data$antigen,data_summary$antigen)])/
  data_summary$sd[match(long_data$antigen,data_summary$antigen)]

# long_scaled_data_summary <- long_data %>%
#   group_by(antigen)%>%
#   summarise("min"=min(scaled_concentration),
#             "max"=max(scaled_concentration),
#             "median"=median(scaled_concentration),
#             "mean"=mean(scaled_concentration),
#             "sd"=sd(scaled_concentration))


ggplot(long_scaled_data_summary, aes(x=median, y=sd))+
  geom_point(alpha=0.4)+
  geom_smooth(method="lm")+
  theme_minimal()+
  theme(legend.position="none")


wide_data <- long_data %>%
  pivot_wider(names_from = antigen, values_from = concentration, id_cols = id)

big_pca <-  prcomp(wide_data[,-1], center = T)

# pca_plot_data <- as.data.frame(big_pca$x)
pca_plot_data <- as.data.frame(cbind(wide_data, big_pca$x))


pca_plot_data <- cbind(pca_plot_data, demographic_data[match(wide_data$id, demographic_data$id),-1])

pca_theme <- theme(legend.title = element_text(size=15),
                   legend.text = element_text(size=13),
                   axis.title = element_text(size=15),
                   axis.text = element_text(size=12))

ggplot(pca_plot_data, aes(x=PC1, y=PC2, fill=factor(motherid)))+
  geom_point(shape=21, size=2)+
  xlab(paste("PC1 ", data.frame(summary(big_pca)[6])[2,1]*100, "%", sep = ""))+
  ylab(paste("PC2 ", data.frame(summary(big_pca)[6])[2,2]*100, "%", sep = ""))+
  # geom_density_2d(aes(color=Txarm))+
  theme_minimal()+
  pca_theme+
  theme(legend.position = "none")



long_combo_data <- cbind(long_data, demographic_data[match(long_data$id, demographic_data$id),-1])

purrrf <- long_combo_data %>%
  group_by(antigen) %>%
  nest() %>%
  mutate(model=map(data, ~lm(concentration~Txarm+gender, data=.))) %>%
  # mutate(summary=map(model, ~summary(.))) %>%
  mutate(raw_p=map_dbl(model, ~summary(.)$coefficients[11]))%>%
  ungroup()%>%
  mutate(padj=p.adjust(raw_p))

sigs <- filter(purrrf, raw_p<0.01)


long_combo_data %>%
  filter(antigen %in% sigs$antigen)%>%
  ggplot(., aes(x=gender, y=concentration))+
  geom_boxplot(aes(fill=gender))+
  # geom_violin(aes(fill=Txarm),  draw_quantiles = c(0.25, 0.5, 0.75))+
  facet_wrap(~antigen)+
  theme_minimal()
  

pca_plot_data
ggplot(pca_plot_data, aes(x=PC1, y=PC2, fill=PF3D7_0530100.2o2))+
  geom_point(shape=21, size=2)+
  xlab(paste("PC1 ", data.frame(summary(big_pca)[6])[2,1]*100, "%", sep = ""))+
  ylab(paste("PC2 ", data.frame(summary(big_pca)[6])[2,2]*100, "%", sep = ""))+
  # geom_density_2d(aes(color=Txarm))+
  theme_minimal()+
  pca_theme+
  theme(legend.position = "none")+
  viridis::scale_fill_viridis(option = "B")



# split PCA by treatment arm
wide_combo_data <- cbind(wide_data, demographic_data[match(wide_data$id, demographic_data$id),-1])

dp_data <- filter(wide_combo_data, Txarm=="DP")
sp_data <- filter(wide_combo_data, Txarm=="SP")

dp_pca <-  prcomp(dp_data[,2:991], center = T)
sp_pca <-  prcomp(sp_data[,2:991], center = T)

# pca_plot_data <- as.data.frame(big_pca$x)
dp_pca_plot_data <- as.data.frame(cbind(dp_data, dp_pca$x))
sp_pca_plot_data <- as.data.frame(cbind(sp_data, sp_pca$x))

# 
# dp_pca_plot_data <- cbind(dp_pca_plot_data, demographic_data[match(dp_data$id, demographic_data$id),-1])
# sp_pca_plot_data <- cbind(sp_pca_plot_data, demographic_data[match(sp_data$id, demographic_data$id),-1])


ggplot(dp_pca_plot_data, aes(x=PC1, y=PC2, fill=PF3D7_0530100.2o2))+
  geom_point(shape=21, size=2)+
  xlab(paste("PC1 ", data.frame(summary(dp_pca)[6])[2,1]*100, "%", sep = ""))+
  ylab(paste("PC2 ", data.frame(summary(dp_pca)[6])[2,2]*100, "%", sep = ""))+
  # geom_density_2d(aes(color=Txarm))+
  theme_minimal()+
  pca_theme+
  theme(legend.position = "none")


ggplot(sp_pca_plot_data, aes(x=PC1, y=PC2, fill=PF3D7_0530100.2o2))+
  geom_point(shape=21, size=2)+
  xlab(paste("PC1 ", data.frame(summary(sp_pca)[6])[2,1]*100, "%", sep = ""))+
  ylab(paste("PC2 ", data.frame(summary(sp_pca)[6])[2,2]*100, "%", sep = ""))+
  # geom_density_2d(aes(color=Txarm))+
  theme_minimal()+
  pca_theme+
  theme(legend.position = "none")


# look for patterns in other covariates
num_cols <- colnames(demographic_data)[lapply(demographic_data, typeof)=="double"]
num_cols <- num_cols[-c(1,2, 15)]
char_cols <- colnames(demographic_data)[lapply(demographic_data, typeof)=="character"]
char_cols <- char_cols[-c(1,3)]

long_demo_data <- demographic_data %>%
  pivot_longer(cols=num_cols, names_to = "num_covariate", values_to = "num_value")%>%
  pivot_longer(cols=char_cols, names_to = "char_covariate", values_to = "char_value")


ggplot(long_demo_data, aes(x=Txarm, num_value))+
  geom_point()+
  geom_violin(aes(fill=Txarm))+
  facet_wrap(~num_covariate, scales = "free")+
  theme_minimal()

