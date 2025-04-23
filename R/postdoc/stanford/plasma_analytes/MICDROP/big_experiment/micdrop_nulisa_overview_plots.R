library(tidyr)
library(dplyr)
library(ggplot2)
library(purrr)

`%notin%` <- Negate(`%in%`)

clean_data <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/clean_data_with_meta.csv")

clean_data%>%
  distinct(id, timepoint, total_n_malaria_12)%>%
  group_by(total_n_malaria_12)%>%
  summarise(n())

# PCA ####
big_palette_hex <- c("#000000", "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
                     "#FFDBE5", "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
                     "#5A0007", "#809693", "#FEFFE6", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
                     "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92", "#FF90C9", "#B903AA", "#D16100",
                     "#DDEFFF", "#000035", "#7B4F4B", "#A1C299", "#300018", "#0AA6D8", "#013349", "#00846F",
                     "#372101", "#FFB500", "#C2FFED", "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09",
                     "#00489C", "#6F0062", "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66",
                     "#885578", "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
                     "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F", "#938A81",
                     "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757", "#C8A1A1", "#1E6E00",
                     "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C", "#772600", "#D790FF", "#9B9700",
                     "#549E79", "#FFF69F", "#201625", "#72418F", "#BC23FF", "#99ADC0", "#3A2465", "#922329",
                     "#5B4534", "#FDE8DC", "#404E55", "#0089A3", "#CB7E98", "#A4E804", "#324E72", "#6A3A4C")

# big_palette <- c("240,163,255","0,117,220","153,63,0","76,0,92","25,25,25","0,92,49","43,206,72","255,204,153","128,128,128","148,255,181","143,124,0","157,204,0","194,0,136","0,51,128","255,164,5","255,168,187","66,102,0","255,0,16","94,241,242","0,153,143","224,255,102","116,10,255","153,0,0","255,255,128","255,255,0","255,80,5")
# 
# big_palette_hex <- sapply(strsplit(big_palette, ","), function(x)
#   rgb(x[1], x[2], x[3], maxColorValue=255))


id_columns <- c("sample", "qc_failed",  "plate", "id", "gender", "ageinwks", "timepoint", "log_qpcr", "total_n_malaria_12", "mstatus", "febrile", "fever", "rogerson", "SGA")

wide_df2 <- clean_data %>%
  filter(targetName %notin% c("IFNA2", "CTSS", "LTA|LTB"))%>%
  # filter(sample_id %notin% c("X384_S.1_D1V93U", "X744_A.1_D1KT5Z", "X323_A14_D1DZ2D", "X496_NM7_D1CAYS", "X316_S.1_D12FNR", "X667_NM7_D1GBYA", "X 176 S_t14 D1FXRJ", "X164_NM0_D1WA6Q"))%>%
  pivot_wider(names_from = targetName, values_from = conc, id_cols = all_of(id_columns))


rownames(wide_df2) <- wide_df2$sample

# each row = measurement; each column = feature
big_pca <-  prcomp(wide_df2[,(length(id_columns)+1):ncol(wide_df2)], center = T)
pca_plot_data <- as.data.frame(cbind(wide_df2, big_pca$x))

# [1] "sample"             "qc_failed"          "plate"              "id"                 "gender"             "ageinwks"           "timepoint"         
# [8] "log_qpcr"           "total_n_malaria_12" "mstatus"            "febrile"            "fever"


set.seed(123)

large_pca12 <- pca_plot_data %>%
  mutate(pcr_cat=round(log_qpcr))%>%
  # mutate(qc_failed=as.factor(qc_failed))%>%
  pivot_longer(cols=c(id_columns[-c(1, 4, 7, 8)], pcr_cat), names_to = "meta_column", values_to = "meta_value", values_transform = as.character)%>%
  arrange(meta_value)%>%
  ggplot(., aes(x=PC1, y=PC2, color=factor(meta_value)))+
  geom_point()+
  stat_ellipse()+
  facet_wrap(~meta_column)+
  xlab(paste("PC1 ", data.frame(summary(big_pca)[6])[2,1]*100, "%", sep = ""))+
  ylab(paste("PC2 ", data.frame(summary(big_pca)[6])[2,2]*100, "%", sep = ""))+
  theme_minimal()+
  scale_color_manual(values=(sample(big_palette_hex, size = 48)))+
  # viridis::scale_color_viridis(discrete = T)+
  theme(legend.position="none",
        axis.text = element_blank())

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/figures/large_covariate_pc1_pc2.png", large_pca12, height=10, width = 10, bg="white")



medium_pca12 <- pca_plot_data %>%
  mutate(pcr_cat=round(log_qpcr))%>%
  # mutate(qc_failed=as.factor(qc_failed))%>%
  pivot_longer(cols=c("febrile", "mstatus", "gender", "qc_failed"), names_to = "meta_column", values_to = "meta_value", values_transform = as.character)%>%
  arrange(meta_value)%>%
  ggplot(., aes(x=PC1, y=PC2, color=factor(meta_value)))+
  geom_point()+
  stat_ellipse()+
  facet_wrap(~meta_column)+
  xlab(paste("PC1 ", data.frame(summary(big_pca)[6])[2,1]*100, "%", sep = ""))+
  ylab(paste("PC2 ", data.frame(summary(big_pca)[6])[2,2]*100, "%", sep = ""))+
  theme_minimal()+
  scale_color_manual(values=c("#FFB500", "#B903AA", "#8FB0FF", "grey", "darkred"))+
  # viridis::scale_color_viridis(discrete = T)+
  theme(legend.position="none",
        axis.text = element_blank())

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/figures/medium_covariate_pc1_pc2.png", medium_pca12, height=6, width = 6, bg="white")


plate_pca12 <- pca_plot_data %>%
  mutate(pcr_cat=round(log_qpcr))%>%
  # mutate(qc_failed=as.factor(qc_failed))%>%
  ggplot(., aes(x=PC1, y=PC2, color=factor(plate)))+
  geom_point()+
  stat_ellipse()+
  xlab(paste("PC1 ", data.frame(summary(big_pca)[6])[2,1]*100, "%", sep = ""))+
  ylab(paste("PC2 ", data.frame(summary(big_pca)[6])[2,2]*100, "%", sep = ""))+
  theme_minimal()+
  viridis::scale_color_viridis(discrete = T)+
  # viridis::scale_color_viridis(discrete = T)+
  theme(axis.text = element_blank(),
        legend.title = element_blank())

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/figures/plate_pca12.png", plate_pca12, height=4, width = 6, bg="white")


age_pca12 <- pca_plot_data %>%
  mutate(pcr_cat=round(log_qpcr))%>%
  filter(timepoint!="68 weeks")%>%
  # mutate(qc_failed=as.factor(qc_failed))%>%
  ggplot(., aes(x=PC1, y=PC2, color=factor(timepoint, levels=c("8 weeks", "24 weeks", "52 weeks", "68 weeks"))))+
  geom_point()+
  stat_ellipse()+
  xlab(paste("PC1 ", data.frame(summary(big_pca)[6])[2,1]*100, "%", sep = ""))+
  ylab(paste("PC2 ", data.frame(summary(big_pca)[6])[2,2]*100, "%", sep = ""))+
  theme_minimal()+
  viridis::scale_color_viridis(discrete = T, direction = -1)+
  # viridis::scale_color_viridis(discrete = T)+
  theme(axis.text = element_blank(),
        legend.title = element_blank())

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/figures/age_pca12.png", age_pca12, height=4, width = 6, bg="white")



large_pca34 <- pca_plot_data %>%
  mutate(pcr_cat=round(log_qpcr))%>%
  # mutate(qc_failed=as.factor(qc_failed))%>%
  pivot_longer(cols=c(id_columns[-c(1, 4, 7, 8)], pcr_cat), names_to = "meta_column", values_to = "meta_value", values_transform = as.character)%>%
  arrange(meta_value)%>%
  ggplot(., aes(x=PC3, y=PC4, color=factor(meta_value)))+
  geom_point()+
  stat_ellipse()+
  facet_wrap(~meta_column)+
  xlab(paste("PC3 ", data.frame(summary(big_pca)[6])[2,3]*100, "%", sep = ""))+
  ylab(paste("PC4 ", data.frame(summary(big_pca)[6])[2,4]*100, "%", sep = ""))+
  theme_minimal()+
  scale_color_manual(values=(big_palette_hex))+
  # viridis::scale_color_viridis(discrete = T)+
  theme(legend.position="none",
        axis.text = element_blank())

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/figures/large_covariate_pc3_pc4.png", large_pca34, height=10, width = 10, bg="white")


large_pca56 <- pca_plot_data %>%
  mutate(pcr_cat=round(log_qpcr))%>%
  # mutate(qc_failed=as.factor(qc_failed))%>%
  pivot_longer(cols=c(id_columns[-c(1, 4, 7, 8)], pcr_cat), names_to = "meta_column", values_to = "meta_value", values_transform = as.character)%>%
  arrange(meta_value)%>%
  ggplot(., aes(x=PC5, y=PC6, color=factor(meta_value)))+
  geom_point()+
  stat_ellipse()+
  facet_wrap(~meta_column)+
  xlab(paste("PC5 ", data.frame(summary(big_pca)[6])[2,5]*100, "%", sep = ""))+
  ylab(paste("PC6 ", data.frame(summary(big_pca)[6])[2,6]*100, "%", sep = ""))+
  theme_minimal()+
  scale_color_manual(values=(big_palette_hex))+
  # viridis::scale_color_viridis(discrete = T)+
  theme(legend.position="none",
        axis.text = element_blank())

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/figures/large_covariate_pc5_pc6.png", large_pca56, height=10, width = 10, bg="white")
# loadings_df <- data.frame(big_pca$rotation)
# loadings_df$targetName <- rownames(loadings_df)
# loadings_df$targetName <- factor(loadings_df$targetName, levels = loadings_df$targetName[order(loadings_df$PC1)])
# 
# pc2_cols <- colorspace::sequential_hcl(nrow(loadings_df), palette = "Purple Yellow")
# names(pc2_cols) <- loadings_df$targetName[order(loadings_df$PC2)]
# 
# pc2_plot <- ggplot(loadings_df, aes(x=factor(targetName, levels = targetName[order(loadings_df$PC2)]), y=PC2, fill=targetName))+
#   geom_bar(stat = "identity")+
#   scale_fill_manual(values = pc2_cols)+
#   theme_minimal()+
#   theme(axis.title.x = element_blank(),
#         axis.text.x = element_text(angle = 90, hjust=1),
#         legend.position = "none")

# parasitemia plots


clean_data %>%
  filter(mstatus!=1, timepoint=="52 weeks", targetName %in% c("IL10", "IL15", "TLR3"))%>%
  ggplot(., aes(x=log_qpcr, y=conc, color=factor(total_n_malaria_12)))+
  geom_point()+
  theme_minimal()+
  facet_wrap(~targetName, scales="free")+
  scale_color_manual(values=c("black", "darkgrey", "darkred", "red"))


n_para_qpcr_plot <- clean_data %>%
  filter(mstatus!=1, timepoint!="68 weeks", targetName %in% c("IL10", "IL15", "TLR3"))%>%
  mutate(timepoint=factor(timepoint, levels=c("8 weeks", "24 weeks", "52 weeks")))%>%
  arrange(total_n_para_12)%>%
  ggplot(., aes(x=log_qpcr, y=conc, color=factor(total_n_para_12)))+
  geom_point()+
  theme_minimal()+
  facet_wrap(~targetName+rev(factor(timepoint)), scales="fixed")+
  viridis::scale_color_viridis(discrete = T)+
  theme(legend.title = element_blank())

ggsave("~/Downloads/n_para_qpcr_plot.png", n_para_qpcr_plot, width=12, height=6, dpi=444, bg="white")
