library(dplyr)
library(patchwork)
bc1<- haven::read_dta("~/postdoc/stanford/clinical_data/BC1/MergedAntibodyData_ChildClinical.dta")

slim_bc <- data.frame(bc1 %>%
  select(id, wealthcat, totalmalariapreg, anymalariapreg, gender, gestage, anyparasitemia, timepoint, grep("log", colnames(bc1), value = TRUE)))

# NAs are present. counts NAs per column
count_na <- function(x){table(is.na(x))}

na_tables<- lapply(slim_bc, count_na)
na_counts <- do.call(rbind, na_tables)

# NAs mostly come from logGST(146) and logpd(448); get rid of those columns, NA.omit and let's roll

slimmer_bc <- slim_bc %>%
  select(-logpd, logGST) %>%
  na.omit()

big_pca <-  prcomp(slimmer_bc[,grep("log", colnames(slimmer_bc), value = TRUE)], center = T)
pca_plot_data <- as.data.frame(cbind(slimmer_bc[, 1:8], big_pca$x))

# pca loadings figure

loadings_df <- data.frame(big_pca$rotation)
loadings_df$antibody <- rownames(loadings_df)
loadings_df$antibody <- factor(loadings_df$antibody, levels = loadings_df$antibody[order(loadings_df$PC1)])

pc1_cols <- colorspace::sequential_hcl(nrow(loadings_df), palette = "Purple Yellow")
names(pc1_cols) <- loadings_df$antibody[order(loadings_df$PC1)]
                                        
PC1_plot <- ggplot(loadings_df, aes(x=factor(antibody, levels = antibody[order(loadings_df$PC1)]), y=PC1, fill=antibody))+
  geom_bar(stat = "identity")+
  scale_fill_manual(values = pc1_cols)+
  theme_minimal()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust=1),
        legend.position = "none")


pc2_cols <- colorspace::sequential_hcl(nrow(loadings_df), palette = "Purple Yellow")
names(pc2_cols) <- loadings_df$antibody[order(loadings_df$PC2)]

loadings_df$antibody <- factor(loadings_df$antibody, levels = loadings_df$antibody[order(loadings_df$PC2)])
PC2_plot <- ggplot(loadings_df, aes(x=factor(antibody, levels = antibody[order(loadings_df$PC2)]), y=PC2, fill=antibody))+
  geom_bar(stat = "identity")+
  scale_fill_manual(values = pc2_cols)+
  theme_minimal()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, hjust=1),
        legend.position = "none")


combo_plot <- PC1_plot+PC2_plot


# ginormous_df <- pivot_longer(pca_plot_data, cols = 1:8, names_to = "demo_data", values_to = "demo_value")

# tl;dr this is the only one where there's seperation
ggplot(pca_plot_data, aes(x=PC1, y=PC2, color=factor(timepoint)))+
  geom_point()+
  xlab(paste("PC1 ", data.frame(summary(big_pca)[6])[2,1]*100, "%", sep = ""))+
  ylab(paste("PC2 ", data.frame(summary(big_pca)[6])[2,2]*100, "%", sep = ""))+
  theme_minimal()+
  #ggrepel::geom_label_repel(aes_string(label = "id"), show.legend = FALSE)+ 
  ggtitle("Antibody Titre PCA")

# other PCAs ####
wealthcat_plot <- ggplot(pca_plot_data, aes(x=PC1, y=PC2, color=factor(wealthcat)))+
  geom_point()+
  xlab(paste("PC1 ", data.frame(summary(big_pca)[6])[2,1]*100, "%", sep = ""))+
  ylab(paste("PC2 ", data.frame(summary(big_pca)[6])[2,2]*100, "%", sep = ""))+
  theme_minimal()+
  #ggrepel::geom_label_repel(aes_string(label = "id"), show.legend = FALSE)+ 
  ggtitle("Antibody Titre PCA")

  
totalmalariapreg_plot <- ggplot(pca_plot_data, aes(x=PC1, y=PC2, color=factor(totalmalariapreg)))+
    geom_point()+
    xlab(paste("PC1 ", data.frame(summary(big_pca)[6])[2,1]*100, "%", sep = ""))+
    ylab(paste("PC2 ", data.frame(summary(big_pca)[6])[2,2]*100, "%", sep = ""))+
    theme_minimal()+
    #ggrepel::geom_label_repel(aes_string(label = "id"), show.legend = FALSE)+ 
    ggtitle("Antibody Titre PCA")
  
anymalariapreg_plot <-  ggplot(pca_plot_data, aes(x=PC1, y=PC2, color=factor(anymalariapreg)))+
    geom_point()+
    xlab(paste("PC1 ", data.frame(summary(big_pca)[6])[2,1]*100, "%", sep = ""))+
    ylab(paste("PC2 ", data.frame(summary(big_pca)[6])[2,2]*100, "%", sep = ""))+
    theme_minimal()+
    #ggrepel::geom_label_repel(aes_string(label = "id"), show.legend = FALSE)+ 
    ggtitle("Antibody Titre PCA")
  
gender_plot <-  ggplot(pca_plot_data, aes(x=PC1, y=PC2, color=factor(gender)))+
    geom_point()+
    xlab(paste("PC1 ", data.frame(summary(big_pca)[6])[2,1]*100, "%", sep = ""))+
    ylab(paste("PC2 ", data.frame(summary(big_pca)[6])[2,2]*100, "%", sep = ""))+
    theme_minimal()+
    #ggrepel::geom_label_repel(aes_string(label = "id"), show.legend = FALSE)+ 
    ggtitle("Antibody Titre PCA")
  
  
gestage_plot <-  ggplot(pca_plot_data, aes(x=PC1, y=PC2, color=gestage))+
    geom_point()+
    xlab(paste("PC1 ", data.frame(summary(big_pca)[6])[2,1]*100, "%", sep = ""))+
    ylab(paste("PC2 ", data.frame(summary(big_pca)[6])[2,2]*100, "%", sep = ""))+
    theme_minimal()+
    scale_color_gradient(low="#073376", high="#FFAE42")+
    #ggrepel::geom_label_repel(aes_string(label = "id"), show.legend = FALSE)+ 
    ggtitle("Antibody Titre PCA")
  
  
anyparasitemia_plot <-  ggplot(pca_plot_data, aes(x=PC1, y=PC2, color=anyparasitemia))+
    geom_point()+
    xlab(paste("PC1 ", data.frame(summary(big_pca)[6])[2,1]*100, "%", sep = ""))+
    ylab(paste("PC2 ", data.frame(summary(big_pca)[6])[2,2]*100, "%", sep = ""))+
    theme_minimal()+
    scale_color_gradient(low="#073376", high="#FFAE42")+
    ggtitle("Antibody Titre PCA")
  
  
timepoint_plot <-  ggplot(pca_plot_data, aes(x=PC1, y=PC2, color=factor(timepoint)))+
    geom_point()+
    xlab(paste("PC1 ", data.frame(summary(big_pca)[6])[2,1]*100, "%", sep = ""))+
    ylab(paste("PC2 ", data.frame(summary(big_pca)[6])[2,2]*100, "%", sep = ""))+
    theme_minimal()+
    #ggrepel::geom_label_repel(aes_string(label = "id"), show.legend = FALSE)+ 
    ggtitle("Antibody Titre PCA")
  
# stats ####  

