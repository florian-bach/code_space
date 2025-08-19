##########################################
# 1. Load libraries
##########################################
library(tidyr)
library(dplyr)
library(FactoMineR)   # PCA
library(kml3d)        # 3D longitudinal clustering
library(kml)
library(ggplot2)

`%notin%`=Negate(`%in%`)
##########################################
# 2. Load your data (long format)
##########################################
mic_drop_key <- haven::read_dta("~/Downloads/MIC-DROP treatment assignments.dta")

clean_data <- read.csv("~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/clean_data_with_meta.csv")

df <- clean_data%>%
  filter(timepoint%in%c("8 weeks", "24 weeks", "52 weeks"),
         sample%notin%c("11660_tp24_Repeat"))%>%
  mutate(timepoint=factor(timepoint, levels=c("8 weeks", "24 weeks", "52 weeks", "68 weeks")))%>%
  filter(targetName %notin% c("CTSS", "LTA|LTB", "IFNA2"))%>%
  pivot_wider(names_from = targetName, values_from = conc, names_prefix = "Feature_")


##########################################
# 3. PCA for dimensionality reduction
##########################################

feature_cols <- grep("^Feature", names(df))  # feature columns

ncp=10
pca_res <- PCA(df[, feature_cols], scale.unit = TRUE, ncp = ncp, graph = FALSE)

# Extract first 3 PCs for trajectory clustering
pc_scores <- as.data.frame(pca_res$ind$coord[, 1:ncp])
names(pc_scores) <- paste("PC", 1:ncp, sep="")

##########################################
# 4. Merge PCA scores back with SampleID and Time
##########################################
df_pca <- bind_cols(df %>% select(id, timepoint), pc_scores)

##########################################
# 5. Reshape to wide format for kml3d
##########################################
# kml3d expects: rows = subjects, time along one dimension, variables = PC coords
df_wide <- df_pca %>%
  pivot_wider(names_from = timepoint, values_from = paste("PC", 1:ncp, sep=""))

# Drop SampleID for numeric matrix
mat <- df_wide %>% select(-id)

##########################################
# 6. Create kml3d data object
##########################################
traj_data <- clusterLongData3d(
  array(as.matrix(mat),
        dim = c(nrow(df_wide), length(unique(df_pca$timepoint)), ncp),
        dimnames = list(NULL, paste0("Time", sort(unique(df_pca$timepoint))),
                        paste("PC", 1:ncp, sep="")))
)

##########################################
# 7. Run kml3d clustering
##########################################
# Try 2 to 6 clusters, pick the best using Calinski–Harabasz
kml3d(traj_data, nbClusters = 2:10, nbRedrawing = 20)

##########################################
# 8. Get cluster assignments
##########################################
      # getClusters() lives here
# after you've run: kml3d(cld3d, nbClusters = 2:8, ...)

# 1) read the active criterion matrix (rows = c2, c3, ..., columns = redraws)
# crit_mat <- traj_data["criterionValuesAsMatrix", traj_data["criterionActif"]]
# 
# # 2) pick the K whose best redraw gives the highest value
# 2. Create a data matrix of trajectories (flatten over time & PCs)
# This must match the shape used in kml3d
k_values <- 2:8
traj_matrix <- mat # from your earlier wide-format step

# 3. Compute CH index for each K
ch_scores <- c()

for(k in k_values){
  # Get cluster labels for first redraw
  labs <- getClusters(traj_data, nbCluster = k, asInteger = T)
  
  # Compute distance matrix
  dist_mat <- daisy(traj_matrix)
  
  # Compute Calinski–Harabasz
  stats <- fpc::cluster.stats(d = dist_mat, clustering = labs)
  ch_scores = c(ch_scores, stats$ch)
  
  bestK <- k_values[which.max(ch_scores)]
  message("Best K according to CH index: ", bestK)
  
}

# 4. Pick bestK


# 3) extract per-sample cluster labels (as integers)
clusters <- kml::getClusters(traj_data, nbCluster = bestK, clusterRank = 1, asInteger = TRUE)

# 4) join back to IDs
df_clusters <- tibble(
  id = df_wide$id,
  Cluster  = paste("cluster ", as.integer(clusters), sep="")
)

##########################################
# 9. Visualize mean trajectories per cluster
##########################################
df_plot <- df_pca %>%
  left_join(df_clusters, by = "id") %>%
  group_by(Cluster, timepoint) %>%
  summarise(across(starts_with("PC"), mean), .groups = "drop")

(trajectory_plot1 <- ggplot(df_plot, aes(x = timepoint, y = PC1, color = factor(Cluster), group = Cluster)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  theme_minimal() +
  facet_wrap(~Cluster)+
  labs(title = "Cluster Mean Trajectories (PC1)",
       color = "Cluster"))

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/figures/trajectory_plot1.png",trajectory_plot1, width=8, height=4 , bg="white")


# inspect loadings ####

loadings <- data.frame(pca_res$var$coord)   # matrix: variables x PCs

top_analytes1 <- loadings%>%
  mutate(feature=rownames(.))%>%
  arrange(desc(abs(Dim.1)))%>%
  select(feature)

top_analytes5 <- loadings%>%
  mutate(feature=rownames(.))%>%
  arrange(desc(abs(Dim.5)))%>%
  select(feature)

top_9_genes1 <- substr(start=9, stop=20, x=top_analytes1$feature)[1:9]
next_top_9_genes1 <- substr(start=9, stop=20, x=top_analytes1$feature)[10:18]

top_9_genes5 <- substr(start=9, stop=20, x=top_analytes5$feature)[1:9]
next_top_9_genes5 <- substr(start=9, stop=20, x=top_analytes5$feature)[10:18]

clean_data2 <- clean_data %>%
  left_join(., df_clusters, by="id", )%>%
  left_join(., mic_drop_key, by="id")%>%
  mutate(treatmentarm=case_match(treatmentarm,
                          1~"Placebo",
                          2~"DP 1 year",
                          3~"DP 2 years"))
  

# individual analytes 

top9_pc1_loadings <- clean_data2%>%
  filter(timepoint%in%c("8 weeks", "24 weeks", "52 weeks"),
         sample%notin%c("11660_tp24_Repeat"))%>%
  mutate(timepoint=factor(timepoint, levels=c("8 weeks", "24 weeks", "52 weeks")))%>%
  filter(targetName %in% top_9_genes1)%>%
  # filter(Cluster %in% c("cluster 2", "cluster 4"))%>%
  ggplot(., aes(x=timepoint, y=conc, fill = factor(Cluster)))+
  geom_boxplot(outliers = F)+
  facet_wrap(~targetName,nrow=3, scales="free")+
  theme_minimal()

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/figures/top9_pc1_loadings.png",top9_pc1_loadings, width=8, height=6, bg="white")


top9_pc5_loadings <- clean_data2%>%
  filter(timepoint%in%c("8 weeks", "24 weeks", "52 weeks"),
         sample%notin%c("11660_tp24_Repeat"))%>%
  mutate(timepoint=factor(timepoint, levels=c("8 weeks", "24 weeks", "52 weeks")))%>%
  filter(targetName %in% top_9_genes5)%>%
  # filter(Cluster %in% c("cluster 2", "cluster 4"))%>%
  ggplot(., aes(x=timepoint, y=conc, fill = factor(Cluster)))+
  geom_boxplot(outliers = F)+
  facet_wrap(~targetName,nrow=3, scales="free")+
  theme_minimal()

ggsave("~/postdoc/stanford/plasma_analytes/MICDROP/big_experiment/figures/top9_pc1_loadings.png",top9_pc1_loadings, width=8, height=6, bg="white")



clean_data2%>%
  filter(timepoint%in%c("8 weeks", "24 weeks", "52 weeks"),
         sample%notin%c("11660_tp24_Repeat"))%>%
  mutate(timepoint=factor(timepoint, levels=c("8 weeks", "24 weeks", "52 weeks")))%>%
  filter(targetName %in% next_top_9_genes)%>%
  filter(Cluster %in% c("cluster 2", "cluster 4"))%>%
  ggplot(., aes(x=timepoint, y=conc, fill = factor(gender_categorical)))+
  geom_boxplot(outliers = F)+
  facet_wrap(~targetName)+
  theme_minimal()

hi = clean_data2%>%
  filter(treatmentarm!="DP 2 years")%>%
  mutate(any_mala=if_else(total_n_malaria==0, "none", "some"))%>%
  mutate(any_para=if_else(total_n_para==0, "none", "some"))%>%
  distinct(id, Cluster, gender, treatmentarm, any_mala, any_para)


table(hi$Cluster, hi$treatmentarm)
table(hi$Cluster, hi$any_mala)
table(hi$Cluster, hi$any_para)
