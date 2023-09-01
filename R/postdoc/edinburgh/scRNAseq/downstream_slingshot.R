library(Seurat)
library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(tidyr)
library(slingshot)

red_half <- colorspace::sequential_hcl(8, palette = "Lajolla")
blue_half <- colorspace::sequential_hcl(6, palette = "Mako")[-1]

cluster_palette <- c(red_half, blue_half)





data=NULL
gc()
data <- readRDS("~/postdoc/edinburgh/scRNAseq/seurat_objecitons/harmont_slingshot_end11.sce")
# data <- readRDS("~/postdoc/edinburgh/scRNAseq/seurat_objecitons/curvy_harmont_slingshot_end11.sce")

colnames(colData(data))

try <- data.frame("umap1"=reducedDim(data, "UMAP")[,1],
                  "umap2"=reducedDim(data, "UMAP")[,2],
                  "seurat_clusters"=data[["seurat_clusters"]],
                  "sample_ID"=data[["sample_ID"]],
                  "slingPseudotime_1"=data[["slingPseudotime_1"]],
                  "slingPseudotime_2"=data[["slingPseudotime_2"]],
                  "slingPseudotime_3"=data[["slingPseudotime_3"]],
                  "slingPseudotime_4"=data[["slingPseudotime_4"]],
                  "slingPseudotime_5"=data[["slingPseudotime_5"]],
                  "slingPseudotime_6"=data[["slingPseudotime_6"]]
)

try <- pivot_longer(try, cols=c(slingPseudotime_1,
                                slingPseudotime_2,
                                slingPseudotime_3,
                                slingPseudotime_4,
                                slingPseudotime_5,
                                slingPseudotime_6
),
names_to="slinger", values_to="sling_value")

try %>%
  # filter(slinger %in% c("slingPseudotime_5", "slingPseudotime_6"))%>%
  filter(!is.na(sling_value))%>%
  # mutate(alpha=if_else(is.na(sling_value), -10, sling_value))%>%
  # arrange(alpha)%>%
  arrange(sling_value)%>%
  ggplot(., aes(x=umap1, y=umap2, color=sling_value))+
  geom_point(na.rm = TRUE)+
  #facet_wrap(~slinger)+
  theme_minimal()+
  scale_color_viridis_c(option = "B", na.value = "white")

na_grey <- try %>%
  filter(slinger %in% c("slingPseudotime_5", "slingPseudotime_6"))%>%
  #filter(!is.na(sling_value))%>%
  mutate(pseudotime=if_else(is.na(sling_value), NA, sling_value))%>%
  arrange(!is.na(pseudotime), pseudotime)%>%
  ggplot(., aes(x=umap1, y=umap2, color=pseudotime))+
  geom_point(na.rm = TRUE, size=1.2)+
  facet_wrap(~slinger)+
  theme_minimal()+
  scale_color_viridis_c(option = "B")

ggsave("~/postdoc/edinburgh/scRNAseq/scg/vac63c/final_revision/final_figures/grey_slingshot.pdf", na_grey, width=8, height=4)


na_black <- try %>%
  filter(slinger %in% c("slingPseudotime_5", "slingPseudotime_6"))%>%
  #filter(!is.na(sling_value))%>%
  mutate(pseudotime=if_else(is.na(sling_value), NA, sling_value))%>%
  arrange(!is.na(pseudotime), pseudotime)%>%
  ggplot(., aes(x=umap1, y=umap2, color=pseudotime))+
  geom_point(na.rm = TRUE, size=1.2)+
  facet_wrap(~slinger)+
  theme_minimal()+
  scale_color_viridis_c(option = "B", na.value = "black")

ggsave("~/postdoc/edinburgh/scRNAseq/scg/vac63c/final_revision/final_figures/black_slingshot.pdf", na_black, width=8, height=4)



na_black_all <- try %>%
  # filter(slinger %in% c("slingPseudotime_5", "slingPseudotime_6"))%>%
  #filter(!is.na(sling_value))%>%
  mutate(pseudotime=if_else(is.na(sling_value), NA, sling_value))%>%
  arrange(!is.na(pseudotime), pseudotime)%>%
  ggplot(., aes(x=umap1, y=umap2, color=pseudotime))+
  geom_point(na.rm = TRUE, size=1.2)+
  facet_wrap(~slinger)+
  theme_minimal()+
  scale_color_viridis_c(option = "B", na.value = "black")

ggsave("~/postdoc/edinburgh/scRNAseq/scg/vac63c/final_revision/final_figures/black_all_slingshot.pdf", na_black_all, width=12, height=8)




na_grey2 <- try %>%
  filter(slinger %in% c("slingPseudotime_5", "slingPseudotime_6"))%>%
  # filter(seurat_clusters %in% c("10", "11"))%>%
  #filter(!is.na(sling_value))%>%
  mutate(pseudotime=if_else(is.na(sling_value), NA, sling_value))%>%
  arrange(!is.na(pseudotime), pseudotime)%>%
  ggplot(., aes(x=umap1, y=umap2, color=pseudotime))+
  geom_point(na.rm = TRUE, size=1.2)+
  facet_grid(slinger~seurat_clusters)+
  theme_minimal()+
  scale_color_viridis_c(option = "B")

ggsave("~/postdoc/edinburgh/scRNAseq/scg/vac63c/final_revision/final_figures/grey_slingshot2.pdf", na_grey2, width=8, height=4)


# alternative way of slingthooting
rd <- reducedDims(data)$UMAP
cl <- data$seurat_clusters

colors <- colorRampPalette(RColorBrewer::brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(data$slingPseudotime_6, breaks=100)]

plot(rd, col = plotcol, pch=16, asp = 1)
lines(slingshot::SlingshotDataSet(crv1), lwd = 3, col = 'black')
