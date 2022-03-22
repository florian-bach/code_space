library(readr)
library(Rmagic)
library(phateR)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

test_data <- Read10X(data.dir = "~/postdoc/scRNAseq/cell_ranger_adventure_7/sample_feature_bc_matrix/")
vac63c <- CreateSeuratObject(counts = test_data, project = "vac63c_cd4_tcells", min.cells = 3, min.features = 200)


vac63c[["percent.mt"]] <- PercentageFeatureSet(vac63c, pattern = "^MT-")

#filter this bad boy
pbmc <- subset(vac63c, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

# scaling data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)

#determine 2k most variable genes
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
variable_genes <- VariableFeatures(pbmc)

#subset data to only include those genes
small_pbmc <- subset(pbmc, features =  variable_genes)

# PHATE wants cells on rows, features on columns
small_pbmc_matrix <- data.frame(t(small_pbmc@assays$RNA@scale.data))

# randomly subsample 5555 cells to speed up plotting
set.seed(1234); small_pbmc_matrix <- small_pbmc_matrix[floor(runif(5555, min=1, max=nrow(small_pbmc_matrix))),]


phate_object <- phate(small_pbmc_matrix, knn = 7, decay = 70, t=7)






feature1 <- "CXCL2"
plot1 <- ggplot(phate_object)+
  geom_point(aes(PHATE1, PHATE2, color=small_pbmc_matrix[,feature1]))+
  labs(color=feature1)+
  ggtitle(paste("t =", i))+
  ggtitle(feature1)+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_blank())+
  viridis::scale_color_viridis(option="B")

feature2 <- "LMNA"
plot2 <- ggplot(phate_object)+
  geom_point(aes(PHATE1, PHATE2, color=small_pbmc_matrix[,feature2]))+
  labs(color=feature2)+
  ggtitle(feature2)+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_blank())+
  viridis::scale_color_viridis(option="B")

feature3 <- "IFNG"
plot3 <- ggplot(phate_object)+
  geom_point(aes(PHATE1, PHATE2, color=small_pbmc_matrix[,feature3]))+
  labs(color=feature3)+
  ggtitle(feature3)+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_blank())+
  viridis::scale_color_viridis(option="B")

feature4 <- "TRAV38.2DV8"
plot4 <- ggplot(phate_object)+
  geom_point(aes(PHATE1, PHATE2, color=small_pbmc_matrix[,feature4]))+
  labs(color=feature4)+
  ggtitle(feature4)+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_blank())+
  viridis::scale_color_viridis(option="B")

combo_plot <- cowplot::plot_grid(plot1, plot2, plot3, plot4, nrow = 2)

ggsave(paste("~/postdoc/scRNAseq/exploratory_plots/phate", feature1, feature2, feature3, paste(feature4, ".png", sep=''), sep='_'), combo_plot, height=8, width=10, bg="white")


