library(Seurat)
library(reticulate)
use_python("/Users/fbach/.virtualenvs/r-reticulate/bin/python", required = TRUE)
leidenalg <- reticulate::import("leidenalg")
# reticulate::import("numpy")
reticulate::py_install(packages = "pandas")
reticulate::py_install(packages = "leidenalg")

inferno <- colorspace::sequential_hcl("inferno", n=8)


hdfs <- data.frame(matrix(0,6))
hdfs$path <- list.files("~/Downloads/GSE182536_RAW/", pattern = ".h5", full.names = TRUE)
hdfs$file <- list.files("~/Downloads/GSE182536_RAW/", pattern = ".h5", full.names = FALSE)
hdfs$name <- substr(hdfs$file, 12, 18)

t_cells <- Seurat::Read10X_h5(hdfs$path[1])
t_cells <- CreateSeuratObject(counts = t_cells)
t_cells$sample <- hdfs$name[1]

for(i in 2:length(hdfs$path)){
  count_matrix <- Seurat::Read10X_h5(hdfs$path[i])
  count_seurat <- CreateSeuratObject(counts = count_matrix)
  count_seurat$sample <- hdfs$name[i]
  t_cells <- merge(t_cells, count_seurat)
}


t_cells$memory <- ifelse(
  is.na(stringr::str_match(t_cells$sample, "N[0-9]*")),
        "memory", "naive")


t_cells$participant <- dplyr::recode(t_cells$sample,
                              "Mem104-"="v104",
                              "N104-5G"="v104",
                              "MEM79-5"="v79",
                              "N79-5GE"="v79",
                              "mem83-5"="v83",
                              "N83-5GE"="v83"
                              )


t_cells[["percent.mt"]] <- PercentageFeatureSet(t_cells, pattern = "^MT-")
VlnPlot(t_cells, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

t_cells <- subset(t_cells, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

t_cells <- NormalizeData(t_cells, normalization.method = "LogNormalize", scale.factor = 10000)
t_cells <- FindVariableFeatures(t_cells, selection.method = "vst", nfeatures = 2500)
all.genes <- rownames(t_cells)
t_cells <- ScaleData(t_cells, features = all.genes)



t_cells <- RunPCA(t_cells, features = VariableFeatures(object = t_cells))

t_cells <- harmony::RunHarmony(t_cells, "participant")
t_cells <- RunUMAP(t_cells, dims = 1:15, reduction = "harmony")
t_cells <- FindNeighbors(t_cells, dims = 1:15, reduction = "harmony")
t_cells <- FindClusters(t_cells, resolution = 0.8, algorithm = 4)

# t_cells <- subset(t_cells, seurat_clusters != 11)
DimPlot(t_cells, label=TRUE)

# 
FeaturePlot(t_cells, features = c("NKG7", "GNLY", "GZMA", "GZMK", "GZMB", "PRF1", "GZMH", "CD38", "ICOS", "CXCR5", "HLA-DRB1", "TCF7"), order=TRUE, cols=inferno)
FeaturePlot(t_cells, features = c("LAG3", "TNFRSF4", "TNFRSF18"), order=TRUE, cols=inferno)
# FeaturePlot(t_cells, features = c("TCF7", "MKI67", "GZMB", "GZMK", "GNLY", "NKG7", "PRF1", "GZMH", "KLRB1", "CD27"), order=TRUE, cols=inferno)
# FeaturePlot(t_cells, features = c("TCF7", "MKI67", "GZMB", "GZMK", "GNLY", "NKG7", "PRF1", "GZMH", "KLRB1", "CD27", "TNF","TNFRSF1B"), order=TRUE, cols=inferno)
# 
# DimPlot(t_cells, group.by = "participant")
DotPlot(t_cells, features = c("TCF7", "MKI67", "GZMB", "GZMK", "GNLY", "NKG7", "PRF1", "GZMH", "KLRB1", "CD27", "TNF","TNFRSF1B"))


cluster_counts <- data.frame(table(t_cells@active.ident, t_cells@meta.data$participant, t_cells@meta.data$memory))
colnames(cluster_counts) <- c("Cluster_ID", "Participant", "Memory", "Count")

cluster_counts <- cluster_counts %>%
  group_by(Participant, Memory) %>%
  mutate("Percentage"=Count/sum(Count))


cluster_counts %>%
  filter(Memory=="memory")%>%
  ggplot(., aes(x=Cluster_ID, y=Percentage, fill=Participant))+
  geom_boxplot()+
  geom_point(position = position_dodge(width = 0.75), aes(color=Participant, group=Participant))+
  scale_y_continuous(labels = scales::label_percent())+
  scale_fill_manual(values=c("#AA3377", "#4477AA", "deepskyblue"))+
  #scale_color_manual(values = volunteer_colors)+
  facet_wrap(~Cluster_ID, scales="free")+
  theme_minimal()+
  theme(strip.text = element_text(hjust = 0))


FeaturePlot(t_cells, features = c("GZMB", "CCL3", "CCL4", "ADRB2", "FGFBP2", "SPON2", "ADGRG1"), order=TRUE, cols=inferno)


tr1_markers <- c("CD38", "LAG3","HAVCR2", "PDCD1", "KLRB1", "IL21", "IL10", "IFNG", "GZMA", "GZMK", "CCL5", "CXCR3", "CXCR4", "CXCR6", "CCR2", "CCR5", "CCR7", "CCR9")
#no icos, no hla-dr, no pdl1
tr1_plot <- FeaturePlot(t_cells, features = tr1_markers, order=TRUE, cols = rev(inferno))
tr1_dot_plot <- DotPlot(t_cells, features = tr1_markers)+theme(axis.text.x = element_text(angle=90, hjust=1))

ggsave("~/postdoc/stanford/misc/furtado_tr1_markers.png", tr1_dot_plot, width = 8, height = 4, bg="white")

cxcr <- grep("CXCR[1-9]", rownames(t_cells), value = TRUE)
cxcl <- grep("CXCL[1-9]", rownames(t_cells), value = TRUE)
ccr <- grep("CCR[1-9]", rownames(t_cells), value = TRUE)
ccl <- grep("CCL[1-9]", rownames(t_cells), value = TRUE)
tbs <- grep("TBX[1-9]", rownames(t_cells), value = TRUE)
gzm <- c("PRF1", "GNLY", "NKG7", grep("GZM[A-Z]", rownames(t_cells), value = TRUE))
CX3C <-grep("^CX3C", rownames(t_cells), value = TRUE)
FeaturePlot(t_cells, features = CX3C, order=TRUE)
#tcf3 negative
tcf <- grep("TCF[1-9]", rownames(t_cells), value = TRUE)
foxp <- grep("FOXP[1-9]", rownames(t_cells), value = TRUE)
foxo <- grep("FOXO[1-9]", rownames(t_cells), value = TRUE)
ifna <- grep("^IFNA[1-9]", rownames(t_cells), value = TRUE)
ifnb <- grep("^IFNB[1-9]", rownames(t_cells), value = TRUE)
ils <- grep("^IL[1-9]", rownames(t_cells), value = TRUE)
ils <- grep("^IL[1-9]", rownames(t_cells), value = TRUE)
notch <- grep("^NOTCH[1-9]", rownames(t_cells), value = TRUE)
stat <- grep("^STAT[1-9]", rownames(t_cells), value = TRUE)
fcr <- grep("^FCGR[1-9]", rownames(t_cells), value = TRUE)
NCAM <- grep("^NCAM[1-9]", rownames(t_cells), value = TRUE)
KLR <- grep("^KLR[A-Z][1-9]", rownames(t_cells), value = TRUE)

FeaturePlot(t_cells, KLR, order=TRUE)

grep("LAG[1-9]", rownames(t_cells), value = TRUE)

isgs <- c("STAT1", "STAT2",  "IRF1", "IRF2", "IRF7", "IRF9", "MYD88",
                               "TLR4", "IDO1", "IDO2", "ACOD1", "GBP1",
                               "GBP2", "GBP3", "GBP4", "GBP5", "GBP6", "SOD1", "SOD2", "SOD3",
                               "S100A8", "S100A9", "HIF1A", "HMOX1", "HMOX2", "ICAM1",
                               "CD40", "PDCD1", "CD274", "PDCD1LG2", "CXCL11", "CXCL10", "CCL2",
                               "CCL25", "IL27", "CCL23", "TNFSF13B", "IL1RN", "TNF", "IL15",
                               "IL1B", "CSF1", "TNFSF13", "TGFB1", "IL1A", "IL18", "IL18", "IL7",
                               "IL12B", "CSF2", "LTA", "IFNB1", "IFNA1", "IL10", "IL6", "IL12A",
                               "CXCL8", "AIM2", "OASL", "MX1", "MX2", "DDX58", "IL18", "CXCL9", "IL6", "IL21")

replication <- readxl::read_excel("~/Downloads/GO_term_summary_20240213_170759.xlsx")
replication_genes <- toupper(unique(replication$Symbol))

# no CCDC88A, GMNN, ORC6, PURA; DONSON? ZPR2?
FeaturePlot(t_cells, features = c("LAG3", "CCDC88A", "GMNN", "PURA", "DONSON", "ZPR1"), order=TRUE, cols = rev(inferno))

rep_plot <- FeaturePlot(t_cells, features = c("LAG3", replication_genes[1:90]), order=TRUE, cols = rev(inferno))
ggsave("mega_rep1_90.png", rep_plot, width=10, height=60, bg="white", limitsize = F)


rep_plot <- FeaturePlot(t_cells, features = c("LAG3", replication_genes[91:190]), order=TRUE, cols = rev(inferno))
ggsave("mega_rep_91_190.png", rep_plot, width=10, height=60, bg="white", limitsize = F)


rep_plot <- FeaturePlot(t_cells, features = c("LAG3", replication_genes[191:length(replication_genes)]), order=TRUE, cols = rev(inferno))
ggsave("mega_rep_191_290.png", rep_plot, width=10, height=60, bg="white", limitsize = F)


rep_plot <- FeaturePlot(t_cells, features = c("LAG3", replication_genes[291:390]), order=TRUE, cols = rev(inferno))
ggsave("mega_rep_291_390.png", rep_plot, width=10, height=60, bg="white", limitsize = F)

rep_plot <- FeaturePlot(t_cells, features = c("LAG3", replication_genes[391:487]), order=TRUE, cols = rev(inferno))
ggsave("mega_rep_391_487.png", rep_plot, width=10, height=60, bg="white", limitsize = F)

