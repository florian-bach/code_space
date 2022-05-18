tertiary_dimplot <- DimPlot(pbmc, reduction = "umap", label = TRUE)
ggsave("~/postdoc/scRNAseq/tertiary/figures/tertiary_cluster_umap.png", tertiary_dimplot, height=4, width=5, bg="white")

cluster_freqs <- data.frame("cluster_ID"=names(table(pbmc@meta.data$seurat_clusters)),
                            "percentage"=unname(table(pbmc@meta.data$seurat_clusters))/ncol(pbmc)
)

tertiary_freq_plot <- ggplot(cluster_freqs, aes(x=factor(cluster_ID, levels = 0:13), y=percentage.Freq))+
  geom_bar(stat="identity", aes(fill=cluster_ID))+
  theme_minimal()+
  ggtitle("tertiary")+
  scale_y_continuous(label=scales::percent)+
  geom_text(aes(label= round(percentage.Freq, digits = 3)*100), vjust= -0.2, size=3.5)+
  theme(legend.position = "none",
        axis.title.x = element_blank())

ggsave("~/postdoc/scRNAseq/tertiary/figures/tertiary_freq_plot.png", tertiary_freq_plot, height=4, width=5, bg="white")


cyto_chemo <- FeaturePlot(pbmc, c("NKG7", "CCL5", "GZMB", "TNF", "GNLY", "CCL4"))
ggsave("~/postdoc/scRNAseq/tertiary/figures/cyto_chemo_umap.png",cyto_chemo, height=12, width=9, bg="white")

memory_plot <- FeaturePlot(pbmc, c("SELL","CCR7", "IL7R", "LMNA"))
ggsave("~/postdoc/scRNAseq/tertiary/figures/memory_umap.png",memory_plot, height=8, width=9, bg="white")















primary_dimplot <- DimPlot(pbmc, reduction = "umap", label = TRUE)
ggsave("~/postdoc/scRNAseq/primary/figures/primary_cluster_umap.png", primary_dimplot, height=4, width=5, bg="white")

cluster_freqs <- data.frame("cluster_ID"=names(table(pbmc@meta.data$seurat_clusters)),
                            "percentage"=unname(table(pbmc@meta.data$seurat_clusters))/ncol(pbmc)
)

primary_freq_plot <- ggplot(cluster_freqs, aes(x=factor(cluster_ID, levels = 0:17), y=percentage.Freq))+
  geom_bar(stat="identity", aes(fill=cluster_ID))+
  theme_minimal()+
  ggtitle("primary")+
  scale_y_continuous(label=scales::percent)+
  geom_text(aes(label= round(percentage.Freq, digits = 3)*100), vjust= -0.2, size=3.5)+
  theme(legend.position = "none",
        axis.title.x = element_blank())

ggsave("~/postdoc/scRNAseq/primary/figures/primary_freq_plot.png", primary_freq_plot, height=4, width=5, bg="white")


cyto_chemo <- FeaturePlot(pbmc, c("NKG7", "CCL5", "GZMB", "TNF", "GNLY", "CCL4"))
ggsave("~/postdoc/scRNAseq/primary/figures/cyto_chemo_umap.png",cyto_chemo, height=12, width=9, bg="white")

memory_plot <- FeaturePlot(pbmc, c("SELL","CCR7", "IL7R", "LMNA"))
ggsave("~/postdoc/scRNAseq/primary/figures/memory_umap.png",memory_plot, height=8, width=9, bg="white")
