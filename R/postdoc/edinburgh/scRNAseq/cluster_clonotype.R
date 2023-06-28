
library(scRepertoire)
library(Seurat)
library(ggplot2)
clonetype_palette <- colorspace::sequential_hcl(n=7, palette="Lajolla")

# previous shizniz ####
big_data <-readRDS("~/postdoc/edinburgh/scRNAseq/seurat_objecitons/no_12_harmony.seurat")

# first_seurat <- subset(big_data, subset = n_infection=="First")
# third_seurat <- subset(big_data, subset = n_infection=="Third")

first_contig <- read.csv("~/postdoc/edinburgh/scRNAseq/scg/vac63c/first_all_contig_annotations.csv")
third_contig <- read.csv("~/postdoc/edinburgh/scRNAseq/scg/vac63c/third_all_contig_annotations.csv")

# fix the barcodes that had to be modified before merging back in the day
first_contig$barcode <- paste("First_", substr(first_contig$barcode, 1, 16), sep="")
third_contig$barcode <- paste("Third_", substr(third_contig$barcode, 1, 16), sep="")

combo_contig <- rbind(first_contig, third_contig)
contig.list <- createHTOContigList(combo_contig, big_data, group.by = "seurat_clusters")


combo_combined <- combineTCR(contig.list,
                             samples = paste("cluster", names(contig.list), sep=""),
                             cells ="T-AB")

# saveRDS(combo_combined, "~/postdoc/edinburgh/scRNAseq/scg/vac63c/combo_clusters_first_contigs_combined")

# meat ####

#combo_combined <- readRDS("~/postdoc/edinburgh/scRNAseq/seurat_objecitons/combo_contigs_combined")

combo_homo <- clonalHomeostasis(combo_combined, cloneCall = "aa")

ggsave("~/vac63c/final_revision/combo_homo.png", combo_homo, height=4, width=16, bg="white")

combo_prop_plot <- clonalProportion(combo_combined, cloneCall = "aa") 

ggsave("~/vac63c/final_revision/combo_prop_plot.png", combo_prop_plot, height=4, width=16, bg="white")


cluster_list <- paste("cluster", seq(0,12), sep="")

amino_cluster_alluvial <- compareClonotypes(combo_combined,
                                      numbers = 50,
                                      #clonotypes = ten_eleven_shared_tcrs,
                                      samples = cluster_list,
                                      cloneCall="aa",
                                      graph = "alluvial",
                                      #exportTable = TRUE
                                      )
  
amino_cluster_alluvial <- amino_cluster_alluvial+theme(legend.position = "none")
ggsave("~/postdoc/edinburgh/scRNAseq/final_revision/amino_cluster_alluvial.png", amino_cluster_alluvial, height=4, width=8, bg="white", dpi=444, limitsize = FALSE)



gene_cluster_alluvial <- compareClonotypes(combo_combined,
                                      numbers = 50,
                                      #samples = cluster_list,
                                      cloneCall="gene",
                                      graph = "alluvial",
                                      #exportTable = TRUE
)

gene_cluster_alluvial <- gene_cluster_alluvial+theme(legend.position = "none")
ggsave("~/postdoc/edinburgh/scRNAseq/final_revision/gene_cluster_alluvial.png", gene_cluster_alluvial, height=4, width=8, bg="white", dpi=444, limitsize = FALSE)



fresh_barcodes <- lapply(combo_combined, function(x)substr(x$barcode, regexpr("_", x$barcode)+1, nchar(x$barcode)))
combo_combined2 <- do.call(rbind, combo_combined)
combo_combined2$barcode <- substr(combo_combined2$barcode, regexpr("_", combo_combined2$barcode)+1, nchar(combo_combined2$barcode))

# 
clono_bins <- c(one=1, two=2, three=3, four=4, five=5, six=6, seven=7, eight=8, ten=10, twenty=20, fifty=50)
split_combo_combined2 <- split(combo_combined2, combo_combined2$sample)
names(split_combo_combined2) <- names(combo_combined)
big_data <- combineExpression(df = split_combo_combined2,
                              sc = big_data,
                              cloneCall="aa",
                              proportion = FALSE,
                              cloneTypes=clono_bins
                              )

slot(big_data, "meta.data")$cloneType <- factor(slot(big_data, "meta.data")$cloneType,
                                              levels = unique(slot(big_data, "meta.data")$cloneType))

clonal_umap <- DimPlot(big_data, group.by = "cloneType", split.by = "volunteer", cols = clonetype_palette, order = TRUE, pt.size = 2, ncol = 3)
ggsave("~/postdoc/edinburgh/scRNAseq/final_revision/clonal_umap.png", clonal_umap, height=8, width=12, bg="white", dpi=444, limitsize = FALSE)




clono_freq <- amino_cluster_alluvial$data

repeat_clonos <- clono_freq%>%
  group_by(Clonotypes)%>%
  summarise("how_often"=n())%>%
  filter(how_often>1)

clono_count <- clono_freq %>%
  filter(Clonotypes %in% repeat_clonos$Clonotypes, Sample %in% c("cluster10", "cluster11"))%>%
  group_by(Clonotypes, Sample)%>%
  summarise("how_often"=n())

ten_eleven_shared_tcrs <- c("CAMSAGEKLTF_CSASPGMGGYTF", "CAVTGNQFYF_CAWSASSRDTQYF")

big_data <- highlightClonotypes(big_data, 
                                cloneCall= "aa", 
                                sequence = as.character(clono_count$Clonotypes))

shared_tcr_umap <- DimPlot(big_data, group.by = "highlight", split.by = "volunteer", pt.size = 2, order=TRUE, ncol=3) + 
  theme(plot.title = element_blank())
ggsave("~/postdoc/edinburgh/scRNAseq/final_revision/shared_tcr_umap.png", shared_tcr_umap, height=8, width=12, bg="white", dpi=444, limitsize = FALSE)

shared_tcr_umap2 <- DimPlot(big_data, group.by = "highlight", split.by = "seurat_clusters", pt.size = 2, order=TRUE, ncol=3) + 
  guides(color=guide_legend(ncol = 1))+
  theme(plot.title = element_blank(), 
        legend.position = "bottom")
ggsave("~/postdoc/edinburgh/scRNAseq/final_revision/cluster_wise_shared_tcr_umap.png", shared_tcr_umap2, height=20, width=8, bg="white", dpi=444, limitsize = FALSE)

shared_tcr_umap3 <- DimPlot(big_data, group.by = "highlight", split.by = "volunteer", pt.size = 2, order=TRUE, ncol=3) + 
  theme(plot.title = element_blank())
ggsave("~/postdoc/edinburgh/scRNAseq/final_revision/volunteer_wise_shared_tcr_umap.png", shared_tcr_umap3, height=8, width=12, bg="white", dpi=444, limitsize = FALSE)

multi_clones <- subset(big_data, seurat_clusters%in%c(7,8,11))

dimmer <- DimPlot(multi_clones, label=TRUE)
shared_tcr_umap4 <- DimPlot(object = multi_clones, group.by = "highlight", split.by = "volunteer", pt.size = 2, order=TRUE, ncol=3)

left_half <- shared_tcr_umap4$data
right_half <- dimmer$data
left_half$barcode <- rownames(left_half)
right_half$barcode <- rownames(right_half)

both_halves <- left_join(left_half, right_half, by="barcode")

volunteer_wise_shared_tcr_umap2 <- ggplot(both_halves, aes(x=UMAP_1.x, y=UMAP_2.x, color=highlight, label=ident))+
  geom_point()+
  ggrepel::geom_label_repel(data=na.omit(both_halves))+
  facet_wrap(~volunteer)+
  theme_minimal()

ggsave("~/postdoc/edinburgh/scRNAseq/final_revision/volunteer_wise_shared_tcr_umap2.png", volunteer_wise_shared_tcr_umap2, height=8, width=18, bg="white", dpi=444, limitsize = FALSE)


library(ggraph)

#No Identity filter
clonalNetwork(big_data, 
              reduction = "umap", 
              identity = "ident",
              filter.clones = clono_count$Clonotypes,
              filter.identity = c(7, 8, 11),
              cloneCall = "aa")

saveRDS(big_data, "~/postdoc/edinburgh/scRNAseq/seurat_objecitons/no_12_harmony_with_tcr.seurat")
