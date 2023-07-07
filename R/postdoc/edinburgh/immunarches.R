library(immunarch)
library(Seurat)


red_half <- colorspace::sequential_hcl(8, palette = "Lajolla")
blue_half <- colorspace::sequential_hcl(6, palette = "Mako")[-1]

cluster_palette <- c(red_half, blue_half)

# load, clean, annotate TCR data from first infection ####
first_file_path = "~/postdoc/edinburgh/scRNAseq/scg/vac63c/first_all_contig_annotations.csv"

# Load 10x data, fix barcode to match seurat
first_clonos <- repLoad(first_file_path)
first_clonos$data$first_all_contig_annotations$Barcode <- paste("First_", substr(first_clonos$data$first_all_contig_annotations$Barcode, 1, 16), sep="")

#load the clonotype table we made with scRepertoir
first_combined <- readRDS("~/postdoc/edinburgh/scRNAseq/seurat_objecitons/first_contigs_combined")
# glue it together
first_long <- do.call(rbind, first_combined)
# head(first_long$barcode)
# [1] "v320_Baseline_First_AAACCTGAGAATCTCC"
# [2] "v320_Baseline_First_AAACCTGAGAGCTATA"
# [3] "v320_Baseline_First_AAACCTGAGCCCAGCT"
# [4] "v320_Baseline_First_AAACCTGTCTGCGGCA"
# [5] "v320_Baseline_First_AAACGGGAGTTGAGAT"
# [6] "v320_Baseline_First_AAACGGGGTAGAGTGC"

# table(first_clonos$data$first_all_contig_annotations$Barcode %in% substr(first_long$barcode, nchar(first_long$barcode)-21, nchar(first_long$barcode)) )
# 
# substr(first_long$barcode, nchar(first_long$barcode)-21, nchar(first_long$barcode))

#grab repertoire df
first_clono_data <- first_clonos$data$first_all_contig_annotations
# grab sample id from scRepertoire object, matching the barcodes
first_clono_data$sample_id <- first_long$sample[match(first_clono_data$Barcode, substr(first_long$barcode, nchar(first_long$barcode)-21, nchar(first_long$barcode)))]
first_clono_data <- first_clono_data[!is.na(first_clono_data$sample_id),]
# create cell_id which is essentially the complicated barcode from scRepertoire
first_clono_data$cell_id <- paste(first_clono_data$sample_id, first_clono_data$Barcode, sep="")

# do it all again for third

third_file_path = "~/postdoc/edinburgh/scRNAseq/scg/vac63c/third_all_contig_annotations.csv"

#Load 10x data with repLoad
third_clonos <- repLoad(third_file_path)
third_clonos$data$third_all_contig_annotations$Barcode <- paste("Third_", substr(third_clonos$data$third_all_contig_annotations$Barcode, 1, 16), sep="")
#need to split third_clonos$data into n dfs according to sample_id

#load the clonotype table we made with scRepertoir

third_combined <- readRDS("~/postdoc/edinburgh/scRNAseq/seurat_objecitons/third_contigs_combined")
# glue it together
third_long <- do.call(rbind, third_combined)
# head(third_long$barcode)

third_clono_data <- third_clonos$data$third_all_contig_annotations
# grab sample id from scRepertoire object, matching the barcodes
third_clono_data$sample_id <- third_long$sample[match(third_clono_data$Barcode, substr(third_long$barcode, nchar(third_long$barcode)-21, nchar(third_long$barcode)))]
third_clono_data <- third_clono_data[!is.na(third_clono_data$sample_id),]
# create cell_id which is essentially the complicated barcode from scRepertoire
third_clono_data$cell_id <- paste(third_clono_data$sample_id, third_clono_data$Barcode, sep="")

#combo <- readRDS("/Volumes/lab_prasj/BIG_Flo/vac63c/final_revision/no_12_harmony.seurat")
first_seurat <- readRDS("~/postdoc/edinburgh/scRNAseq/seurat_objecitons/first_no_12_harmony.seurat")
third_seurat <- readRDS("~/postdoc/edinburgh/scRNAseq/seurat_objecitons/third_no_12_harmony.seurat")
# saveRDS(first_seurat, "/labs/prasj/BIG_Flo/vac63c/final_revision/first_no_12_harmony.seurat")
# saveRDS(third_seurat, "/labs/prasj/BIG_Flo/vac63c/final_revision/third_no_12_harmony.seurat")


# cluster & time vdj business ####
# 
# 
# #combo <- readRDS("/labs/prasj/BIG_Flo/vac63c/final_revision/final_10k.seurat")
# combo <- readRDS("/labs/prasj/BIG_Flo/vac63c/final_revision/no_12_harmony.seurat")
# first_seurat <- subset(combo, subset = n_infection=="First")
# 
# #get cluster idents frome seurat
# first_idents <- Idents(first_seurat)
# first_clono_data$seurat_cluster <- first_idents[match(first_clono_data$Barcode, names(first_idents))]
# first_clono_data$timepoint <- substr(first_clono_data$sample_id, 6, nchar(first_clono_data$sample_id))
# first_clono_data$cluster_time <- paste(first_clono_data$seurat_cluster, first_clono_data$timepoint, sep="_")
# 
# #split by cluster time, add metadata
# first_clono_split <- split(first_clono_data, first_clono_data$cluster_time)
# first_clonos$data <- first_clono_split
# first_clonos$meta <- as_tibble(data.frame(Sample = names(first_clono_split)))
# 
# first_clonos$meta$Timepoint <- substr(first_clonos$meta$Sample, regexpr("_", first_clonos$meta$Sample)+1, nchar(first_clonos$meta$Sample))
# first_clonos$meta$Seurat_Cluster <- substr(first_clonos$meta$Sample, 1, regexpr("_", first_clonos$meta$Sample))
# 
# # do it all again for third #
# 
# # WE NEED TO REPEAT THE ABOVE THREE STEPS WITH THE SEURAT IDENTS TO CREATE A LIST OF CLUSTER-LEVEL REPERTOIRES
# 
# #combo <- readRDS("/labs/prasj/BIG_Flo/vac63c/final_revision/final_10k.seurat")
# third_seurat <- subset(combo, subset = n_infection=="Third")
# #get cluster idents frome seurat
# third_idents <- Idents(third_seurat)
# third_clono_data$seurat_cluster <- third_idents[match(third_clono_data$Barcode, names(third_idents))]
# third_clono_data$timepoint <- substr(third_clono_data$sample_id, 6, nchar(third_clono_data$sample_id))
# third_clono_data$cluster_time <- paste(third_clono_data$seurat_cluster, third_clono_data$timepoint, sep="_")
# 
# #split by cluster time, add metadata
# third_clono_split <- split(third_clono_data, third_clono_data$cluster_time)
# third_clonos$data <- third_clono_split
# third_clonos$meta <- as_tibble(data.frame(Sample = names(third_clono_split)))
# 
# third_clonos$meta$Timepoint <- substr(third_clonos$meta$Sample, regexpr("_", third_clonos$meta$Sample)+1, nchar(third_clonos$meta$Sample))
# third_clonos$meta$Seurat_Cluster <- substr(third_clonos$meta$Sample, 1, regexpr("_", third_clonos$meta$Sample))
# 
# # put the two immunarches together ##
# 
# combo_clonos <- first_clonos
# 
# #iterate through list, combining each df
# combo_clonos$data <- lapply(1:length(first_clonos$data), function(x) rbind(first_clonos$data[[x]], third_clonos$data[[x]]))
# # meta is not a list of dfs but simply a tibble so this is not necessary
# #combo_clonos$meta <- lapply(1:length(first_clonos$meta), function(x) rbind(first_clonos$meta[[x]], third_clonos$meta[[x]]))
# 
# names(combo_clonos$data) <- names(first_clonos$data)
# # combo_clonos$meta$Timepoint <- substr(first_sample_clonos$meta$Sample, regexpr("_", first_sample_clonos$meta$Sample)+1, nchar(first_sample_clonos$meta$Sample))
# 
# # calculate gene usage statistics for human beta chain, exclude ambiguous assignments with more than one chain
# imm_gu <- geneUsage(combo_clonos$data, "hs.trbv", .norm=TRUE, .ambig="exc")
# 
# imm_gu_js <- geneUsageAnalysis(imm_gu, .method = "js", .verbose = TRUE)
# imm_gu_cor <- geneUsageAnalysis(imm_gu, .method = "cor", .verbose = TRUE)
# 
# p1 <- vis(imm_gu_js, .title = "Gene usage JS-divergence", .leg.title = "JS", .text.size = 1.5)
# p2 <- vis(imm_gu_cor, .title = "Gene usage correlation", .leg.title = "Cor", .text.size = 1.5)
# 
# p <- p1 + p2
# ggsave("~/postdoc/edinburgh/scRNAseq/final_revision/combo_cluster_time_wise_heats.png", p, height=8, width=16)
# 
# 
# # 
# # 
# # imm_raref <- repDiversity(combo_clonos$data, "raref", .verbose = F)
# # p1 <- vis(imm_raref, .do.norm=TRUE)
# # p2 <- vis(imm_raref, .by = "Timepoint", .meta = combo_clonos$meta)
# # 
# # p <- p1 + p2
# # ggsave("~/postdoc/edinburgh/scRNAseq/final_revision/normed_combo_cluster_time_wise_raref.png", p, height=8, width=16)
# # 
# # 
# # p1 <- vis(imm_raref)
# # p2 <- vis(imm_raref, .by = "Timepoint", .meta = combo_clonos$meta)
# # 
# # p <- p1 + p2
# # ggsave("~/postdoc/edinburgh/scRNAseq/final_revision/not_normed_combo_cluster_time_wise_raref.png", p, height=8, width=16)
# 
# # multidimensional scaling
# imm_cl_pca <- geneUsageAnalysis(imm_gu, "js+pca+kmeans", .verbose = TRUE)
# imm_cl_mds <- geneUsageAnalysis(imm_gu, "js+mds+kmeans", .verbose = TRUE)
# imm_cl_tsne <- geneUsageAnalysis(imm_gu, "js+tsne+kmeans", .perp = .01, .verbose = TRUE,)
# 
# p1 <- vis(imm_cl_pca, .plot = "clust")
# p2 <- vis(imm_cl_mds, .plot = "clust")
# p3 <- vis(imm_cl_tsne, .plot = "clust")
# 
# p <- p1 + p2 + p3
# ggsave("~/postdoc/edinburgh/scRNAseq/final_revision/combo_cluster_time_wise_mds.png", p, height=8, width=16)
# 
# 
# 
# time_vis <- vis(imm_gu, .by = "Timepoint", .meta = combo_clonos$meta)
# time_vis2 <- vis(imm_gu, .by = "Timepoint", .meta = combo_clonos$meta, .plot = "box")
# grid_vis <- vis(imm_gu, .grid = T)
# 
# vdj_freqs <- time_vis$data
# write.csv(vdj_freqs, file = "~/postdoc/edinburgh/scRNAseq/final_revision/final_figures/cluster_time_wise_vdj_freqs.csv", row.names = FALSE)
# 
# 
# ggsave("~/postdoc/edinburgh/scRNAseq/final_revision/beta_combo_cluster_time_wise_time_vis.png", time_vis, height=6, width=15, limitsize = FALSE)
# ggsave("~/postdoc/edinburgh/scRNAseq/final_revision/beta_combo_cluster_time_wise_time_vis2.png", time_vis2, height=15, width=15, limitsize = FALSE)
# ggsave("~/postdoc/edinburgh/scRNAseq/final_revision/beta_combo_cluster_time_wise_grid_vis.png", grid_vis, height=24, width=32, limitsize = FALSE)
# 
# 
# 
# # cluster-wise vdj business ####
first_clono_split <- split(first_clono_data, first_clono_data$seurat_cluster)
first_clonos <- list()
first_clonos$data <- first_clono_split
first_clonos$meta <- as_tibble(data.frame(Sample = names(first_clono_split)))

first_clonos$meta$Timepoint <- substr(first_clonos$meta$Sample, regexpr("_", first_clonos$meta$Sample)+1, nchar(first_clonos$meta$Sample))
first_clonos$meta$Seurat_Cluster <- substr(first_clonos$meta$Sample, 1, regexpr("_", first_clonos$meta$Sample))


third_clono_split <- split(third_clono_data, third_clono_data$seurat_cluster)
third_clonos$data <- third_clono_split
third_clonos$meta <- as_tibble(data.frame(Sample = names(third_clono_split)))

third_clonos$meta$Timepoint <- substr(third_clonos$meta$Sample, regexpr("_", third_clonos$meta$Sample)+1, nchar(third_clonos$meta$Sample))
third_clonos$meta$Seurat_Cluster <- substr(third_clonos$meta$Sample, 1, regexpr("_", third_clonos$meta$Sample))


combo_clonos <- first_clonos

combo_clonos$data <- lapply(1:length(first_clonos$data), function(x) rbind(first_clonos$data[[x]], third_clonos$data[[x]]))
names(combo_clonos$data) <- names(first_clonos$data)

imm_gu <- geneUsage(combo_clonos$data, "hs.trbv", .norm=TRUE, .ambig="exc")

# 
# time_vis <- vis(imm_gu, .by = "Seurat_Cluster", .meta = combo_clonos$meta)
# time_vis2 <- vis(imm_gu, .by = "Seurat_Cluster", .meta = combo_clonos$meta)
# grid_vis <- vis(imm_gu, .grid = T)
# 
# vdj_freqs <- time_vis$data
# write.csv(vdj_freqs, file = "~/postdoc/edinburgh/scRNAseq/final_revision/final_figures/cluster_wise_vdj_freqs.csv", row.names = FALSE)
# 
# ggsave("~/postdoc/edinburgh/scRNAseq/final_revision/beta_combo_cluster_wise_time_vis.png", time_vis, height=6, width=24, limitsize = FALSE)
# ggsave("~/postdoc/edinburgh/scRNAseq/final_revision/beta_combo_cluster_wise_time_vis2.png", time_vis2, height=6, width=24, limitsize = FALSE)
# ggsave("~/postdoc/edinburgh/scRNAseq/final_revision/beta_combo_cluster_wise_grid_vis.png", grid_vis, height=24, width=32, limitsize = FALSE)
# 
# 
# imm_gu_js <- geneUsageAnalysis(imm_gu, .method = "js", .verbose = TRUE)
# imm_gu_cor <- geneUsageAnalysis(imm_gu, .method = "cor", .verbose = TRUE)
# 
# p1 <- vis(imm_gu_js, .title = "Gene usage JS-divergence", .leg.title = "JS", .text.size = 1.5)
# p2 <- vis(imm_gu_cor, .title = "Gene usage correlation", .leg.title = "Cor", .text.size = 1.5)
# 
# p <- p1 + p2
# ggsave("~/postdoc/edinburgh/scRNAseq/final_revision/combo_cluster_wise_heats.png", p, height=8, width=16)
# 
# 
# # multidimensional scaling
# imm_cl_pca <- geneUsageAnalysis(imm_gu, "js+pca+kmeans", .verbose = TRUE)
# imm_cl_mds <- geneUsageAnalysis(imm_gu, "js+mds+kmeans", .verbose = TRUE)
# imm_cl_tsne <- geneUsageAnalysis(imm_gu, "js+tsne+kmeans", .perp = .01, .verbose = TRUE)
# 
# p1 <- vis(imm_cl_pca, .plot = "clust")
# p2 <- vis(imm_cl_mds, .plot = "clust")
# p3 <- vis(imm_cl_tsne, .plot = "clust")
# 
# p <- p1 + p2 + p3
# ggsave("~/postdoc/edinburgh/scRNAseq/final_revision/combo_cluster_wise_mds.png", p, height=8, width=16)

# imm_gu <- geneUsage(combo_clonos$data, "hs.trav", .norm=TRUE, .ambig="exc")

imm_gu_df <- as.data.frame(t(do.call(rbind, imm_gu)))

imm_gu_df <- imm_gu_df %>%
  pivot_longer(cols=colnames(imm_gu_df)[2:14], names_to="seurat_cluster", values_to="perc", values_transform = list("perc" = as.numeric))%>%
  mutate("num_perc"=if_else(is.na(perc), 0, perc))

cluster_sizes <- imm_gu_df %>%
  group_by(seurat_cluster)%>%
  summarise("cluster_size"=sum(num_perc, na.rm = TRUE))

stacked_bar <- imm_gu_df%>%
  ggplot(aes(x=factor(seurat_cluster, levels=0:12), y=num_perc, fill=Names))+
  geom_bar(stat = "identity")+
  theme_minimal()+
  scale_fill_manual(values=colorspace::sequential_hcl(n=47, palette="Lajolla"))

ggsave("~/postdoc/edinburgh/scRNAseq/final_revision/v_alpha_cluster_barstack.png", stacked_bar, height = 4, width=6, bg="white")


point_label <- imm_gu_df%>%
  ggplot(aes(x=factor(seurat_cluster, levels=0:12), y=num_perc, fill=Names, label=Names))+
  geom_point(shape=21)+
  theme_minimal()+
  ggrepel::geom_label_repel()+
  scale_fill_manual(values=colorspace::sequential_hcl(n=47, palette="Lajolla"))

ggsave("~/postdoc/edinburgh/scRNAseq/final_revision/v_alpha_cluster_point_label.png", point_label, height = 6, width=8, bg="white")

line_plot <- imm_gu_df%>%
  ggplot(aes(x=factor(seurat_cluster, levels=0:12), y=num_perc, fill=Names, label=Names))+
  geom_point(shape=21)+
  geom_line(aes(group=Names))+
  facet_wrap(~Names, scales="fixed")+
  theme_minimal()+
  scale_fill_manual(values=colorspace::sequential_hcl(n=47, palette="Lajolla"))+
  theme(legend.position="none",
        axis.text.x = element_blank())

ggsave("~/postdoc/edinburgh/scRNAseq/final_revision/v_alpha_cluster_line_plot.pdf", line_plot, height = 16, width=8, bg="white")




# hill numbers ####
div_hill <- repDiversity(combo_clonos$data, "hill")
p1 <- vis(div_hill, .by = c("Sample"), .meta = combo_clonos$meta, .col="v")
ggsave("~/postdoc/edinburgh/scRNAseq/final_revision/combo_cluster_wise_hill.png", p1, height=6, width=18)

div_gini <- repDiversity(combo_clonos$data, "gini")
gini_df <- data.frame("gini"=div_gini,
                     "seurat_cluster"=as.factor(seq(0, length(div_gini)-1)))

gini_plot <- ggplot(gini_df, aes(x=seurat_cluster, y=gini, fill=seurat_cluster))+
  geom_bar(stat="identity")+
  theme_minimal()+
  scale_fill_manual(values=cluster_palette)+
  theme(axis.text.x = element_blank(),
        legend.position = "none")
ggsave("~/postdoc/edinburgh/scRNAseq/final_revision/combo_cluster_wise_gini.png", gini_plot, height=6, width=6, bg="white")
ggsave("~/postdoc/edinburgh/scRNAseq/final_revision/combo_cluster_wise_gini.pdf", gini_plot, height=6, width=6, bg="white")

#p4 <- vis(div_d50)
# p6 <- vis(div_div)


# p <- p1 + p2
# ggsave("~/postdoc/edinburgh/scRNAseq/final_revision/combo_cluster_wise_diversity2.png", p, height=6, width=18)

# p <- p3 + p6
# ggsave("~/postdoc/edinburgh/scRNAseq/final_revision/combo_cluster_wise_hill.png", p, height=6, width=18)
#
# p3 <- vis(div_hill, .by = c("Seurat_Cluster", "N_Infection"), .meta = Seurat_Cluster$meta)
# p <- p3 + p6
# ggsave("~/postdoc/edinburgh/scRNAseq/final_revision/combo_cluster_wise_n_infection_hill.png", p, height=6, width=18)


# p <- p4
# ggsave("~/postdoc/edinburgh/scRNAseq/final_revision/combo_cluster_wise_diversity4.png", p, height=6, width=9)

# sample-wise vdj business ####

#re-split by sample_id
first_sample_split <- split(first_clono_data, first_clono_data$sample_id)
first_sample_clonos <- list()
first_sample_clonos$data <- first_sample_split
first_sample_clonos$meta <- as_tibble(data.frame(Sample = names(first_sample_split)))

first_sample_clonos$meta$Timepoint <- substr(first_sample_clonos$meta$Sample, regexpr("_", first_sample_clonos$meta$Sample)+1, nchar(first_sample_clonos$meta$Sample))
#first_sample_clonos$meta$Seurat_Cluster <- substr(first_sample_clonos$meta$Sample, 1, regexpr("_", first_sample_clonos$meta$Sample))


third_clono_split <- split(third_clono_data, third_clono_data$sample_id)
third_sample_clonos <- list()
third_sample_clonos$data <- third_clono_split
third_sample_clonos$meta <- as_tibble(data.frame(Sample = names(third_clono_split)))

third_sample_clonos$meta$Timepoint <- substr(third_sample_clonos$meta$Sample, regexpr("_", third_sample_clonos$meta$Sample)+1, nchar(third_sample_clonos$meta$Sample))
#third_sample_clonos$meta$Seurat_Cluster <- substr(third_sample_clonos$meta$Sample, 1, regexpr("_", third_sample_clonos$meta$Sample))

# put the immunarches together

combo_clonos <- first_sample_clonos


combo_clonos$data <- c(first_sample_clonos$data, third_sample_clonos$data)
combo_clonos$meta <- rbind(first_sample_clonos$meta, third_sample_clonos$meta)
names(combo_clonos$data) <- c(names(first_sample_clonos$data), names(third_sample_clonos$data))

imm_gu <- geneUsage(combo_clonos$data, "hs.trbv", .norm = T, .ambig="exc")

time_vis <- vis(imm_gu, .by = "Timepoint" , .meta = combo_clonos$meta)
time_vis2 <- vis(imm_gu, .by = "Timepoint", .meta = combo_clonos$meta, .plot = "box")
grid_vis <- vis(imm_gu, .grid = T)

vdj_freqs <- time_vis$data
write.csv(vdj_freqs, file = "~/postdoc/edinburgh/scRNAseq/final_revision/final_figures/beta_sample_wise_vdj_freqs.csv", row.names = FALSE)


ggsave("~/postdoc/edinburgh/scRNAseq/final_revision/beta_combo_sample_wise_time_vis.png", time_vis, height=6, width=15, limitsize = FALSE)
ggsave("~/postdoc/edinburgh/scRNAseq/final_revision/beta_combo_sample_wise_time_vis2.png", time_vis2, height=15, width=15, limitsize = FALSE)
ggsave("~/postdoc/edinburgh/scRNAseq/final_revision/beta_combo_sample_wise_grid_vis.png", grid_vis, height=24, width=32, limitsize = FALSE)




imm_gu <- geneUsage(combo_clonos$data, "hs.trav", .norm = T, .ambig="exc")

time_vis <- vis(imm_gu, .by = "Timepoint", .meta = combo_clonos$meta)
time_vis2 <- vis(imm_gu, .by = "Timepoint", .meta = combo_clonos$meta, .plot = "box")
grid_vis <- vis(imm_gu, .grid = T)

vdj_freqs <- time_vis$data
write.csv(vdj_freqs, file = "~/postdoc/edinburgh/scRNAseq/final_revision/final_figures/alpha_sample_wise_vdj_freqs.csv", row.names = FALSE)

ggsave("~/postdoc/edinburgh/scRNAseq/final_revision/alpha_combo_sample_wise_time_vis.png", time_vis, height=6, width=15, limitsize = FALSE)
ggsave("~/postdoc/edinburgh/scRNAseq/final_revision/alpha_combo_sample_wise_time_vis2.png", time_vis2, height=15, width=15, limitsize = FALSE)
ggsave("~/postdoc/edinburgh/scRNAseq/final_revision/alpha_combo_sample_wise_grid_vis.png", grid_vis, height=24, width=32, limitsize = FALSE)


div_gini <- repDiversity(combo_clonos$data, "gini")
gini_df <- data.frame("gini"=div_gini, 
                      "sample"=names(combo_clonos$data))

gini_plot <- ggplot(gini_df, aes(x=sample, y=gini, fill=seurat_cluster))+
  geom_bar(stat="identity")+
  theme_minimal()
ggsave("~/postdoc/edinburgh/scRNAseq/final_revision/combo_sample_wise_gini.png", gini_plot, height=6, width=6, bg="white")

imm_gu_df <- as.data.frame(t(do.call(rbind, imm_gu)))

imm_gu_df <- imm_gu_df %>%
  pivot_longer(cols=colnames(imm_gu_df)[2:13], names_to="sample", values_to="perc", values_transform = list("perc" = as.numeric))%>%
  mutate("num_perc"=if_else(is.na(perc), 0, perc))

stacked_bar <- imm_gu_df%>%
  ggplot(aes(x=factor(sample), y=num_perc, fill=Names))+
  geom_bar(stat = "identity")+
  theme_minimal()+
  scale_fill_manual(values=colorspace::sequential_hcl(n=47, palette="Lajolla"))+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90))

ggsave("~/postdoc/edinburgh/scRNAseq/final_revision/v_beta_sample_barstack.png", stacked_bar, height = 4, width=6, bg="white")


point_label <- imm_gu_df%>%
  ggplot(aes(x=factor(sample), y=num_perc, fill=Names, label=Names))+
  geom_point(shape=21)+
  theme_minimal()+
  ggrepel::geom_label_repel()+
  scale_fill_manual(values=colorspace::sequential_hcl(n=47, palette="Lajolla"))

ggsave("~/postdoc/edinburgh/scRNAseq/final_revision/v_beta_sample_point_label.png", point_label, height = 6, width=8, bg="white")

line_plot <- imm_gu_df%>%
  ggplot(aes(x=factor(sample), y=num_perc, fill=Names, label=Names))+
  geom_point(shape=21)+
  geom_line(aes(group=Names))+
  facet_wrap(~Names, scales="fixed")+
  theme_minimal()+
  scale_fill_manual(values=colorspace::sequential_hcl(n=47, palette="Lajolla"))+
  theme(legend.position="none",
        axis.text.x = element_blank())

ggsave("~/postdoc/edinburgh/scRNAseq/final_revision/v_beta_sample_line_plot.png", line_plot, height = 16, width=8, bg="white")


# cluster & n_infection vdj business ####


#combo <- readRDS("/labs/prasj/BIG_Flo/vac63c/final_revision/final_10k.seurat")
combo <- readRDS("/labs/prasj/BIG_Flo/vac63c/final_revision/no_12_harmony.seurat")
first_seurat <- subset(combo, subset = n_infection=="First")

#get cluster idents frome seurat
first_idents <- Idents(first_seurat)
first_clono_data$seurat_cluster <- first_idents[match(first_clono_data$Barcode, names(first_idents))]
first_clono_data$timepoint <- substr(first_clono_data$sample_id, 6, nchar(first_clono_data$sample_id))
first_clono_data$cluster_time <- paste(first_clono_data$seurat_cluster, first_clono_data$timepoint, sep="_")
first_clono_data$n_infection <- "First"
first_clono_data$cluster_n <- paste(first_clono_data$n_infection, first_clono_data$seurat_cluster, sep="_")

#split by cluster time, add metadata
first_clono_split <- split(first_clono_data, first_clono_data$cluster_n)
first_clonos <- list()
first_clonos$data <- first_clono_split
first_clonos$meta <- as_tibble(data.frame(Sample = names(first_clono_split)))

first_clonos$meta$Timepoint <- substr(first_clonos$meta$Sample, regexpr("_", first_clonos$meta$Sample)+1, nchar(first_clonos$meta$Sample))
first_clonos$meta$Seurat_Cluster <- substr(first_clonos$meta$Sample, 1, regexpr("_", first_clonos$meta$Sample))
first_clonos$meta$N_Infection <- "First"

# do it all again for third #

# WE NEED TO REPEAT THE ABOVE THREE STEPS WITH THE SEURAT IDENTS TO CREATE A LIST OF CLUSTER-LEVEL REPERTOIRES

#combo <- readRDS("/labs/prasj/BIG_Flo/vac63c/final_revision/final_10k.seurat")
# third_seurat <- subset(combo, subset = n_infection=="Third")
#get cluster idents frome seurat
third_idents <- Idents(third_seurat)
third_clono_data$seurat_cluster <- third_idents[match(third_clono_data$Barcode, names(third_idents))]
third_clono_data$timepoint <- substr(third_clono_data$sample_id, 6, nchar(third_clono_data$sample_id))
third_clono_data$cluster_time <- paste(third_clono_data$seurat_cluster, third_clono_data$timepoint, sep="_")
third_clono_data$n_infection <- "Third"
third_clono_data$cluster_n <- paste(third_clono_data$n_infection, third_clono_data$seurat_cluster, sep="_")

#split by cluster time, add metadata
third_clono_split <- split(third_clono_data, third_clono_data$cluster_n)
third_clonos <- list()
third_clonos$data <- third_clono_split
third_clonos$meta <- as_tibble(data.frame(Sample = names(third_clono_split)))

third_clonos$meta$Timepoint <- substr(third_clonos$meta$Sample, regexpr("_", third_clonos$meta$Sample)+1, nchar(third_clonos$meta$Sample))
third_clonos$meta$Seurat_Cluster <- substr(third_clonos$meta$Sample, 1, regexpr("_", third_clonos$meta$Sample))
third_clonos$meta$N_Infection <- "Third"

# put the two immunarches together ##

combo_clonos <- first_clonos

#iterate through list, combining each df

combo_clonos$data <- c(first_clonos$data, third_clonos$data)
combo_clonos$meta <- rbind(first_clonos$meta, third_clonos$meta)
names(combo_clonos$data) <- c(names(first_clonos$data), names(third_clonos$data))

# calculate gene usage statistics for human beta chain, exclude ambiguous assignments with more than one chain
imm_gu <- geneUsage(combo_clonos$data, "hs.trbv", .norm=TRUE, .ambig="exc")

imm_gu_js <- geneUsageAnalysis(imm_gu, .method = "js", .verbose = TRUE)
imm_gu_cor <- geneUsageAnalysis(imm_gu, .method = "cor", .verbose = TRUE)

p1 <- vis(imm_gu_js, .title = "Gene usage JS-divergence", .leg.title = "JS", .text.size = 1.5)
p2 <- vis(imm_gu_cor, .title = "Gene usage correlation", .leg.title = "Cor", .text.size = 1.5)

p <- p1 + p2
ggsave("~/postdoc/edinburgh/scRNAseq/final_revision/combo_cluster_n_infection_wise_heats.png", p, height=8, width=16)

time_vis <- vis(imm_gu, .by = "Timepoint", .meta = combo_clonos$meta)
time_vis2 <- vis(imm_gu, .by = "Timepoint", .meta = combo_clonos$meta, .plot = "box")
grid_vis <- vis(imm_gu, .grid = T)

vdj_freqs <- time_vis$data
write.csv(vdj_freqs, file = "~/postdoc/edinburgh/scRNAseq/final_revision/final_figures/cluster_n_infection_wise_vdj_freqs.csv", row.names = FALSE)

ggsave("~/postdoc/edinburgh/scRNAseq/final_revision/beta_combo_cluster_n_infection_vis.png", time_vis, height=6, width=15, limitsize = FALSE)
ggsave("~/postdoc/edinburgh/scRNAseq/final_revision/beta_combo_cluster_n_infection_vis2.png", time_vis2, height=15, width=15, limitsize = FALSE)
ggsave("~/postdoc/edinburgh/scRNAseq/final_revision/beta_combo_ccluster_n_infection_grid_vis.png", grid_vis, height=24, width=32, limitsize = FALSE)



div_gini <- repDiversity(combo_clonos$data, "gini")
gini_df <- data.frame("gini"=div_gini, 
                      "cluster"=substr(names(combo_clonos$data), 7,8),
                      "n_infection"=substr(names(combo_clonos$data), 1,5))

gini_plot <- ggplot(gini_df, aes(x=factor(cluster, levels=seq(0,12)), y=gini, fill=cluster))+
  geom_bar(stat="identity")+
  facet_wrap(~n_infection)+
  theme_minimal()+
  scale_fill_manual(values=rep(cluster_palette, 2))+
  theme(legend.position = "none",
        axis.title.x = element_blank())
  
ggsave("~/postdoc/edinburgh/scRNAseq/final_revision/combo_cluster_n_infection_gini.png", gini_plot, height=4, width=6, bg="white")
ggsave("~/postdoc/edinburgh/scRNAseq/final_revision/combo_cluster_n_infection_gini.pdf", gini_plot, height=4, width=6, bg="white")

#bespoke heatmaps ####

library(ComplexHeatmap)

data <- read.csv("~/postdoc/edinburgh/scRNAseq/scg/vac63c/final_revision/final_figures/cluster_wise_vdj_freqs.csv")
# data <- read.csv("~/postdoc/edinburgh/scRNAseq/scg/vac63c/final_revision/final_figures/cluster_time_wise_vdj_freqs.csv")

data2 <- data %>%
  select(Gene, Freq, Sample)%>%
  mutate(Freq=if_else(is.na(Freq), 0, Freq))%>%
  pivot_wider(names_from = Sample, values_from = Freq, id_cols = Gene)%>%
  select(-Gene)

spearman <- cor(data2, method = "spearman")

baseline_dist <- dist(spearman, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)
baseline_hclust <- hclust(baseline_dist)

#check.names=FALSE here makes sure that the +/- symbols parse and spaces aren't dots
# colnames(spearman) <- gsub(".", " ", colnames(spearman), fixed=T)
# rownames(spearman) <- gsub(".", " ", rownames(spearman), fixed=T)

baseline_spearman_df  <- data.frame(spearman, check.names = FALSE)
colnames(baseline_spearman_df) <- gsub(".", " ", colnames(baseline_spearman_df), fixed=T)
spearman_matrix <- as.matrix(baseline_spearman_df[rownames(baseline_spearman_df)[rev(baseline_hclust$order)],colnames(baseline_spearman_df)[rev(baseline_hclust$order)]])

floor <- min(spearman_matrix)/2
col_fun_pearson <- circlize::colorRamp2(
  c(floor, floor+(1-floor)/2,1),
  c("#0859C6", "black", "#FFA500"))
# col_fun_pearson <- circlize::colorRamp2(c(0, 0.5, 1), c("#0859C6", "black", "#FFA500"))

#doesnt work for pdf..
# colnames(pearson_matrix) <- gsub("IFNy", "IFNγ", colnames(pearson_matrix))
# rownames(pearson_matrix) <-  gsub("IFNy", "IFNγ",rownames(pearson_matrix))


spearman_heatmap <- Heatmap(matrix = spearman_matrix,
                            cluster_rows = TRUE,
                            cluster_columns=TRUE,
                            show_row_dend = FALSE,
                            show_column_dend = TRUE,
                            show_heatmap_legend = TRUE,
                            name = "Spearman Rho",
                            #cluster_columns = FALSE,
                            column_names_gp = gpar(fontsize = 6),
                            row_names_gp = gpar(fontsize = 6),
                            row_names_side = "left",
                            col = col_fun_pearson,
                            column_names_rot = 45)

pdf("~/postdoc/edinburgh/scRNAseq/scg/vac63c/final_revision/final_figures/cluster_vbeta_spearman.pdf", height=4.5, width = 5)
draw(spearman_heatmap, padding=unit(c(2,2,2,2), "mm"))
dev.off()

# 
# pdf("~/postdoc/edinburgh/scRNAseq/scg/vac63c/final_revision/final_figures/cluster_time_vbeta_spearman.pdf", height=4.5, width = 5)
# draw(spearman_heatmap, padding=unit(c(2,2,2,2), "mm"))
# dev.off()








data <- read.csv("~/postdoc/edinburgh/scRNAseq/scg/vac63c/final_revision/final_figures/cluster_time_wise_vdj_freqs.csv")

data2 <- data %>%
  select(Gene, Freq, Sample)%>%
  mutate(Freq=if_else(is.na(Freq), 0, Freq))%>%
  pivot_wider(names_from = Sample, values_from = Freq, id_cols = Gene)%>%
  select(-Gene)

spearman <- cor(data2, method = "spearman")

baseline_dist <- dist(spearman, method = "euclidean", diag = FALSE, upper = FALSE, p = 2)
baseline_hclust <- hclust(baseline_dist)

#check.names=FALSE here makes sure that the +/- symbols parse and spaces aren't dots
# colnames(spearman) <- gsub(".", " ", colnames(spearman), fixed=T)
# rownames(spearman) <- gsub(".", " ", rownames(spearman), fixed=T)

baseline_spearman_df  <- data.frame(spearman, check.names = FALSE)
colnames(baseline_spearman_df) <- gsub(".", " ", colnames(baseline_spearman_df), fixed=T)
spearman_matrix <- as.matrix(baseline_spearman_df[rownames(baseline_spearman_df)[rev(baseline_hclust$order)],colnames(baseline_spearman_df)[rev(baseline_hclust$order)]])

col_fun_pearson <- circlize::colorRamp2(c(min(spearman_matrix), (1-min(spearman_matrix))+min(spearman_matrix)/2,1), c("#0859C6", "black", "#FFA500"))
col_fun_pearson <- circlize::colorRamp2(c(0, (1-min(spearman_matrix))+min(spearman_matrix)/2,1), c("#0859C6", "black", "#FFA500"))

#doesnt work for pdf..
# colnames(pearson_matrix) <- gsub("IFNy", "IFNγ", colnames(pearson_matrix))
# rownames(pearson_matrix) <-  gsub("IFNy", "IFNγ",rownames(pearson_matrix))


spearman_heatmap <- Heatmap(matrix = spearman_matrix,
                            cluster_rows = TRUE,
                            cluster_columns=TRUE,
                            show_row_dend = FALSE,
                            show_column_dend = TRUE,
                            show_heatmap_legend = TRUE,
                            name = "Pearson r",
                            #cluster_columns = FALSE,
                            column_names_gp = gpar(fontsize = 6),
                            row_names_gp = gpar(fontsize = 6),
                            row_names_side = "left",
                            col = col_fun_pearson,
                            column_names_rot = 45)

pdf("~/postdoc/edinburgh/scRNAseq/scg/vac63c/final_revision/final_figures/cluster_vbeta_spearman.pdf", height=4.5, width = 5)
draw(spearman_heatmap, padding=unit(c(2,2,2,2), "mm"))
dev.off()

