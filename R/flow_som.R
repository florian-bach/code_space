library(FlowSOM) 
library(readxl) 
library(flowCore) 
library(premessa)
library(matrixStats)
library(ggplot2)
library(reshape2)
library(ConsensusClusterPlus)
library(dplyr)
library(RColorBrewer)
library(pheatmap)
library(limma)
library(ggrepel)



setwd("/Users/s1249052/PhD/flow_data/vac69a/t_cells_only/better_gating/")

files_list <- list.files(path=".", pattern="*.fcs")

# concatenate_fcs_files(files_list, "concat02.fcs")
# 
# fcs_raw <- list()
# for (i in files_list){
#   print(i)
#   j <- read.FCS(i, transformation = FALSE, truncate_max_range = FALSE)
#   assign(paste(substr(i, nchar(i)-10, nchar(i)-4)), j)
#   fcs_raw[[paste(substr(i, nchar(i)-10, nchar(i)-4))]] <- j}
# 
# head(fcs_raw)
# 
df_lol <- as.matrix(fcs_raw[[1]]@parameters@data[["desc"]])
channel_names <- df_lol[1:69]

# timecourse for only Volunteer 03
flo_set <- read.flowSet(files_list[21:25], transformation = FALSE, truncate_max_range = FALSE)
name_frame <- as.matrix(flo_set[[1]]@parameters@data[["desc"]])

channel_names <- substr(name_frame, 7, nchar(name_frame))


## arcsinh transformation and column subsetting
fcs <- fsApply(flo_set, function(x, cofactor = 5){
  colnames(x) <- channel_names
  expr <- exprs(x)
  expr <- asinh(expr[, channel_names[c(3, 14:15, 23:57, 63, 65)]] / cofactor)
  exprs(x) <- expr
  x
})

fcs

expr <- fsApply(fcs, exprs)
dim(expr)


rng <- colQuantiles(expr, probs = c(0.01, 0.99))
expr01 <- t((t(expr) - rng[, 1]) / (rng[, 2] - rng[, 1]))
expr01[expr01 < 0] <- 0
expr01[expr01 > 1] <- 1

md <- read_excel("metadata.xlsx", col_names = T)
sample_ids <- rep(md$file_name[21:25], fsApply(flo_set, nrow))

sample_ids <- substr(sample_ids, nchar(sample_ids)-9, nchar(sample_ids)-7)
time_point <- gsub("_1", "C+1", sample_ids)



ggdf <- data.frame(time_point = time_point, expr)
ggdf_long <- gather(ggdf, marker, expression, colnames(ggdf)[2:41])


ggplot(ggdf_long, aes(x = expression, color = time_point, 
                 group = time_point)) +
  geom_density() +
  facet_wrap(~ marker, scales = "free") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1), 
        strip.text = element_text(size = 7), axis.text = element_text(size = 5)) +
  guides(color = guide_legend(ncol = 1))
  





expr_median_sample_tbl <- data.frame(sample_id = sample_ids, expr) %>%
  group_by(sample_id) %>% 
  summarize_all(funs(median))

expr_median_sample <- t(expr_median_sample_tbl[, -1])
colnames(expr_median_sample) <- expr_median_sample_tbl$sample_id


mds <- plotMDS(expr_median_sample, plot = FALSE)

ggdf <- data.frame(MDS1 = mds$x, MDS2 = mds$y, 
                   sample_id = colnames(expr_median_sample))

ggdf$volunteer <- as.character(ggdf$sample_id)
ggdf$volunteer <- substr(ggdf$volunteer, nchar(ggdf$volunteer)-1, nchar(ggdf$volunteer))

mm <- match(ggdf$sample_id, md$file_name)
ggdf$volunteer <- md$volunteer[21:25]

ggplot(ggdf, aes(x = MDS1, y = MDS2, color = volunteer)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_label_repel(aes(label = sample_id)) +
  theme_bw()







fsom <- ReadInput(fcs, transform = FALSE, scale = FALSE)
set.seed(1234)
som <- BuildSOM(fsom, colsToUse = fcs@colnames[c(3:31, 33:40)])


#### clustering on som

codes <- som$map$codes
plot_outdir <- "consensus_plots"
nmc <- 42
mc <- ConsensusClusterPlus(t(codes), maxK = nmc, reps = 100, 
                           pItem = 0.9, pFeature = 1, title = plot_outdir, plot = "png", 
                           clusterAlg = "hc", innerLinkage = "average", finalLinkage = "average", 
                           distance = "euclidean", seed = 1234)





## Get cluster ids for each cell
code_clustering1 <- mc[[nmc]]$consensusClass
cell_clustering1 <- code_clustering1[som$map$mapping[,1]]

for (i in seq(1, 5)){
  assign(paste("d", i, sep='_'), som$metaData[[i]][[2]]-som$metaData[[i]][[1]]+1)
}


time_points <- rep(c("C-1", "C+8", "C+10", "DoD", "T+6"), times=c(d_1, d_2, d_3, d_4, d_5))

dr <- data.frame(cluster_id = code_clustering1[som$map$mapping[,1]], time_point = time_points, 
                 expr[, marker_levels])

mm <- match(dr$sample_id, md$sample_id)
dr$condition <- md$condition[mm]



color_clusters <- rep(c("#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72", 
                    "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3", 
                    "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D", 
                    "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999", 
                    "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", "#808000", 
                    "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00"), times=2)




plot_clustering_heatmap_wrapper <- function(expr, expr01, 
                                            cell_clustering, color_clusters, cluster_merging = NULL){
  
  # Calculate the median expression
  expr_median <- data.frame(expr, cell_clustering = cell_clustering1) %>%
    group_by(cell_clustering) %>% 
    summarize_all(funs(median))
  colnames(expr_median) <- c("cell_clustering", colnames(expr))  
  expr01_median <- data.frame(expr01, cell_clustering = cell_clustering1) %>%
    group_by(cell_clustering) %>% 
    summarize_all(funs(median))
  colnames(expr01_median) <- c("cell_clustering", colnames(expr01)) 
  # Calculate cluster frequencies
  clustering_table <- as.numeric(table(cell_clustering1))
  
  # This clustering is based on the markers that were used for the main clustering
  d <- dist(expr_median[, colnames(expr)], method = "euclidean")
  cluster_rows <- hclust(d, method = "average")
  
  expr_heat <- as.matrix(expr01_median[, colnames(expr01)])
  rownames(expr_heat) <- expr01_median$cell_clustering
  
  labels_row <- paste0(rownames(expr_heat), " (", 
                       round(clustering_table / sum(clustering_table) * 100, 2), "%)")
  labels_col <- colnames(expr_heat)
  
  # Row annotation for the heatmap
  annotation_row <- data.frame(cluster = factor(expr01_median$cell_clustering))
  rownames(annotation_row) <- rownames(expr_heat)
  
  color_clusters <- color_clusters[1:nlevels(annotation_row$cluster)]
  names(color_clusters) <- levels(annotation_row$cluster)
  annotation_colors <- list(cluster = color_clusters)
  annotation_legend <- FALSE
  
  if(!is.null(cluster_merging)){
    cluster_merging$new_cluster <- factor(cluster_merging$new_cluster)
    annotation_row$cluster_merging <- cluster_merging$new_cluster 
    color_clusters <- color_clusters[1:nlevels(cluster_merging$new_cluster)]
    names(color_clusters) <- levels(cluster_merging$new_cluster)
    annotation_colors$cluster_merging <- color_clusters
    annotation_legend <- TRUE
  }
  
  # Colors for the heatmap
  color <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(100)
  
  pheatmap(expr_heat, color = color, 
           cluster_cols = FALSE, cluster_rows = cluster_rows, 
           labels_col = labels_col, labels_row = labels_row, 
           display_numbers = TRUE, number_color = "black", 
           fontsize = 8, fontsize_number = 4,
           annotation_row = annotation_row, annotation_colors = annotation_colors, 
           annotation_legend = annotation_legend)
  
}

marker_levels <- c("CD4",
                   "CD8",
                   "Vd2",
                   "TCRgd",
                   "Va7.2",
                   "CD38",
                   "HLA-DR",
                   "ICOS",
                   "CD28",
                   "PD1",
                   "TIM-3",
                   "CD95",
                   "BCL-2",
                   "CD27",
                   "Perforin",
                   "GZB",
                   "Tbet",
                   "Ki-67",
                   "CD127",
                   "IntegrinB7",
                   "CD56",
                   "CD16",
                   "CD161",
                   "CD49d",
                   "CD103",
                   "CD25",
                   "FoxP3",
                   "CD39",
                   "CTLA4",
                   "CLA",
                   "CXCR5",
                   "CD57",
                   "CD45RA",
                   "CD45RO",
                   "CCR7")



plot_clustering_heatmap_wrapper(expr = expr[,marker_levels],
                                expr01 = expr01[, marker_levels],
                                cell_clustering = cell_clustering1, color_clusters = color_clusters)