library(ggplot2)
library(tidyr)
library(dplyr)
library(ComplexHeatmap)



# # expression per cluster ####
# shplit_cells <- function(x, by) {
#   stopifnot(is.character(by), by %in% colnames(colData(x)))
#   cd <- data.frame(colData(x))
#   cd$cluster_id <- cluster_ids(x, k="flo_merge")
#   dt <- data.table::data.table(cd, i = seq_len(ncol(x)))
#   dt_split <- split(dt, by = by, sorted = TRUE, flatten = FALSE)
#   map_depth(dt_split, length(by), "i")
# }
# 
# 
# ahgg <- function(x, by, fun = c("median", "mean", "sum")) {
#   fun <- switch(match.arg(fun), 
#                 median = rowMedians, mean = rowMeans, sum = rowSums)
#   cs <- shplit_cells(x, by)
#   pb <- map_depth(cs, -1, function(i) {
#     if (length(i) == 0) return(numeric(nrow(x)))
#     fun(assay(x, "exprs")[, i, drop = FALSE])
#   })
#   map_depth(pb, -2, function(u) as.matrix(data.frame(
#     u, row.names = rownames(x), check.names = FALSE)))
# }
# 
# 
# # make maxtrix of median expression ####
# merged_daf <- read_full("/home/flobuntu/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/")
# 
# 
# #get rid of resting Vd cluster

# 
# meta_up_t6 <- filterSCE(merged_daf, k = "flo_merge", cluster_id %in% sig_clusters)
# 
# clustering_markers <- read.csv("/home/flobuntu/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/refined_markers.csv", stringsAsFactors = F)
# #merged_daf$cluster_id <- cluster_ids(merged_daf, k="flo_merge")
# 
# ms <- ahgg(merged_daf[clustering_markers$refined_markers,], by = c("cluster_id", "timepoint"), fun="median")
# 
# ms2 <- lapply(ms, function(x)data.frame(x))
# ms2 <- Map(cbind, ms2, cluster_id = names(ms))
# 
# ms3 <- data.table::rbindlist(ms2)
# ms3[,Marker := unlist(lapply(ms2, rownames))]
# 
# ms3 <- select(ms3, T6, cluster_id, Marker)
# ms4 <- as.matrix(tidyr::spread(ms3, cluster_id, T6))
# 
# rownames(ms4) <- ms4[,1]
# ms5 <- ms4[,2:ncol(ms4)]
# 
# scaled_mat <- apply(apply(ms5, c(1,2), as.numeric), MARGIN = 1, function(x)scales::rescale(x, to=c(0, 1)))
# 
# sig_scaled_mat <- scaled_mat[rownames(scaled_mat)%in%sig_clusters,]
# 
# reordered_sig_scaled_mat <- sig_scaled_mat[,match(clustering_markers$refined_markers, colnames(sig_scaled_mat))]
# 

inferno <- colorspace::sequential_hcl("inferno", n=8)
col_inferno <- circlize::colorRamp2(seq(0,1, by=1/{length(inferno)-1}), inferno)

reordered_sig_scaled_mat <- as.matrix(read.csv("/home/flobuntu/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/cluster_medians_heatmap_t6.csv", header = T, row.names = 1))

sig_clusters <- read.csv("/home/flobuntu/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/sig_t6_clusters.csv", header = TRUE, stringsAsFactors = FALSE)
sig_clusters <- sig_clusters <- sig_clusters[-9,2]

all_t6_data <- read.csv("/home/flobuntu/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/all_t6_data.csv", stringsAsFactors = FALSE, header = T)

t6_map_data <- all_t6_data%>%
  group_by(cluster_id, timepoint) %>%
  mutate(mean_freq=mean(frequency)) %>%
  ungroup() %>%
  select(cluster_id, mean_freq, timepoint)
  
t6_map_data <- filter(t6_map_data, timepoint=="T6")
t6_map_data <- t6_map_data[!duplicated(t6_map_data), ]


t6_map_data$cluster_id <- ifelse(substr(t6_map_data$cluster_id, 1, 1)==" ", substr(t6_map_data$cluster_id, 2, nchar(t6_map_data$cluster_id)), t6_map_data$cluster_id)
t6_map_data$cluster_id <- ifelse(substr(t6_map_data$cluster_id, 1, 1)==" ", substr(t6_map_data$cluster_id, 2, nchar(t6_map_data$cluster_id)), t6_map_data$cluster_id)

t6_map_data$cluster_id <- ifelse(substr(t6_map_data$cluster_id, nchar(t6_map_data$cluster_id), nchar(t6_map_data$cluster_id))==" ", substr(t6_map_data$cluster_id, 1, nchar(t6_map_data$cluster_id)-1), t6_map_data$cluster_id)
t6_map_data$cluster_id <- ifelse(substr(t6_map_data$cluster_id, nchar(t6_map_data$cluster_id), nchar(t6_map_data$cluster_id))==" ", substr(t6_map_data$cluster_id, 1, nchar(t6_map_data$cluster_id)-1), t6_map_data$cluster_id)



t6_map_data <- filter(t6_map_data, cluster_id %in% sig_clusters)
t6_map_data <- t6_map_data[order(t6_map_data$mean_freq, decreasing = T),]

rereordered_sig_scaled_mat <- reordered_sig_scaled_mat[match(t6_map_data$cluster_id, rownames(reordered_sig_scaled_mat)),]


cd3_right_anno <- rowAnnotation(gap = unit(2, "mm"),
                                "Mean Frequency at T6" = anno_barplot(t6_map_data$mean_freq, which="row", axis = TRUE, ylim = c(0, 4)),
                                width = unit(4, "cm")
)
                                


# only run this portion when making frequency boxplots ####

box_map_data <- all_t6_data%>%
  group_by(cluster_id, timepoint) %>%
  mutate(mean_freq=mean(frequency)) %>%
  ungroup() %>%
  #select(cluster_id, mean_freq, timepoint)
  select(cluster_id, frequency, timepoint, volunteer)

box_map_data <- filter(box_map_data, timepoint=="T6")

box_map_data$cluster_id <- ifelse(substr(box_map_data$cluster_id, 1, 1)==" ", substr(box_map_data$cluster_id, 2, nchar(box_map_data$cluster_id)), box_map_data$cluster_id)
box_map_data$cluster_id <- ifelse(substr(box_map_data$cluster_id, 1, 1)==" ", substr(box_map_data$cluster_id, 2, nchar(box_map_data$cluster_id)), box_map_data$cluster_id)

box_map_data$cluster_id <- ifelse(substr(box_map_data$cluster_id, nchar(box_map_data$cluster_id), nchar(box_map_data$cluster_id))==" ", substr(box_map_data$cluster_id, 1, nchar(box_map_data$cluster_id)-1), box_map_data$cluster_id)
box_map_data$cluster_id <- ifelse(substr(box_map_data$cluster_id, nchar(box_map_data$cluster_id), nchar(box_map_data$cluster_id))==" ", substr(box_map_data$cluster_id, 1, nchar(box_map_data$cluster_id)-1), box_map_data$cluster_id)


box_map_data <- filter(box_map_data, cluster_id %in% sig_clusters)
wide_t6_map_data <- spread(box_map_data, volunteer, frequency)
slim_wide_t6_map_data <- as.matrix(select(wide_t6_map_data, -timepoint, -cluster_id))
rownames(slim_wide_t6_map_data) <- wide_t6_map_data$cluster_id

slim_wide_t6_map_data <- slim_wide_t6_map_data[match(t6_map_data$cluster_id, rownames(reordered_sig_scaled_mat)),]

lineage_palette <- c("#AA3377", "#EE6677", "#4477AA", "#66CCEE", "#228833", "#FFA500", "#BBBBBB")
names(lineage_palette) <-c("CD4", "Treg", "CD8", "MAIT", "gd", "DN", "Resting")

short_lineage_palette <- lineage_palette[c(1,3,4,5)]


colcsv <- read.csv("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/cluster_palette.csv", header=T, stringsAsFactors = F)

col_pal <- colcsv$x
names(col_pal) <- colcsv$X

colz <- unname(col_pal[match(rownames(rereordered_sig_scaled_mat), names(col_pal))])
breakz <- names(col_pal[match(rownames(rereordered_sig_scaled_mat), names(col_pal))])

cd3_right_anno_var <- rowAnnotation(gap = unit(2, "mm"),
                                    # "Mean Frequency at T6" = anno_boxplot(slim_wide_t6_map_data, which="row", axis = TRUE, gp=gpar(fill=c("#228833","#AA3377", "#AA3377", "#66CCEE", "#AA3377", "#AA3377", "#AA3377", "#4477AA", "#4477AA" ))),
                                    "Mean Percentage\nof CD3+ at T6" = anno_barplot(t6_map_data$mean_freq, which="row", axis = TRUE, ylim = c(0, 4), gp=gpar(fill=colz)),

                                    width = unit(4, "cm"),
                                    annotation_legend_param = list("Mean Percentage\nof CD3+ at T6" = list(title ="Lineage",
                                                                                                 at = breakz),
                                                                                                 legend_gp = gpar(fill = colz),
                                                                                                 title_position = "topleft")
                                    )
                                                                   



draw(cd3_right_anno_var)


box_lgd <- Legend(labels =  breakz,
                  title = "Cluster_ID",
                  type = "grid",
                  legend_gp = gpar(fill = colz)
)


  
  median_cluster_heat <- Heatmap(matrix = rereordered_sig_scaled_mat,
                       cluster_rows = FALSE,
                       name = "Median Marker Expression",
                       cluster_columns = FALSE,
                       row_names_side = "left",
                       col = col_inferno,
                       rect_gp = gpar(col = "white"),
                       #top_annotation = combo_top_anno,
                       right_annotation = cd3_right_anno_var,
                       show_heatmap_legend = FALSE,
                       column_names_rot = 60,
                       #heatmap_legend_param = list(col = col_fun4, title = "Normalised Frequency", title_position = "topleft"),
                       width = unit(16, "cm"),
                       height = unit(16*9/28, "cm")
  )
  
  
  png("/home/flobuntu/PhD/cytof/vac69a/figures_for_paper/diffcyt/edgeR/sig_cluster_t6_phenotype_heat_var.png", width=13, height=4, units = "in", res=400)
  draw(median_cluster_heat,
       #annotation_legend_list = list(box_lgd),
       merge_legends = FALSE,
       #padding = unit(c(200, 200, 200, 200), "mm")
  )
  dev.off()

# pie chart of relative frequency of activated clusters

pie_data <- data.frame(t6_map_data)


pie_data$lineage <- as.factor(c("gd", "CD4", "CD4", "MAIT", "CD4", "CD4", "CD4", "CD8", "CD8"))
pie_data <- pie_data[order(pie_data$lineage),]

label_pos <- pie_data$mean_freq[1]/2

for(i in 2:nrow(pie_data)){
  new_entry <- sum(sum(pie_data$mean_freq[seq(1, i-1)]), pie_data$mean_freq[i]/2)
  label_pos <- c(label_pos, new_entry)
}

pie_data$label_position <- label_pos



floor <- 0

for(i in 2:nrow(pie_data)){
  new_entry <- sum(pie_data$mean_freq[seq(1, i-1)])
  floor <- c(floor, new_entry)
}

pie_data$floor <- floor


ceiling <- pie_data$mean_freq[1]

for(i in 2:(nrow(pie_data))+1){
  new_entry <- sum(pie_data$mean_freq[seq(1, i-1)])
  ceiling <- c(ceiling, new_entry)
}



pie_data$ceiling <- ceiling



# pie_palette <- c(lineage_palette[c(1,3,4,5)], colorspace::qualitative_hcl("dark3", n=9))
# names(pie_palette)[5:length(pie_palette)] <- pie_data$cluster_id




myAng <-  seq(-30,-360,length.out = 9)
  
ggplot(pie_data)+
  geom_rect(aes(xmin=1, xmax=3, ymin=floor, ymax=ceiling, fill=lineage))+
  geom_rect(aes(xmin=3.2, xmax=3.5, ymin=floor, ymax=ceiling, fill=cluster_id), colour="black")+

  #geom_col(aes(x=1, y=mean_freq, fill=factor(cluster_id, levels=cluster_id)), color="black", )+
  #geom_bar(stat="identity", position = "stack", width = 1, aes(x=1, y=0.1502437))+
  #geom_bar(stat="identity", position = "stack", width=3, aes(x=3, y=mean_freq/sum(mean_freq), fill=lineage))+
  #geom_bar(stat="identity", position = "stack", width=2, color="grey", aes(x=5, y=mean_freq/13.52193, fill=lineage))+
  
  #geom_bar(stat="identity", position = "stack", width=1,  aes(x=4.5, y=mean_freq/sum(mean_freq), fill=lineage))+
  #geom_bar(stat="identity", position = "stack", width=1,  aes(x=1, y=mean_freq/sum(mean_freq)), fill="white", color="white")+
  
  #geom_text(aes(label = cluster_id, x=5, y=label_position, angle=360), size = 3.5, hjust = 0)+
  circlize::circos.text()+
  scale_fill_manual(values=pie_palette)+
  coord_polar(theta = "y", start = 0)+
    theme_minimal()+
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          panel.grid = element_blank(),
          legend.title = element_blank()
          )
  
