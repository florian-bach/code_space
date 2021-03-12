library(CATALYST)
library(diffcyt)
# library(dplyr)
# library(tidyr)
#library(ggplot2)

`%notin%` <- Negate(`%in%`)


setwd("~/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only")
#setwd("~/PhD/cytof/vac63c/normalised_renamed_comped/debarcoded/")

fcs <- list.files(pattern = "fcs")
#fcs <- subset(fcs, !grepl(pattern = "C45", fcs))
#fcs <- subset(fcs, grepl(pattern = "Baseline", fcs))
#fcs <- subset(fcs, !grepl(pattern = "DoD", fcs))
fcs <- subset(fcs, !grepl(pattern = "ctrl", fcs))
fcs <- subset(fcs, !grepl(pattern = "307", fcs))
fcs <- subset(fcs, !grepl(pattern = "302", fcs))


vac63_flowset <- flowCore::read.flowSet(fcs)

# how md was made 
# md <- data.frame("file_name"=fcs,
#                  "volunteer"=paste("v", substr(fcs, 1,3), sep=''),
#                  "timepoint"=substr(fcs, 5,7),
#                  "batch"=rep(cytonorm_metadata$Batch[4:14], each=3)
# )
# 
# md$timepoint <- ifelse(md$timepoint=="Bas", "Baseline", md$timepoint)
# md$timepoint <- gsub("_", "", md$timepoint, fixed=T)
# 
# md$sample_id <- paste(md$volunteer, "_", md$timepoint, sep='')
# 
# write.csv(md, "vac63c_metadata.csv", row.names = FALSE)




#only run this when exclusively looking at ctrl samples, they're not in the metadata file so you need something else

# md <- cbind("file_name"=fcs, "batch"=c("batch_1", "batch_2", "batch_3"), "sample_id"=c("batch_1", "batch_2", "batch_3"))
# 
# panel <- read.csv("~/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/vac63c_panel.csv", header=T)
# 
# sce <- prepData(vac63_flowset, panel, md, md_cols =
#                   list(file = "file_name", id = "sample_id", factors = "batch"))
# 
# 



md <- read.csv("vac63c_metadata.csv", header = T)

md <- subset(md, md$file_name %in% fcs)
md <- md[order(md$timepoint),]


#  
panel <- read.csv("~/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/vac63c_panel.csv", header=T)

sce <- prepData(vac63_flowset, panel, md, md_cols =
                  list(file = "file_name", id = "sample_id", factors = c("volunteer", "timepoint", "n_infection", "batch")))



# p <- plotExprs(sce, color_by = "batch")
# p$facet$params$ncol <- 7                   
# ggplot2::ggsave("marker_expression.png", p, width=12, height=12)                                      


refined_markers <- c("CD4",
                     "CD8",
                     "Vd2",
                     "Va72",
                     #"CXCR5",
                     "CD38",
                     #"CD69",
                     "HLADR",
                     "ICOS",
                     "CD28",
                     "PD1",
                     #"TIM3",
                     "CD95",
                     "BCL2",
                     "CD27",
                     "Perforin",
                     "GZB",
                     #"TCRgd",
                     "Tbet",
                     "Eomes",
                     #"RORgt",
                     #"GATA3",
                     "CTLA4",
                     "Ki67",
                     "CD127",
                     "CD56",
                     #"CD16",
                     "CD161",
                     "CD49d",
                     "CD25",
                     "FoxP3",
                     "CD39",
                     "CX3CR1",
                     "CD57",
                     "CD45RA",
                     "CD45RO",
                     "CCR7")

# all_markers <- c("CD45", "CD3", "CD3", "CD14", "CD16", refined_markers)

set.seed(123);sce <- CATALYST::cluster(sce, features = refined_markers, xdim = 10, ydim = 10, maxK = 50)


# vac63c_control_tcell_cluster_heatmap <- plotExprHeatmap(x = sce,
#                                                         by = "cluster",
#                                                         row_clust = FALSE,
#                                                         col_clust = FALSE,
#                                                         #m = "flo_merge",
#                                                         k= "meta50",
#                                                         bars = TRUE,
#                                                         features = refined_markers)
# 
# pdf("./figures/vac63c_tcell_cluster_heatmap_meta50_flo_names.pdf", height = 10, width = 9)
# vac63c_control_tcell_cluster_heatmap
# dev.off()
# 
# 
# ki67_cd38_plot <- plotScatter(sce, chs = c("Ki67", "CD38"))
# cxcr5_cd4_plot <- plotScatter(sce, chs = c("CXCR5", "CD4"))
# cd56_cd161_plot <- plotScatter(sce, chs = c("CD56", "CD161"))
# 
# 
# pdf("./figures/ki67_cd38_plot.pdf", height = 4, width = 4)
# ki67_cd38_plot
# dev.off()
# 
# 
# pdf("./figures/cxcr5_cd4_plot.pdf", height = 4, width = 4)
# cxcr5_cd4_plot
# dev.off()
# 
# 
# pdf("./figures/cd56_cd161_plot.pdf", height = 4, width = 4)
# cd56_cd161_plot
# dev.off()
# 



meta45_table <- read.csv("/home/flobuntu/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/prim_ter_50_merge.csv")

sce <- CATALYST::mergeClusters(sce, k = "meta50", table = meta45_table, id = "flo_merge", overwrite = TRUE)

t6_edger <- read.csv("~/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/differential_abundance/edgeR/t6_edgeR.csv", header = TRUE, stringsAsFactors = FALSE)
sig_t6_clusters <- subset(t6_edger, t6_edger$p_adj<0.05 & abs(t6_edger$logFC)>1)$cluster_id


primaries <- filterSCE(sce, volunteer %in% c("v313", "v315", "v320"))
tertiaries <- filterSCE(sce, volunteer %in% c("v301", "v304", "v305", "v306", "v308", "v310"))

prim_ei <- metadata(primaries)$experiment_info
prim_design <- createDesignMatrix(prim_ei, c("timepoint", "volunteer"))
colnames(prim_design)

prim_dod_contrast <- createContrast(c(c(0,0,1,0), rep(0, 2)))
ter_dod_contrast <- createContrast(c(c(0,0,1,0), rep(0, 5)))


ter_ei <- metadata(tertiaries)$experiment_info
ter_design <- createDesignMatrix(ter_ei, c("timepoint", "volunteer"))
colnames(ter_design)


prim_t6_contrast <- createContrast(c(c(0,0,0,1), rep(0, 2)))
ter_t6_contrast <- createContrast(c(c(0,0,0,1), rep(0, 5)))



da_dod_prim <- diffcyt(primaries,
                      design = prim_design,
                      contrast = prim_dod_contrast,
                      analysis_type = "DA",
                      method_DA = "diffcyt-DA-edgeR",
                      clustering_to_use = "flo_merge",
                      verbose = T)

da_dod_ter <- diffcyt(tertiaries,
                      design = ter_design,
                      contrast =ter_dod_contrast,
                      analysis_type = "DA",
                      method_DA = "diffcyt-DA-edgeR",
                      clustering_to_use = "flo_merge",
                      verbose = T)

table(rowData(da_dod_prim$res)$p_adj < 0.05)
# FALSE  TRUE 
# 48     1 

table(rowData(da_dod_ter$res)$p_adj < 0.05)
# FALSE  TRUE 
# 18    31 


prim_dod_df <- data.frame(rowData(da_dod_prim$res))
prim_dod_df <- subset(prim_dod_df, prim_dod_df$p_adj<0.05 & abs(prim_dod_df$logFC)>1)


ter_dod_df <- data.frame(rowData(da_dod_ter$res))
ter_dod_df <- subset(ter_dod_df, ter_dod_df$p_adj<0.05 & abs(ter_dod_df$logFC)>1)




dod_prim_diffy <- vac69a.cytof::vac63_diffcyt_boxplot(da_dod_prim, sce, FDR = 0.05, logFC = log2(2))
dod_ter_diffy <- vac69a.cytof::vac63_diffcyt_boxplot(da_dod_ter, sce, FDR = 0.05, logFC = log2(2))

diffy <- dod_ter_diffy

sig_cluster_boxplot_data <- diffy$data

sig_cluster_boxplot_data$batch <- md$batch[match(sig_cluster_boxplot_data$sample_id, md$sample_id)]
sig_cluster_boxplot_data$n_infection <- md$n_infection[match(sig_cluster_boxplot_data$sample_id, md$sample_id)]
sig_cluster_boxplot_data$direction <- ifelse(sig_cluster_boxplot_data)

library(ggplot2)

time_col <- colorspace::sequential_hcl(5, palette = "Purple Yellow")



sig_t6_all_plot <- ggplot(sig_cluster_boxplot_data, aes(x=factor(timepoint, levels=c("Baseline", "DoD", "T6", "C45")), y=frequency))+
  geom_boxplot(aes(fill=n_infection))+
  geom_point(aes(group=n_infection, colour=volunteer), position = position_dodge(width = 0.75))+
  facet_wrap(~cluster_id, scales = "free", ncol = 5, labeller = label_wrap_gen(width = 8))+
  theme_minimal()+
  scale_y_continuous()+
  ylab("% of all CD3+ cells")+
  # scale_fill_manual(values = c("Baseline"=time_col[4],
  #                              "DoD"=time_col[2],
  #                              "T6"=time_col[1],
  #                              "C45"=time_col[5]))+
  scale_fill_manual(values = c("First" = time_col[4],
                               "Third" = time_col[1]))+    
  # scale_fill_manual(values = c("First"="red",
  #                                "Second"="darkblue",
  #                                "Third"="darkgreen"))+
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        axis.title.x = element_blank(),
        strip.text = element_text(size=7))

ggsave("./figures/sig_dod_ter_boxplot.png", sig_t6_all_plot, height=7, width=11)
#ggsave("./figures/sig_dod_boxplot.png", sig_t6_all_plot, height=7.5, width=12)
#ggsave("./figures/sig_c45_boxplot.png", sig_t6_all_plot, height=7.5, width=12)







da_t6_prim <- diffcyt(primaries,
                     design = prim_design,
                     contrast = prim_t6_contrast,
                     #formula = glm_formula,
                     #contrast = t6_glm_contrast,
                     analysis_type = "DA",
                     method_DA = "diffcyt-DA-edgeR",
                     clustering_to_use = "flo_merge",
                     min_cells = 0,
                     verbose = T)


da_t6_ter <- diffcyt(tertiaries,
                      design = ter_design,
                      contrast = ter_t6_contrast,
                      #formula = glm_formula,
                      #contrast = t6_glm_contrast,
                      analysis_type = "DA",
                      method_DA = "diffcyt-DA-edgeR",
                      clustering_to_use = "flo_merge",
                      verbose = T,min_cells = 0
                      )
#View(data.frame(rowData(all_da_t6$res)))


table(rowData(da_t6_prim$res)$p_adj < 0.05)
# FALSE  TRUE 
# 21    26 
table(rowData(da_t6_ter$res)$p_adj < 0.05)
# FALSE  TRUE 
# 33    13 

prim_t6_df <- data.frame(rowData(da_t6_prim$res))


ter_t6_df <- data.frame(rowData(da_t6_ter$res))

write.csv(ter_t6_df, "./differential_abundance/edgeR/ter_t6_df_edger.csv", row.names = FALSE)
write.csv(prim_t6_df, "./differential_abundance/edgeR/prim_t6_df_edger.csv", row.names = FALSE)

sig_prim_t6_df <- subset(prim_t6_df, prim_t6_df$p_adj<0.05 & prim_t6_df$logFC>1)
sig_ter_t6_df <- subset(ter_t6_df, ter_t6_df$p_adj<0.05 & ter_t6_df$logFC>1)

write.table(x = sig_prim_t6_df$cluster_id, "/home/flobuntu/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/prim_ter_sig_clusters.csv", row.names = FALSE, col.names = FALSE)

# sig_t6_clusters[sig_t6_clusters %notin% sig_ter_t6_df$cluster_id]
# [1] "activated CD27- ICOS+HLADR+ EM CD4" "activated HLADR- EM CD4"
# all(sig_t6_clusters %in% sig_prim_t6_df$cluster_id)
# [1] TRUE (17)


t6_prim_diffy <- vac69a.cytof::vac63_diffcyt_boxplot(da_t6_prim, sce, FDR = 0.05, logFC = log2(2))
t6_ter_diffy <- vac69a.cytof::vac63_diffcyt_boxplot(da_t6_ter, sce, FDR = 0.05, logFC = log2(2))

diffy <- t6_prim_diffy

sig_cluster_boxplot_data <- diffy$data

sig_cluster_boxplot_data$batch <- md$batch[match(sig_cluster_boxplot_data$sample_id, md$sample_id)]
sig_cluster_boxplot_data$n_infection <- md$n_infection[match(sig_cluster_boxplot_data$sample_id, md$sample_id)]
sig_cluster_boxplot_data$direction <- ifelse(sig_cluster_boxplot_data)



sig_t6_all_plot <- ggplot(sig_cluster_boxplot_data, aes(x=factor(timepoint, levels=c("Baseline", "DoD", "T6", "C45")), y=frequency))+
  geom_boxplot(aes(fill=n_infection))+
  geom_point(aes(group=n_infection, colour=volunteer), position = position_dodge(width = 0.75))+
  facet_wrap(~cluster_id, scales = "free", ncol = 5, labeller = label_wrap_gen(width = 8))+
  theme_minimal()+
  scale_y_continuous()+
  ylab("% of all CD3+ cells")+
  # scale_fill_manual(values = c("Baseline"=time_col[4],
  #                              "t6"=time_col[2],
  #                              "T6"=time_col[1],
  #                              "C45"=time_col[5]))+
  scale_fill_manual(values = c("First" = time_col[4],
                               "Third" = time_col[1]))+    
  # scale_fill_manual(values = c("First"="red",
  #                                "Second"="darkblue",
  #                                "Third"="darkgreen"))+
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        axis.title.x = element_blank(),
        strip.text = element_text(size=7))

ggsave("./figures/sig_t6_prim_boxplot.png", sig_t6_all_plot, height=13, width=11)




t6_prim_diffy <- vac69a.cytof::vac63_diffcyt_boxplot(da_t6_prim, sce, FDR =1, logFC = log2(1))
t6_ter_diffy <- vac69a.cytof::vac63_diffcyt_boxplot(da_t6_ter, sce, FDR = 1, logFC = log2(1))


# t6_prim_diffy <- vac69a.cytof::vac63_diffcyt_boxplot(da_t6_prim, sce, FDR = 0.05, logFC = log2(2))
# t6_ter_diffy <- vac69a.cytof::vac63_diffcyt_boxplot(da_t6_ter, sce, FDR = 0.05, logFC = log2(2))



t6_prim_diffy <- t6_prim_diffy$data
t6_prim_diffy$batch <- md$batch[match(t6_prim_diffy$sample_id, md$sample_id)]
t6_prim_diffy$n_infection <- md$n_infection[match(t6_prim_diffy$sample_id, md$sample_id)]





t6_ter_diffy <- t6_ter_diffy$data
t6_ter_diffy$batch <- md$batch[match(t6_ter_diffy$sample_id, md$sample_id)]
t6_ter_diffy$n_infection <- md$n_infection[match(t6_ter_diffy$sample_id, md$sample_id)]


combo <- rbind(t6_ter_diffy, t6_prim_diffy)

combo$cluster_id <- as.character(combo$cluster_id)

#combo$clusterf <- ifelse(combo$cluster_id %in% t6_ter_diffy$cluster_id, paste("_", combo$cluster_id, sep=''), combo$cluster_id)

sig_t6_all_plot <- ggplot(combo, aes(x=factor(timepoint, levels=c("Baseline", "DoD", "T6", "C45")), y=frequency))+
  geom_boxplot(aes(fill=n_infection))+
  geom_point(aes(group=n_infection, colour=volunteer), position = position_dodge(width = 0.75))+
  facet_wrap(~cluster_id, scales = "free", ncol = 7, labeller = label_wrap_gen(width = 15))+
  theme_minimal()+
  scale_y_continuous()+
  ylab("% of all CD3+ cells")+
  scale_fill_manual(values = c("First" = NA,
                               "Third" = NA))+    
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        axis.title.x = element_blank(),
        strip.text = element_text(size=7))

# ggsave("./figures/sig_t6_combo_boxplot.png", sig_t6_all_plot, height=14, width=13)
ggsave("./figures/figures_for_paper/supp_all_cluster_freqs.pdf", sig_t6_all_plot, height=14, width=13)



# sig cluster median expression heatmap

shplit_cells <- function(x, by) {
  stopifnot(is.character(by), by %in% colnames(colData(x)))
  cd <- data.frame(colData(x))
  cd$cluster_id <- cluster_ids(x, k="flo_merge")
  dt <- data.table::data.table(cd, i = seq_len(ncol(x)))
  dt_split <- split(dt, by = by, sorted = TRUE, flatten = FALSE)
  purrr::map_depth(dt_split, length(by), "i")
}


ahgg <- function(x, by, fun = c("median", "mean", "sum")) {
  fun <- switch(match.arg(fun),
                median = rowMedians, mean = rowMeans, sum = rowSums)
  cs <- shplit_cells(x, by)
  pb <- purrr::map_depth(cs, -1, function(i) {
    if (length(i) == 0) return(numeric(nrow(x)))
    fun(assay(x, "exprs")[, i, drop = FALSE])
  })
  purrr::map_depth(pb, -2, function(u) as.matrix(data.frame(
    u, row.names = rownames(x), check.names = FALSE)))
}






#read in sig clusters, get rid of resting Vd cluster
sig_clusters <- read.delim("/home/flobuntu/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/prim_ter_sig_clusters.csv", header=FALSE)$V1


ms <- ahgg(sce[refined_markers,], by = c("cluster_id", "timepoint"), fun="median")

ms2 <- lapply(ms, function(x)data.frame(x))
ms2 <- Map(cbind, ms2, cluster_id = names(ms))

ms3 <- data.table::rbindlist(ms2)
ms3[,Marker := unlist(lapply(ms2, rownames))]
ms_mean <- mutate(ms3, "mean_expression"= apply(ms3[,1:4], MARGIN = 1, mean ))

#write.csv(ms_mean, "/home/flobuntu/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/cluster_medians_prim_ter.csv")

ms3 <- select(ms3, T6, cluster_id, Marker)

ms4 <- as.matrix(tidyr::spread(ms3, cluster_id, T6))

rownames(ms4) <- ms4[,1]
ms5 <- ms4[,2:ncol(ms4)]

scaled_mat <- apply(apply(ms5, c(1,2), as.numeric), MARGIN = 1, function(x)scales::rescale(x, to=c(0, 1)))
write.csv(scaled_mat, "/home/flobuntu/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/cluster_medians_prim_ter_T6.csv")

  
scaled_mat <- read.csv("/home/flobuntu/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/cluster_medians_prim_ter_T6.csv", row.names = 1)
reordered_scaled_mat <- scaled_mat[,match(refined_markers, colnames(scaled_mat))]








ter_t6_df <- read.csv("/home/flobuntu/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/differential_abundance/edgeR/ter_t6_df_edger.csv")

ter_sig_clusters <- ter_t6_df %>%
  filter(p_adj <= 0.05, logFC >=1)%>%
  select(cluster_id)



prim_sig_clusters <- read.delim("/home/flobuntu/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/prim_ter_sig_clusters.csv", header=FALSE)$V1

colcsv <- read.csv("~/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/cluster_colours.csv")
cluster_order <- read.csv("/home/flobuntu/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/cluster_order.csv")$x

sig_colcsv <- subset(colcsv, colcsv$cluster_id %in% prim_sig_clusters)

sig_colcsv <- sig_colcsv[rev(match(cluster_order, sig_colcsv$cluster_id)), ]

col_pal <- sig_colcsv$colour
names(col_pal) <- sig_colcsv$cluster_id


both <- ifelse(rownames(reordered_scaled_mat) %in% ter_sig_clusters$cluster_id, "both", "first only")


combo_left_anno_var <-  rowAnnotation(annotation_name_gp = gpar(fontsize=10),
                                      annotation_name_rot = 45,
                                      gap = unit(1.5, "mm"),
                                      "cluster_id"=sig_colcsv$cluster_id,
                                      "n_infection"=both,
                                      show_legend = c(FALSE, TRUE),
                                      
                                      show_annotation_name = TRUE,
                                      simple_anno_size = unit(8, "mm"), # width of the significance bar
                                      col=list("n_infection" = c("both"="darkgrey", "first only"="#36454f"),
                                               "cluster_id" = col_pal),
                                      annotation_legend_param = list(n_infection = list(title ="n_infection",
                                                                                        at = c("both", "first only"),
                                                                                        title_gp=gpar(angle=45),
                                                                                        legend_gp = gpar(fill = c("darkgrey","#36454f")),
                                                                                        title_position = "topleft")
                                                                     # cluster_id = list(title="cluster_id",
                                                                     #                   at = sig_colcsv$cluster_id,
                                                                     #                   legend_gp = gpar(fill = unname(col_pal)))
                                      )
                                      
)





inferno <- colorspace::sequential_hcl("inferno", n=8)
col_inferno <- circlize::colorRamp2(seq(0,1, by=1/{length(inferno)-1}), inferno)

reordered_scaled_mat <- subset(reordered_scaled_mat, rownames(reordered_scaled_mat) %in% prim_sig_clusters)

rereordered_scaled_mat <- reordered_scaled_mat[rev(cluster_order),]


# ht_opt("ROW_ANNO_PADDING" = unit(3,"mm"))
# ht_opt(RESET = TRUE)

all_cluster_heatmap <- Heatmap(matrix = as.matrix(rereordered_scaled_mat),
                               #cluster_rows = FALSE,
                               cluster_rows = TRUE,
                               show_row_dend = TRUE,
                               show_heatmap_legend = FALSE,
                               name = "Median Marker Expression",
                               cluster_columns = FALSE,
                               row_names_side = "left",
                               row_dend_side = "right",
                               col = col_inferno,
                               rect_gp = gpar(col = "white"),
                               left_annotation = combo_left_anno_var,
                               column_names_rot = 45
                               #heatmap_legend_param = list(col = col_fun4, title = "Normalised Frequency", title_position = "topleft"),
                               # width = unit(16, "cm"),
                               # height = unit(16*34/28, "cm")
)

  png("/home/flobuntu/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/figures/sig_cluster_heatmap_var.png", width=16, height=8, units = "in", res=400)
  draw(all_cluster_heatmap, padding = unit(c(2, 25, 2, 15), "mm"))
  dev.off()



all_cluster_heatmap <- Heatmap(#matrix = rereordered_scaled_mat,
  matrix = as.matrix(rereordered_scaled_mat),
  cluster_rows = TRUE,
  show_row_dend = FALSE,
  show_heatmap_legend = FALSE,
  name = "Median Marker Expression",
  cluster_columns = FALSE,
  row_names_side = "left",
  col = col_inferno,
  rect_gp = gpar(col = "white"),
  left_annotation = combo_left_anno_var,
  column_names_rot = 45
  #heatmap_legend_param = list(col = col_fun4, title = "Normalised Frequency", title_position = "topleft"),
  # width = unit(16, "cm"),
  # height = unit(16*34/28, "cm")
)




png("/home/flobuntu/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/figures/sig_t6_cluster_heatmap.png", width=12, height=5, units = "in", res=400)
draw(all_cluster_heatmap, padding = unit(c(2, 25, 2, 15), "mm"))
dev.off()

