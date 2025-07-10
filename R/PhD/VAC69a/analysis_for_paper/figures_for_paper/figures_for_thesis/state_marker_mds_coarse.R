


daf <- read_full("~/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/")

coarse_table <- read.csv("/home/flobuntu/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/merging_tables/most_coarse_merge.csv", header=T, stringsAsFactors = F)



#get rid of spaces at beginning of string

coarse_table$new_cluster <- factor(coarse_table$new_cluster)

merged_daf<- mergeClusters(daf, k = "mod_meta45", table = coarse_table, id = "coarse_merge")


plotClusterHeatmap(merged_daf, hm2=NULL,
                   k = "mod_meta45",
                   m = "coarse_merge",
                   cluster_anno = FALSE,
                   draw_freqs = TRUE,
                   scale = TRUE
)




shplit_cells <- function(x, by) {
  stopifnot(is.character(by), by %in% colnames(colData(x)))
  cd <- data.frame(colData(x))
  cd$cluster_id <- cluster_ids(x, "coarse_merge")
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

refined_markers <- read.csv("~/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/refined_markers.csv")

limma_markers <- c("CD38","ICOS","CD27","Tbet","Ki67","FoxP3","CD127","PD1","BCL2","GZB","CTLA4","CD25","CD95","HLADR","CD28")


ms <- ahgg(merged_daf[limma_markers, ], by = c("cluster_id", "sample_id"))
ms <- lapply(ms, reshape2::melt, varnames = c("antigen", "sample_id"))
ms <- bind_rows(ms, .id = "cluster_id")

ms$timepoint <- substr(ms$sample_id, 5, nchar(as.character(ms$sample_id)))



slimmer_ms <- ms %>%
  mutate("contrast"= paste(antigen, " on ", cluster_id))%>%
  select(contrast, sample_id, value) %>%
  pivot_wider(names_from = sample_id, values_from = value)

slimmer_ms <- data.frame(slimmer_ms)
rownames(slimmer_ms) <- slimmer_ms$contrast

slimmer_ms$contrast <- NULL

mds <- limma::plotMDS(slimmer_ms)
df <- data.frame(MDS1 = mds$x, MDS2 = mds$y)
md <- S4Vectors::metadata(merged_daf)$experiment_info
m <- match(rownames(df), md$sample_id)
df <- data.frame(df, md[m, ])

#df2 <- filter(df, timepoint %in% c("Baseline", "DoD"))

state_markers_mds_coarse <- ggplot(df, aes(x=MDS1, y=MDS2))+
  geom_point(aes(shape=timepoint, color=volunteer, fill=volunteer))+
  #ggrepel::geom_label_repel(aes_string(label = "sample_id"), show.legend = FALSE)+ 
  theme_minimal()+
  xlab("")+
  ggtitle("Cell State Marker Expression")+
  #scale_shape_manual(values = c("Baseline"=21, "C10"=24, "Diagnosis"=22, "T6"=3))+
  scale_color_manual(values = volunteer_palette)+
  scale_fill_manual(values = volunteer_palette)+
  guides(color=guide_legend(title="Volunteer",override.aes = list(size = 1)),
         shape=guide_legend(title="Timepoint", override.aes = list(size = 1)),
         fill=guide_none())+
  theme(legend.title = element_text(size=8),
        legend.text = element_text(size=6),
        plot.title = element_text(size=10, hjust=0.5),
        legend.position = "none")



indie_var_plot <- cowplot::plot_grid(state_markers_mds_coarse, aitchison_cytof, plasma_pca, indie_var_lgd, nrow = 1, rel_widths = c(3,3,3,1))

ggsave("~/PhD/figures_for_thesis/chapter_2/indie_var_plot_coarse_merge.pdf", indie_var_plot, height = 3.5, width=8)




