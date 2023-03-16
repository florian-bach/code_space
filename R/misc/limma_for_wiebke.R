# Make a merging table that organises cells in broad categories like memory / lineage (i.e. not defined by activation / phenotypic markers), these will be the populations within which you will check for differential marker expression
coarse_table <- read.csv("~merging_tables/most_coarse_merge.csv", header=T, stringsAsFactors = F)

merged_sce<- mergeClusters(sce, k = "mod_meta45", table = coarse_table, id = "coarse_merge")

ei <- metadata(merged_sce)$experiment_info

# same rules apply as for differential abundance, pull out the metadata that you want to include in your modelling e.g. condition, batch etc.
design <- createDesignMatrix(ei, c("timepoint", "volunteer"))

levels(ei$timepoint)

###If a design matrix has been used, the entries of contrast correspond to the columns of the design
#matrix and the length of contrast equals the number of columns in the design matrix. If a model formula
#has been used, the entries correspond to the levels of the fixed effect terms;
#and the length equals the number of levels of the fixed effect terms.

FDR_cutoff <- 0.05

#write your contrast like you would for differential abundance, assuming your design is the same (i found that correcting for batch can be useful, so check that out)
pairwise_contrast_t6 <- createContrast(c(c(0, 0, 0, 1), rep(0,5)))

# here you can redefine what "state" markers are to your liking; first obtain a full list of markers and then trim it to your heart's content
states <- names(marker_classes(merged_sce))
states <- states[-c(1:21,53, 54, 56:length(states),
                    match(c("CD20", "TIM3", "KG", "CD3"), states))]

metadata(merged_sce)$id_state_markers <- states

# honestly i can't remember why markers_to_test has to be like that, my old laptop has some issues which means i can't get into it to quickly check, but if markers_to_test = "state" doesn't work try a logic vector like this (length has to be the same as states above)
ds_t6 <- diffcyt(merged_sce,
                 markers_to_test = rep(TRUE, 28),                         
                 design = design,
                 contrast = pairwise_contrast_t6,                    
                 analysis_type = "DS",
                 method_DS =  "diffcyt-DS-limma",         
                 clustering_to_use = "coarse_merge",
                 verbose = TRUE)

#quickly look at results
topTable(ds_t6, top_n = 55, order_by = "cluster_id",
         show_meds = TRUE, format_vals = TRUE, digits = 3)

#store results in table
res_DA_t6 <- topTable(ds_t6, all = TRUE, show_logFC = T)
table(res_DA_t6$p_adj <= 0.05)

sigs_t6 <- subset(res_DA_t6, p_adj<=0.05)


# we need these guys again

shplit_cells <- function(x, by) {
  stopifnot(is.character(by), by %in% colnames(colData(x)))
  cd <- data.frame(colData(x))
  # makre sure to name the correct clustering, whatever you picked right at the top
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

# this is a list of markers that have good resolution
refined_markers <- read.csv("~/refined_markers.csv")

# the following steps create a big matrix with median expression values for each cluster in each sample
ms <- ahgg(merged_sce[refined_markers$refined_markers, ], by = c("cluster_id", "sample_id"))
ms <- lapply(ms, reshape2::melt, varnames = c("antigen", "sample_id"))
ms <- bind_rows(ms, .id = "cluster_id")

# scale marker expression seperately for each cluster and antigen
scaled_ms <- ms %>%
  group_by(cluster_id, antigen) %>%
  mutate("scaled_value"=scale(value)) %>%
  ungroup()

#subset differential expression data to only include cluster-marker combinations with an FDR<0.05 and a 10% shift in expression. you wanto to play with a different percentage cutoff, it's a bit tricky because you can have meaningful shifts in expression that don't move the mean / median much, so look at your data
sig_ds_contrasts <- subset(rowData(ds_t6$res), rowData(ds_t6$res)$p_adj<0.05 & abs(rowData(ds_t6$res)$logFC) > log2(1.1))

# add brackets to the significant cluster-marker combinations for plotting
sig_ds_contrasts <- paste(sig_ds_contrasts$cluster_id, " (", sig_ds_contrasts$marker_id, ")", sep="")

# reshape table to be in the right format for heatmapping
wide_scaled_ms <- scaled_ms %>%
  mutate("contrast"=paste(cluster_id, " (", antigen, ")", sep="")) %>%
  filter(contrast %in% sig_ds_contrasts) %>%
  select(scaled_value, contrast, sample_id, antigen) %>%
  pivot_wider(names_from = sample_id, values_from = scaled_value)

# make data frame (was tibble before)
wide_scaled_ms <- data.frame(wide_scaled_ms)
rownames(wide_scaled_ms) <- wide_scaled_ms$contrast

# ditch contrast column and convert to numerical matrix
wide_scaled_ms$contrast <- NULL
num_wide_scaled_ms <- as.matrix(wide_scaled_ms[, 2:ncol(wide_scaled_ms)])

# this is a plotting trick to make heatmaps more legible by forcing the values to be symmetrical; i thought this was cheeky at first but CATALYST implements this and it does make heatmaps way more legible
num_wide_scaled_ms <- ifelse(num_wide_scaled_ms>abs(min(num_wide_scaled_ms)), abs(min(num_wide_scaled_ms)), num_wide_scaled_ms)




library(ComplexHeatmap)


inferno <- colorspace::sequential_hcl("inferno", n=8)
# make a color palette from blue to yellow via black, adjust the limits (i.e. c(-4, 0, 4)) according to the intensity in your data
col_fun_ds_limma <- circlize::colorRamp2(c(-4, 0, 4), c("#0859C6", "black", "#FFA500"))

# my heatmaps were split vertically by volunteer and horizontally by lineage, in your case this may be more or less complicated; these need to be categorical vectors. the factor levels are the names of the groups and where factor X changes to facotr Y in the vector, there the division will be drawn
# see https://jokergoo.github.io/ComplexHeatmap-reference/book/a-single-heatmap.html#heatmap-split
limma_col_split <- factor(substr(colnames(wide_scaled_ms)[2:ncol(wide_scaled_ms)], 5, 20))
limma_row_split_lineage <- factor(substr(rownames(wide_scaled_ms),0, 4))


median_cluster_heat <- Heatmap(matrix = num_wide_scaled_ms,
                               cluster_rows = FALSE,
                               name = "Scaled Marker Expression",
                               cluster_columns = FALSE,
                               row_names_side = "left",
                               col = col_fun_ds_limma,
                               column_names_gp = gpar(fontsize = 11),
                               column_split = limma_col_split,
                               split = limma_row_split_lineage,
                               rect_gp = gpar(col = "white"),
                               show_heatmap_legend = TRUE,
                               column_names_rot = 45,
                               heatmap_legend_param = list(legend_position = "top",
                                                           col=col_fun_ds_limma,
                                                           title = "Normalised Marker Expression",
                                                           legend_direction = "horizontal",
                                                           title_position = "topcenter",
                                                           legend_width = unit(6.2, "cm"),
                                                           border = FALSE)
                               #width = unit(16, "cm"),
                               #height = unit(16*9/28, "cm")
)

pdf("~/whatever_ds_limma.pdf", width = 8, height=2.6)
draw(median_cluster_heat, heatmap_legend_side = "bottom")
dev.off()
