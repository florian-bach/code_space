library(CATALYST)
library(diffcyt)

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



md <- read.csv("vac63c_metadata.csv", header = T)

md <- subset(md, md$file_name %in% fcs)
md <- md[order(md$timepoint),]


#
panel <- read.csv("~/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/vac63c_panel.csv", header=T)

merged_daf <- prepData(vac63_flowset, panel, md, md_cols =
                  list(file = "file_name", id = "sample_id", factors = c("volunteer", "timepoint", "n_infection", "batch")))

# all_markers <- c("CD45", "CD3", "CD3", "CD14", "CD16", refined_markers)

set.seed(123);merged_daf <- CATALYST::cluster(sce, features = refined_markers, xdim = 10, ydim = 10, maxK = 50)


#merge cells first normally then coarsely
meta45_table <- read.csv("/home/flobuntu/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/prim_ter_50_merge.csv")
merged_daf <- CATALYST::mergeClusters(sce, k = "meta50", table = meta45_table, id = "flo_merge", overwrite = TRUE)


coarse_table <- read.csv("/home/flobuntu/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/vac63c_coarse_merge_table.csv", header=T, stringsAsFactors = F)
merged_daf <- CATALYST::mergeClusters(sce, k = "meta50", table = coarse_table, id = "coarse_merge", overwrite = TRUE)



ei <- metadata(merged_daf)$experiment_info

design <- createDesignMatrix(ei, c("timepoint", "volunteer"))

levels(ei$timepoint)

###If a design matrix has been used, the entries of contrast correspond to the columns of the design
#matrix and the length of contrast equals the number of columns in the design matrix. If a model formula
#has been used, the entries correspond to the levels of the fixed effect terms;
#and the length equals the number of levels of the fixed effect terms.

FDR_cutoff <- 0.05

# limma models with all coefficient for volunteer####

pairwise_contrast_t6 <- createContrast(c(c(0, 0, 0, 1), rep(0,8)))
pairwise_contrast_dod <- createContrast(c(c(0, 0, 1, 0), rep(0,8)))
pairwise_contrast_c45 <- createContrast(c(c(0, 1, 0, 0), rep(0,8)))

# how to test ALLL markers
states <- state_markers(merged_daf)
states_plus <- c(states, "CXCR5", "CD161", "CD27", "FoxP3", "CD56", "CX3CR1", "CD16")


logic <- names(marker_classes(merged_daf)) %in% states

metadata(merged_daf)$id_state_markers <- states_plus



ds_dod <- diffcyt(merged_daf,
                  design = design, contrast = pairwise_contrast_dod,
                  analysis_type = "DS", method_DS = "diffcyt-DS-limma",
                  clustering_to_use = "coarse_merge", verbose = TRUE)


ds_t6 <- diffcyt(merged_daf,
                 design = design, contrast = pairwise_contrast_t6,
                 analysis_type = "DS",  method_DS =  "diffcyt-DS-limma",
                 clustering_to_use = "coarse_merge", verbose = TRUE)

ds_c45 <- diffcyt(merged_daf,
                 design = design, contrast = pairwise_contrast_c45,
                 analysis_type = "DS",  method_DS =  "diffcyt-DS-limma",
                 clustering_to_use = "coarse_merge", verbose = TRUE)





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



ms <- ahgg(merged_daf[refined_markers, ], by = c("cluster_id", "sample_id"))
ms <- lapply(ms, reshape2::melt, varnames = c("antigen", "sample_id"))
ms <- bind_rows(ms, .id = "cluster_id")

ms$timepoint <- substr(ms$sample_id, 6, nchar(as.character(ms$sample_id)))

scaled_ms <- ms %>%
  group_by(cluster_id, antigen) %>%
  mutate("scaled_value"=scale(value)) %>%
  ungroup()



sig_ds_contrasts <- subset(rowData(ds_t6$res), rowData(ds_t6$res)$p_adj<0.05 & abs(rowData(ds_t6$res)$logFC) > log2(1.1))

sig_ds_contrasts <- paste(sig_ds_contrasts$cluster_id, " (", sig_ds_contrasts$marker_id, ")", sep="")



scaled_ms$volunteer <- substr(scaled_ms$sample_id, 0, 4)


wide_scaled_ms <- scaled_ms %>%
  mutate("contrast"=paste(cluster_id, " (", antigen, ")", sep="")) %>%
  filter(contrast %in% sig_ds_contrasts) %>%
  select(scaled_value, contrast, sample_id, antigen) %>%
  pivot_wider(names_from = sample_id, values_from = scaled_value)




wide_scaled_ms <- data.frame(wide_scaled_ms)
rownames(wide_scaled_ms) <- wide_scaled_ms$contrast



wide_scaled_ms$contrast <- NULL
num_wide_scaled_ms <- as.matrix(wide_scaled_ms[, 2:ncol(wide_scaled_ms)])

write.csv(num_wide_scaled_ms, "/home/flobuntu/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/sample_wise_cluster_expression_matrix.csv")
