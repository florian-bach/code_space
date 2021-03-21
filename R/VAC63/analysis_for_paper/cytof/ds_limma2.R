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

coarse_table <- read.csv("/home/flobuntu/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/vac63c_coarse_merge_table.csv", header=T, stringsAsFactors = F)
sce <- CATALYST::mergeClusters(sce, k = "meta50", table = coarse_table, id = "coarse_merge", overwrite = TRUE)



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
                       analysis_type = "DS", method_DS = "diffcyt-DS-limma",
                       clustering_to_use = "coarse_merge", verbose = TRUE)

da_dod_ter <- diffcyt(tertiaries,
                      design = ter_design,
                      contrast =ter_dod_contrast,
                      analysis_type = "DS", method_DS = "diffcyt-DS-limma",
                      clustering_to_use = "coarse_merge", verbose = TRUE)

table(rowData(da_dod_prim$res)$p_adj < 0.05)
# FALSE  TRUE 
# 48     1 

table(rowData(da_dod_ter$res)$p_adj < 0.05)
# FALSE  TRUE 
# 18    31 
ter_dod_df <- data.frame(rowData(da_dod_ter$res))


prim_dod_df <- data.frame(rowData(da_dod_prim$res))
prim_dod_df <- subset(prim_dod_df, prim_dod_df$p_adj<0.05 & abs(prim_dod_df$logFC)>1)
prim_dod_df$n_infection <- "first"

ter_dod_df <- data.frame(rowData(da_dod_ter$res))
#write.csv(rbind(prim_dod_df, ter_dod_df), "./differential_abundance/edgeR/sig_dod_edgeR.csv", row.names = FALSE)

# 
# dod_prim_diffy <- vac69a.cytof::vac63_diffcyt_boxplot(da_dod_prim, sce, FDR = 0.05, logFC = log2(2))
# dod_ter_diffy <- vac69a.cytof::vac63_diffcyt_boxplot(da_dod_ter, sce, FDR = 0.05, logFC = log2(2))

diffy <- dod_ter_diffy

sig_cluster_boxplot_data <- diffy$data

sig_cluster_boxplot_data$batch <- md$batch[match(sig_cluster_boxplot_data$sample_id, md$sample_id)]
sig_cluster_boxplot_data$n_infection <- md$n_infection[match(sig_cluster_boxplot_data$sample_id, md$sample_id)]
sig_cluster_boxplot_data$direction <- ifelse(sig_cluster_boxplot_data)


da_t6_prim <- diffcyt(primaries,
                      design = prim_design,
                      contrast = prim_t6_contrast,
                      analysis_type = "DS", method_DS = "diffcyt-DS-limma",
                      clustering_to_use = "coarse_merge", verbose = TRUE)


da_t6_ter <- diffcyt(tertiaries,
                     design = ter_design,
                     contrast = ter_t6_contrast,
                     analysis_type = "DS", method_DS = "diffcyt-DS-limma",
                     clustering_to_use = "coarse_merge", verbose = TRUE)
#View(data.frame(rowData(all_da_t6$res)))


table(rowData(da_t6_prim$res)$p_adj < 0.05)
# FALSE  TRUE 
# 21    26 
table(rowData(da_t6_ter$res)$p_adj < 0.05)
# FALSE  TRUE 
# 33    13 

prim_c45_contrast <- createContrast(c(c(0,1,0,0), rep(0, 2)))
ter_c45_contrast <- createContrast(c(c(0,1,0,0), rep(0, 5)))




da_c45_prim <- diffcyt(primaries,
                      design = prim_design,
                      contrast = prim_c45_contrast,
                      analysis_type = "DS", method_DS = "diffcyt-DS-limma",
                      clustering_to_use = "coarse_merge", verbose = TRUE)


da_c45_ter <- diffcyt(tertiaries,
                     design = ter_design,
                     contrast = ter_c45_contrast,
                     analysis_type = "DS", method_DS = "diffcyt-DS-limma",
                     clustering_to_use = "coarse_merge", verbose = TRUE)



table(rowData(da_c45_prim$res)$p_adj < 0.05)
# FALSE  TRUE 
# 21    26 
table(rowData(da_c45_ter$res)$p_adj < 0.05)
# FALSE  TRUE 
# 33    13 

c45_ter_df <- data.frame(rowData(da_c45_ter$res))
filter(c45_ter_df, p_adj<0.05)

prim_t6_df <- data.frame(rowData(da_t6_prim$res))


ter_t6_df <- data.frame(rowData(da_t6_ter$res))

# write.csv(ter_t6_df, "./differential_abundance/edgeR/ter_t6_df_edger.csv", row.names = FALSE)
# write.csv(prim_t6_df, "./differential_abundance/edgeR/prim_t6_df_edger.csv", row.names = FALSE)

prim_t6_contrasts <- subset(prim_t6_df, prim_t6_df$p_adj<0.05 & abs(prim_t6_df$logFC) > log2(1.1))
ter_t6_contrasts <- subset(ter_t6_df, ter_t6_df$p_adj<0.05 & abs(ter_t6_df$logFC) > log2(1.1))
ter_dod_contrasts <- subset(ter_dod_df, ter_dod_df$p_adj<0.05 & abs(ter_dod_df$logFC) > log2(1.1))

prim_t6_contrasts <- paste(prim_t6_contrasts$cluster_id, " (", prim_t6_contrasts$marker_id, ")", sep="")
ter_t6_contrasts <- paste(ter_t6_contrasts$cluster_id, " (", ter_t6_contrasts$marker_id, ")", sep="")
ter_dod_contrasts <- paste(ter_dod_contrasts$cluster_id, " (", ter_dod_contrasts$marker_id, ")", sep="")


#write out matrices


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



ter_dod_ms <- ahgg(sce[refined_markers, ], by = c("cluster_id", "sample_id"))
ter_dod_ms <- lapply(ter_dod_ms, reshape2::melt, varnames = c("antigen", "sample_id"))
ter_dod_ms <- bind_rows(ter_dod_ms, .id = "cluster_id")

ter_dod_ms$timepoint <- substr(ter_dod_ms$sample_id, 6, nchar(as.character(ter_dod_ms$sample_id)))

scaled_ms_ter_dod <- ter_dod_ms %>%
  group_by(cluster_id, antigen) %>%
  mutate("scaled_value"=scale(value)) %>%
  ungroup()


wide_scaled_ms <- scaled_ms_ter_dod %>%
  mutate("contrast"=paste(cluster_id, " (", antigen, ")", sep="")) %>%
  filter(contrast %in% ter_t6_contrasts) %>%
  select(scaled_value, contrast, sample_id, antigen) %>%
  pivot_wider(names_from = sample_id, values_from = scaled_value)


wide_scaled_ms <- data.frame(wide_scaled_ms)
rownames(wide_scaled_ms) <- wide_scaled_ms$contrast



wide_scaled_ms$contrast <- NULL
num_wide_scaled_ms <- as.matrix(wide_scaled_ms[, 2:ncol(wide_scaled_ms)])

write.csv(num_wide_scaled_ms, "/home/flobuntu/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/differential_abundance/ds_limma/ter_dod_ms.csv")
write.csv(num_wide_scaled_ms, "/home/flobuntu/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/differential_abundance/ds_limma/prim_t6_ms.csv")
write.csv(num_wide_scaled_ms, "/home/flobuntu/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/differential_abundance/ds_limma/ter_t6_ms.csv")



try1 <- read.csv("/home/flobuntu/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/differential_abundance/ds_limma/ter_dod_ms.csv")
try2 <- read.csv("/home/flobuntu/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/differential_abundance/ds_limma/prim_t6_ms.csv")
try3 <- read.csv("/home/flobuntu/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/differential_abundance/ds_limma/ter_t6_ms.csv")



# 
# write.table(x = sig_prim_t6_df$cluster_id, "/home/flobuntu/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/prim_ter_sig_clusters.csv", row.names = FALSE, col.names = FALSE)
# data.table:::fwrite(x = data.frame(rbind(all_sig_ter_t6_df, all_sig_prim_t6_df)), "/home/flobuntu/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/differential_abundance/edgeR/sig_t6_edgeR.csv")
