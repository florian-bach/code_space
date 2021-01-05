library(CATALYST)
library(diffcyt)
# library(dplyr)
# library(tidyr)
#library(ggplot2)


setwd("~/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only")
fcs <- list.files(pattern = "fcs")
#fcs <- subset(fcs, !grepl(pattern = "C45", fcs))
fcs <- subset(fcs, !grepl(pattern = "ctrl", fcs))
#fcs <- subset(fcs, !grepl(pattern = "DoD", fcs))

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

md <- read.csv("vac63c_metadata.csv", header = T)

md <- subset(md, md$file_name %in% fcs)
md <- md[order(md$timepoint),]

panel <- read.csv("vac63c_panel.csv", header=T)

sce <- prepData(vac63_flowset, panel, md, md_cols =
                  list(file = "file_name", id = "sample_id", factors = c("timepoint", "batch", "volunteer", "n_infection")))

# p <- plotExprs(sce, color_by = "batch")
# p$facet$params$ncol <- 7                   
# ggplot2::ggsave("marker_expression.png", p, width=12, height=12)                                      


refined_markers <- c("CD4",
                   "CD8",
                   "Vd2",
                   "Va72",
                   "CD38",
                   "CD69",
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
                   "CXCR5",
                   "CX3CR1",
                   "CD57",
                   "CD45RA",
                   "CD45RO",
                   "CCR7")

set.seed(123);sce <- CATALYST::cluster(sce, features = refined_markers, xdim = 10, ydim = 10, maxK = 50)

# plotExprHeatmap(x = tertiaries, by = "cluster", row_clust = TRUE, k = "meta45")
plotAbundances(x = tertiaries, by = "cluster_id", k = "meta45", group_by = "batch")

primaries <- filterSCE(sce, volunteer %in% c("v313", "v315", "v320"))
tertiaries <- filterSCE(sce, volunteer %in% c("v301", "v304", "v305", "v306", "v308", "v310"))

# sig_ter <- filterSCE(sce, cluster_id %in% paste(subset(rowData(ter_da_t6$res), rowData(ter_da_t6$res)$p_adj < 0.05)$cluster_id), k = "meta45")
# plotAbundances(x = sig_ter, by = "cluster_id", k = "meta45", group_by = "batch")


sig_prime <- filterSCE(sce, cluster_id %in% paste(subset(rowData(prime_da_t6$res), rowData(prime_da_t6$res)$p_adj < 0.05)$cluster_id), k = "meta45")
plotAbundances(x = sig_prime, by = "cluster_id", k = "meta45", group_by = "timepoint")

prime_ei <- metadata(primaries)$experiment_info
prime_design <- createDesignMatrix(prime_ei, c("timepoint", "volunteer"))

ter_ei <- metadata(tertiaries)$experiment_info
ter_design <- createDesignMatrix(ter_ei, c("timepoint", "volunteer"))

prime_t6_contrast <- createContrast(c(c(0,0,0,1), rep(0, 2)))


prime_da_t6 <- diffcyt(primaries,
                  design = prime_design,
                  contrast = prime_t6_contrast,
                  analysis_type = "DA",
                  method_DA = "diffcyt-DA-edgeR",
                  clustering_to_use = "meta45",
                  verbose = T)

table(rowData(prime_da_t6$res)$p_adj < 0.05)

# FALSE  TRUE 
# 17    28 

plotDiffHeatmap(x=primaries,
                y=rowData(prime_da_t6$res),
                k="meta45",
                assay = "counts",
                lfc=1,
                sort_by = "padj",
                normalize=T,
                all = T)




ter_t6_contrast <- createContrast(c(c(0,0,0,1), rep(0, 5)))

ter_da_t6 <- diffcyt(tertiaries,
                 design = ter_design,
                 contrast = ter_t6_contrast,
                 analysis_type = "DA",
                 method_DA = "diffcyt-DA-edgeR",
                 clustering_to_use = "meta45",
                 verbose = T)

table(rowData(ter_da_t6$res)$p_adj < 0.05)

# FALSE  TRUE 
# 32    13 

plotDiffHeatmap(x=tertiaries,
                y=rowData(ter_da_t6$res),
                k="meta45",
                lfc=2,
                top_n = 15,
                sort_by = "padj",
                normalize=T,
                all = T)
