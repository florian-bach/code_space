library(CATALYST)
library(diffcyt)
# library(dplyr)
# library(tidyr)
#library(ggplot2)


setwd("~/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only")
fcs <- list.files(pattern = "fcs")
#fcs <- subset(fcs, !grepl(pattern = "C45", fcs))
#fcs <- subset(fcs, grepl(pattern = "Baseline", fcs))
#fcs <- subset(fcs, !grepl(pattern = "DoD", fcs))
fcs <- subset(fcs, !grepl(pattern = "ctrl", fcs))

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

plotExprHeatmap(x = sce, by = "cluster", row_clust = FALSE, k = "meta45", bars = TRUE)
plotAbundances(x = sce, by = "cluster_id", k = "meta45", group_by = "batch")

primaries <- filterSCE(sce, volunteer %in% c("v313", "v315", "v320"))
tertiaries <- filterSCE(sce, volunteer %in% c("v301", "v304", "v305", "v306", "v308", "v310"))

# sig_ter <- filterSCE(sce, cluster_id %in% paste(subset(rowData(ter_da_t6$res), rowData(ter_da_t6$res)$p_adj < 0.05)$cluster_id), k = "meta45")
# plotAbundances(x = sig_ter, by = "cluster_id", k = "meta45", group_by = "batch")


# sig_prime <- filterSCE(sce, cluster_id %in% paste(subset(rowData(prime_da_t6$res), rowData(prime_da_t6$res)$p_adj < 0.05)$cluster_id), k = "meta45")
# plotAbundances(x = sig_prime, by = "cluster_id", k = "meta45", group_by = "timepoint")


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




all_ei <- metadata(sce)$experiment_info
all_design <- createDesignMatrix(all_ei, c("timepoint", "volunteer"))
colnames(all_design)
all_t6_contrast <- createContrast(c(c(0,0,0,1), rep(0, 10)))


all_da_t6 <- diffcyt(sce,
                       design = all_design,
                       contrast = all_t6_contrast,
                       analysis_type = "DA",
                       method_DA = "diffcyt-DA-edgeR",
                       clustering_to_use = "meta45",
                       verbose = T)

table(rowData(all_da_t6$res)$p_adj < 0.05)



da_formula <- createFormula(all_ei,
                             cols_fixed = "timepoint",
                             cols_random = "volunteer")


voom_t6_contrast <- createContrast(c(c(0,0,0,1), rep(0, 10)))

voom_all_da_t6 <- diffcyt(sce,
                     contrast = all_t6_contrast,
                     #formula = da_formula,
                     design = all_design,
                     method_DA = "diffcyt-DA-voom",
                     clustering_to_use = "meta45",
                     verbose = T)

voom_dod_contrast <- createContrast(c(c(0,1,0,0), rep(0, 10)))

voom_all_da_dod <- diffcyt(sce,
                          contrast = voom_dod_contrast,
                          #formula = da_formula,
                          design = all_design,
                          method_DA = "diffcyt-DA-voom",
                          clustering_to_use = "meta45",
                          verbose = T)

table(rowData(voom_all_da_dod$res)$p_adj<0.05 & abs(rowData(voom_all_da_dod$res)$logFC) > 1)

paste(subset(rowData(voom_all_da_dod$res), rowData(voom_all_da_dod$res)$p_adj < 0.05 & abs(rowData(voom_all_da_dod$res)$logFC)>1)$cluster_id)

# voom, timepoint+volunteer
# FALSE  TRUE 
# 29    16

# edgeR, timepoint+volunteer
# FALSE  TRUE 
# 30    15


diffy <- vac69a.cytof::vac63_diffcyt_boxplot(all_da_t6, sce, FDR = 0.05, logFC = 1)
  
sig_cluster_boxplot_data <- diffy$data

sig_cluster_boxplot_data$batch <- md$batch[match(sig_cluster_boxplot_data$sample_id, md$sample_id)]
sig_cluster_boxplot_data$n_infection <- md$n_infection[match(sig_cluster_boxplot_data$sample_id, md$sample_id)]

library(ggplot2)

time_col <- colorspace::sequential_hcl(5, palette = "Purple Yellow")



sig_t6_all_plot <- ggplot(sig_cluster_boxplot_data, aes(x=timepoint, y=frequency))+
  geom_boxplot(aes(fill=timepoint))+
  geom_point(aes(colour=n_infection))+
  facet_wrap(~cluster_id, scales = "free", ncol = 5)+
  theme_minimal()+
  ylab("% of all CD3+ cells")+
  scale_fill_manual(values = c("Baseline"=time_col[4],
                               "DoD"=time_col[2],
                               "T6"=time_col[1],
                               "C45"=time_col[5]))+
  scale_colour_manual(values = c("First"="red",
                                 "Second"="darkblue",
                                 "Third"="darkgreen"))+
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        axis.title.x = element_blank())
  
ggsave("/home/flobuntu/PhD/cytof/vac63c/figures/sig_t6_boxplot.png", sig_t6_all_plot, height=5, width=8)

sig_all_t6 <- filterSCE(sce, cluster_id %in% paste(subset(rowData(all_da_t6$res), rowData(all_da_t6$res)$p_adj < 0.05 & abs(rowData(all_da_t6$res)$logFC)>1)$cluster_id), k = "meta45")
#plotAbundances(x = sig_ter, by = "cluster_id", k = "meta45", group_by = "batch")
sig_t6_cluster_phenotype <- plotExprHeatmap(x = sig_all_t6,
                features = refined_markers,
                by = "cluster",
                row_clust = FALSE,
                col_clust = FALSE,
                k = "meta45",
                bars = TRUE,
                perc=TRUE,
                hm_pal = colorspace::sequential_hcl("inferno", n=8))

png("/home/flobuntu/PhD/cytof/vac63c/figures/sig_t6_cluster_phenotype.png", height=6, width=8, units = "in", res = 600)
ComplexHeatmap::draw(sig_t6_cluster_phenotype)
dev.off()

