# gating by 2d density ####


library(openCyto)
library(ggcyto)
# library(dplyr)
# library(tidyr)
#library(ggplot2)
# 
`%notin%` <- Negate(`%in%`)
# 
# 
# .scatter <- function(gs, chs, gate_id = NULL,
#                      subset = ifelse(is.null(gate_id), "root", "_parent_")) {
#   p <- ggcyto(gs, max_nrow_to_plot = 1e5,
#               aes_string(chs[1], chs[2]), subset) +
#     geom_hex(bins = 100) + facet_wrap(~ name, ncol = 5) +
#     (if (is.null(gate_id)) list() else geom_gate(gate_id)) +
#     ggtitle(NULL) + theme_bw(base_size = 8) + theme(
#       aspect.ratio = 1,
#       legend.position = "none",
#       panel.grid.minor = element_blank(),
#       strip.background = element_rect(fill = NA),
#       axis.text = element_text(color = "black"),
#       axis.text.x = element_text(angle = 45, hjust = 1))
#   suppressMessages(p + coord_equal(expand = FALSE,
#                                    xlim = c(-1, 11), ylim = c(-1, 11)))
# }



setwd("~/PhD/cytof/vac63c/normalised_renamed_comped/debarcoded/")

fcs <- list.files(pattern = "fcs")
#fcs <- subset(fcs, !grepl(pattern = "C45", fcs))
#fcs <- subset(fcs, grepl(pattern = "Baseline", fcs))
fcs <- subset(fcs, grepl(pattern = "C45", fcs))
#fcs <- subset(fcs, grepl(pattern = "Baseline", fcs))

vac63_flowset <- flowCore::read.flowSet(fcs)

vac63_gs <- GatingSet(vac63_flowset)


#define DNA channels
panel <- read.csv("~/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/vac63c_panel.csv", header=T)
dna <- grep("^Ir", panel$fcs_colname, value = TRUE)



gs_add_gating_method(vac63_gs,
                     alias = "cells",
                     pop = "+", parent = "root",
                     dims = "Ir191Di,Ce140Di",
                     gating_method = "flowClust.2d",
                     gating_args = "K=1,quantile=0.97,target=c(5,5)")



gs_add_gating_method(vac63_gs,
                     alias = "singlets",
                     pop = "+", parent = "cells",
                     dims = paste(dna, collapse = ","),
                              gating_method = "flowClust.2d",
                              gating_args = "K=1,quantile=0.97,target=c(5,5)")
                     
df <- gs_pop_get_stats(vac63_gs,
                       type = "percent",
                       nodes = c("cells", "singlets"))

df


fs <- gs_pop_get_data(vac63_gs, "/cells/singlets") # get data from ’GatingSet’
write.flowSet(x = fs, outdir = "./whole_blood_single_cells/")


# gating by flowSOM cluster exclusion ####




library(CATALYST)
#library(diffcyt)
# library(dplyr)
# library(tidyr)
#library(ggplot2)

`%notin%` <- Negate(`%in%`)


setwd("~/PhD/cytof/vac63c/normalised_renamed_comped/debarcoded/whole_blood_single_cells/pbmcs/")
#setwd("~/PhD/cytof/vac63c/normalised_renamed_comped/debarcoded/")


fcs <- list.files(pattern = "fcs")
fcs <- subset(fcs, !grepl(pattern = "ctrl", fcs))

#comment next line
#fcs <- subset(fcs, grepl(pattern = "C45", fcs))


md <- read.csv("~/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/vac63c_metadata.csv", header = T)

#comment next line
#md <- subset(md, grepl(pattern = "DoD", md$timepoint))


wb_md <- cbind("file_name"=fcs, md[,2:6])



#fcs <- subset(fcs, grepl(pattern = "Baseline", fcs))
#fcs <- subset(fcs, grepl(pattern = "Baseline", fcs))
#fcs <- subset(fcs, grepl(pattern = "DoD", fcs))


vac63_flowset <- flowCore::read.flowSet(fcs)

# 
# sampling_ceiling <- 45000
# Being reproducible is a plus
set.seed(1234)

# sample.int takes a sample of the specified size from the elements of x using either with or without replacement.
# smaller_vac63_flowset <- flowCore::fsApply(vac63_flowset, function(ff) {
#   idx <- sample.int(nrow(ff), min(sampling_ceiling, nrow(ff)))
#   ff[idx,]  # alt. ff[order(idx),]
# })
# 


smaller_vac63_flowset <- flowCore::fsApply(vac63_flowset, function(ff){
    idx <- sample.int(nrow(ff), nrow(ff)/min(flowCore::fsApply(vac63_flowset, nrow))*8000)
    ff[idx,]
  })

cbind(flowCore::fsApply(vac63_flowset, nrow), flowCore::fsApply(smaller_vac63_flowset, nrow)[,1])

vac63_flowset <- NULL

wb_md <- subset(wb_md, wb_md$file_name %in% fcs)
# wb_md <- wb_md[order(wb_md$timepoint),]


#  
panel <- read.csv("~/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/vac63c_panel.csv", header=T)

#uncomment for pbmc check
panel <- panel[-c(1,2,65:68),]

sce <- prepData(smaller_vac63_flowset, panel, wb_md, md_cols =
                  list(file = "file_name", id = "sample_id", factors = c("volunteer", "timepoint", "n_infection", "batch")))



# p <- plotExprs(sce, color_by = "batch")
# p$facet$params$ncol <- 7                   
# ggplot2::ggsave("marker_expression.png", p, width=12, height=12)                                      


refined_markers <- {c("CD4",
                     "CD8",
                     "Vd2",
                     "Va72",
                     "CXCR5",
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
                     "TCRgd",
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
                     "CCR7")}

all_markers <- c("CD45", "CD3", "CD3", "CD14", "CD16", refined_markers)

set.seed(123);sce <- CATALYST::cluster(sce, features = all_markers, xdim = 12, ydim = 12, maxK = 45)





# 
vac63c_control_tcell_cluster_heatmap <- plotExprHeatmap(x = sce, by = "cluster", row_clust = FALSE, col_clust = FALSE,  k = "flo_merge", bars = TRUE, features = all_markers)
# pdf("./figures/flo_merge_vac63c_pbmc_cluster_heatmap.pdf", height = 8, width = 11)
png("./figures/flo_merge_vac63c_pbmc_cluster_heatmap.png", height = 8, width = 11, units = "in", res = 250)
vac63c_control_tcell_cluster_heatmap
dev.off()
# 
# 
#pbmc <- filterSCE(sce, cluster_id %notin% c(21,22,28:30), k="meta30")

# fs <- sce2fcs(pbmc, split_by="sample_id")
# 
# flowCore::write.flowSet(x = fs, outdir = "./pbmcs/")


# edger ####


meta45_table <- read.csv("/home/flobuntu/PhD/cytof/vac63c/normalised_renamed_comped/debarcoded/whole_blood_single_cells/pbmcs/prelim_meta45_flo_merge.csv")

sce <- CATALYST::mergeClusters(sce, k = "meta45", table = meta45_table, id = "flo_merge")



library(diffcyt)

all_ei <- metadata(sce)$experiment_info
all_design <- createDesignMatrix(all_ei, c("timepoint", "volunteer"))
colnames(all_design)



all_dod_contrast <- createContrast(c(c(0,0,1,0), rep(0, 10)))

all_da_dod <- diffcyt(sce,
                      design = all_design,
                      contrast = all_dod_contrast,
                      analysis_type = "DA",
                      method_DA = "diffcyt-DA-edgeR",
                      clustering_to_use = "flo_merge",
                      verbose = T)

table(rowData(all_da_dod$res)$p_adj < 0.05)
# FALSE  TRUE 
# 1    29



all_t6_contrast <- createContrast(c(c(0,0,0,1), rep(0, 10)))

all_da_t6 <- diffcyt(sce,
                     design = all_design,
                     contrast = all_t6_contrast,
                     analysis_type = "DA",
                     method_DA = "diffcyt-DA-edgeR",
                     clustering_to_use = "flo_merge",
                     verbose = T)

table(rowData(all_da_t6$res)$p_adj < 0.05)
# FALSE  TRUE 
# 17    13 


all_c45_contrast <- createContrast(c(c(0,1,0,0), rep(0, 10)))

all_da_c45 <- diffcyt(sce,
                      design = all_design,
                      contrast = all_c45_contrast,
                      analysis_type = "DA",
                      method_DA = "diffcyt-DA-edgeR",
                      clustering_to_use = "flo_merge",
                      verbose = T)

table(rowData(all_da_c45$res)$p_adj < 0.05)
# FALSE  TRUE 
# 29     1


da_c45 <- rowData(all_da_c45$res)
da_t6<- rowData(all_da_t6$res)
da_dod <- rowData(all_da_dod$res)


plotDiffHeatmap(x=sce, y=da_dod, fdr=0.05, k="meta45", all = TRUE, lfc = 1)






#set log2FC to 1.5 for DoD because there are just too many...
diffy <- vac69a.cytof::vac63_diffcyt_boxplot(all_da_c45, sce, FDR = 0.05, logFC = log2(2))

sig_cluster_boxplot_data <- diffy$data

sig_cluster_boxplot_data$batch <- md$batch[match(sig_cluster_boxplot_data$sample_id, md$sample_id)]
sig_cluster_boxplot_data$n_infection <- md$n_infection[match(sig_cluster_boxplot_data$sample_id, md$sample_id)]

library(ggplot2)

time_col <- colorspace::sequential_hcl(5, palette = "Purple Yellow")



sig_t6_all_plot <- ggplot(sig_cluster_boxplot_data, aes(x=factor(timepoint, levels=c("Baseline", "DoD", "T6", "C45")), y=frequency))+
  geom_boxplot(aes(fill=n_infection))+
  geom_point(aes(group=n_infection, colour=volunteer), position = position_dodge(width = 0.75))+
  facet_wrap(~cluster_id, scales = "free", ncol = 5, labeller = label_wrap_gen())+
  theme_minimal()+
  scale_y_continuous(trans = "sqrt")+
  ylab("% of all CD3+ cells")+
  # scale_fill_manual(values = c("Baseline"=time_col[4],
  #                              "DoD"=time_col[2],
  #                              "T6"=time_col[1],
  #                              "C45"=time_col[5]))+
  scale_fill_manual(values = c("First" = time_col[4],
                               "Second" = time_col[2],
                               "Third" = time_col[1]))+    
  # scale_fill_manual(values = c("First"="red",
  #                                "Second"="darkblue",
  #                                "Third"="darkgreen"))+
  theme(axis.text.x = element_text(angle = 45, hjust=1),
        axis.title.x = element_blank())

#ggsave("./figures/sig_t6_boxplot.png", sig_t6_all_plot, height=7.5, width=12)
#ggsave("./figures/sig_dod_boxplot.png", sig_t6_all_plot, height=7.5, width=12)
ggsave("./figures/sig_c45_boxplot.png", sig_t6_all_plot, height=7.5, width=12)












