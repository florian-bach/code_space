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
  
  
  vac63c_control_tcell_cluster_heatmap <- plotExprHeatmap(x = sce,
                                                          by = "cluster",
                                                          row_clust = FALSE,
                                                          col_clust = FALSE,
                                                          #m = "flo_merge",
                                                          k= "meta50",
                                                          bars = TRUE,
                                                          features = refined_markers)
  
  pdf("./figures/vac63c_tcell_cluster_heatmap_meta50_flo_names.pdf", height = 10, width = 9)
  vac63c_control_tcell_cluster_heatmap
  dev.off()
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



vac63c_control_tcell_cluster_heatmap <- plotExprHeatmap(x = sce, by = "cluster", row_clust = FALSE, col_clust = FALSE,  k = "flo_merge", bars = TRUE, features = refined_markers)

pdf("./figures/labeled_vac63c_tcell_cluster_heatmap.pdf", height = 8, width = 12)
vac63c_control_tcell_cluster_heatmap
dev.off()




#plotAbundances(x = sce, by = "cluster_id", k = "meta35", group_by = "batch")

primaries <- filterSCE(sce, volunteer %in% c("v313", "v315", "v320"))
tertiaries <- filterSCE(sce, volunteer %in% c("v301", "v304", "v305", "v306", "v308", "v310"))

all_ei <- metadata(sce)$experiment_info
all_design <- createDesignMatrix(all_ei, c("timepoint", "volunteer", "n_infection"))
colnames(all_design)

glm_formula <- createFormula(all_ei, cols_fixed = c("timepoint", "n_infection"), cols_random = "volunteer")


all_dod_contrast <- createContrast(c(c(0,0,1,0), rep(0, 8)))

all_da_dod <- diffcyt(sce,
                      design = all_design,
                      contrast = createContrast(c(0,0,1,0,0)),
                      analysis_type = "DA",
                      method_DA = "diffcyt-DA-GLMM",
                      clustering_to_use = "meta45",
                      verbose = T)

table(rowData(all_da_dod$res)$p_adj < 0.05)
# FALSE  TRUE 
# 1    29

dod_df <- data.frame(rowData(all_da_dod$res))
dod_df <- subset(dod_df, dod_df$p_adj<0.05 & abs(dod_df$logFC)>1)

all_t6_contrast <- createContrast(c(c(0,0,0,1), rep(0, 8)))
t6_glm_contrast <- createContrast(c(1,0))

all_da_t6 <- diffcyt(sce,
                     design = all_design,
                     contrast = all_t6_contrast,
                     #formula = glm_formula,
                     #contrast = t6_glm_contrast,
                     analysis_type = "DA",
                     method_DA = "diffcyt-DA-voom",
                     #clustering_to_use = "meta45",
                     verbose = T)
#View(data.frame(rowData(all_da_t6$res)))


table(rowData(all_da_t6$res)$p_adj < 0.05)
# FALSE  TRUE 
# 17    13 

t6_df <- data.frame(rowData(all_da_t6$res))
t6_df <- subset(t6_df, t6_df$p_adj<0.05 & abs(t6_df$logFC)>1)


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


#plotDiffHeatmap(x=sce, y=da_t6, fdr=0.05, k="meta30")





#set log2FC to 1.5 for DoD because there are just too many...
diffy <- vac69a.cytof::vac63_diffcyt_boxplot(all_da_t6, sce, FDR = 0.05, logFC = log2(2))
  
sig_cluster_boxplot_data <- diffy$data

sig_cluster_boxplot_data$batch <- md$batch[match(sig_cluster_boxplot_data$sample_id, md$sample_id)]
sig_cluster_boxplot_data$n_infection <- md$n_infection[match(sig_cluster_boxplot_data$sample_id, md$sample_id)]
sig_cluster_boxplot_data$direction <- ifelse(sig_cluster_boxplot_data)

library(ggplot2)

time_col <- colorspace::sequential_hcl(5, palette = "Purple Yellow")



sig_t6_all_plot <- ggplot(sig_cluster_boxplot_data, aes(x=factor(timepoint, levels=c("Baseline", "DoD", "T6", "C45")), y=frequency))+
  geom_boxplot(aes(fill=n_infection))+
  geom_point(aes(group=n_infection, colour=volunteer), position = position_dodge(width = 0.75))+
  facet_wrap(~cluster_id, scales = "free", ncol = 5, labeller = label_wrap_gen(width = 10, multi_line = TRUE))+
  theme_minimal()+
  scale_y_continuous()+
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
  
ggsave("./figures/sig_t6_boxplot.png", sig_t6_all_plot, height=7.5, width=12)
#ggsave("./figures/sig_dod_boxplot.png", sig_t6_all_plot, height=7.5, width=12)
#ggsave("./figures/sig_c45_boxplot.png", sig_t6_all_plot, height=7.5, width=12)




sig_all_t6 <- filterSCE(sce, cluster_id %in% paste(subset(rowData(all_da_t6$res), rowData(all_da_t6$res)$p_adj < 0.05 & abs(rowData(all_da_t6$res)$logFC)>1)$cluster_id), k = "flo_merge")
#plotAbundances(x = sig_ter, by = "cluster_id", k = "meta45", group_by = "batch")
sig_t6_cluster_phenotype <- plotExprHeatmap(x = sig_all_t6,
                features = refined_markers,
                by = "cluster",
                row_clust = FALSE,
                col_clust = FALSE,
                k = "flo_merge",
                bars = TRUE,
                perc=TRUE,
                hm_pal = colorspace::sequential_hcl("inferno", n=8))

png("./figures/sig_t6_cluster_phenotype.png", height=6, width=8, units = "in", res = 600)
ComplexHeatmap::draw(sig_t6_cluster_phenotype, show_heatmap_legend = FALSE)
dev.off()





sig_all_dod <- filterSCE(sce, cluster_id %in% paste(subset(rowData(all_da_dod$res), rowData(all_da_dod$res)$p_adj < 0.05 & abs(rowData(all_da_dod$res)$logFC)>1.5)$cluster_id), k = "flo_merge")
#plotAbundances(x = sig_ter, by = "cluster_id", k = "meta45", group_by = "batch")
sig_dod_cluster_phenotype <- plotExprHeatmap(x = sig_all_dod,
                                            features = refined_markers,
                                            by = "cluster",
                                            row_clust = FALSE,
                                            col_clust = FALSE,
                                            k = "flo_merge",
                                            bars = TRUE,
                                            perc=TRUE,
                                            hm_pal = colorspace::sequential_hcl("inferno", n=8))

png("./figures/sig_dod_cluster_phenotype.png", height=6, width=8, units = "in", res = 600)
ComplexHeatmap::draw(sig_dod_cluster_phenotype, show_heatmap_legend = FALSE)
dev.off()




sig_all_c45 <- filterSCE(sce, cluster_id %in% paste(subset(rowData(all_da_c45$res), rowData(all_da_c45$res)$p_adj < 0.05 & abs(rowData(all_da_c45$res)$logFC)>1)$cluster_id), k = "flo_merge")
#plotAbundances(x = sig_ter, by = "cluster_id", k = "meta45", group_by = "batch")
sig_c45_cluster_phenotype <- plotExprHeatmap(x = sig_all_c45,
                                            features = refined_markers,
                                            by = "cluster",
                                            row_clust = FALSE,
                                            col_clust = FALSE,
                                            k = "flo_merge",
                                            bars = TRUE,
                                            perc=TRUE,
                                            hm_pal = colorspace::sequential_hcl("inferno", n=8))

png("./figures/sig_c45_cluster_phenotype.png", height=6, width=8, units = "in", res = 600)
ComplexHeatmap::draw(sig_c45_cluster_phenotype, show_heatmap_legend = FALSE)
dev.off()


# system.time(sce <- scater::runUMAP(sce,
#                                     subset_row=refined_markers,
#                                     exprs_values = "exprs",
#                                     scale=T))
# 
# 
# big_table <- vac69a.cytof::prep_sce_for_ggplot(sce)
# 
# data.table::fwrite(big_table, "vac63c_all_Tcells_with_UMAP.csv")


library(vac69a.cytof)
library(ggplot2)

big_table <- data.table::fread("vac63c_all_Tcells_with_UMAP.csv")



short_big_table_t6 <- subset(big_table, big_table$timepoint=="T6")

short_big_table_t6 <- short_big_table_t6 %>%
  group_by(n_infection) %>%
  sample_n(5*10^4) %>%
  ungroup()
  

# short_big_table_t6 <- short_big_table_t6[seq(1,nrow(short_big_table_t6), by=3), ]
# 

short_big_table_t6  %>%
  group_by(n_infection) %>%
  summarise(n=n())


# n_infection       n
# <chr>         <int>
# 1 First        663311
# 2 Second       313358
# 3 Third       1016788

for(channel in refined_markers){
  plotplot <- flo_umap(short_big_table_t6, color_by = channel, facet_by = "n_infection")
  ggsave(paste("./figures/umaps/t6_n_infection_umap_", channel, ".png", sep = ""), plotplot,  height=4, width=8)
  print(channel)
}