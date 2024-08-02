library(Libra)

`%notin%` <- Negate(`%in%`)
sce <- readRDS("/Volumes/lab_prasj/BIG_Flo/kassie_bd_rhapsody/combo_sce")
sce$cell_type=ifelse(sce$seurat_clusters %in% c(5, 9, 12), "B Cell",
                           ifelse(sce$seurat_clusters %in% c(11, 17), "Red Cell", "T Cell"))

sce <- sce[,sce$stim!="PBMC"]
d0_irbc <- sce[,sce$stim%in%c("d0", "iRBC")]
d0_urbc <- sce[,sce$stim%in%c("d0", "uRBC")]
irbc_urbc <- sce[,sce$stim%in%c("iRBC", "uRBC")]

d0_irbc_DE <- Libra::run_de(d0_irbc,
             cell_type_col = "cell_type",
             label_col = "stim",
             replicate_col="donor",
             # de_family = "singlecell",
             de_method = "edgeR",
             de_type = "LRT",
             # de_met hod = "wilcox",
             min_reps = 3,
             min_cells = 15)

d0_urbc_DE <- run_de(d0_urbc,
                     cell_type_col = "cell_type",
                     label_col = "stim",
                     replicate_col="donor",
                     # de_family = "singlecell",
                     de_method = "edgeR",
                     de_type = "LRT",
                     # de_method = "wilcox",
                     min_reps = 3,
                     min_cells = 15)

irbc_urbc$cell_type=ifelse(irbc_urbc$seurat_clusters %in% c(5, 9, 12), "B Cell",
                           ifelse(irbc_urbc$seurat_clusters %in% c(11, 17), "Red Cell", "T Cell"))

irbc_urbc_DE <- Libra::run_de(irbc_urbc,
                     cell_type_col = "cell_type",
                     label_col = "stim",
                     replicate_col="sample_id",
                     # de_family = "singlecell",
                     de_method = "edgeR",
                     de_type = "LRT",
                     # de_method = "wilcox",
                     # min_reps = 3,
                     # min_cells = 15
)

sig_irbc_urbc_T <- irbc_urbc_DE %>%
  filter(cell_type=="T Cell", p_val_adj<0.05)

sig_irbc_urbc_B <- irbc_urbc_DE %>%
  filter(cell_type=="B Cell", p_val_adj<0.05)

list_of_DE <- list(d0_irbc_DE, d0_urbc_DE, irbc_urbc_DE)
list_of_sig <- lapply(list_of_DE, function(x) subset(x, x$p_val_adj<0.05 & abs(x$avg_logFC)>2))

"%notin%" <- Negate("%in%")
urbc_markers <- d0_urbc_DE[d0_urbc_DE$gene%notin%d0_irbc_DE,]
sig_urbc_markers <- subset(urbc_markers, p_val_adj<0.05 & cell_type %in% c(1,3))



# scran pseudobulk ####

sum_by <- c("sample_id")
summed <- scuttle::aggregateAcrossCells(irbc_urbc, id=SummarizedExperiment::colData(irbc_urbc)[,sum_by])
colnames(summed) <- apply(SummarizedExperiment::colData(summed)[,sum_by], 1, function(x)paste(x, collapse="_"))

out <- scran::pseudoBulkDGE(summed, 
                            label = summed$sample_id,
                            # vector or factor of length equal to ncol(x), specifying the experimental condition for each column
                            # col.data = SummarizedExperiment::colData(summed),
                            condition = summed$stim,
                            # A formula to be used to construct a design matrix from variables in col.data
                            design = ~stim,
                            coef="irbc"
                            # sig at t6 third infection
                            # contrast = matrix(c(0,1,0,0,0,0,0))
)

library(scran)

info <- DataFrame(sample=sce$sample_id, cluster=sce$seurat_clusters, stim=sce$stim)
pseudo <- sumCountsAcrossCells(sce, info)
# pseudo$stim <- sce$stim[match(pseudo$sample, sce$sample_id)]

out <- pseudoBulkDGE(pseudo, 
                     label=pseudo$cluster,
                     condition=pseudo$stim,
                     design=~stim,
                     contrast=c(0,1),
                     coef="iRBC"
)



# genes that are DEG from d0 to urbc but not DEG from d0 to irbc

sig_d0_urbc_DE <- d0_urbc_DE%>%
  filter(cell_type=="T Cell", p_val_adj<0.05)

sig_d0_irbc_DE <- d0_irbc_DE%>%
  filter(cell_type=="T Cell", p_val_adj<0.05)

non_overlap <- sig_d0_irbc_DE[sig_d0_irbc_DE$gene %notin% sig_d0_urbc_DE$gene,]
non_overlap2 <- sig_d0_urbc_DE[sig_d0_urbc_DE$gene %notin% sig_d0_irbc_DE$gene,]
