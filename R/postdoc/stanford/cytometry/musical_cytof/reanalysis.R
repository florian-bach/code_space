# CATALYST ####
library(CATALYST)
library(tidyr)
library(dplyr)
library(ggplot2)
library(diffcyt)
# data generation ####
`%notin%` <- Negate(`%in%`)

kids_palette <- colorspace::sequential_hcl(8, palette = "Emrld")
names(kids_palette) <- as.character(c(268, 324, 137, 176, 353, 161, 363, 571))

adults_palette <- colorspace::sequential_hcl(5, palette = "OrRd")[1:4]
names(adults_palette) <- as.character(c(575, 627, 563, 643))


hard_downsample <- function(fs, event_number){
  flowCore::fsApply(fs, function(ff){
    idx <- sample.int(nrow(ff), min(event_number, nrow(ff)))
    ff[idx,]
  })
}

musical_panel <- read.csv("~/postdoc/stanford/cytometry/CyTOF/MUSICAL/pilot75/musical_panel_edit.csv")[,2:4]

musical_metadata <- read.csv("~/postdoc/stanford/cytometry/CyTOF/MUSICAL/pilot75/single_cell_metadata.csv")
musical_metadata$batch <- as.character(musical_metadata$batch)

# kids only
metadata_to_read <- filter(musical_metadata, subject_id %in% c(268, 324, 137, 176, 353, 161, 363, 571))

musical_flowset <- ncdfFlow::read.ncdfFlowSet(metadata_to_read$file_path)
set.seed(1234)
musical_flowset <- hard_downsample(musical_flowset, event_number = 10000)
# 
# 

sce <- prepData(musical_flowset,
                musical_panel,
                metadata_to_read,
                transform = FALSE,
                md_cols =
                  list(file = "file_name",
                       id = "sample_id",
                       factors = c("subject_id", "timepoint", "class", "batch")),
                panel_cols =
                  list(channel = "fcs_colname",
                       antigen = "antigen",
                       class = "class"))


musical_flowset <- NULL
gc(full = TRUE)

assay(sce, "exprs") <- assay(sce, "counts")

# p <- plotExprs(sce, color_by = "batch")
# p$facet$params$ncol <- 6
# p

cluster_markers <- musical_panel$antigen[-c(1,2,4:5,40:42, 44:48)]
type_markers <- musical_panel$antigen[musical_panel$class=="type"&!is.na(musical_panel$class)]


cluster_markers <- c("CD3",
                     "CD4",
                     "CD8",
                     "TCRgd",
                     "Vd2TCR",
                     "CD161",
                     "CD56",
                     "CD94",
                     "CD19",
                     "CD20",
                     "CD21",
                     "IgD",
                     "CD38",
                     "CD27",
                     "CD25",
                     "CD127",
                     "CD39",
                     "PD-1",
                     "CXCR5",
                     "ICOS",
                     "CCR6",
                     "CXCR3",
                     "CD49b",
                     "HLADR",
                     "CD14",
                     "CD16",
                     "CD85j",
                     "CD11b",
                     "CD11c",
                     "CD123",
                     "CD33",
                     "CD86",
                     "LAG3",
                     "CD57",
                     "CD45RA",
                     "CCR7")

# sce <- filterSCE(sce, batch!=41824)

sce <- CATALYST::cluster(sce,
                         features = cluster_markers,
                         xdim = 12, ydim = 12, maxK = 50, seed = 1234)

# phenotypic heatmaps & UMAPs####

unclean_heatmap <- plotExprHeatmap(sce,
                by="cluster_id",
                assay = "exprs",
                fun="median",
                k="meta50",
                row_anno = TRUE,
                bars = T,
                features=cluster_markers,
                row_clust = FALSE,
                col_clust = FALSE)


png("~/postdoc/stanford/cytometry/CyTOF/MUSICAL/figures/reanalysis/unordered_unclean_cluster_heat.png", width=10, height=10, units = "in", res=400)
ComplexHeatmap::draw(unclean_heatmap)
dev.off()

# mixed_pops <- filterSCE(sce, cluster_id %in% c(30, 19), k="meta50")
# 
# mixed_pops_heatmap <- plotExprHeatmap(mixed_pops,
#                                    by="cluster_id",
#                                    assay = "exprs",
#                                    fun="median",
#                                    k="som144",
#                                    m="meta50",
#                                    row_anno = TRUE,
#                                    features=cluster_markers,
#                                    row_clust = FALSE,
#                                    col_clust = FALSE)
# 
# 
# png("~/postdoc/stanford/cytometry/CyTOF/MUSICAL/figures/reanalysis/mixed_pops_heat.png", width=10, height=2, units = "in", res=400)
# ComplexHeatmap::draw(mixed_pops_heatmap)
# dev.off()

## relabel clusters ####
rough_merge <- read.csv("~/postdoc/stanford/cytometry/CyTOF/MUSICAL/pilot75/rough_activated_50_merge.csv", header = TRUE)

rough_merge$new <- gsub("  ", " ", rough_merge$new)
rough_merge$new <- gsub("^ ", "", rough_merge$new)
rough_merge$new <- gsub(" $", "", rough_merge$new)

sce <- mergeClusters(sce, k="meta50", rough_merge, "rough_merge", overwrite = T)
saveRDS(sce, "~/postdoc/stanford/cytometry/CyTOF/MUSICAL/pilot75/new_sce.RDS") 
sort(unique(rough_merge$new))
rough_merge_palette <- c(
  colorspace::sequential_hcl(n=4, "Peach")[1:2],
  colorspace::sequential_hcl(n=11, "BluGrn"),
  colorspace::sequential_hcl(n=4, "Peach")[3],
  "red",
  colorspace::sequential_hcl(n=6, "Purp")[1:5],
  "darkred",
  "yellow",
  colorspace::sequential_hcl(n=4, "Peach")[4],
  "grey",
  "#380000",
  "orange")

  
# names(rough_merge_palette) <- na.omit(lineage_df$new)
  
clean_heat <- plotExprHeatmap(sce,
                            by="cluster_id",
                            assay = "exprs",
                            fun="median",
                            k="rough_merge",
                            # row_anno = TRUE,
                            features=cluster_markers,
                            row_clust = TRUE,
                            col_clust = FALSE)

png("~/postdoc/stanford/cytometry/CyTOF/MUSICAL/figures/reanalysis/clean_heat.png", width=10, height=10, units = "in", res=400)
ComplexHeatmap::draw(clean_heat)
dev.off()

system.time(sce <- runDR(sce,
                         seed=1234,
                         dr="UMAP",
                         cells=10000,
                         features=cluster_markers)
)

subset_umap <- plotDR(sce,  color_by = "rough_merge", scale = TRUE, k_pal = rough_merge_palette)
# plotDR(sce,  color_by = "Dead", scale = TRUE)

marker_umap <- plotDR(sce,  color_by = sort(type_markers), scale = TRUE)
marker_umap2 <- plotDR(sce,  color_by = sort(cluster_markers), scale = TRUE)

ggplot2::ggsave("~/postdoc/stanford/cytometry/CyTOF/MUSICAL/figures/reanalysis/type_markers_umap.png", marker_umap, height = 12, width=12, bg="white")
ggplot2::ggsave("~/postdoc/stanford/cytometry/CyTOF/MUSICAL/figures/reanalysis/all_markers_umap.png", marker_umap2, height = 12, width=12, bg="white")

ggplot2::ggsave("~/postdoc/stanford/cytometry/CyTOF/MUSICAL/figures/reanalysis/overview_umap.png", subset_umap, height = 12, width=12, bg="white")


# plot cluster abundances ####

cell_matrix <- data.frame(colData(sce))
hash_map <- data.frame("som144"=cluster_ids(sce, k="som144"),
                       "rough_merge"=cluster_ids(sce, k="rough_merge"))
hash_map <- hash_map[!duplicated(hash_map$som144),]

cell_matrix$cluster_id <- hash_map$rough_merge[match(cell_matrix$cluster_id, hash_map$som144)]

cell_freqs <- cell_matrix %>%
  reframe(cell_count = n(), .by=c(subject_id, timepoint, class, cluster_id))%>%
  mutate(cell_freq=cell_count/10000,
         timepoint=factor(timepoint, levels=c("baseline", "day 0", "day 7", "day 14")))

big_box <- cell_freqs%>%
  # filter(grepl("activated", cluster_id))%>%
  ggplot(., aes(x=timepoint, y=cell_freq, fill=class))+
  geom_boxplot()+
  scale_y_continuous(labels = scales::label_percent())+
  facet_wrap(~cluster_id, scales = "free")+
  theme_minimal()+
  scale_fill_manual(values = viridis::magma(n=3))

ggplot2::ggsave("~/postdoc/stanford/cytometry/CyTOF/MUSICAL/figures/reanalysis/big_box.png", big_box, height = 12, width=18, bg="white")

big_line <- cell_freqs%>%
  # filter(grepl("activated", cluster_id))%>%
  ggplot(., aes(x=timepoint, y=cell_freq, color=class, group=interaction(subject_id, class)))+
  geom_point()+
  geom_line()+
  scale_y_continuous(labels = scales::label_percent())+
  theme_minimal()+
  facet_wrap(~cluster_id, scales = "free")+
  scale_color_manual(values = viridis::magma(n=3))

ggplot2::ggsave("~/postdoc/stanford/cytometry/CyTOF/MUSICAL/figures/reanalysis/big_line.png", big_line, height = 12, width=18, bg="white")


base_zero_contrast <- t(matrix(c(rep(0,12), 1, 0)))
# base_seven_contrast <- t(matrix(c(rep(0,12), 0, 1,0)))
library(purrr)
cell_freq_purrf <- cell_freqs %>%
  group_by(cluster_id)%>%
  mutate(subject_id=as.character(subject_id))%>%
  nest()%>%
  mutate(model=map(data, ~lm(cell_freq~timepoint*class+subject_id, data=.)))%>%
  mutate(summary=map(model, ~summary(.))) %>%
  mutate(base_zero_p=map_dbl(summary, ~coef(.)[55])) %>%
  mutate(base_14_p=map_dbl(summary, ~coef(.)[56])) %>%
  
  ungroup()%>%
  mutate(base_zero_padj=p.adjust(base_zero_p, method="BH"),
         base_14_padj=p.adjust(base_14_p, method="BH"))
  
sig_freq <- cell_freq_purrf %>%
  filter(base_zero_p < 0.1)


sig_box <- cell_freqs%>%
  filter(cluster_id %in% sig_freq$cluster_id)%>%
  filter(timepoint %in% c("baseline", "day 0", "day 7", "day 14"))%>%
  ggplot(., aes(x=factor(timepoint, levels=c("baseline", "day 0", "day 7", "day 14")), y=cell_freq, fill=class))+
  geom_boxplot(outliers = F)+
  scale_y_continuous(labels = scales::label_percent())+
  facet_wrap(~cluster_id, scales = "free")+
  theme_minimal()+
  scale_fill_manual(values = viridis::magma(n=3))+
  theme(axis.title.x = element_blank())

ggplot2::ggsave("~/postdoc/stanford/cytometry/CyTOF/MUSICAL/figures/reanalysis/sig_box.png", sig_box, height = 8, width=8, bg="white")

