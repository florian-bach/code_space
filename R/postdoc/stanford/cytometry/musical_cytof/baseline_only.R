library(CATALYST)
library(dplyr)
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

musical_panel <- read.csv("~/postdoc/stanford/cytometry/CyTOF/MUSICAL/pilot75/musical_panel_edit.csv")

musical_metadata <- read.csv("~/postdoc/stanford/cytometry/CyTOF/MUSICAL/pilot75/single_cell_metadata.csv")
musical_metadata$batch <- as.character(musical_metadata$batch)

# long_file_list <- list.files("/Users/fbach/postdoc/stanford/cytometry/CyTOF/MUSICAL/pilot75/peacoQC/PeacoQC_results/fcs_files/", full.names = TRUE)
# long_file_list <- list.files("/Users/fbach/postdoc/stanford/cytometry/CyTOF/MUSICAL/pilot75/big_fcs/", full.names = TRUE)
# files_to_read <- subset(musical_metadata, subject_id==1217, select=file_path)

metadata_to_read <- musical_metadata %>%
  group_by(subject_id) %>%
  filter(timepoint=="baseline")%>%
  filter(n()==2)

musical_flowset <- ncdfFlow::read.ncdfFlowSet(metadata_to_read$file_path)

# set.seed(1234)
# musical_flowset <- hard_downsample(musical_flowset, event_number = 150000)
# 


sce <- prepData(musical_flowset,
                musical_panel,
                metadata_to_read,
                transform = FALSE,
                md_cols =
                  list(file = "file_name",
                       id = "sample_id",
                       factors = c("subject_id", "timepoint", "class", "batch")))


musical_flowset <- NULL
gc(full = TRUE)

assay(sce, "exprs") <- assay(sce, "counts")

cluster_markers <- musical_panel$antigen[-c(1,2,4:5,40:42, 44:48)]
type_markers <- musical_panel$antigen[musical_panel$class=="type"&!is.na(musical_panel$class)]


cluster_markers <- factor(cluster_markers, levels = c("CD3",
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
                                                      "CCR7"))


sce <- CATALYST::cluster(sce,
                         features = cluster_markers,
                         xdim = 12, ydim = 12, maxK = 50, seed = 1234)

big_heat <- plotExprHeatmap(sce,
                              by="cluster_id",
                              assay = "exprs",
                              fun="median",
                              k="meta45",
                              features=levels(cluster_markers),
                              row_clust = FALSE,
                              col_clust = FALSE)


png("~/postdoc/stanford/cytometry/CyTOF/MUSICAL/figures/baselineonly/da_base_class.png",  width=10, height=5, units = "in", res=400)
ComplexHeatmap::draw(big_heat)
dev.off()

library(diffcyt)

ei <- metadata(sce)$experiment_info

design <- createDesignMatrix(ei, c("class", "subject_id"))


colnames(design)
# [1] "(Intercept)"      "classsymptomatic" "subject_id161"    "subject_id176"    "subject_id268"   
# [6] "subject_id324"    "subject_id353"    "subject_id363"
#timepoint=baseline is dummy; 

base_symp <- createContrast(c(c(0, 1), rep(0,6)))


da_base_symp <- diffcyt(sce,
                     design = design,
                     contrast = base_symp,
                     analysis_type = "DA",
                     method_DA = "diffcyt-DA-edgeR",
                     clustering_to_use = "meta45",
                     verbose = T)

table(rowData(da_base_symp$res)$p_adj < 0.1)

da_base_symp_res <- rowData(da_base_symp$res)
plotDiffHeatmap(sce, da_base_symp_res, all=TRUE)


library(tidyr)
library(dplyr)
library(ggplot2)

## base vs day0####
results_table <- data.frame(diffcyt::topTable(da_base_symp, all=T, show_counts = TRUE, show_props = TRUE))

long_results <- results_table %>%
  pivot_longer(cols = starts_with("props"), names_to = "sample_prop", values_to = "freq")%>%
  pivot_longer(cols = starts_with("counts"), names_to = "sample_counts", values_to = "count")%>%
  mutate("subject_id"=substr(sample_prop, 14, 16),
         "timepoint"=substr(sample_prop, nchar(sample_prop)-5, nchar(sample_prop)),
         "class"=substr(sample_prop, 27, 38))%>%
  mutate(timepoint=case_when(timepoint=="seline"~"baseline",
                             timepoint=="_day.0"~ "day 0",
                             timepoint=="day.14"~ "day 14",
                             timepoint=="_day.7"~ "day 7"))


symp_vs_time <- long_results %>%
  # filter(p_adj<0.1)%>%
  mutate(age_class=if_else(subject_id %in% names(adults_palette), "adult", "child"))%>%
  ggplot(aes(x=class, y=freq/100, fill=class))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(color=subject_id))+
  # geom_line(aes(color=subject_id, group=subject_id))+
  ggtitle("Differntially abundant baseline vs. day 0, all symptomatic individuals")+
  # geom_boxplot(aes(fill=age_class), outlier.shape = NA)+
  facet_wrap(~cluster_id, scales="free", nrow=9)+
  scale_y_continuous(labels = scales::label_percent())+
  xlab("")+
  ylab("")+
  theme_minimal()+
  scale_color_manual(values=c(kids_palette, adults_palette))+
  scale_fill_manual(values=c("darkred", "darkblue"))+
  theme(legend.position = "none")

ggsave("~/postdoc/stanford/cytometry/CyTOF/MUSICAL/figures/baselineonly/da_base_class.png", symp_vs_time, height=15, width=12, dpi=400, bg="white")
