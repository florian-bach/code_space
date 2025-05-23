# CATALYST ####
library(CATALYST)

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

metadata_to_read <- subset(musical_metadata, class == "asymptomatic")
metadata_to_read$num_timepoint <- as.numeric(factor(metadata_to_read$timepoint))

musical_flowset <- ncdfFlow::read.ncdfFlowSet(metadata_to_read$file_path)

musical_flowset <- hard_downsample(musical_flowset, event_number = 150000)



sce <- prepData(musical_flowset,
                musical_panel,
                metadata_to_read,
                transform = FALSE,
                md_cols =
                  list(file = "file_name",
                       id = "sample_id",
                       factors = c("subject_id", "timepoint", "class", "batch", "num_timepoint")))


musical_flowset <- NULL
gc(full = TRUE)

assay(sce, "exprs") <- assay(sce, "counts")
# p <- plotExprs(sce, color_by = "batch", features = "Dead")
# p$facet$params$ncol <- 6
# p

cluster_markers <- musical_panel$antigen[-c(1,2,4:5,40:42, 44:48)]
type_markers <- musical_panel$antigen[musical_panel$class=="type"&!is.na(musical_panel$class)]

# sce <- filterSCE(sce, batch!=41824)
set.seed(1234)
sce <- CATALYST::cluster(sce,
                         features = cluster_markers,
                         xdim = 12, ydim = 12, maxK = 50)
#remove putative doublets
# sce <- filterSCE(sce, cluster_id %notin% c(98, 109, 110, 143))
# plotCounts(sce, group_by = "sample_id", color_by = "class")



big_heat <- plotMultiHeatmap(sce,
                             k = "meta45",
                             # m="meta30",
                             hm1=cluster_markers,
                             hm2 = "abundances",
                             row_clust = TRUE,
                             col_clust = TRUE)

png("~/postdoc/stanford/cytometry/CyTOF/MUSICAL/redownload_pilot75/redownload_big_fcs/single_cells/figures/asymptomatic_only_big_heat.png", width=20, height=12, units = "in", res=400)
ComplexHeatmap::draw(big_heat)
dev.off()

big_heat <- plotMultiHeatmap(sce,
                             k = "meta45",
                             # m= "meta30",
                             hm1=cluster_markers,
                             bars = TRUE,
                             hm2 = "abundances",
                             row_clust = FALSE,
                             col_clust = TRUE)

png("~/postdoc/stanford/cytometry/CyTOF/MUSICAL/redownload_pilot75/redownload_big_fcs/single_cells/unordered_asymptomatic_only_peacoq_big_heat.png", width=20, height=12, units = "in", res=400)
ComplexHeatmap::draw(big_heat)
dev.off()

# big_heat <- plotMultiHeatmap(sce,
#                              m = "meta45",
#                              k= "som144",
#                              hm1=cluster_markers,
#                              hm2 = "abundances",
#                              row_clust = FALSE,
#                              col_clust = TRUE)
# 
# png("~/postdoc/stanford/cytometry/CyTOF/MUSICAL/redownload_pilot75/sandbox_out/peacoQC/huge_unordered_asymptomatic_only_peacoq_big_heat.png", width=20, height=20, units = "in", res=400)
# ComplexHeatmap::draw(big_heat)
# dev.off()

# 73s for 5000 cells / sample
# 162s for 10000 cells / sample
system.time(sce <- runDR(sce,
                         dr="UMAP",
                         cells=10000,
                         features=cluster_markers)
)

time_umap <- plotDR(sce,  color_by = "timepoint", facet_by = "subject_id", scale = TRUE)
# plotDR(sce,  color_by = "Dead", scale = TRUE)

marker_umap <- plotDR(sce,  color_by = type_markers, scale = TRUE)

ggplot2::ggsave("~/postdoc/stanford/cytometry/CyTOF/MUSICAL/figures/asymptomatic_only/marker_umap.png", marker_umap, height = 12, width=12, bg="white")
ggplot2::ggsave("~/postdoc/stanford/cytometry/CyTOF/MUSICAL/figures/asymptomatic_only/time_individual_umap.png", time_umap, height = 12, width=12, bg="white")

# filter dead ####
# sce <- filterSCE(sce, cluster_id %notin% c(1, 2, 4, 9, 13, 22, 28, 31, 37), k="meta45")
# sce <- CATALYST::cluster(sce, features = cluster_markers, xdim = 12, ydim = 12, maxK = 50)

# differential abundance, tracking cluster percentages through time  ####
library(diffcyt)

ei <- metadata(sce)$experiment_info
ei$age_class <- ifelse(ei$subject_id %in% c(268, 324, 137, 176, 353, 161, 363, 571), "child", "adult")
ei$num_timepoint <- as.numeric(ei$num_timepoint)

FDR_cutoff <- 0.05

design <- createDesignMatrix(ei, c("num_timepoint", "subject_id"))


colnames(design)
# [1] "(Intercept)"   "num_timepoint" "subject_id161" "subject_id176" "subject_id268" "subject_id324" "subject_id353"
# [8] "subject_id363" "subject_id571"



base_0_age <-    createContrast(c(c(0, 1), rep(0,7)))
# base_14_age <- createContrast(c(c(0, 0, 1, 0), rep(0,6)))
# zero_14_o_age <-  createContrast(c(c(0, 1, -1, 0), rep(0,6)))

da_base_0_age <- diffcyt(sce,
                         design = design,
                         contrast = base_0_age,
                         analysis_type = "DA",
                         method_DA = "diffcyt-DA-edgeR",
                         clustering_to_use = "meta45",
                         verbose = T)

da_base_14_age <- diffcyt(sce,
                          design = design,
                          contrast = zero_14_o_age,
                          analysis_type = "DA",
                          method_DA = "diffcyt-DA-edgeR",
                          clustering_to_use = "meta45",
                          verbose = T)


da_zero_14_o_age <- diffcyt(sce,
                          design = design,
                          contrast = zero_14_o_age,
                          analysis_type = "DA",
                          method_DA = "diffcyt-DA-edgeR",
                          clustering_to_use = "meta45",
                          verbose = T)






table(rowData(da_base_0_age$res)$p_adj < FDR_cutoff)
# FALSE 
# 45 


table(rowData(da_base_14_age$res)$p_adj < FDR_cutoff)
# FALSE
#    45
table(rowData(da_zero_14_o_age$res)$p_adj < FDR_cutoff)
# FALSE 
# 45

## visualise diffcyt results ####
library(tidyr)
library(dplyr)
library(ggplot2)

## base vs day0####
results_table <- data.frame(diffcyt::topTable(da_base_0_age, all=T, show_counts = TRUE, show_props = TRUE))

long_results <- results_table %>%
  pivot_longer(cols = starts_with("props"), names_to = "sample_prop", values_to = "freq")%>%
  pivot_longer(cols = starts_with("counts"), names_to = "sample_counts", values_to = "count")%>%
  mutate("subject_id"=substr(sample_prop, 14, 16),
         "timepoint"=substr(sample_prop, nchar(sample_prop)-5, nchar(sample_prop)),
         "class"=substr(sample_prop, 24, 30))%>%
  mutate(timepoint=case_when(timepoint=="seline"~"baseline",
                             timepoint=="_day.0"~ "day 0",
                             timepoint=="day.14"~ "day 14",
                             timepoint=="_day.7"~ "day 7"))


symp_vs_time <- long_results %>%
  filter(p_adj<0.5, timepoint != "day 7")%>%
  mutate(age_class=if_else(subject_id %in% names(adults_palette), "adult", "child"))%>%
  ggplot(aes(x=factor(timepoint, levels=c("baseline", "day 0", "day 7", "day 14")), y=freq/100))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(color=subject_id))+
  geom_line(aes(color=subject_id, group=subject_id))+
  ggtitle("Differntially abundant baseline vs. day 0, all asymptomatic individuals")+
  facet_wrap(~cluster_id, scales="free")+
  scale_y_continuous(labels = scales::label_percent())+
  xlab("")+
  ylab("")+
  theme_minimal()+
  scale_color_manual(values=c(kids_palette, adults_palette))

ggsave("~/postdoc/stanford/cytometry/CyTOF/MUSICAL/figures/asymptomatic_only/da_base_0_all_asymp.png", symp_vs_time, height=8, width=12, dpi=400, bg="white")



## day 0 v day 14 ####
results_table <- data.frame(diffcyt::topTable(da_zero_14_o_age, all=T, show_counts = TRUE, show_props = TRUE))

long_results <- results_table %>%
  pivot_longer(cols = starts_with("props"), names_to = "sample_prop", values_to = "freq")%>%
  pivot_longer(cols = starts_with("counts"), names_to = "sample_counts", values_to = "count")%>%
  mutate("subject_id"=substr(sample_prop, 14, 16),
         "timepoint"=substr(sample_prop, nchar(sample_prop)-5, nchar(sample_prop)),
         "class"=substr(sample_prop, 24, 30))%>%
  mutate(timepoint=case_when(timepoint=="seline"~"baseline",
                             timepoint=="_day.0"~ "day 0",
                             timepoint=="day.14"~ "day 14",
                             timepoint=="_day.7"~ "day 7"))

symp_vs_time <- long_results %>%
  filter(p_adj<0.05, timepoint != "day 7")%>%
  mutate(age_class=if_else(subject_id %in% names(adults_palette), "adult", "child"))%>%
  ggplot(aes(x=factor(timepoint, levels=c("baseline", "day 0", "day 7", "day 14")), y=freq/100))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(color=subject_id))+
  geom_line(aes(color=subject_id, group=subject_id))+
  ggtitle("Differentially abundant day 0 vs. day 14, all asymptomatic individuals")+
  facet_wrap(~cluster_id, scales="free")+
  scale_y_continuous(labels = scales::label_percent())+
  xlab("")+
  ylab("")+
  theme_minimal()+
  scale_color_manual(values=c(kids_palette, adults_palette))


ggsave("~/postdoc/stanford/cytometry/CyTOF/MUSICAL/figures/asymptomatic_only/da_zero_14_all_asymp.png", symp_vs_time,  height=8, width=12, dpi=400, bg="white")




## day baseline vs individual ####

# how accesory files were made #### 
## make new metadata ####
file_list <- list.files("/Users/fbach/postdoc/stanford/cytometry/CyTOF/MUSICAL/redownload_pilot75/redownload_big_fcs/")
long_file_list <- list.files("/Users/fbach/postdoc/stanford/cytometry/CyTOF/MUSICAL/redownload_pilot75/redownload_big_fcs/", full.names = TRUE)

batch <- paste("batch_", substr(file_list, 0, 6), sep="")


short_file_names <- substr(file_list, 39, nchar(file_list))
short_file_names <- gsub("final-", "", short_file_names)

subject_id <- substr(short_file_names, 1, 4)
subject_id <- gsub("_", "", subject_id)

timepoint <- substr(short_file_names, 7, 11)
timepoint <- dplyr::case_when(timepoint=="14_01"~"day 14",
                              timepoint=="7_01_"~"day 7",
                              timepoint=="0_01_"~"day 0",
                              timepoint=="basel"~"baseline",
                              timepoint=="mune_"~"immune",
                              timepoint=="IMCct"~"control")

class <- substr(short_file_names, 4, 6)

class <- dplyr::case_when(class=="7_H"~"control",
                          class=="_A_"~"asympmtomatic",
                          class=="_S_"~"symptomatic",
                          class=="_im"~"immune")

cytof_metadata <- data.frame("batch"=as.character(batch),
                             "subject_id"=subject_id,
                             "timepoint"=timepoint,
                             "class"=class,
                             "file_name"=file_list,
                             "file_path"=long_file_list,
                             "sample_id"=paste(batch, subject_id, class, timepoint, sep="_"))

write.csv(cytof_metadata, "~/postdoc/stanford/cytometry/CyTOF/MUSICAL/redownload_pilot75/redownload_raw_metadata.csv", row.names = FALSE)

## make old metadata ####
file_list <- list.files("~/postdoc/stanford/cytometry/CyTOF/MUSICAL/redownload_pilot75/redownload_big_fcs/single_cells/", pattern = ".fcs$")
long_file_list <- list.files("~/postdoc/stanford/cytometry/CyTOF/MUSICAL/redownload_pilot75/redownload_big_fcs/single_cells/", pattern = ".fcs$", full.names = TRUE)

batch <- substr(file_list, 0, 6)


short_file_names <- substr(file_list, 39, nchar(file_list))
short_file_names <- gsub("final-", "", short_file_names)

subject_id <- substr(short_file_names, 1, 4)
subject_id <- gsub("_", "", subject_id)

timepoint <- substr(short_file_names, 7, 11)
timepoint <- dplyr::case_when(timepoint=="14_01"~"day 14",
                              timepoint=="7_01_"~"day 7",
                              timepoint=="0_01_"~"day 0",
                              timepoint=="basel"~"baseline",
                              timepoint=="mune_"~"immune",
                              timepoint=="IMCct"~"control")

class <- substr(short_file_names, 4, 6)

class <- dplyr::case_when(class=="7_H"~"control",
                          class=="_A_"~"asympmtomatic",
                          class=="_S_"~"symptomatic",
                          class=="_im"~"immune")

cytof_metadata <- data.frame("batch"=as.character(batch),
                             "subject_id"=subject_id,
                             "timepoint"=timepoint,
                             "class"=class,
                             "file_name"=file_list,
                             "file_path"=long_file_list,
                             "sample_id"=paste(batch, subject_id, timepoint, sep="_"))

write.csv(cytof_metadata, "~/postdoc/stanford/cytometry/CyTOF/MUSICAL/pilot75/single_cell_metadata.csv", row.names = FALSE)


## get panel ####
raw_musical_panel <- data.frame(premessa::read_parameters(long_file_list[1]))

musical_panel <- data.frame("fcs_colname"=rownames(raw_musical_panel),
                            "antigen"= substr(raw_musical_panel[,1], 7, 12))

musical_panel$antigen <- ifelse(nchar(musical_panel$antigen)<3, musical_panel$fcs_colname, musical_panel$antigen)
musical_panel$antigen <- gsub("length", "Event_length", musical_panel$antigen)

musical_panel$class <- "type"

write.csv(musical_panel, "~/postdoc/stanford/cytometry/CyTOF/MUSICAL/pilot75/musical_panel.csv", row.names = FALSE)

# scatterplots ####


b_t_doubles_all <- CATALYST::plotScatter(sce, chs = c("CD3", "IgD"),facet_by = "cluster_id",  k="meta45")
b_t_doubles_select <- CATALYST::plotScatter(doublets, chs = c("CD3", "IgD"), facet_by = "cluster_id", k="meta45")

ggsave("~/postdoc/stanford/cytometry/CyTOF/MUSICAL/figures/asymptomatic_only/b_t_doubles_all.png", b_t_doubles_all,  height=32, width=8, dpi=400, bg="white")
ggsave("~/postdoc/stanford/cytometry/CyTOF/MUSICAL/figures/asymptomatic_only/b_t_doubles_select.png", b_t_doubles_select,  height=6, width=6, dpi=400, bg="white")

b_mono_doubles_all <- CATALYST::plotScatter(sce, chs = c("CD14", "CD19"), facet_by = "cluster_id", k="meta45")
b_mono_doubles_select <- CATALYST::plotScatter(doublets, chs = c("CD14", "CD19"), facet_by = "cluster_id", color_by = "subject_id")

ggsave("~/postdoc/stanford/cytometry/CyTOF/MUSICAL/figures/asymptomatic_only/b_mono_doubles_all.png", b_mono_doubles_all,  height=32, width=8, dpi=400, bg="white")
ggsave("~/postdoc/stanford/cytometry/CyTOF/MUSICAL/figures/asymptomatic_only/b_mono_doubles_select.png", b_mono_doubles_select,  height=6, width=6, dpi=400, bg="white")
