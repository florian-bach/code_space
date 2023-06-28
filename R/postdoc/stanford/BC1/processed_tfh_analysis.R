library(CATALYST)
library(dplyr)
library(tidyr)


fcs_dir <- "~/postdoc/stanford/clinical_data/BC1/tfh_data/processed_data/"
#fcs_dir <- "~/postdoc/stanford/clinical_data/BC1/tfh_data/raw_tfh_batch1/Preprocessed/"

fcs_file_paths <- list.files(fcs_dir, pattern="*.fcs", full.names = TRUE, recursive = FALSE)[-54]
fcs_file_names <- list.files(fcs_dir, pattern="*.fcs", full.names = FALSE, recursive = FALSE)[-54]

md <- data.frame("file_paths"=fcs_file_paths,
                 "file_names"=fcs_file_names,
                 "sample_id"=substr(fcs_file_names, regexpr("B1 ", fcs_file_names)+3, regexpr("B1 ", fcs_file_names)+6))

cleaned_tfh <- read.csv("~/postdoc/stanford/clinical_data/BC1/tfh_data/cTfh_all.csv")
colnames(cleaned_tfh) <- c("idt", "tfh_q4", "tfh_q3","tfh_q2","tfh_q1","tfh_of_cd4", "parent2t", "tfh_drop")
drop_removed_tfh <- subset(cleaned_tfh, tfh_drop=="No")

md <- subset(md, file_names %in% drop_removed_tfh$idt)


combo_data <- read.csv("~/postdoc/stanford/clinical_data/BC1/antibody_infs.csv")

one_year <- combo_data %>% filter(timepoint==3)%>%
  pivot_wider(names_from = antigen, values_from =  conc)%>%
  mutate(sample_id = substr(id, 2,6))%>%
  filter(sample_id %in% md$sample_id)%>%
  select(-id)

md <- subset(md, sample_id %in% one_year$sample_id)

md <- left_join(md, one_year, by="sample_id")

tfh_markers <- c("CD3", "CD4", "CD45RA", "CXCR5", "PD1", "CXCR3", "CCR4", "CCR6")


premessa_table <- premessa::read_parameters(list.files(fcs_dir, pattern="*.fcs", full.names = TRUE, recursive = FALSE))
panel <- data.frame("fcs_colname"=rownames(premessa_table[1]),
                    "marker_name"=premessa_table[[1]])

panel$marker_class <- ifelse(panel$marker_name %in% tfh_markers, "type", "state")
panel$antigen <- ifelse(panel$marker_name=="", panel$fcs_colname, panel$marker_name)
panel$marker_name <- ifelse(panel$marker_name=="", panel$fcs_colname, panel$marker_name)


processed_tfhs <- flowCore::read.flowSet(md$file_paths)

sce <- CATALYST::prepData(processed_tfhs, panel, md, FACS = TRUE,
                          fix_chs="common",
                          #features=tfh_markers,
                          md_cols = list(file = "file_names",
                                         id = "sample_id",
                                         factors = colnames(md)[-c(2,3)]
                          ),
                          panel_cols = list(channel = "fcs_colname",
                                            antigen = "marker_name",
                                            class = "marker_class")
)

# clustering ####


set.seed(1234);sce <- CATALYST::cluster(sce, features = tfh_markers, xdim = 10, ydim = 10, maxK = 30)

sce_pheno <-CATALYST::plotExprHeatmap(sce, by="cluster_id",
                                       k="meta24",
                                       features = tfh_markers,
                                       perc = TRUE,
                                       bars = TRUE,
                                       # row_anno = TRUE,
                                       # #col_anno = TRUE,
                                       row_clust = FALSE,
                                       col_clust = FALSE,
                                       # row_dend = TRUE,
                                       # #col_dend = TRUE,
                                       #hm_pal = colorspace::sequential_hcl("inferno", n=8)
                                       )

png("~/postdoc/stanford/clinical_data/BC1/tfh_data/figures/batch1_peacoqd_small_heat_dropped.png", width=8, height=8, units = "in", res=400)
ComplexHeatmap::draw(sce_pheno)
dev.off()



cd3_cd19_cells_only <- CATALYST::filterSCE(sce, k="meta24", cluster_id %in% c(1, 15, 18:20))
set.seed(1234);cd3_cd19_cells_only <- CATALYST::cluster(cd3_cd19_cells_only, features = tfh_markers, xdim = 10, ydim = 10, maxK = 30)



cd3_cd19_pheno <-CATALYST::plotExprHeatmap(cd3_cd19_cells_only, by="cluster_id",
                                      k="meta24",
                                      features = tfh_markers,
                                      perc = TRUE,
                                      bars = TRUE,
                                      scale = "last",
                                      # row_anno = TRUE,
                                      # #col_anno = TRUE,
                                      row_clust = FALSE,
                                      col_clust = FALSE,
                                      # row_dend = TRUE,
                                      # #col_dend = TRUE,
                                      #hm_pal = colorspace::sequential_hcl("inferno", n=8)
)

png("~/postdoc/stanford/clinical_data/BC1/tfh_data/figures/cd3_cd19_pheno.png", width=8, height=8, units = "in", res=400)
ComplexHeatmap::draw(cd3_cd19_pheno)
dev.off()





cd3_cd19_multi <-  CATALYST::plotMultiHeatmap(cd3_cd19_cells_only,
                                        k="meta24",
                                        hm2 = "abundances",
                                        hm1 = "type",
                                        #features = tfh_markers,
                                        perc = TRUE,
                                        bars = TRUE,
                                        #row_anno = TRUE,
                                        col_anno = FALSE,
                                        row_clust = FALSE,
                                        col_clust = TRUE,
                                        #row_dend = TRUE,
                                        #col_dend = TRUE,
                                        #hm1_pal = colorspace::sequential_hcl("inferno", n=8)
                                        )

png("~/postdoc/stanford/clinical_data/BC1/tfh_data/figures/batch1_peacoqd_cd4t_cells_only_reclustered.png", width=24, height=8, units = "in", res=400)
ComplexHeatmap::draw(cd3_cd19_multi)
dev.off()

all_cd4_only <- CATALYST::filterSCE(cd3_cd19_cells_only, k="meta15", cluster_id %in% c(11:24))
set.seed(1234);all_cd4_only <- CATALYST::cluster(all_cd4_only, features = tfh_markers[-c(1:2)], xdim = 10, ydim = 10, maxK = 30)


all_cd4_pheno <-CATALYST::plotExprHeatmap(all_cd4_only, by="cluster_id",
                                           #m="meta10",
                                           k="meta10",
                                           features = tfh_markers,
                                           perc = TRUE,
                                           bars = TRUE,
                                           #scale = "never",
                                           # row_anno = TRUE,
                                           # #col_anno = TRUE,
                                           row_clust = FALSE,
                                           col_clust = FALSE,
                                           # row_dend = TRUE,
                                           # #col_dend = TRUE,
                                           #hm_pal = colorspace::sequential_hcl("inferno", n=8)
)

png("~/postdoc/stanford/clinical_data/BC1/tfh_data/figures/all_cd4_pheno.png", width=8, height=8, units = "in", res=400)
ComplexHeatmap::draw(all_cd4_pheno)
dev.off()

saveRDS(all_cd4_only, "~/postdoc/stanford/clinical_data/BC1/tfh_data/all_cd4_only.RDS")


# tfh subsetting ####
# 
# 
# 
# 
# 
# clean_tfh_only <- CATALYST::filterSCE(cd3_cd19_cells_only, k="meta15", cluster_id %in% c(9, 13))
# set.seed(1234);clean_tfh_only <- CATALYST::cluster(clean_tfh_only, features = tfh_markers[-c(1:4)], xdim = 10, ydim = 10, maxK = 30)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# tfh_multi <-  CATALYST::plotMultiHeatmap(clean_tfh_only,
#                                         k="meta6",
#                                         hm2 = "abundances",
#                                         hm1="type",
#                                         #features = tfh_markers,
#                                         perc = TRUE,
#                                         bars = TRUE,
#                                         #row_anno = TRUE,
#                                         col_anno = FALSE,
#                                         row_clust = FALSE,
#                                         col_clust = TRUE,
#                                         #row_dend = TRUE,
#                                         #col_dend = TRUE,
#                                         #hm1_pal = colorspace::sequential_hcl("inferno", n=8)
#                                         )
# 
# png("~/postdoc/stanford/clinical_data/BC1/tfh_data/figures/batch1_peacoqd_clean_tfh_only.png", width=24, height=8, units = "in", res=400)
# ComplexHeatmap::draw(tfh_multi)
# dev.off()
# 
# 
# 
# 
# tfh_pheno <-  CATALYST::plotExprHeatmap(clean_tfh_only,
#                                         by="cluster_id",
#                                         k="meta8",
#                                         features = tfh_markers[-c(1:4)],
#                                         perc = TRUE,
#                                         bars = TRUE,
#                                         row_anno = FALSE,
#                                         col_anno = FALSE,
#                                         row_clust = FALSE,
#                                         col_clust = TRUE,
#                                         #row_dend = TRUE,
#                                         #col_dend = TRUE,
#                                         #hm1_pal = colorspace::sequential_hcl("inferno", n=8)
# )
# 
# png("~/postdoc/stanford/clinical_data/BC1/tfh_data/figures/tfh_only_pheno_heat.png", width=24, height=8, units = "in", res=400)
# ComplexHeatmap::draw(tfh_pheno)
# dev.off()
# 





# saveRDS(sce, file="~/postdoc/stanford/clinical_data/BC1/tfh_data/raw_tfh_batch1/all.RDS")
# saveRDS(t_cells_only, file="~/postdoc/stanford/clinical_data/BC1/tfh_data/raw_tfh_batch1/t_cells_only.RDS")
# saveRDS(clean_cd4_t_cells_only, file="~/postdoc/stanford/clinical_data/BC1/tfh_data/raw_tfh_batch1/clean_cd4_t_cells_only.RDS")

#set.seed(1234);t_cells_only <- CATALYST::cluster(t_cells_only, features = tfh_markers, xdim = 10, ydim = 10, maxK = 30)

# diffcyt ####
library(diffcyt)


ei <- metadata(all_cd4_only)$experiment_info
gender_design <- createDesignMatrix(ei, c("gender"))
# 

#design <- model.matrix(~ei$time+ei$volunteer:ei$alt)
# batch_design <- createDesignMatrix(ei, c("timepoint", "t"))

FDR_cutoff <- 0.1

binary_contrast <- createContrast(c(0, 1))

gender_diff <- diffcyt(all_cd4_only,
                  design = gender_design,
                  contrast = binary_contrast,
                  analysis_type = "DA",
                  method_DA = "diffcyt-DA-edgeR",
                  clustering_to_use = "meta12",
                  verbose = T)

table(rowData(gender_diff$res)$p_adj < FDR_cutoff)


# function to calculate percentages of clusters, keep p values
make_ggplotable <- function(da){
  
  res = data.frame(diffcyt::topTable(da, all=T, show_counts = T))
  
  res %>%
    pivot_longer(cols=starts_with("counts_"), names_to = "sample_id", values_to = "counts")%>%
    mutate(sample_id = gsub("counts_", "", sample_id))%>%
    group_by(sample_id)%>%
    mutate(perc=counts/sum(counts)) -> df
  
  return(df)
  
}


cluster_by_gender <- make_ggplotable(gender_diff)

cluster_by_gender <- left_join(cluster_by_gender, ei, by="sample_id")
cluster_by_gender$gender <- ifelse(cluster_by_gender$gender==1, "male", ifelse(cluster_by_gender$gender==2, "female", NA))

ggplot(cluster_by_gender, aes(x=gender, y=perc, fill=factor(gender)))+
  geom_point()+
  geom_boxplot()+
  facet_wrap(~cluster_id, scales = "free")+
  scale_y_continuous(labels=scales::label_percent())+
  theme_minimal()+
  theme(legend.position = "none")



hp_design <- createDesignMatrix(ei, c("anyHPfinalx"))

hp_diff <- diffcyt(all_cd4_only,
                       design = hp_design,
                       contrast = binary_contrast,
                       analysis_type = "DA",
                       method_DA = "diffcyt-DA-edgeR",
                       clustering_to_use = "meta12",
                       verbose = T)

table(rowData(hp_diff$res)$p_adj < FDR_cutoff)


cluster_by_hp <- make_ggplotable(hp_diff)

cluster_by_hp <- left_join(cluster_by_hp, ei, by="sample_id")
cluster_by_hp$hp <- ifelse(cluster_by_hp$anyHPfinal==1, "placental malaria", ifelse(cluster_by_hp$anyHPfinal==0, "no pathology", "no data"))

dev.new()
ggplot(cluster_by_hp, aes(x=hp, y=perc, fill=factor(hp)))+
  geom_point(aes(color=p_adj<0.1))+
  geom_boxplot()+
  facet_wrap(~cluster_id, scales = "free")+
  theme_minimal()+
  scale_color_manual(values = c("Black", "Red"))+
  scale_y_continuous(labels=scales::label_percent())+
  theme(legend.position = "none")




inf_0_12_design <- createDesignMatrix(ei, c("anyHPfinalx"))

inf_0_12_diff <- diffcyt(all_cd4_only,
                   design = inf_0_12_design,
                   contrast = binary_contrast,
                   analysis_type = "DA",
                   method_DA = "diffcyt-DA-edgeR",
                   clustering_to_use = "meta12",
                   verbose = T)

table(rowData(inf_0_12_diff$res)$p_adj < FDR_cutoff)


cluster_by_hp <- make_ggplotable(inf_0_12_diff)

cluster_by_hp <- left_join(cluster_by_hp, ei, by="sample_id")
cluster_by_hp$hp <- ifelse(cluster_by_hp$anyHPfinal==1, "placental malaria", ifelse(cluster_by_hp$anyHPfinal==0, "no pathology", "no data"))

dev.new()
ggplot(cluster_by_hp, aes(x=hp, y=perc, fill=factor(hp)))+
  geom_point(aes(color=p_adj<0.1))+
  geom_boxplot()+
  facet_wrap(~cluster_id, scales = "free")+
  theme_minimal()+
  scale_color_manual(values = c("Black", "Red"))+
  scale_y_continuous(labels=scales::label_percent())+
  theme(legend.position = "none")


# ggcyto ####

colData(all_cd4_only)$meta10 <- cluster_ids(all_cd4_only, "meta10")

plotScatter(
  all_cd4_only,
  chs=c("CD45RA", "CXCR3"),
  #facet_by = "meta10"
)


plotScatter(
  all_cd4_only,
  chs=c("CCR4", "CCR6"),
  #facet_by = "meta10"
)

library(ggcyto)

ggcyto(processed_tfhs, aes(x = `CD45RA`, y = `CXCR3`))+
  scale_x_logicle()+
  scale_y_logicle()+
  geom_hex(bins=256)
