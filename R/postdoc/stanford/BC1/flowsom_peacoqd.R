
fcs_dir <- "~/postdoc/stanford/clinical_data/BC1/tfh_data/raw_tfh_batch1/Preprocessed"

fcs_file_paths <- list.files(fcs_dir, pattern="*.fcs", full.names = TRUE, recursive = FALSE)[-54]
fcs_file_names <- list.files(fcs_dir, pattern="*.fcs", full.names = FALSE, recursive = FALSE)[-54]

md <- data.frame("file_paths"=fcs_file_paths,
                 "file_names"=fcs_file_names,
                 "sample_id"=substr(fcs_file_names, regexpr("B1 ", fcs_file_names)+3, regexpr("B1 ", fcs_file_names)+6))


cleaned_tfh <- read.csv("~/postdoc/stanford/clinical_data/BC1/tfh_data/cTfh_all.csv")
colnames(cleaned_tfh) <- c("idt", "tfh_q4", "tfh_q3","tfh_q2","tfh_q1","tfh_of_cd4", "parent2t", "tfh_drop")
drop_removed_tfh <- subset(cleaned_tfh, tfh_drop=="No")

md <- subset(md, file_names %in% drop_removed_tfh$idt)

md <- md[1:24,]

premessa_table <- premessa::read_parameters(list.files(fcs_dir, pattern="*.fcs", full.names = TRUE, recursive = FALSE))
panel <- data.frame("fcs_colname"=rownames(premessa_table[1]),
                    "marker_name"=premessa_table[[1]],
                    "marker_class"="type")
panel$antigen <- ifelse(panel$marker_name=="", panel$fcs_colname, panel$marker_name)
panel$marker_name <- ifelse(panel$marker_name=="", panel$fcs_colname, panel$marker_name)


processed_tfhs <- flowCore::read.flowSet(as.character(md$file_paths), truncate_max_range = FALSE, )

sce <- CATALYST::prepData(processed_tfhs, panel, md, FACS = TRUE,
                          fix_chs="common",
                          features=panel$fcs_colname[panel$antigen%in%tfh_markers],
                          md_cols = list(file = "file_names",
                                         id = "sample_id",
                                         factors = c("file_paths")
                          ),
                          panel_cols = list(channel = "fcs_colname",
                                            antigen = "marker_name",
                                            class = "marker_class")
)

# clustering ####

tfh_markers <- c("CD3", "CD4", "CD45RA", "CXCR5", "PD1", "CXCR3", "CCR4", "CCR6")

sce <- CATALYST::cluster(sce, features = tfh_markers, xdim = 10, ydim = 10, maxK = 30, verbose = TRUE, seed=1234)

system.time(sce <-  CATALYST::runDR(sce,
                                    dr="UMAP",
                                    cells=10000,
                                    features=tfh_markers,
                                    scale=T)
)
# HDF5Array::saveHDF5SummarizedExperiment(sce, dir="~/postdoc/stanford/clinical_data/BC1/tfh_data/raw_tfh_batch1/Preprocessed/sce", prefix="", replace=FALSE,
#                                         chunkdim=NULL, level=NULL, as.sparse=NA,
#                                         verbose=NA)


# sce <- HDF5Array::loadHDF5SummarizedExperiment("~/postdoc/stanford/clinical_data/BC1/tfh_data/raw_tfh_batch1/Preprocessed/sce/", prefix = "")

small_heat <-CATALYST::plotExprHeatmap(sce, by="cluster_id",
                                       k="meta24",
                                       features = tfh_markers,
                                       perc = TRUE,
                                       bars = TRUE,
                                       # row_anno = TRUE,
                                       # #col_anno = TRUE,
                                       # row_clust = TRUE,
                                       # #col_clust = TRUE,
                                       # row_dend = TRUE,
                                       # #col_dend = TRUE,
                                       hm_pal = colorspace::sequential_hcl("inferno", n=8))

png("~/postdoc/stanford/clinical_data/BC1/tfh_data/figures/peacoqd_small_heat_dropped.png", width=8, height=8, units = "in", res=400)
ComplexHeatmap::draw(small_heat)
dev.off()

t_cells_only <- CATALYST::filterSCE(sce, k="meta24", cluster_id %in% c(1234))


big_heat <-  CATALYST::plotMultiHeatmap(sce,
                                        k = "meta20",
                                        m="meta10",
                                        hm1="type",
                                        hm2 = NULL,
                                        row_clust = TRUE,
                                        col_clust = TRUE,
                                        col_dend = TRUE,
                                        scale = "first")


png("~/postdoc/stanford/clinical_data/BC1/tfh_data/figures/peacoqd_big_heat_dropped.png", width=30, height=4, units = "in", res=400)
ComplexHeatmap::draw(big_heat)
dev.off()



marker_umap <-  CATALYST::plotDR(sce,  color_by = tfh_markers, scale = TRUE)
ggplot2::ggsave("~/postdoc/stanford/clinical_data/BC1/tfh_data/figures/peacoqd_comp_marker_umap.png", marker_umap, height = 10, width=10, bg="white")

# 
