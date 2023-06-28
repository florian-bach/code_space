library(CATALYST)

#make metadata table

# file_paths<- list.files("~/postdoc/stanford/clinical_data/BC1/tfh_data/SingleCell_TFH/", pattern="*.fcs", full.names = TRUE)
# file_names <- list.files("~/postdoc/stanford/clinical_data/BC1/tfh_data/SingleCell_TFH/", pattern="*.fcs")
# 
# md <- data.frame("file_paths"=file_paths,
#                  "file_names"=file_names,
#                  "sample_id"=substr(file_paths, 95, 98),
#                  "id"=paste(1, substr(file_paths, 95, 98), sep=""))
# 




# md$batch <- ""
# md$sex <- ""



# write.csv(md, "~/postdoc/stanford/clinical_data/BC1/tfh_data/metadata.csv", quote = FALSE, row.names = FALSE)
md <- read.csv("~/postdoc/stanford/clinical_data/BC1/tfh_data/metadata.csv", header=T, stringsAsFactors = TRUE)

cleaned_tfh <- read.csv("~/postdoc/stanford/clinical_data/BC1/tfh_data/cTfh_all.csv")
colnames(cleaned_tfh) <- c("idt", "tfh_q4", "tfh_q3","tfh_q2","tfh_q1","tfh_of_cd4", "parent2t", "tfh_drop")
cleaned_tfh$idt <- paste("1", substr(cleaned_tfh$idt, 15, 18), sep='')
drop_removed_tfh <- subset(cleaned_tfh, tfh_drop=="No")

md <- subset(md, id %in% drop_removed_tfh$idt)


#make panel
# premessa_table <- premessa::read_parameters(list.files("~/postdoc/stanford/clinical_data/BC1/tfh_data/SingleCell_TFH/", pattern="*.fcs", full.names = TRUE))
# panel <- data.frame("fcs_colname"=rownames(premessa_table[1]),
#                     "marker_name"=premessa_table[1]$export_T_FH_PANEL_B1_1001_Single_Cells.fcs,
#                     "marker_class"="type")
# panel$antigen <- ifelse(panel$marker_name=="", panel$fcs_colname, panel$marker_name)
# write.csv(panel, "~/postdoc/stanford/clinical_data/BC1/tfh_data/tfh_panel", quote = FALSE, row.names = FALSE)
panel <- read.csv("~/postdoc/stanford/clinical_data/BC1/tfh_data/tfh_panel", header = T, stringsAsFactors = F)

### read in flowfiles using flowCore
tfhs <- flowCore::read.flowSet(as.character(md$file_paths), truncate_max_range = FALSE)
channel_names <- markernames(tfhs[[1]])
markernames(tfhs) <- channel_names


### CATALYST ####
#construct daFrame #
## md has to have particular properties: file_name=NAMED LIST (chr), ID and everything else in factors
tfh_markers <- panel$marker_name[1:8]

sce <- CATALYST::prepData(tfhs, panel, md, FACS = TRUE, fix_chs="common",
                          transform = FALSE,
                          features=panel$fcs_colname[1:8],
                          md_cols = list(file = "file_names",
                                         id = "sample_id",
                                         factors = c("id", "file_paths")
                          ),
                          panel_cols = list(channel = "fcs_colname",
                                            antigen = "marker_name",
                                            class = "marker_class")
)

# clustering ####

set.seed(123);sce <- CATALYST::cluster(sce, features = tfh_markers, xdim = 10, ydim = 10, maxK = 50)

big_heat <- plotMultiHeatmap(sce,
                 k = "meta20",
                 m="meta10",
                 hm1="type",
                 hm2 = "abundances",
                 row_clust = TRUE,
                 col_clust = TRUE,
                 col_dend = TRUE,
                 scale = "first")


png("~/postdoc/stanford/clinical_data/BC1/tfh_data/figures/big_heat_dropped.png", width=30, height=4, units = "in", res=400)
ComplexHeatmap::draw(big_heat)
dev.off()


system.time(sce <- runDR(sce,
                         dr="UMAP",
                         cells=10000,
                         features=tfh_markers,
                         scale=T)
)

sce <- HDF5Array::loadHDF5SummarizedExperiment(dir ="~/postdoc/stanford/clinical_data/BC1/tfhsummarizedexperiment/")

marker_umap <- plotDR(sce,  color_by = tfh_markers, scale = TRUE)
ggplot2::ggsave("~/postdoc/stanford/clinical_data/BC1/tfh_data/figures/marker_umap.png", marker_umap, height = 10, width=10, bg="white")

cluster_umap <- plotDR(sce,  color_by = "meta16")
ggplot2::ggsave("~/postdoc/stanford/clinical_data/BC1/tfh_data/figures/cluster_umap.png", cluster_umap, height = 4, width=4, bg="white")

HDF5Array::saveHDF5SummarizedExperiment(sce, dir="~/postdoc/stanford/clinical_data/BC1/tfhsummarizedexperiment/", prefix="", replace=FALSE,
                             chunkdim=NULL, level=NULL, as.sparse=NA,
                             verbose=NA)
#ggcyto ####
