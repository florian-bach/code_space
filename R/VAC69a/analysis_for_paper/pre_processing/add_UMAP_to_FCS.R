library(ggcyto)
library(flowCore)
library(scater)
library(SingleCellExperiment)
library(CATALYST)

ggcyto(smaller_vac69a, aes(x=`BCL-2`, y=CD38))+
  geom_hex(bins=64)+
  scale_x_flowjo_fasinh()+
  scale_y_flowjo_fasinh()


fcs_files <- grep("fcs", list.files(), value = T)

primaries <- c()

timepoints <- c("C-1", "C_10", "DoD", "T_6")

for(i in timepoints){
  one_moment <- grep(i, fcs_files, value=T)
  primaries <- c(timepoints, one_moment)
}

### read in flowfiles using flowCore
vac69a <- read.flowSet(md$file_name)




daf <- prepData(vac69a, panel, md, md_cols =
                  list(file = "file_name", id = "sample_id", factors = c("timepoint", "batch", "volunteer")),
                panel_cols = list(channel = "fcs_colname", antigen = "marker_name", class =
                                    "marker_class"))

daf <- scater::runUMAP(daf,
                            subset_row=refined_markers,
                            exprs_values = "exprs",
                            scale=T)

#this is potentially inefficent... maybe investigate!
big_table <- data.frame(t(data.frame(assays(daf)$exprs)))
big_table <- data.frame(cbind(big_table, colData(daf)))

slim_umap <- data.frame(reducedDim(daf, "UMAP"))
colnames(slim_umap) <- c("UMAP1", "UMAP2")

big_table <- data.frame(cbind(big_table, slim_umap), stringsAsFactors = F)

data.table::fwrite(big_table, "big_table.csv")


for (i in md$file_name){
  
  post_dim_red <- subset(big_table, sample_id==as.character(subset(md, md$file_name==i)$sample_id))
  dim_red <- cbind(post_dim_red$UMAP1, post_dim_red$UMAP2)
  colnames(dim_red) <- c("UMAP1", "UMAP2")
  
  old_frame <- read.FCS(i)
  new_frame <- fr_append_cols(old_frame, dim_red)
  
  write.FCS(filename=paste("./post_umap/",substr(i, 0, nchar(i)-4),"_post_umap.fcs", sep=''), new_frame)
  print(paste(unique(post_dim_red$sample_id), "done"))
  }
