library(flowCore)
library(vac69a.cytof)

daf <- read_full("~/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/")

refined_markers <- read.csv("~/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/refined_markers.csv")

daf <- scater::runUMAP(daf,
                            subset_row=refined_markers[,1],
                            exprs_values = "exprs",
                            scale=T)

#this is potentially inefficent... maybe investigate!
# big_table <- data.frame(t(data.frame(assays(daf)$exprs)))
# big_table <- data.frame(cbind(big_table, colData(daf)))
# 
# slim_umap <- data.frame(reducedDim(daf, "UMAP"))
# colnames(slim_umap) <- c("UMAP1", "UMAP2")
# 
# big_table <- data.frame(cbind(big_table, slim_umap), stringsAsFactors = F)

big_table <- prep_sce_for_ggplot(daf)

data.table::fwrite(big_table, "~/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/all_cells_with_UMAP.csv")


setwd("~/PhD/cytof/vac69a/reprocessed/reprocessed_relabeled_comped/T_cells_only/")
md <- read.csv("meta_data.csv", stringsAsFactors = F)

for (i in md$file_name){
  
  post_dim_red <- subset(big_table, sample_id==as.character(subset(md, md$file_name==i)$sample_id))
  dim_red <- cbind(post_dim_red$UMAP1, post_dim_red$UMAP2)
  colnames(dim_red) <- c("UMAP1", "UMAP2")
  
  old_frame <- read.FCS(i)
  new_frame <- fr_append_cols(old_frame, dim_red)
  
  write.FCS(filename=paste("./post_umap/",substr(i, 0, nchar(i)-4),"_post_umap.fcs", sep=''), new_frame)
  print(paste(unique(post_dim_red$sample_id), "done"))
  }
