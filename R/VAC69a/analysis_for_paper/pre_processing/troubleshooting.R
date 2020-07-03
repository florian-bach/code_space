library(readxl)    
library(HDCytoData)                                                               

RNGversion("3.5.3")

url <- "http://imlspenticton.uzh.ch/robinson_lab/cytofWorkflow"
md <- "PBMC8_metadata.xlsx"                                    
download.file(file.path(url, md), destfile = md, mode = "wb")  
md <- read_excel(md)                                           
head(data.frame(md))    

fs <- Bodenmiller_BCR_XL_flowSet() 





panel <- "PBMC8_panel_v3.xlsx"                                     
download.file(file.path(url, panel), destfile = panel, mode = "wb")
panel <- read_excel(panel) 


md$condition <- factor(md$condition, levels = c("Ref", "BCRXL"))        
md$sample_id <- factor(md$sample_id,                                    
                       levels = md$sample_id[order(md$condition)])                         

# construct SingleCellExperiment                                        
sce <- prepData(fs, panel, md, features = panel$fcs_colname)

set.seed(1234)                                                                      
sce <- runDR(sce, dr = "UMAP", cells = 1e3, features = "type") 

plotDR(sce, "UMAP", color_by = "CD4")





scaled <- scater::runUMAP(sce,
                subset_row=type_markers(sce),
                exprs_values = "exprs",
                scale=T)


non_scaled <- scater::runUMAP(sce,
                          subset_row=type_markers(sce),
                          exprs_values = "exprs",
                          scale=F)

scaled_plot <- plotDR(scaled, color_by="CD4")

non_scaled_plot <- plotDR(non_scaled, color_by="CD4")


cowplot::plot_grid(scaled_plot, non_scaled_plot, ncol=1)

