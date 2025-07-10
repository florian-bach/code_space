library(CytoNorm)
library(FlowSOM)
library(dplyr)
library(ggcyto)
library(flowCore)
library(ggplot2)


# cytonorm ####

setwd("~/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only")

fcs <- list.files(pattern = "fcs")
ctrl_fcs <- subset(fcs, grepl(pattern = "ctrl", x = fcs))
base_fcs <- subset(fcs, grepl(pattern = "Baseline", x = fcs))

cytonorm_metadata <- data.frame(File=c(ctrl_fcs, base_fcs),
                                Path=paste("~/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/", c(ctrl_fcs, base_fcs)),
                                Type=c(rep("ctrl", length(ctrl_fcs)), rep("baseline", length(base_fcs))),
                                Batch=paste("Batch_", c(1,2,3,1,1,1,2,2,2,2,3,3,3,3), sep='')
)

train_data <- dplyr::filter(cytonorm_metadata, Type == "ctrl")
validation_data <- dplyr::filter(cytonorm_metadata, Type == "baseline")



ff <- flowCore::read.flowSet(ctrl_fcs)
channels <- flowCore::colnames(ff)[c(14, 22:56, 62)]
transformList <- flowCore::transformList(channels,
                                         cytofTransform)

transformList.reverse <- flowCore::transformList(channels,
                                                 cytofTransform.reverse)
# 
# 
# 
# 
# fsom <- prepareFlowSOM(ctrl_fcs,
#                        channels,
#                        nCells = 25000,
#                        FlowSOM.params = list(xdim = 10,
#                                              ydim = 10,
#                                              nClus = 50,
#                                              scale = FALSE),
#                        transformList = transformList,
#                        seed = 1)
# 
# cvs <- testCV(fsom,
#               cluster_values = c(30, 35, 40)) 

# tis appropriate yo.. let's train the model!

model <- CytoNorm.train(files = train_data$File,
                        labels = train_data$Batch,
                        channels = channels,
                        transformList = transformList,
                        FlowSOM.params = list(nCells = 25000, 
                                              xdim = 10,
                                              ydim = 10,
                                              nClus = 45,
                                              scale = FALSE),
                        normMethod.train = QuantileNorm.train,
                        normParams = list(nQ = 101,
                                          goal = "mean"),
                        seed = 1,
                        verbose = TRUE)



# batch normalise data


CytoNorm.normalize(model = model,
                   files = validation_data$File,
                   labels = validation_data$Batch,
                   transformList = transformList,
                   transformList.reverse = transformList.reverse,
                   normMethod.normalize = QuantileNorm.normalize,
                   outputDir = "Normalized",
                   prefix = "Norm_",
                   clean = TRUE,
                   verbose = TRUE)



#inspect results by ggcyto ####


setwd("~/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only")

fcs <- list.files(pattern = "fcs")
ctrl_fcs <- subset(fcs, grepl(pattern = "ctrl", x = fcs))
  base_fcs <- subset(fcs, grepl(pattern = "T6", x = fcs))
  
  pre_norm <- read.flowSet(base_fcs[c(1:2, 4:5, 10:11)])
  post_norm <- read.flowSet(paste("./Normalized/Norm_", base_fcs[c(1:2, 4:5, 10:11)], sep=''))
  
  panel <- read.csv("vac63c_panel.csv", header=T)
  
  get_channel <- function(antigen){
    panel$fcs_colname[match(antigen, panel$antigen)]
  }
  
  norm_aes <- aes(x=!!get_channel("GATA3"), 
                  y=!!get_channel("ICOS"))
  
  prenorm_plots <- ggcyto(pre_norm, norm_aes)+
    #geom_hex(bins=45)+
    scale_x_flowjo_fasinh()+
    scale_y_flowjo_fasinh()+
    facet_wrap(~name, labeller = labeller(name = label_wrap_gen(10, multi_line = TRUE)))
  
  post_norm_plots <- ggcyto(post_norm, norm_aes)+
    #geom_hex(bins=45)+
    scale_x_flowjo_fasinh()+
    scale_y_flowjo_fasinh()+
    facet_wrap(~name, labeller = labeller(name = label_wrap_gen(10, multi_line = TRUE)))
  
  
  
  gg_prenorm_plots <- as.ggplot(prenorm_plots)
  
  gg_prenorm_plots <- gg_prenorm_plots+
    geom_point(shape=".")+
    stat_density_2d(contour = TRUE, bins=64, color="red")+
    coord_cartesian(xlim = c(0, 550),
                    ylim = c(0, 550))
  
  ggsave("gg_prenorm_plots.png", gg_prenorm_plots, height=4, width=8)
  
  gg_post_norm_plots <- as.ggplot(post_norm_plots)
  
  gg_post_norm_plots <- gg_post_norm_plots+
    geom_point(shape=".")+
    stat_density_2d(contour = TRUE, bins=64, color="red")+
    coord_cartesian(xlim = c(0, 550),
                    ylim = c(0, 550))
  
  ggsave("gg_post_norm_plots.png", gg_post_norm_plots, height=4, width = 8)



#inspect results by comparing flowsom results ####

library(CATALYST)
library(flowCore)

setwd("~/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only")

fcs <- list.files(path = "~/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only", pattern = "fcs", full.names = TRUE)

prenorm_baseline_fcs <- subset(fcs, grepl(pattern = "Baseline", x = fcs))
prenorm_t6_fcs <- subset(fcs, grepl(pattern = "T6", x = fcs))

pre_norm_sample <- c(prenorm_baseline_fcs[c(1:2, 4:5, 10:11)],prenorm_t6_fcs[c(1:2, 4:5, 10:11)]) 

post_fcs <- list.files(path = "~/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/Normalized", pattern = "fcs", full.names = TRUE)
post_norm_baseline_fcs <- subset(post_fcs, grepl(pattern = "Baseline", x = post_fcs))
post_norm_t6_fcs <- subset(post_fcs, grepl(pattern = "T6", x = post_fcs))

post_norm_sample <- c(post_norm_baseline_fcs[c(1:2, 4:5, 10:11)],post_norm_t6_fcs[c(1:2, 4:5, 10:11)]) 

norm_check <- c(pre_norm_sample, post_norm_sample)

post_norm <- read.flowSet(norm_check)

panel <- read.csv("vac63c_panel.csv", header=T)

new_file_names <- substr(norm_check[1:12], start = 72, stop = nchar(norm_check[1:12]))
new_file_names <- c(new_file_names, substr(norm_check[13:24], start = 83, stop = nchar(norm_check[13:24])))
md <- data.frame("file_name"=new_file_names,
                 "batch"=rep(c("batch_1", "batch_1", "batch_2", "batch_2", "batch_3", "batch_3"), 4),
                 "timepoint"=rep(c("Baseline", "T6", "Baseline", "T6"), each=6),
                 "volunteer" = rep(c("301", "302", "305", "306", "315", "320"), 4),
                 "normalised"=c(rep("raw", 12), rep("normalised", 12)))

md$sample_id <- paste(md$volunteer, md$timepoint, md$normalised, sep="_")           

sce <- prepData(post_norm, panel, md, md_cols =
                  list(file = "file_name", id = "sample_id", factors = c("batch", "timepoint", "volunteer", "normalised")))



refined_markers <- c("CD4",
                     "CD8",
                     "Vd2",
                     "Va72",
                     "CD38",
                     "CD69",
                     "HLADR",
                     "ICOS",
                     "CD28",
                     "PD1",
                     #"TIM3",
                     "CD95",
                     "BCL2",
                     "CD27",
                     "Perforin",
                     "GZB",
                     "Tbet",
                     "Eomes",
                     #"RORgt",
                     #"GATA3",
                     "CTLA4",
                     "Ki67",
                     "CD127",
                     "CD56",
                     #"CD16",
                     "CD161",
                     "CD49d",
                     "CD25",
                     "FoxP3",
                     "CD39",
                     "CXCR5",
                     "CX3CR1",
                     "CD57",
                     "CD45RA",
                     "CD45RO",
                     "CCR7")

set.seed(123);sce <- CATALYST::cluster(sce, features = "type", xdim = 10, ydim = 10, maxK = 50)

#plotExprHeatmap(x = sce, by = "cluster", row_clust = FALSE, k = "meta45", bars = TRUE)
plotAbundances(x = sce, by = "cluster_id", k = "meta45", group_by = "normalised", shape_by = "volunteer")
ggplot2::ggsave("norm_vs_raw_flowsom.png", height=11, width=8)





plotExprHeatmap(x = sce, by = "cluster", row_clust = FALSE, k = "meta45", bars = FALSE, features = refined_markers)







# library(CATALYST)
# library(dplyr)
# library(flowCore)
# library(flowWorkspace)
# library(ggcyto)
# library(ggplot2)
# library(mvtnorm)
# library(openCyto)
# 
# qc_theme <- list(
#   theme_bw(base_size = 8), theme(
#     panel.grid.minor = element_blank(),
#     panel.grid.major.x = element_blank(),
#     plot.title = element_text(face = "bold"),
#     axis.text = element_text(color = "black"),
#     axis.text.x = element_text(angle = 45, hjust = 1)))
# 
# 
# fcs <- list.files("~/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/", "fcs", full.names = TRUE)
# base_fcs <- subset(fcs, grepl(pattern = "Baseline", x = fcs))
# 
# (sce <- prepData(base_fcs, transform = TRUE, cofactor = 5))
# 
# ctrl_fcs <- subset(fcs, grepl(pattern = "ctrl", x = fcs))
# 
# 
# 
# # compute 98th count quantiles via back-transformation
# # (using same cofactor as always) & average across replicates
# cf <- int_metadata(sce)$cofactor
# ref <- data.frame(
#   target = colnames(ref),
#   mean_count = colMeans(sinh(ref)*cf))
# 
# # initialize correction factor of 1 for all channels
# cfs <- setNames(rep(1, nrow(sce)), rownames(sce))
# 
# # compute batch correction factors for relevant channels
# cs <- assay(sce, "compcounts")
# csR <- cs[ref$target, sce$group == "R"]
# run <- rowQuantiles(csR, probs = 0.98)
# cfs[ref$target] <- run / ref$mean_count
# 
# # apply marker-specific batch correction (bc)
# cs <- sweep(cs, 1, cfs, "/")
# assay(sce, "bccounts") <- cs
# 
# # apply arcsinh-transformation
# assay(sce, "bcexprs") <- asinh(cs/cf)
# 
# 
# 
