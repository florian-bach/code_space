war

#inspect results

pre_norm <- read.flowSet(base_fcs[c(1:2, 4:5, 10:11)])
post_norm <- read.flowSet(paste("./Normalized/Norm_", base_fcs[c(1:2, 4:5, 10:11)], sep=''))

norm_aes<-aes(x=Sm152Di, y=Nd146Di)

prenorm_plots <- ggcyto(pre_norm, norm_aes)+
  #geom_hex(bins=45)+
  scale_x_flowjo_fasinh()+
  scale_y_flowjo_fasinh()

post_norm_plots <- ggcyto(post_norm, norm_aes)+
  #geom_hex(bins=45)+
  scale_x_flowjo_fasinh()+
  scale_y_flowjo_fasinh()


gg_prenorm_plots <- as.ggplot(prenorm_plots)

gg_prenorm_plots+
  geom_point(shape=".")+
  stat_density_2d(contour = TRUE, bins=64, color="red")+
  coord_cartesian(xlim = c(0, 550),
                  ylim = c(0, 550))

ggsave("gg_prenorm_plots.png", height=4, width=8)

gg_post_norm_plots <- as.ggplot(post_norm_plots)

gg_post_norm_plots+
  geom_point(shape=".")+
  stat_density_2d(contour = TRUE, bins=64, color="red")+
  coord_cartesian(xlim = c(0, 550),
                  ylim = c(0, 550))

ggsave("gg_post_norm_plots.png", height=4, width = 8)


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
