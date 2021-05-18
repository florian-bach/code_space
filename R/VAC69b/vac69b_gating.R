library(openCyto)
library(ggcyto)
# library(dplyr)
# library(tidyr)
#library(ggplot2)
# 
`%notin%` <- Negate(`%in%`)
# 
# 
# .scatter <- function(gs, chs, gate_id = NULL,
#                      subset = ifelse(is.null(gate_id), "root", "_parent_")) {
#   p <- ggcyto(gs, max_nrow_to_plot = 1e5,
#               aes_string(chs[1], chs[2]), subset) +
#     geom_hex(bins = 100) + facet_wrap(~ name, ncol = 5) +
#     (if (is.null(gate_id)) list() else geom_gate(gate_id)) +
#     ggtitle(NULL) + theme_bw(base_size = 8) + theme(
#       aspect.ratio = 1,
#       legend.position = "none",
#       panel.grid.minor = element_blank(),
#       strip.background = element_rect(fill = NA),
#       axis.text = element_text(color = "black"),
#       axis.text.x = element_text(angle = 45, hjust = 1))
#   suppressMessages(p + coord_equal(expand = FALSE,
#                                    xlim = c(-1, 11), ylim = c(-1, 11)))
# }



setwd("~/PhD/cytof/vac69b/")

md <- read.csv("./T_cells_only/metadata.csv")
md <- subset(md, md$volunteer %in% c("v07"))


vac69b_flowset <- flowCore::read.flowSet(md$file_name)

vac69b_gs <- GatingSet(vac69b_flowset)


#define DNA channels
panel <- read.csv("~/PhD/cytof/vac69b/T_cells_only/vac69b_panel.csv", header=T)
dna <- grep("^Ir", panel$fcs_colname, value = TRUE)



gs_add_gating_method(vac69b_gs,
                     alias = "cells",
                     pop = "+", parent = "root",
                     dims = "Ir191Di,Ce140Di",
                     gating_method = "flowClust.2d",
                     gating_args = "K=1,quantile=0.97,target=c(5,5)")



gs_add_gating_method(vac69b_gs,
                     alias = "singlets",
                     pop = "+", parent = "cells",
                     dims = paste(dna, collapse = ","),
                     gating_method = "flowClust.2d",
                     gating_args = "K=1,quantile=0.97,target=c(5,5)")

df <- gs_pop_get_stats(vac69b_gs,
                       type = "percent",
                       nodes = c("cells", "singlets"))

df


fs <- gs_pop_get_data(vac69b_gs, "/cells/singlets") # get data from ’GatingSet’
write.flowSet(x = fs, outdir = "./whole_blood_single_cells/")
