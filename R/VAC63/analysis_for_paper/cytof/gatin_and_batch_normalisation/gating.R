library(CATALYST)
library(diffcyt)
library(openCyto)
library(ggcyto)
# library(dplyr)
# library(tidyr)
#library(ggplot2)

`%notin%` <- Negate(`%in%`)


.scatter <- function(gs, chs, gate_id = NULL,
                     subset = ifelse(is.null(gate_id), "root", "_parent_")) {
  p <- ggcyto(gs, max_nrow_to_plot = 1e5,
              aes_string(chs[1], chs[2]), subset) +
    geom_hex(bins = 100) + facet_wrap(~ name, ncol = 5) +
    (if (is.null(gate_id)) list() else geom_gate(gate_id)) +
    ggtitle(NULL) + theme_bw(base_size = 8) + theme(
      aspect.ratio = 1,
      legend.position = "none",
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = NA),
      axis.text = element_text(color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1))
  suppressMessages(p + coord_equal(expand = FALSE,
                                   xlim = c(-1, 11), ylim = c(-1, 11)))
}



setwd("~/PhD/cytof/vac63c/normalised_renamed_comped/debarcoded/")

fcs <- list.files(pattern = "fcs")
#fcs <- subset(fcs, !grepl(pattern = "C45", fcs))
#fcs <- subset(fcs, grepl(pattern = "Baseline", fcs))
#fcs <- subset(fcs, !grepl(pattern = "DoD", fcs))
fcs <- subset(fcs, grepl(pattern = "ctrl", fcs))

vac63_flowset <- flowCore::read.flowSet(fcs)

vac63_gs <- GatingSet(vac63_flowset)


#define DNA channels
panel <- read.csv("~/PhD/cytof/vac63c/normalised_renamed_comped/T_cells_only/vac63c_panel.csv", header=T)
dna <- grep("^Ir", panel$fcs_colname, value = TRUE)



gs_add_gating_method(vac63_gs,
                     alias = "cells",
                     pop = "+", parent = "root",
                     dims = paste(dna[1],"Ce140Di", sep = ", "),
                     gating_method = "flowClust.2d",
                     gating_args = "K=1,quantile=0.97,target=c(5,5)")



gs_add_gating_method(vac63_gs,
                     alias = "singlets",
                     pop = "+", parent = "cells",
                     dims = paste(dna, collapse = ","),
                              gating_method = "flowClust.2d",
                              gating_args = "K=1,quantile=0.97,target=c(5,5)")
                     
df <- gs_pop_get_stats(vac63_gs,
                       type = "percent",
                       nodes = c("cells", "singlets"))

df

ggcyto(vac63_gs, aes(x="Ir191Di", y="Ce140Di"), subset = "root")+
  geom_hex(bins=100)+
  geom_gate("cells")+
  scale_x_flowCore_fasinh()+
  scale_y_flowCore_fasinh()



ggcyto(vac63_gs, aes(x="Ir191Di", y="Ir193Di"), subset = "root")+
  geom_hex(bins=100)+
  geom_gate("singlets")+
  scale_x_flowCore_fasinh()+
  scale_y_flowCore_fasinh()



fs <- gs_pop_get_data(vac63_gs, "/cells/singlets") # get data from ’GatingSet’
write.flowSet(x = fs, outdir = "./whole_blood_single_cells/")



