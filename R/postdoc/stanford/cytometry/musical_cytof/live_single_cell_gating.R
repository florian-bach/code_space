library(flowGate)
library(flowCore)
library(flowWorkspace)
library(openCyto)
library(ncdfFlow)

musical_panel <- read.csv("~/postdoc/stanford/cytometry/CyTOF/MUSICAL/pilot75/musical_panel_edit.csv")
mass_channels <- musical_panel$fcs_colname[-c(1,2,44:48)]

path_to_fcs <- "/Users/fbach/postdoc/stanford/cytometry/CyTOF/MUSICAL/redownload_pilot75/redownload_big_fcs/"
files_to_read <- list.files(path_to_fcs, full.names = TRUE, pattern = "*.fcs$")

fs <- read.ncdfFlowSet(files = files_to_read)

asinhTrans <- flowCore::arcsinhTransform(transformationId="defaultArcsinhTransform", a=0, b=1/5, c=0)


translist <- transformList(mass_channels, asinhTrans)

ff_t <- flowCore::transform(fs, translist)

gs <- GatingSet(ff_t)

gs_add_gating_method(gs,
                      alias = "cells",
                      parent = "root",
                      dims = "Ir191Di,Ce140Di",
                      gating_method = "boundary",
                      gating_args = "min = c(5, -1), max=c(8,2)")
 
# gs_remove_gating_method(gs)

cell_gate_plot <- autoplot(gs, "cells", bins=150)
ggsave("/Users/fbach/postdoc/stanford/cytometry/CyTOF/MUSICAL/redownload_pilot75/redownload_big_fcs/figures/ir191_vs_cd140_all.png", height = 20, width=20, bg="white", dpi=444)


gs_add_gating_method(gs,
                     alias = "cells2",
                     parent = "cells",
                     dims = "Ir191Di,Ir193Di",
                     gating_method = "flowClust.2d",
                     gating_args = "K=1, quantile=0.90")

# gs_remove_gating_method(gs)

cell_gate_plot <- autoplot(gs, "cells2", bins=150)
ggsave("/Users/fbach/postdoc/stanford/cytometry/CyTOF/MUSICAL/redownload_pilot75/redownload_big_fcs/figures/ir191_vs_ir193.png", height = 20, width=20, bg="white", dpi=444)



gs_add_gating_method(gs,
                     alias = "singlets",
                     parent = "cells2",
                     dims = "Ir191Di,Event_length",
                     gating_method = "flowClust.2d",
                     gating_args = "K=1, quantile=0.90")



singlet_gate_plot <- autoplot(gs, "singlets", bins=150)
ggsave("/Users/fbach/postdoc/stanford/cytometry/CyTOF/MUSICAL/redownload_pilot75/redownload_big_fcs/figures/singlet_gate_all.png", height = 20, width=20, bg="white", dpi=444)



# gs_add_gating_method(gs,
#                      alias = "singlets2",
#                      parent = "singlets",
#                      dims = "Residual,Event_length",
#                      gating_method = "flowClust.2d",
#                      gating_args = "K=1, quantile=0.90")
# 
# 
# 
# singlet_gate_plot2 <- autoplot(gs, "singlets2", bins=150)
# ggsave("/Users/fbach/postdoc/stanford/cytometry/CyTOF/MUSICAL/redownload_pilot75/redownload_big_fcs/figures/singlet_gate_all2.png", height = 20, width=20, bg="white", dpi=444)




gs_add_gating_method(gs,
                     alias = "live_singlets",
                     pop = "+", parent = "singlets",
                     dims = "Ir193Di,Dead",
                     gating_method = "flowClust.2d",
                     gating_args = "K=1 , quantile=0.90")

live_singlets_gate_plot <- autoplot(gs, "live_singlets", bins=150)
ggsave("/Users/fbach/postdoc/stanford/cytometry/CyTOF/MUSICAL/redownload_pilot75/redownload_big_fcs/figures/live_singlets.png", height = 20, width=20, bg="white", dpi=444)


cleaned_fs <- gs_pop_get_data(gs, y="live_singlets")

write.flowSet(cleaned_fs,
              # filename = "single_cell_",
              outdir = "/Users/fbach/postdoc/stanford/cytometry/CyTOF/MUSICAL/redownload_pilot75/sandbox/")



# check on the gated populations ####
path_to_fcs <- "/Users/fbach/postdoc/stanford/cytometry/CyTOF/MUSICAL/redownload_pilot75/sandbox/"
files_to_read <- list.files(path_to_fcs, full.names = TRUE, pattern = "*.fcs$")

fs <- ncdfFlow::read.ncdfFlowSet(files = files_to_read)
flowCore::fsApply(fs, nrow)
