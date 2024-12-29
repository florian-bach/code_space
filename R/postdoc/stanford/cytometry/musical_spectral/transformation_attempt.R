library(flowGate)
library(flowCore)
library(flowWorkspace)
library(openCyto)
library(ggcyto)


# fcs_files <- list.files("~/Downloads/FCSFilesLiveSinglets-selected", pattern = ".fcs")
# full_paths <-  list.files("~/Downloads/FCSFilesLiveSinglets-selected", pattern = ".fcs", full.names = TRUE)

fcs_files <- list.files("~/Downloads/fcs-selected", pattern = ".fcs")
full_paths <-  list.files("~/Downloads/fcs-selected", pattern = ".fcs", full.names = TRUE)

csv_to_fcs <- openCyto:::.read.FCS.csv(file = "~/Downloads/csv scale values-selected/livecells_mus4_D2 Well_038_Live Cells.csv")

musical_panel <- readxl::read_excel("~/Downloads/MUSICALSpectralFlowPanel.xlsx")
musical_panel$Fluorophore <- paste(musical_panel$Fluorophore, "A", sep="-")

musical_flowset <- read.flowSet(files = full_paths)

# fs_non_neg = fs[rowMeans(exprs(fs[,1:2])>0) == 1,]

flowjo_biex <- flowjo_biexp(
  # channelRange = 5*10^5,
  # maxValue = 2^22,
  pos = 4.5,
  neg = 0
  # widthBasis = -20,
  # inverse = FALSE
)


translist <- transformList(musical_panel$Fluorophore, flowjo_biex)

ff_t <- flowCore::transform(musical_flowset, translist)

autoplot(try, x="CD19", y="CD20", bins=256)+
  scale_x_flowjo_biexp()+
  scale_y_flowjo_biexp()

ggplot(try, aes(x="CD3", y="CD8"))+
  geom_hex(bins = 128)
  axis_x_inverse_trans()+
  axis_y_inverse_trans()

  
autoplot(fs, "CD3")+scale_x_flowJo_biexp()

# cytotree####

fcs_files <- list.files("~/Downloads/FCSFilesLiveSinglets-selected", pattern = ".fcs")
full_paths <-  list.files("~/Downloads/FCSFilesLiveSinglets-selected", pattern = ".fcs", full.names = TRUE)

musical_panel <- readxl::read_excel("~/Downloads/MUSICALSpectralFlowPanel.xlsx")
musical_panel$Fluorophore <- paste(musical_panel$Fluorophore, "A", sep="-")

fs <- ncdfFlow::read.ncdfFlowSet(files = full_paths)

trans <- estimateLogicle(fs, colnames(fcsfile_comp_clean[,7:38]))
fcsfile_comp_clean_trans <- transform(fcsfile_comp_clean, trans)
# See information
cyt



#cytoexplorer ####

library(CytoExploreR)


full_paths <-  list.files("~/Downloads/FCSFilesLiveSinglets-selected", pattern = ".fcs", full.names = TRUE)

gs <- cyto_setup("~/Downloads/FCSFilesLiveSinglets-selected", clean=T, )

# trans_biex <- cyto_transformer_biex(gs,
#                                     channelRange = 5*10^5,
#                                     maxValue = 2^22,
#                                     pos = 4.5,
#                                     neg = 0,
#                                     widthBasis = -200,
#                                     inverse = FALSE)

trans_logicle <- cyto_transformer_logicle(x = gs,
                                          t = 2^22,
                                          w = 3.3,
                                          m = 6.62,
                                          a = 0)



trans_gs <- cyto_transform(gs, trans=trans_logicle)

cyto_plot(gs, channels = c("CD3", "CD8"))


