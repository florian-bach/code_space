library(ggplot2)
library(tidyr)
library(RColorBrewer)
library(gridExtra)
library(cowplot)

setwd("/Users/s1249052//PhD/cytof/better_gating/double_flowsoms/FlowSOM_big_timecourse_06b_results/results/cluster_medians/")
data <- read.csv("aggregate_cluster_medians.csv")

data_clean <- select(data, colnames(data)[1:36])

colnames(data_clean) <- c('ClusterID', 'MetaclusterID', '115In_CD57', '141Pr_HLA-DR', '142Nd_BCL-2', '143Nd_CD45RA', '144Nd_GZB', '145Nd_CD4', '146Nd_Vd2',
                          '148Nd_ICOS', '149Sm_CXCR5', '150Nd_CD95', '151Eu_CD103', '153Eu_Va7.2', '154Sm_TIM-3', '155Gd_PD1',
                          '156Gd_CD161', '158Gd_CD27', '159Tb_FoxP3', '160Gd_CTLA4', '161Dy_Tbet', '162Dy_IntegrinB7', '163Dy_CD28', '164Dy_Ki-67',
                          '165Ho_CD45RO', '166Er_CD56', '167Er_CCR7', '168Er_CD127', '169Tm_CD38', '171Yb_CD49d', '172Yb_CD25', '173Yb_CD39',
                          '174Yb_CLA', '175Lu_Perforin', '198Pt_CD8', '209Bi_CD16')

colnames(data_clean)[3:36] <- substr(colnames(data_clean)[3:36], 7, 30)
head(data_clean)

vol_03_t6 <- gather(data_clean, Channel, Median, c(colnames(data_clean)[3:36]))
vol_03_t6$MetaclusterID <- NULL

vol_03_t6 <- vol_03_t6[1:200,]
cluster_list <- split(vol_03_t6, vol_03_t6$ClusterID)



facetted <- ggplot(vol_03_t6, aes(x = Channel, y = Median, group = ClusterID))+
  geom_line(color="black")+
  coord_polar(theta = "x", direction = -1)+
  facet_wrap(~ ClusterID, ncol=8)+
  theme(axis.text.x = element_text(vjust = 3, size=12, color="black"),
        axis.text.y = element_blank(),
        axis.ticks.y  = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        aspect.ratio = 1.5
  )

ggsave("facetted.pdf", plot=facetted, width=15, height=60, device = "pdf",  dpi=350, limitsize = FALSE)  




list_of_graphics <- vector("list", length=36)

for (i in cluster_list){
  
  graphic <- ggplot(i, aes(x = Channel, y = Median, group = ClusterID)) +
    geom_line(color = "black") +
    coord_polar(theta = "x", direction = -1)+
    ggtitle(paste("Cluster_", as.character(unique(i$ClusterID)), sep=''))+
    theme(axis.text.x = element_text(size=16),
          axis.text.y = element_blank(),
          axis.ticks.y  = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank()
          )
   
    
    #print(graphic)
  
  list_of_graphics[[as.numeric(unique(i$ClusterID))]] <- graphic
  print(as.numeric(unique(i$ClusterID)))
  }


grid.arrange(grobs=grob, widths=c(100,100,100,100,100,100), ncol=6)

grob <- marrangeGrob(list_of_graphics, ncol=6, nrow=6)

multiplot <- align_plots(list_of_graphics[1:36], align = "hv")

######## sandbox


ggplot(cluster_list[[1]], aes(x = Channel, y = Median, group = ClusterID)) +
  geom_line(color = "black") +
  coord_polar(theta = "x", direction = -1)+
  ggtitle(paste("Cluster_", as.character(unique(cluster_list[[1]][[1]])), sep=''))+
  theme(axis.text.x = element_text(size=20),
        axis.text.y = element_blank(),
        axis.ticks.y  = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())


ggtitle(paste("Cluster_", as.character(unique(i$ClusterID)), sep=''))


mtcars$model<-rownames(mtcars)


mtcarsscaled <- as.data.frame(lapply(mtcars, ggplot2:::rescale01))
mtcarsscaled$model <- rownames(mtcars)
mtcarsmelted <- reshape2::melt(mtcarsscaled)


ggplot(mtcarsmelted, aes(x=variable, y=value))+
  geom_polygon(aes(group = model, color = model), fill=NA, show.legend = F)+
  geom_point()+
  coord_radar()
  

coord_radar <- function (theta = "x", start = 0, direction = 1) 
{
  theta <- match.arg(theta, c("x", "y"))
  r <- if (theta == "x") 
    "y"
  else "x"
  ggproto("CordRadar", CoordPolar, theta = theta, r = r, start = start, 
          direction = sign(direction),
          is_linear = function(coord) TRUE)
}




