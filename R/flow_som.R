library(FlowSOM) 
library(readxl) 
library(flowCore) 
library(premessa)
library(matrixStats)
library(ggplot2)
library(reshape2)
library(ConsensusClusterPlus)
library(dplyr)
library(RColorBrewer)
library(pheatmap)
library(limma)
library(ggrepel)
library(tidyr)



setwd("/Users/s1249052/PhD/flow_data/vac69a/t_cells_only/better_gating/")

files_list <- list.files(path=".", pattern="*.fcs")

# concatenate_fcs_files(files_list, "concat02.fcs")
# 
fcs_raw <- list()
for (i in files_list){
  print(i)
  j <- read.FCS(i, transformation = FALSE, truncate_max_range = FALSE)
  assign(paste(substr(i, nchar(i)-10, nchar(i)-4)), j)
  fcs_raw[[paste(substr(i, nchar(i)-10, nchar(i)-4))]] <- j}

head(fcs_raw)
# 
df_lol <- as.matrix(fcs_raw[[1]]@parameters@data[["desc"]])
channel_names <- df_lol[1:69]

# timecourse for only Volunteer 03
flo_set <- read.flowSet(files_list[21:25], transformation = FALSE, truncate_max_range = FALSE)

timepoint <- c(1, 10, 8, 12, 18)
#1 = baseline, 2 = C+10, 3= C+8, 4=DoD, 5=T+6)

new_flo_set <- new('flowSet', colnames=colnames(tmp)[1:68])

for (i in seq(1,5)){
  tmp <- flo_set[[i]]
  col <- rep(timepoint[i], nrow(tmp))
  cols <- matrix(col, dimnames = list(NULL, "Timepoint"))
  tmp <- fr_append_cols(tmp, cols)
  assign(paste("frame_", i, sep=''), tmp)
}

new_flo_set <- new('flowSet', c(frame_1, frame_2), colnames=colnames(tmp))

tmp <- flo_set[[1]]
col <- rep(1, nrow(tmp))
cols <- matrix(col, dimnames = list(NULL, "Timepoint"))
tmp <- fr_append_cols(tmp, cols)




name_frame <- as.matrix(flo_set[[1]]@parameters@data[["desc"]])

channel_names <- substr(name_frame, 7, nchar(name_frame))


## arcsinh transformation and column subsetting
fcs <- fsApply(flo_set, function(x, cofactor = 5){
  colnames(x) <- channel_names
  expr <- exprs(x)
  expr <- asinh(expr[, channel_names[c(3, 14:15, 23:57, 63, 65)]] / cofactor)
  exprs(x) <- expr
  x
})

fcs

# extract expression matrix

expr <- fsApply(fcs, exprs)
dim(expr)

#do 0_1 transform
rng <- colQuantiles(expr, probs = c(0.01, 0.99))
expr01 <- t((t(expr) - rng[, 1]) / (rng[, 2] - rng[, 1]))
expr01[expr01 < 0] <- 0
expr01[expr01 > 1] <- 1

md <- read_excel("metadata.xlsx", col_names = T)
sample_ids <- rep(md$file_name[21:25], fsApply(flo_set, nrow))

sample_ids <- substr(sample_ids, nchar(sample_ids)-9, nchar(sample_ids)-7)
time_point <- gsub("_1", "C+1", sample_ids)



ggdf <- data.frame(time_point = time_point, expr)
ggdf_long <- gather(ggdf, marker, expression, colnames(ggdf)[2:41])


# this makes a series of overlaid histograms showing expression of each marker at different timepoints.. couple of interesting plots but otherwise
# not very informative

# ggplot(ggdf_long, aes(x = expression, color = time_point, 
#                  group = time_point)) +
#   geom_density() +
#   facet_wrap(~ marker, scales = "free") +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1), 
#         strip.text = element_text(size = 7), axis.text = element_text(size = 5)) +
#   guides(color = guide_legend(ncol = 1))
  




# make a table of median expression values per timepoint for multidimensional scaling
expr_median_sample_tbl <- data.frame(sample_id = sample_ids, expr) %>%
  group_by(sample_id) %>% 
  summarize_all(list(median)) %>%
  ungroup()

#transpose expression matrix so rows are markers, kinda like long format i guess
expr_median_sample <- t(expr_median_sample_tbl[, -1])
colnames(expr_median_sample) <- expr_median_sample_tbl$sample_id

# creates object called mds that contains X and Y variables, some other infor
mds <- plotMDS(expr_median_sample, plot = FALSE)

# mut mds object as dataframe with sample_id column
ggdf <- data.frame(MDS1 = mds$x, MDS2 = mds$y, 
                   sample_id = colnames(expr_median_sample))

ggdf$volunteer <- "V03"

# this would be for when you have multiple people in the same dataset

# ggdf$volunteer <- as.character(ggdf$sample_id)
# ggdf$volunteer <- substr(ggdf$volunteer, nchar(ggdf$volunteer)-1, nchar(ggdf$volunteer))

#mm <- match(ggdf$sample_id, md$file_name)
#ggdf$volunteer <- md$volunteer[21:25]

ggplot(ggdf, aes(x = MDS1, y = MDS2, color = volunteer)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_label_repel(aes(label = sample_id)) +
  theme_bw()






# create fsom object that contains fcs (the 5 fcs files from Volunteer 03), and some metadata like column names, number of cells per timepoint etc
# set seed
rm(som,fsom)
set.seed(1234)
fsom <- ReadInput(fcs, transform = FALSE, scale = FALSE)


#actually build the SOM, choosing which channels to use for mapping
som <- FlowSOM(fcs, colsToUse = fcs@colnames[c(3:31, 33:40)], maxMeta=12)
som2 <- FlowSOM(fcs, colsToUse = fcs@colnames[c(3:31, 33:40)], maxMeta=12)

bonjour<- NewData(som, fcs[[1]])
table(som$data[,5] == som2$data[,5])

#### clustering on som; nmc = number of clusters

codes <- som$map$codes
plot_outdir <- "consensus_plots"
nmc <- 42
mc <- ConsensusClusterPlus(t(codes), maxK = nmc, reps = 100, 
                           pItem = 0.9, pFeature = 1, title = plot_outdir, plot = "png", 
                           clusterAlg = "hc", innerLinkage = "average", finalLinkage = "average", 
                           distance = "euclidean", seed = 1234)





## Get cluster ids for each cell
code_clustering1 <- mc[[nmc]]$consensusClass
cell_clustering1 <- code_clustering1[som$map$mapping[,1]]

# create a bunch of variable that contain how many cells are in each timepoint
for (i in seq(1, 5)){
  assign(paste("d", i, sep='_'), som$metaData[[i]][[2]]-som$metaData[[i]][[1]]+1)
}

# shove them into a vector in the right proportions
time_points <- rep(c("C-1", "C+8", "C+10", "DoD", "T+6"), times=c(d_1, d_2, d_3, d_4, d_5))


# this is a more biologically meaningful ordering of the channels contained in out CyTOF panel

marker_levels <- c("CD4",
                   "CD8",
                   "Vd2",
                   "TCRgd",
                   "Va7.2",
                   "CD38",
                   "HLA-DR",
                   "ICOS",
                   "CD28",
                   "PD1",
                   "TIM-3",
                   "CD95",
                   "BCL-2",
                   "CD27",
                   "Perforin",
                   "GZB",
                   "Tbet",
                   "Ki-67",
                   "CD127",
                   "IntegrinB7",
                   "CD56",
                   "CD16",
                   "CD161",
                   "CD49d",
                   "CD103",
                   "CD25",
                   "FoxP3",
                   "CD39",
                   "CTLA4",
                   "CLA",
                   "CXCR5",
                   "CD57",
                   "CD45RA",
                   "CD45RO",
                   "CCR7")




# make a data frame with cluster ID, timepoint and expression matrix for each cell
dr <- data.frame(cluster_id = code_clustering1[som$map$mapping[,1]], time_point = time_points, 
                 expr[, marker_levels])

dr_long <- gather(dr, Marker, Intensity, colnames(dr)[3:37])

#split by timepoint
dr_list <- split(dr, dr$time_point)

#1 = baseline, 2 = C+10, 3= C+8, 4=DoD, 5=T+6

# make a list of tables for cluster abundances at each timepoint
dr_list_tbl <- lapply(dr_list, function(x) table(x$cluster_id))

# sum all the cluster abundances and put them in a list (dictionary); kind of a sanity check
dr_list_counts <- lapply(dr_list_tbl, function(x) sum(x))

# convert absolute cell counts into percentages
dr_list_freq <- list()
for (i in seq(1,5)){
dr_list_freq[[i]] <- dr_list_tbl[[i]]/dr_list_counts[[i]]*100
}

# convert list of tables into list of data frames
dr_list_freq <- lapply(dr_list_freq, function(x) as.data.frame(x))

# make one dataframe that contains cluster abundances through time
edger_input1 <- data.frame(ClusterID= dr_list_freq[[1]]$Var, Baseline_1=dr_list_freq[[1]]$Freq, C_8_1=dr_list_freq[[3]]$Freq, C_10_1=dr_list_freq[[2]]$Freq, DoD_1=dr_list_freq[[4]]$Freq, T_6_1=dr_list_freq[[5]]$Freq)






####
####

#### do it again, do it again: for replicates

####
####
hi<-som$data[5:10,5]
rm(som,fsom)
fsom1 <- ReadInput(fcs, transform = FALSE, scale = FALSE)

#actually build the SOM, choosing which channels to use for mapping
som1 <- BuildSOM(fsom1, colsToUse = fcs@colnames[c(3:31, 33:40)])
hey<-som1$data[5:10,5]
# helo2 <- som1$map$sdValues



head(som1$map$mapping, n=20)
som1 <- BuildMST(som1,tSNE=TRUE)
PlotStars(som1, view='tSNE')
PlotMarker(som1,"HLA-DR")


length(som1$map$mapping)


codes <- som2$map$codes
plot_outdir <- "consensus_plots"
nmc <- 42
mc1 <- ConsensusClusterPlus(t(codes), maxK = nmc, reps = 100, 
                           pItem = 0.9, pFeature = 1, title = plot_outdir, plot = "png", 
                           clusterAlg = "hc", innerLinkage = "average", finalLinkage = "average", 
                           distance = "euclidean", seed = 1234)

# the output here is actually somewhat different!



## Get cluster ids for each cell
code_clustering2 <- mc1[[nmc]]$consensusClass
cell_clustering2 <- code_clustering2[som1$map$mapping[,1]]

for (i in seq(1, 5)){
  assign(paste("d", i, sep='_'), som1$metaData[[i]][[2]]-som1$metaData[[i]][[1]]+1)
}


# time_points <- rep(c("C-1", "C+8", "C+10", "DoD", "T+6"), times=c(d_1, d_2, d_3, d_4, d_5))

dr2 <- data.frame(cluster_id = code_clustering2[som1$map$mapping[,1]], time_point = time_points, 
                 expr[, marker_levels])

dr2_list <- split(dr2, dr2$time_point)

#1 = baseline, 2 = C+10, 3= C+8, 4=DoD, 5=T+6
dr2_list_tbl <- lapply(dr2_list, function(x) table(x$cluster_id))

dr2_list_counts <- lapply(dr2_list_tbl, function(x) sum(x))


dr2_list_freq <- list()
for (i in seq(1,5)){
  dr2_list_freq[[i]] <- dr2_list_tbl[[i]]/dr2_list_counts[[i]]*100
}

dr2_list_freq <- lapply(dr2_list_freq, function(x) as.data.frame(x))
edger_input2 <- data.frame(Baseline_2=dr2_list_freq[[1]]$Freq, C_8_2=dr2_list_freq[[3]]$Freq, C_10_2=dr2_list_freq[[2]]$Freq, DoD_2=dr2_list_freq[[4]]$Freq, T_6_2=dr2_list_freq[[5]]$Freq)


########       time to edge some garbage yo





short <- cbind(edger_input1, edger_input2)
short_long <- gather(short, Timepoint, Frequency, colnames(short)[2:11])

half <- filter(short_long, Timepoint %in% c("Baseline_1", "Baseline_2"))

ggplot(data=half, aes(x= as.factor(ClusterID), y=Frequency, fill=Timepoint))+
  geom_bar(stat="identity", position = "dodge")

# make up group matrix, matching "replicates" to each other. also remove _UMAP.fcs suffix
groups <- rep(substr(colnames(short)[1:12],1, nchar(colnames(short))-9), times=2)

# create design matrix & give it proper column names
design <- model.matrix(~0 + groups)
colnames(design) <- substr(colnames(short)[1:12],1, nchar(colnames(short))-9)

# create DGEList object, estimate disperson, fit GLM
y <- DGEList(short, group=groups, lib.size = colSums(short))
fit <- estimateDisp(y, design)
fit <- glmQLFit(fit, design, robust=TRUE) 

























##### figure code
my_palette <- c("#D53E4F","#D96459","#F2AE72","#588C73","#1A9CC7")

cluster_heatmap_theme <- function(){
  theme(panel.border = element_blank(),
        axis.text.y.left = element_text(size=35),
        axis.line.y.left = element_blank(),
        axis.line.y.right = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 33),
        #axis.text.y.right = element_text(size = 35),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.title = element_blank(),
        legend.position = "none",
        plot.title = element_text(size = 45, hjust = 0.5),
        plot.margin = unit(c(1,0,1,0), "cm"))
}


Volunteer_plots <- dr_long %>%
  group_by(time_point)%>%
  filter(Marker == "CD4") %>%
  arrange(desc(Intensity)) %>%
  cluster_levels <- do(unique(.$cluster_id))
  ungroup() %>%
  group_by(time_point) %>%
              do(ggsave(
              ggplot(data = ., aes(x = factor(cluster_id), y = factor(Marker), group=time_point))+
              geom_tile(aes(fill=Intensity), color="white")+
              scale_fill_gradientn(colors=rev(my_palette))+
              scale_y_discrete(position = "left")+
              xlab(NULL)+
              ggtitle(paste(.$time_point[1], "HeatMap", sep='_')),
              ggtitle(paste(i))+
              cluster_heatmap_theme(),
              filename = paste(.$time_point[1], "HeatMap", ".pdf", sep=''),
              device="pdf"))
            
  
  
   dr_long %>%
       group_by(time_point)%>%
       filter(Marker == "CD4") %>%
       arrange(desc(Intensity)) %>%
       assign("cluster_levels", unique(.$cluster_id)) %>%
       do(print(cluster_levels))
  
  
   dr_long %>%
           group_by(time_point)%>%
           filter(Marker == "CD4") %>%
           arrange(desc(Intensity)) %>%
           select(cluster_id) %>%
           {unique(.) ->> cluster_levels} %>%
           print(cluster_levels) %>%
           ungroup()
  


# for(i in unique(dr_long$time_point)){
#   # specific_levels <- NULL
#   # ifelse(i %in% c("02","06"), assign("result", element_text(size=35)), assign("result", element_blank()))
#   # ifelse(i %in% c("05","09"), assign("result1", "right"), assign("result1", "left"))
#   
#   sub_set <- filter(dr_long, time_point == i)
#   
#   specific_levels <- sub_set %>% 
#     filter(Marker == "CD4") %>%
#     arrange(desc(Intensity))
#   
#   specific_levels <- c(as.character(specific_levels$clusterID))
#   
#   assign(paste(i),
#          ggplot(data = sub_set, aes_(x = factor(sub_set$clusterID, levels = specific_levels), y = factor(sub_set$Marker, levels = rev(marker_levels)), group=sub_set$time_point))+
#            geom_tile(aes(fill=Intensity), color="white")+
#            scale_fill_gradientn(colors=rev(my_palette))+
#            scale_y_discrete(position = "left")+
#            xlab(NULL)+
#            ggtitle(paste(i))+
#            theme(panel.border = element_blank(),
#                  axis.text.y.left = element_text(size=35),
#                  axis.line.y.left = element_blank(),
#                  axis.line.y.right = element_blank(),
#                  axis.ticks.y = element_blank(),
#                  axis.title.y = element_blank(),
#                  axis.text.x = element_text(size = 33),
#                  #axis.text.y.right = element_text(size = 35),
#                  panel.grid.major = element_blank(),
#                  panel.grid.minor = element_blank(),
#                  axis.line = element_line(colour = "black"),
#                  legend.title = element_blank(),
#                  legend.position = "none",
#                  plot.title = element_text(size = 45, hjust = 0.5),
#                  plot.margin = unit(c(1,0,1,0), "cm"))
#   )
# } 
# 
# 
# 
# ggsave("sandbox.pdf", grid.arrange(Volunteer_02, Volunteer_03, Volunteer_05, Volunteer_06, Volunteer_07, Volunteer_09, ncol=3, nrow=2, layout_matrix = rbind(c(1,2,3),
#                                                                                                                                                              c(4,5,6))
# ),  width = 40, height = 40, limitsize = F)















expri <- cbind.data.frame(expr, time_points)
expri01 <- cbind.data.frame(expr01, time_points)

expri$time_points <- as.character(expri$time_points)
expri01$time_points <- as.character(expri01$time_points)

mm <- match(dr$sample_id, md$sample_id)
dr$condition <- md$condition[mm]







color_clusters <- rep(c("#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72", 
                    "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3", 
                    "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D", 
                    "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999", 
                    "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", "#808000", 
                    "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00"), times=2)




plot_clustering_heatmap_wrapper <- function(time_point, hypersine, zero_one, 
                                            cell_clustering, color_clusters, cluster_merging = NULL){
  
  # Calculate the median expression
  expr_median <- data.frame(hypersine, cell_clustering = cell_clustering1[match(time_point, expri$time_points):match(time_point, expri$time_points)+nrow(hypersine)]) %>%
    group_by(cell_clustering) %>% 
    summarize_if(is.numeric, funs(median))
    colnames(expr_median) <- c("cell_clustering", colnames(hypersine))  
    
  expr01_median <- data.frame(zero_one, cell_clustering = cell_clustering1[match(time_point, expri$time_points):match(time_point, expri$time_points)+nrow(hypersine)]) %>%
    group_by(cell_clustering) %>% 
    summarize_if(is.numeric, funs(median))
    colnames(expr01_median) <- c("cell_clustering", colnames(zero_one)) 
    
  # Calculate cluster frequencies
  clustering_table <- as.numeric(table(cell_clustering = cell_clustering1[match(time_point, expri$time_points):match(time_point, expri$time_points)+nrow(hypersine)]))
  
  # This clustering is based on the markers that were used for the main clustering
  d <- dist(expr_median[, colnames(hypersine)], method = "euclidean")
  cluster_rows <- hclust(d, method = "average")
  
  expr_heat <- as.matrix(expr01_median[, colnames(zero_one)])
  rownames(expr_heat) <- expr01_median$cell_clustering
  
  labels_row <- paste0(rownames(expr_heat), " (", 
                       round(clustering_table / sum(clustering_table) * 100, 2), "%)")
  labels_col <- colnames(expr_heat)
  
  # Row annotation for the heatmap
  annotation_row <- data.frame(cluster = factor(expr01_median$cell_clustering))
  rownames(annotation_row) <- rownames(expr_heat)
  
  color_clusters <- color_clusters[1:nlevels(annotation_row$cluster)]
  names(color_clusters) <- levels(annotation_row$cluster)
  annotation_colors <- list(cluster = color_clusters)
  annotation_legend <- FALSE
  
  if(!is.null(cluster_merging)){
    cluster_merging$new_cluster <- factor(cluster_merging$new_cluster)
    annotation_row$cluster_merging <- cluster_merging$new_cluster 
    color_clusters <- color_clusters[1:nlevels(cluster_merging$new_cluster)]
    names(color_clusters) <- levels(cluster_merging$new_cluster)
    annotation_colors$cluster_merging <- color_clusters
    annotation_legend <- TRUE
  }
  
  # Colors for the heatmap
  color <- colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlBu")))(100)
  
  pheatmap(expr_heat, color = color, 
           cluster_cols = FALSE, cluster_rows = cluster_rows, 
           labels_col = labels_col, labels_row = labels_row, 
           display_numbers = TRUE, number_color = "black", 
           fontsize = 8, fontsize_number = 4,
           annotation_row = annotation_row, annotation_colors = annotation_colors, 
           annotation_legend = annotation_legend)
  
}


baseline_heatmap <- plot_clustering_heatmap_wrapper(time_point = "DoD",
                                                    hypersine = expri %>%
                                                      filter(., time_points=="DoD") %>%
                                                      select(-time_points),
                                                    zero_one =  expri01 %>%
                                                      filter(., time_points=="DoD") %>%
                                                      select(-time_points),
                                                    cell_clustering = cell_clustering1[match(time_point, expri$time_points):match(time_point, expri$time_points)+nrow(hypersine)],
                                                    color_clusters = color_clusters)
#33699

DoD_heatmap <- plot_clustering_heatmap_wrapper(hypersine = expri %>%
                                                 filter(., time_points=="DoD") %>%
                                                 select(-time_points),
                                               zero_one =  expri01 %>%
                                                 filter(., time_points=="DoD") %>%
                                                 select(-time_points),
                                               cell_clustering = cell_clustering1[sum(d_1, d_2, d_3):sum(d_1, d_2, d_3, d_4)],
                                               color_clusters = color_clusters)

Treat6_heatmap <- plot_clustering_heatmap_wrapper(hypersine = expri %>%
                                                    filter(., time_points=="T+6") %>%
                                                    select(-time_points),
                                                  zero_one =  expri01 %>%
                                                    filter(., time_points=="T+6") %>%
                                                    select(-time_points),
                                                  cell_clustering = cell_clustering1[sum(d_1, d_2, d_3, d_4):length(cell_clustering1)],
                                                  color_clusters = color_clusters)





















