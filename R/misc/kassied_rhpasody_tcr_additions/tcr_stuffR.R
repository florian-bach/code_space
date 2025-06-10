library(scRepertoire)

gd_only <- readRDS("/Volumes/lab_prasj/BIG_Flo/kassie_bd_rhapsody/gd_only_combined.TCR")
# gd_only <- readRDS("/Volumes/lab_prasj/BIG_Flo/kassie_bd_rhapsody/gd_only.RDS")
gd_present_in_both <- read.csv("/Volumes/lab_prasj/BIG_Flo/kassie_bd_rhapsody/gd_clones_in_both_stims.csv")
# list_repertoire <- readRDS("/Volumes/lab_prasj/BIG_Flo/kassie_bd_rhapsody/list_repertoire.RDS")
# 
gd_only.TCR <- combineExpression(list_repertoire,
                                 gd_only,
                                  cloneCall="aa",
                                  group.by = "sample_id",
                                  proportion = T)
gd_only.TCR <- subset(gd_only.TCR, stim!="PBMC")

alluvial <- clonalCompare(gd_only.TCR, 
              group.by = "stim",
              clones = gd_present_in_both$CTaa,
              cloneCall="aa", 
              graph = "alluvial")+theme(legend.position = "none")

scRepertoire::alluvialClones(gd_only.TCR,
                   y.axes=c('donor'), cloneCall = "CTaa")
