library(fgsea)
library(dplyr)
library(ggplot2)
library(msigdbr)

m_df<- msigdbr(species = "Homo sapiens", category = "C7")
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

#pbmc.genes <- read.csv("/labs/prasj/BIG_Flo/kassie_bd_rhapsody/cluster_markers.csv")
pbmc.genes <- read.csv("/Volumes/lab_prasj/BIG_Flo/kassie_bd_rhapsody/cluster_markers.csv")

irbc = c(0, 1, 3)
urbc = c(2,5,6)
d0 = c(4, 7, 14, 15)

top100 <- pbmc.genes %>%
  mutate("specific"=if_else(cluster %in% irbc, "irbc", if_else(cluster %in% urbc, "urbc", if_else(cluster %in% d0, "d0", "neither"))))%>%
  filter(specific %in% c("irbc", "urbc", "d0"), p_val_adj<0.1)%>%
  group_by(specific)%>%
  arrange(desc(avg_log2FC), desc(p_val_adj))%>%
  filter(!duplicated(gene))%>%
  ungroup()
  

# for(i in unique(top100$specific)){
  
cluster.genes <- top100 %>%
  dplyr::filter(specific == "irbc") %>%
  # arrange(desc(avg_log2FC), p_val_adj)%>%
  dplyr::select(gene, p_val_adj)

ranks <- tibble::deframe(cluster.genes)

system.time(
  fgseaRes <- fgsea(fgsea_sets[1:10],
                  stats = ranks,
                  minSize  = 15,
                  maxSize  = 500)
)

