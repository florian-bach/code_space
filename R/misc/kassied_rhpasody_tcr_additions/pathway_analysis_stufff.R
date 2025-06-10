library(clusterProfiler)
organism = "org.Hs.eg.db"
library(organism, character.only = TRUE)
library(tidyr)
library(dplyr)


#irbc vs day0 ####

cd4_irbc_vs_day0 <- read.csv("/Volumes/lab_prasj/BIG_Flo/kassie_bd_rhapsody/cd4cell_changes_culture3.csv")
cd4_irbc_vs_day0$cell_type <- "CD4" 

cd8_irbc_vs_day0 <- read.csv("/Volumes/lab_prasj//BIG_Flo/kassie_bd_rhapsody/cd8_changes_culture3.csv")
cd8_irbc_vs_day0$cell_type <- "CD8" 

b_irbc_vs_day0 <- read.csv("/Volumes/lab_prasj/BIG_Flo/kassie_bd_rhapsody/bcell_changes_culture3.csv")
b_irbc_vs_day0$cell_type <- "B" 

combo_list <- bind_rows(cd4_irbc_vs_day0, cd8_irbc_vs_day0, b_irbc_vs_day0)%>%
  mutate(gene=X)

# cell_type `sum(p_val_adj < 0.05)`
# <chr>                       <int>
# 1 B                            1028
# 2 CD4                          1021
# 3 CD8                          1270

## urbc vs day0 ####
cd4_urbc_vs_day0 <- read.csv("/Volumes/lab_prasj/BIG_Flo/kassie_bd_rhapsody/cd4cell_changes_culture2.csv")
cd4_urbc_vs_day0$cell_type <- "CD4" 

cd8_urbc_vs_day0 <- read.csv("/Volumes/lab_prasj//BIG_Flo/kassie_bd_rhapsody/cd8_changes_culture2.csv")
cd8_urbc_vs_day0$cell_type <- "CD8" 

b_urbc_vs_day0 <- read.csv("/Volumes/lab_prasj/BIG_Flo/kassie_bd_rhapsody/bcell_changes_culture2.csv")
b_urbc_vs_day0$cell_type <- "B" 


combo_list <- bind_rows(cd4_urbc_vs_day0, cd8_urbc_vs_day0, b_urbc_vs_day0)%>%
  mutate(gene=X)
# cell_type `sum(p_val_adj < 0.05)`
# <chr>                       <int>
# 1 B                            1069
# 2 CD4                          1017
# 3 CD8                          1256
## irbc vs urbc ####
cd4_irbc_vs_urbc <- read.csv("/Volumes/lab_prasj/BIG_Flo/kassie_bd_rhapsody/cd4cell_changes_culture.csv")
cd4_irbc_vs_urbc$cell_type <- "CD4" 

cd8_irbc_vs_urbc <- read.csv("/Volumes/lab_prasj//BIG_Flo/kassie_bd_rhapsody/cd8_changes_culture.csv")
cd8_irbc_vs_urbc$cell_type <- "CD8" 

b_irbc_vs_urbc <- read.csv("/Volumes/lab_prasj/BIG_Flo/kassie_bd_rhapsody/bcell_changes_culture.csv")
b_irbc_vs_urbc$cell_type <- "B" 

combo_list <- bind_rows(cd4_irbc_vs_urbc, cd8_irbc_vs_urbc, b_irbc_vs_urbc)%>%
  mutate(gene=X)

# cell_type `sum(p_val_adj < 0.05)`
# <chr>                       <int>
# 1 B                              77
# 2 CD4                            16
# 3 CD8                            55




# cluster profiler ####

top100 <- combo_list %>%
  filter(p_val_adj<0.05)%>%
  group_by(cell_type)%>%
  mutate(abs_diff=abs(pct.1-pct.2))%>%
  arrange(desc(abs_diff))%>%
  filter(!duplicated(gene))%>%
  ungroup()



list_of_top100 <- split(top100, top100$cell_type)
dedup_list_of_genes <- lapply(list_of_top100, function(x) subset(x, !duplicated(x$gene)))


significant_KEGGs <- vector(length = length(list_of_top100), mode = "list")
names(significant_KEGGs) <- names(dedup_list_of_genes)

significant_GOs <- vector(length = length(list_of_top100), mode = "list")
names(significant_GOs) <- names(dedup_list_of_genes)

for(i in names(significant_GOs)){
  
  print(i)
  
  # translate gene symbold to ENTREZ
  ids <- bitr(dedup_list_of_genes[[i]]$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = organism)
  # remove duplicate IDs
  dedup_ids = ids[!duplicated(ids$SYMBOL),]
  
  # Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
  df2 = dedup_list_of_genes[[i]][dedup_list_of_genes[[i]]$gene %in% dedup_ids$SYMBOL,]
  
  # Create a new column in df2 with the corresponding ENTREZ IDs
  df2$Y = dedup_ids$ENTREZID[match(df2$gene, dedup_ids$SYMBOL)]
  
  # Create a vector of the gene unuiverse
  kegg_gene_list <- df2$avg_log2FC
  
  # Name vector with ENTREZ ids
  names(kegg_gene_list) <- df2$Y
  
  # omit any NA values 
  kegg_gene_list<-na.omit(kegg_gene_list)
  
  # sort the list in decreasing order (required for clusterProfiler)
  kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)

  kegg_organism = "hsa"
  kk2 <- gseKEGG(geneList     = kegg_gene_list,
                 organism     = kegg_organism,
                 # nPerm        = 10000,
                 # minGSSize    = 15,
                 # maxGSSize    = 500,
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH",
                 scoreType = "pos",
                 keyType       = "ncbi-geneid")
  
  significant_KEGGs[[i]] <- kk2@result$Description
  
  kk2_go <- gseGO(
    kegg_gene_list,
    ont = "ALL",
    scoreType = "pos",
    OrgDb=organism,
    exponent = 1,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH"
  )
  
  significant_GOs[[i]] <- kk2_go@result$Description
  
}



# topGO ####

b_irbc_vs_day0 <- read.csv("/Volumes/lab_prasj/BIG_Flo/kassie_bd_rhapsody/nkcell_changes_culture2.csv")

top100 <- b_irbc_vs_day0 %>%
  mutate(gene=X, abs_diff=abs(pct.1-pct.2))%>%
  filter(p_val_adj<0.05, abs_diff>0.2 | abs(avg_log2FC)>0.5)%>%
  filter(!duplicated(gene))%>%
  arrange(desc(abs(avg_log2FC)))%>%
  ungroup()


# gene list should be an array of p values (or fc / diff?) with entrezid names
symbol_entrez <- bitr(top100$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = organism)
symbol_entrez$avg_log2FC <- top100$avg_log2FC[match(symbol_entrez$SYMBOL, top100$gene)]

geneList <- symbol_entrez$avg_log2FC
names(geneList) <-symbol_entrez$ENTREZID
# ordered_geneList <- sort(geneList)

# Create topGOData object
GOdata <- new("topGOdata",
              ontology = "BP",
              allGenes = geneList,
              geneSelectionFun = function(x)x,
              annot = annFUN.org , mapping = "org.Hs.eg.db")

resultKS <- runTest(GOdata, algorithm = "weight01", statistic = "ks")
tab <- GenTable(GOdata, raw.p.value = resultKS, topNodes = length(resultKS@score), numChar = 120)
tab$padj=p.adjust(tab$raw.p.value, method = "fdr")
tab$Term[tab$padj<0.1]
