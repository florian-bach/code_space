library(data.table)
library(dplyr)
library(tidyr)


setwd("/home/flobuntu/PhD/RNAseq/vac69a/all/xls")

all_files <- list.files(pattern = "*xls")
names(all_files) <- c("C14_Baseline","DoD_Baseline","C56_Baseline","T6_Baseline","T6_DoD")

list_of_files <- lapply(all_files, function(x)fread(x, header = T, stringsAsFactors = F))
list_of_sig <- lapply(list_of_files, function(x) filter(x, padj<0.05))
list_of_sig_unique <- lapply(list_of_sig, function(x)
      x %>%
     group_by(Description) %>%
     top_n(n = -1, wt = padj)
     )

list_of_sig_unique_FC <- lapply(list_of_sig_unique, function(x) filter(x, abs(log2FoldChange)>=1))


#DE_counts with FC correction 
# C14_Baseline DoD_Baseline C56_Baseline  T6_Baseline       T6_DoD 
# 141         1318           15          228         1051

# No FC correction
# C14_Baseline DoD_Baseline C56_Baseline  T6_Baseline       T6_DoD 
# 294         3795           19          446         3057 
# 
# $DoD_Baseline
# [1] 3795
# 
# $C56_Baseline
# [1] 19
# 
# $T6_Baseline
# [1] 446
# 
# $T6_DoD
# [1] 3057



C14_Baseline <- list_of_sig_unique_FC[[1]]
DoD_Baseline <- list_of_sig_unique_FC[[2]]
C56_Baseline <- list_of_sig_unique_FC[[3]]
T6_Baseline <- list_of_sig_unique_FC[[4]]
T6_DoD <- list_of_sig_unique_FC[[5]]

#volcano
list_of_named_files <- lapply(names(list_of_files), function(x) {
  list_of_files[[x]] %>% 
    mutate(file_name = x)
})

list_of_named_unique <- lapply(list_of_named_files, function(x)
  x %>%
    group_by(Description) %>%
    top_n(n = -1, wt = padj)
)

big_table <- do.call(rbind, list_of_named_unique)

library(ggplot2)

comp_levels <- c("C14_Baseline", "DoD_Baseline", "T6_DoD", "T6_Baseline", "C56_Baseline")

(all_volcanoes <- ggplot(big_table, aes(x=log2FoldChange, y=-log10(padj)))+
  geom_point(shape="o", aes(color=ifelse(abs(log2FoldChange)>1 & padj<0.05, "red", "black")))+
  theme_minimal()+
  #scale_y_continuous(limits = c(0, 0.1))+
  scale_color_manual(values = c("red"="red", "black"="black"))+
  facet_wrap(~factor(file_name, levels = comp_levels),
             ncol=5, scales="fixed")+
  theme(legend.position = "none"))

ggsave("./figures/all_volcanoes.png", all_volcanoes, width=12, height=5)


# sig_gene_counts <- rbind(data.frame(lapply(list_of_sig_unique_FC, nrow), "FC"="With Correction"),
#   data.frame(lapply(list_of_sig_unique, nrow), "FC"="No Correction")
# )
sig_gene_counts <- data.frame(lapply(list_of_sig_unique_FC, nrow))
sig_gene_counts <- gather(sig_gene_counts, Comparison, DE_Genes, colnames(sig_gene_counts)[1:5])

(sig_gene_count_plots <- ggplot(sig_gene_counts, aes(x=factor(Comparison, levels=comp_levels), y=DE_Genes))+
  geom_bar(stat="identity", aes(fill=Comparison))+
  geom_text(aes(label=DE_Genes), vjust = -0.2)+
  theme_minimal()+
  ylab("Number of Differentially Expressed Genes\n")+
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x = element_text(angle=45, hjust=1))+
  coord_cartesian(ylim = c(0,1500)))

ggsave("./figures/sig_gene_count_plots.png", sig_gene_count_plots, width=4.5, height=4)

# simple heatmaps ####



ils <- subset(DoD_Baseline, grepl("interleukin", DoD_Baseline$Description))
ils <- ils[!is.na(ils$Symbol),]
ils <- ils[order(ils$log2FoldChange),]

cytokine <- subset(DoD_Baseline, grepl("cytokine", DoD_Baseline$Description))
cytokine <- cytokine[!is.na(cytokine$Symbol),]
cytokine <- cytokine[order(cytokine$log2FoldChange),]



chemokines <- subset(DoD_Baseline, grepl("chemokine", DoD_Baseline$Description))
chemokines <- chemokines[!is.na(chemokines$Symbol),]
chemokines <- chemokines[order(chemokines$log2FoldChange),]


T_cells <- subset(DoD_Baseline, grepl("T cell", DoD_Baseline$Description))


ggplot(chemokines, aes(x="", y=factor(Symbol, levels=chemokines$Symbol)))+
  geom_tile(aes(fill=log2FoldChange))+
  scale_fill_gradient2(midpoint = 0, low="#0859C6", high="#FFA500", mid="black")+
  theme_void()+
  ggtitle("Chemokines at DoD Relative to Baseline")+
  theme(axis.text.y = element_text(),
        legend.title = element_text(hjust=0.5),
        plot.title = element_text(hjust=0.5))


ggplot(ils, aes(x="", y=factor(Symbol, levels=ils$Symbol)))+
  geom_tile(aes(fill=log2FoldChange))+
  scale_fill_gradient2(midpoint = 0, low="#0859C6", high="#FFA500", mid="black")+
  theme_void()+
  ggtitle("Interleukins at DoD Relative to Baseline")+
  theme(axis.text.y = element_text(),
        legend.title = element_text(hjust=0.5),
        plot.title = element_text(hjust=0.5))

ggplot(cytokine, aes(x="", y=factor(Symbol, levels=cytokine$Symbol)))+
  geom_tile(aes(fill=log2FoldChange))+
  scale_fill_gradient2(midpoint = 0, low="#0859C6", high="#FFA500", mid="black")+
  theme_void()+
  ggtitle("Interleukins at DoD Relative to Baseline")+
  theme(axis.text.y = element_text(),
        legend.title = element_text(hjust=0.5),
        plot.title = element_text(hjust=0.5))


# venn diagrams ####

library(VennDiagram)


countz <- list(T6_Baseline$FeatureID, DoD_Baseline$FeatureID)



venn.diagram( x = countz,
              category.names = c("T6_Baseline" , "DoD_Baseline"),
              filename = "./figures/T6_DoD_Baseline_Venn.png",
              output = TRUE ,
              imagetype="png" ,
              height = 980 , 
              width = 1280 , 
              resolution = 300,
              cat.cex = 0.6,
              cat.fontface = "bold",
              cat.default.pos = "outer",
              cat.pos = c(1, 1),
              cat.dist = c(0.095, 0.055)
              )

countz <- list(C14_Baseline$FeatureID, DoD_Baseline$FeatureID)

venn.diagram( x = countz,
              category.names = c("C14_Baseline" , "DoD_Baseline"),
              filename = "./figures/C14_DoD_Baseline_Venn.png",
              output = TRUE ,
              imagetype="png" ,
              height = 980 , 
              width = 1280 , 
              resolution = 300,
              cat.cex = 0.6,
              cat.fontface = "bold",
              cat.default.pos = "outer",
              cat.pos = c(1, 1),
              cat.dist = c(0.095, 0.055)
              )


countz <- list(C56_Baseline$FeatureID, DoD_Baseline$FeatureID)

venn.diagram( x = countz,
              category.names = c("C56_Baseline" , "DoD_Baseline"),
              filename = "./figures/C56_DoD_Baseline_Venn.png",
              output = TRUE ,
              imagetype="png" ,
              height = 980 , 
              width = 1280 , 
              resolution = 300,
              cat.cex = 0.6,
              cat.fontface = "bold",
              cat.default.pos = "outer",
              )

countz <- list(C56_Baseline$FeatureID, DoD_Baseline$FeatureID, C14_Baseline$FeatureID)


venn.diagram( x = countz,
              category.names = c("C56_Baseline" , "DoD_Baseline", "C14_Baseline"),
              filename = "./figures/C56_DoD_Baseline_Venn.png",
              output = TRUE ,
              imagetype="png" ,
              height = 980 , 
              width = 1280 , 
              resolution = 300,
              cat.cex = 0.6,
              cat.fontface = "bold",
              cat.default.pos = "outer",
)


