library(data.table)
library(dplyr)
library(tidyr)


setwd("/home/flobuntu/PhD/RNAseq/vac69a/all/xls")

all_files <- list.files(pattern = "*xls")
names(all_files) <- c("C14_Baseline","DoD_Baseline","C56_Baseline","T6_Baseline","T6_DoD")

list_of_files <- lapply(all_files, function(x)fread(x, header = T, stringsAsFactors = F))
list_of_sig <- lapply(list_of_files, function(x) x %>%
                        filter(padj<0.05) %>%
                        mutate("FC"=ifelse(log2FoldChange>0, "up", "down"))
)
                      

list_of_sig_unique <- lapply(list_of_sig, function(x)
      x %>%
     group_by(Description) %>%
     top_n(n = -1, wt = padj)
     )

list_of_sig_unique_FC <- lapply(list_of_sig_unique, function(x) filter(x, abs(log2FoldChange)>=log2(1.5) & !is.na(Symbol)))



list_of_sig_unique_FC_named <- lapply(names(list_of_sig_unique_FC), function(x) {
  list_of_sig_unique_FC[[x]] %>%
    mutate(file_name = x)
});names(list_of_sig_unique_FC_named) <- names(list_of_sig_unique_FC)


setwd("gene_lists")


lapply(names(list_of_sig_unique_FC_named), function(x)
  write.table(list_of_sig_unique_FC_named[[x]]$Symbol, paste(x, "_ALL_significant_symbol_only.txt", sep=''), sep = "\t", quote = F, row.names = F, col.names = F)
)

UP_list_of_sig_unique_FC_named <- lapply(list_of_sig_unique_FC_named, function(x) subset(x, log2FoldChange>log2(1.5)))

lapply(names(UP_list_of_sig_unique_FC_named), function(x)
  write.table(UP_list_of_sig_unique_FC_named[[x]]$Symbol, paste(x, "_UP_significant_symbol_only.txt", sep=''), sep = "\t", quote = F, row.names = F, col.names = F)
)


DOWN_list_of_sig_unique_FC_named <- lapply(list_of_sig_unique_FC_named, function(x) subset(x, log2FoldChange<log2(1.5)))

lapply(names(DOWN_list_of_sig_unique_FC_named), function(x)
  write.table(DOWN_list_of_sig_unique_FC_named[[x]]$Symbol, paste(x, "_DOWN_significant_symbol_only.txt", sep=''), sep = "\t", quote = F, row.names = F, col.names = F)
)





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



C14_Baseline <- list_of_sig_unique_FC_named[[1]]
DoD_Baseline <- list_of_sig_unique_FC_named[[2]]
C56_Baseline <- list_of_sig_unique_FC_named[[3]]
T6_Baseline <- list_of_sig_unique_FC_named[[4]]
T6_DoD <- list_of_sig_unique_FC_named[[5]]

#volcano
list_of_named_files <- lapply(names(list_of_files), function(x) {
  list_of_files[[x]] %>% 
    mutate("file_name" = x) %>%
    mutate("FC"=ifelse(log2FoldChange>0, "up", "down"))
})

list_of_named_unique <- lapply(list_of_named_files, function(x)
  x %>%
    group_by(Symbol) %>%
    top_n(n = -1, wt = padj)
)

list_of_named_unique_no_na <- lapply(list_of_named_unique, function(x)
  filter(x, !is.na(Symbol))
)


big_table <- do.call(rbind, list_of_named_unique_no_na)

fwrite(big_table, "all_unique_genes_cleaned.csv")


library(ggplot2)

comp_levels <- c("C14_Baseline", "DoD_Baseline", "T6_DoD", "T6_Baseline", "C56_Baseline")

(all_volcanoes <- ggplot(big_table, aes(x=log2FoldChange, y=-log10(padj)))+
  geom_point(shape="o", aes(color=ifelse(abs(log2FoldChange)>log2(1.5) & padj<0.05, "red", "black")))+
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

try <- lapply(list_of_sig_unique_FC_named, function(x) split(x, x$FC))

sig_gene_counts <- data.frame(
    lapply(try, function(x)
    lapply(x, function(y) nrow(y))
    )
  )

sig_gene_counts <- gather(sig_gene_counts, Comparison, DE_Genes)
sig_gene_counts$Direction <- ifelse(grepl("down", sig_gene_counts, fixed=T), "down", "up")
sig_gene_counts$Comparison <- substr(sig_gene_counts$Comparison, 1, nchar(sig_gene_counts$Comparison)-ifelse(sig_gene_counts$Direction=="down", 5, 3))
sig_gene_counts$DE_Genes <- ifelse(sig_gene_counts$Direction=="down", -sig_gene_counts$DE_Genes, sig_gene_counts$DE_Genes)


(sig_gene_count_plots <- ggplot(sig_gene_counts, aes(x=factor(Comparison, levels=comp_levels), y=DE_Genes))+
  geom_bar(stat="identity", aes(fill=Comparison))+
  geom_text(aes(label=abs(DE_Genes), vjust= -0.2), data = subset(sig_gene_counts, sig_gene_counts$Direction=="up"))+
  geom_text(aes(label=abs(DE_Genes), vjust= 1.2), data = subset(sig_gene_counts, sig_gene_counts$Direction=="down"))+
  theme_minimal()+
  scale_fill_manual(values=colorspace::sequential_hcl(6, palette = "Purple Yellow"))+
  ylab("Number of Differentially Expressed Genes\n")+
  theme(legend.position = "none",
        axis.title.x=element_blank(),
        axis.text.x = element_text(angle=45, hjust=1)))

ggsave("./figures/sig_gene_count_plots_log2fc058.png", sig_gene_count_plots)


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


