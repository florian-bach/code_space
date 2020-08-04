
# CyTOF MDS ####

# cluster frequency

cytof_mds <- data.frame(cmdscale(robCompositions::aDist(t(wide_cytof))))

colnames(cytof_mds) <- c("MDS1", "MDS2")

cytof_mds$Sample_ID <- rownames(cytof_mds)
cytof_mds$Timepoint <- substr(cytof_mds$Sample_ID, 5, nchar(cytof_mds$Sample_ID))
cytof_mds$Volunteer <- substr(cytof_mds$Sample_ID, 1, 3)


my_paired_palette <- c("#FB9A99","#E31A1C","#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C")
time_col=colorspace::sequential_hcl(5, palette = "Purple Yellow")

time_palette <- c("Baseline"=time_col[4], "C10"=time_col[3], "DoD"=time_col[2], "T6"=time_col[1])

ggplot(cytof_mds, aes(x=MDS1, y=MDS2))+
  geom_point(aes(shape=Timepoint, colour=Volunteer))+
  theme_minimal()+
  scale_color_manual(values = my_paired_palette)
                      
ggsave("~/PhD/multi_omics/Aitchison's_CyTOF.png")     


#median marker expression

cs_by_s <- split(seq_len(ncol(merged_daf)), merged_daf$sample_id)
es <- as.matrix(assay(merged_daf, "exprs"))
ms <- vapply(cs_by_s, function(cs) rowMedians(es[, cs, drop = FALSE]), 
             numeric(nrow(merged_daf)))
rownames(ms) <- rownames(merged_daf)
ms <- subset(ms, rownames(ms)%in%c(state_markers(merged_daf)))

mds <- limma::plotMDS(ms, plot = FALSE)
df <- data.frame(MDS1 = mds$x, MDS2 = mds$y)
md <- metadata(merged_daf)$experiment_info
m <- match(rownames(df), md$sample_id)
df <- data.frame(df, md[m, ])

df2 <- filter(df, timepoint %in% c("Baseline", "DoD"))

# \
  #ggrepel::geom_label_repel(aes_string(label = "sample_id"), show.legend = FALSE)+ 
  ggplot()+
  geom_point(df, aes_(x=df$MDS1, y=df$MDS2, color=df$volunteer))+  
  geom_segment(aes_(x=df2$MDS1[grepl("Baseline", df2$timepoint)],
                   xend=df2$MDS1[grepl("DoD", df2$timepoint)],
                   y=df2$MDS2[grepl("Baseline", df2$timepoint)],
                   yend=df2$MDS2[grepl("DoD", df2$timepoint)]),
                   arrow =arrow(length = unit(0.2, "cm")))+
  theme_minimal()+
  scale_color_manual(values = my_paired_palette)


df2 <- df %>%
  group_by(timepoint) %>%
  mutate(paste(unique(volunteer), unique(timepoint), sep="_"))
  #pivot_wider(values_from = c(MDS1,MDS2), names_from = c(timepoint, volunteer)) %>%
  ungroup() %>%
  gather(Position, Position_Name, colnames(.)[5:ncol(.)])





# plasma mds ####

plasma_mds <- data.frame(cmdscale(robCompositions::aDist(t(ms))))

colnames(plasma_mds) <- c("MDS1", "MDS2")

plasma_mds$Sample_ID <- rownames(plasma_mds)
plasma_mds$Timepoint <- substr(plasma_mds$Sample_ID, 5, nchar(plasma_mds$Sample_ID))
plasma_mds$Volunteer <- substr(plasma_mds$Sample_ID, 1, 3)

