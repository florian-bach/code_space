caro_data <- readRDS("~/Downloads/NULISA_genes_expression.rds")
dim(caro_data$pheno)
dim(caro_data$expr)

long_expr <- data.frame(caro_data$expr)%>%
  mutate(gene=rownames(.))%>%
  pivot_longer(cols = names(data.frame(caro_data$expr)), names_to = "novogene code", values_to = "reads")%>%
  pivot_wider(names_from = "gene", values_from = reads)

expr_frame <- caro_data$pheno%>%
  select(cohortid, infection, timepoint_category, `novogene code`)%>%
  left_join(., long_expr, by="novogene code")%>%
  pivot_longer(cols = names(long_expr)[2:length(names(long_expr))], names_to = "targetName", values_to = "RNA")

nulisa_data <- read.csv("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/clean_musical_combo_with_metadata.csv")

nulisa_data <- nulisa_data %>%
  mutate(timepoint = factor(timepoint, levels=c("bad_baseline", "baseline", "day0", "day7", "day14", "day28")))%>%
  mutate(protein=concentration)%>%
  mutate(cohortid=id,  infection=infectiontype,  timepoint_category=timepoint)


combo_data <- expr_frame%>%
  mutate(timepoint_category=case_when(timepoint_category=="Baseline"~"baseline",
                                      timepoint_category=="Day 0"~"day0",
                                      timepoint_category=="Day 14"~"day14",
                                      timepoint_category=="Day 28"~"day28",
                                      timepoint_category=="Day 7"~"day7"))%>%
  inner_join(., nulisa_data, by=c("cohortid", "infection", "timepoint_category", "targetName"))%>%
  mutate(infection_time=paste(infection, timepoint_category, sep = "_"))%>%
  filter(infection %in% c("A", "S"), timepoint_category %in% c("baseline", "day0", "day14"))


combo_data%>%
  filter(targetName%in%c("TLR3", "SELE", "S100A9"))%>%
  ggplot(., aes(x=RNA, y=protein, color=infection_time))+
  geom_smooth(method="lm", alpha = 0.15)+
  ggpubr::stat_cor(method = "spearman")+
  geom_point()+
  facet_wrap(~targetName, scales="free")+
  theme_minimal()
 

combo_data%>%
  filter(targetName%in%c("CD4", "IFNW1", "IL10", "IL5RA", "LAG3", "SDC1", "TLR3", "VCAM1"))%>%
  ggplot(., aes(x=RNA, y=protein, color=infection_time))+
  geom_smooth(method="lm", alpha = 0.15)+
  # ggpubr::stat_cor(method = "spearman")+
  geom_point()+
  facet_wrap(~targetName, scales="free")+
  scale_color_manual(values=colorspace::diverging_hcl(n = 5, palette = "Berlin")[c(3,2,1,4,5)])+
  theme_minimal()

top9_genes <- coexpression_purf%>%
  ungroup()%>%
  slice_max(n = 9, order_by = padj)


combo_data%>%
  filter(targetName%in%top9_genes$targetName)%>%
  ggplot(., aes(x=RNA, y=protein, color=infection_time))+
  geom_smooth(method="lm", alpha = 0.15)+
  # ggpubr::stat_cor(method = "spearman")+
  geom_point()+
  facet_wrap(~targetName, scales="free")+
  scale_color_manual(values=colorspace::diverging_hcl(n = 5, palette = "Berlin")[c(3,2,1,4,5)])+
  theme_minimal()

combo_data%>%
  filter(targetName%in%para_sigs$targetName)%>%
  ggplot(., aes(x=RNA, y=protein, color=infection_time))+
  geom_smooth(method="lm", alpha = 0.15)+
  # ggpubr::stat_cor(method = "spearman", aes(x=RNA, y=protein), inherit.aes = F)+
  ggpubr::stat_cor(method = "spearman", size=2)+
  geom_point()+
  facet_wrap(~targetName, scales="free")+
  scale_color_manual(values=colorspace::diverging_hcl(n = 5, palette = "Berlin")[c(3,2,1,4,5)])+
  theme_minimal()



coexpression_purf <- combo_data%>%
  group_by(targetName, infection_time)%>%
  nest()%>%
  mutate(spearman=map(data, ~cor.test(x = .$RNA, y = .$protein, method="spearman")))%>%
  mutate(rho=map_dbl(spearman, ~.$estimate))%>%
  mutate(p=map_dbl(spearman, ~.$p.value))%>%
  group_by(infection_time)%>%
  mutate(padj=p.adjust(p, method="fdr"))
  
top9_genes <- coexpression_purf%>%
  ungroup()%>%
  slice_min(n = 9, order_by = padj)




para_sigs <- para_regression%>%
  filter(padj<0.1)%>%
  distinct(targetName)


nulisa_data%>%
  ungroup()%>%
  filter(targetName %in% top9_genes$targetName)%>%
  ggplot(., aes(x=concentration, y=IL10, color=infectiontype))+
  geom_smooth(method="lm", alpha = 0.15)+
  # # ggpubr::stat_cor(method = "spearman", aes(x=RNA, y=protein), inherit.aes = F)+
  # ggpubr::stat_cor(method = "spearman", size=2)+
  geom_point()+
  # scale_color_manual(values=colorspace::diverging_hcl(n = 5, palette = "Berlin")[c(3,2,1,4,5)])+
  theme_minimal()
  
  
