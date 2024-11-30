library(emmeans)
library(purrr)
library(tidyr)
library(dplyr)
library(ggplot2)

`%notin%` <- Negate(`%in%`)

# read clean data ####
clean_data <- read.csv("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/clean_musical_combo_with_metadata.csv")

nulisa_analytes <- unique(clean_data$targetName)

slim_clean_data <- clean_data %>%
  mutate(combined_id=as.numeric(substr(sample_id, 1, 3)))%>%
  mutate(gene=targetName)
  select(combined_id, timepoint_imm, infectiontype, gene, concentration)

rnaseq <- read.csv("~/Downloads/A_S_All_ExpressionMatrix.csv")

rnaseq_metadata <- read.csv("~/Downloads/new_rna-seq_metadata.csv")

combo_rnaseq <- rnaseq %>%
  mutate(gene=X)%>%
  select(-X)%>%
  filter(gene %in% nulisa_analytes)%>%
  pivot_longer(cols=colnames(.)[1:56], names_to = "Novogene.sample.name", values_to = "rna")%>%
  left_join(., rnaseq_metadata, by="Novogene.sample.name")%>%
  select(gene, Novogene.sample.name, rna, combined_id, timepoint_imm, infectiontype)%>%
  left_join(., slim_clean_data, by = c("combined_id", "timepoint_imm", "infectiontype", "gene"))%>%
  distinct(combined_id, timepoint_imm, infectiontype, gene, rna, concentration)


detectable_genes <- combo_rnaseq %>%
  group_by(gene)%>%
  summarise("n"=sum(rna==0))%>%
  filter(n<40)




combo_rnaseq_purf <- combo_rnaseq%>%
  filter(infectiontype %in% c("A", "S"))%>%
  filter(gene %in% detectable_genes$gene)%>%
  group_by(gene)%>%
  nest()%>%
  mutate(correlation=map(data, ~cor.test(.$rna, .$concentration, method = "spearman")))%>%
  mutate(p=map_dbl(correlation, ~.$p.value),
         rho=map_dbl(correlation, ~.$estimate))%>%# do(broom::tidy(cor.test(.$concentration, .$freq, method="spearman")))%>%
  ungroup()%>%
  mutate(padj=p.adjust(p))

sig_cor <- combo_rnaseq_purf%>%
  filter(p<0.05)


# significant correlations

rna_nulisa_cor <- combo_rnaseq %>%
  filter(gene %in% sig_cor$gene)%>%
  filter(timepoint_imm %in% c(-1, 0, 7, 14))%>%
  ggplot(., aes(x=rna, y=concentration))+
  geom_point(aes(color=factor(timepoint_imm, levels = c(-1, 0, 7, 14))))+
  geom_smooth(method="lm")+
  scale_x_log10()+
  ggpubr::stat_cor(method="spearman", color="darkred",r.digits = 2, p.digits = 2, size=2)+
  facet_wrap(~gene, scales = "free", ncol = 5)+
  scale_color_manual(values=viridis::magma(5)[c(1,4,3,2)])+
  theme_minimal()+
  theme(legend.title = element_blank())

ggsave("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/rna_nulisa_cor.png", rna_nulisa_cor, width = 12.5, height = 20, dpi=444, bg="white")

# insignificant correlations 
no_rna_nulisa_cor <- combo_rnaseq %>%
  filter(gene %notin% sig_cor$gene & gene %in% detectable_genes$gene)%>%
  filter(timepoint_imm %in% c(-1, 0, 7, 14))%>%
  ggplot(., aes(x=rna, y=concentration))+
  geom_point(aes(color=factor(timepoint_imm, levels = c(-1, 0, 7, 14))))+
  geom_smooth(method="lm")+
  scale_x_log10()+
  ggpubr::stat_cor(method="spearman", color="darkred",r.digits = 2, p.digits = 2, size=2)+
  facet_wrap(~gene, scales = "free", ncol = 10)+
  scale_color_manual(values=viridis::magma(5)[c(1,4,3,2)])+
  theme_minimal()+
  theme(legend.title = element_blank())

ggsave("~/postdoc/stanford/plasma_analytes/MUSICAL/combo/figures/no_rna_nulisa_cor.png", no_rna_nulisa_cor, width = 25, height = 20, dpi=444, bg="white")

 

# linear regression ####

combo_rnaseq_purfff <- combo_rnaseq%>%
  filter(infectiontype %in% c("A", "S"))%>%
  # filter(gene %in% detectable_genes$gene)%>%
  group_by(gene)%>%
  nest()%>%
  mutate(model=map(data, ~lm(concentration~rna+factor(combined_id), data=.)))%>%
  mutate(summary=map(model, ~summary(.)))
        
 