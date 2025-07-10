old_cluster_counts <- read.csv("~/Downloads/final_cluster_counts.csv")
new_cluster_counts <- read.csv("~/postdoc/edinburgh/scRNAseq/scg/big_data_cluster_counts.csv")
old_activated_counts <- read.csv("~/Downloads/final_activated_cluster_counts.csv")


filter(old_cluster_counts, Cluster_ID==12)
filter(new_cluster_counts, Cluster_ID==10)
#filter(new_cluster_counts, Cluster_ID==9)
filter(old_activated_counts, Cluster_ID==9)

old_purr <- old_cluster_counts %>%
  group_by(Cluster_ID, N_Infection)%>%
  nest()%>%
  mutate(model=map(data, ~glm.nb(Count~Timepoint+Volunteer, data=.))) %>%
  mutate(summary=map(model, ~summary(.))) %>%
  mutate(raw_p=map_dbl(model, ~summary(.)$coefficients[14]))%>%
  ungroup()%>%
  group_by(N_Infection)%>%
  mutate(padj=p.adjust(raw_p))%>%
  filter(padj<0.1)


new_purr <- new_cluster_counts %>%
  group_by(Cluster_ID, N_Infection)%>%
  nest()%>%
  mutate(model=map(data, ~glm.nb(Count~Timepoint+Volunteer, data=.))) %>%
  mutate(summary=map(model, ~summary(.))) %>%
  mutate(raw_p=map_dbl(model, ~summary(.)$coefficients[14]))%>%
  ungroup()%>%
  group_by(N_Infection)%>%
  mutate(padj=p.adjust(raw_p))%>%
  filter(padj<0.1)

activated_old_purr <- old_activated_counts %>%
  group_by(Cluster_ID, N_Infection)%>%
  nest()%>%
  mutate(model=map(data, ~glm.nb(Count~Timepoint+Volunteer, data=.))) %>%
  #mutate(summary=map(model, ~summary.glm(.))) %>%
  mutate(summary=map(model, ~summary(.))) %>%
  mutate(raw_p=map_dbl(model, ~summary(.)$coefficients[14]))%>%
  ungroup()%>%
  group_by(N_Infection)%>%
  mutate(padj=p.adjust(raw_p))

top30 <- read.csv("~/postdoc/edinburgh/scRNAseq/scg/top30_signature_markers.csv")
all_genes <- read.csv("~/postdoc/edinburgh/scRNAseq/scg/signature_markers.csv")


top30%>%
  filter(cluster==10)


old_cluster_counts%>%
    group_by(Volunteer, Timepoint)%>%
    summarise("all_cells"=sum(Count))%>%
    ungroup()%>%
    mutate("frac"=all_cells/sum(all_cells))

old_activated_counts%>%
  group_by(Volunteer, Timepoint)%>%
  summarise("all_cells"=sum(Count))%>%
  ungroup()%>%
  mutate("frac"=all_cells/sum(all_cells))
