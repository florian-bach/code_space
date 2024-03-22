library(ggplot2)
library(dplyr)

id_colors <- c("darkblue", "orange")
comp_pal <- c("asymptomatic"="lightgrey", "uncomplicated"="black", "complicated"="orange", "severe"="darkred")

data <- read.csv("~/postdoc/stanford/cytometry/attune/FB1/tables/cytokine_table-2.csv")
data <- data[,-10]
data <- data[c(-7:-10),]

colnames(data) <- c("Sample",
                    "CD4 IL2+",
                    "CD4 IFNg+",
                    "CD4 IL10+",
                    "CD4 TNF+",
                    "CD4 CD137+",
                    "CD8 CD137+",
                    "CD8 IFNg+",
                    "CD8 TNF+")

data$id <- as.character(substr(data$Sample, 1, 3))
data$condition <- ifelse(
  grepl("PHA", data$Sample), "PHA",
  ifelse(grepl("iRBC", data$Sample), "iRBC", "media"))

long_data <- data %>%
  pivot_longer(colnames(.)[2:9], names_to = "cytokine", values_to = "percentage")%>%
  mutate(percentage=percentage/100, 
         condition=factor(condition, levels=c("media", "iRBC", "PHA")))

quick_plot <- ggplot(long_data, aes(x=condition, y=percentage))+
  geom_point(aes(color=id))+
  facet_grid(~cytokine)+
  scale_y_log10(labels = scales::label_percent())+
  scale_color_manual()+
  theme_minimal()+
  theme(axis.text.x = element_text(angle=90, hjust=1),
        axis.title = element_blank())

ggsave("~/postdoc/stanford/cytometry/attune/FB1/figures/quick_cytokine_plot.png", quick_plot, width=8, height=2, bg="white")






# parasitaemia plots

mic_drop <- haven::read_dta("~/postdoc/stanford/clinical_data/MICDROP/visit_databases/2023_07/MICDROP all visit database through July 31st 2023.dta")

# merge parasitemia data so that qPCR takes precedent when both slide and qPCR are present
mic_drop <- mic_drop %>%
  mutate(mstatus = case_match(mstatus,
                              0~"no malaria",
                              1~"uncomplicated",
                              2~"complicated",
                              3~"quinine for AL failure",
                              4~"Q/AS failure"))%>%
  mutate(visit_id = paste(id, date, sep=""))%>%
  mutate("parasitaemia_method" = if_else(qPCRdich==1, "qPCR", if_else(BSdich==1, "smear", "dunno")))%>%
  mutate(any_parsdens = if_else(is.na(qPCRparsdens) & !is.na(pardens), pardens, qPCRparsdens))%>%
  mutate(parasitaemia_method = if_else(is.na(qPCRparsdens) & !is.na(pardens), "smear", parasitaemia_method))

parasitaemia_plot <- mic_drop %>%
  filter(id %in% c(10460, 10958),
         !is.na(any_parsdens),
         # ageinwks < 25
  )%>%
  ggplot(., aes(x=ageinwks, y=as.numeric(any_parsdens)+0.001))+
  geom_point(aes(color=factor(mstatus, #levels=c("0",
                              #       "1",
                              #      "2",
                              #     "3")
  )))+
  geom_line(alpha=0.3, aes(group=id))+
  facet_wrap(~ id)+
  ylab("qPCR parasites / μl\n")+
  xlab("age in weeks")+
  scale_y_log10(breaks=c(1/100, 1, 10^2, 10^4, 10^6))+
  theme_minimal()+
  geom_vline(xintercept = 24, linetype="dashed")+
  scale_color_manual(values=comp_pal)+
  # scale_shape_manual(values=c(16,15))+
  guides(color=guide_legend(title=""))+
  theme(axis.text.x = element_text(size=5, angle=90, vjust=0.5))

ggsave("~/postdoc/stanford/cytometry/attune/FB1/figures/parasitaemia_plot.png", parasitaemia_plot, width=8, height=4, bg="white")

