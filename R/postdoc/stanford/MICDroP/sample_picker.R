library(tidyr)
library(dplyr)
library(xlsx)
library(ggplot2)

`%notin%` <- Negate(`%in%`)
is.blank <- function(x){sapply(x, function(y) {ifelse(y=="", TRUE, FALSE)})}

raw_data <- haven::read_dta("~/postdoc/stanford/clinical_data/MICDROP/specimen_QC/2023_12/MICDSpecimenBoxDec23_withclinical.dta")

# for finding putative sampling visits we'll filter the database to only include visits around the proper sampling timepoint dates, with a week plus/minus
sample_ages <- c(8, 24, 52, 68, 84, 104, 120)
sample_ages_minus <- sample_ages-1
sample_ages_plus <- sample_ages+1

sample_ranges <- sort(c(sample_ages, sample_ages_minus, sample_ages_plus))

#turn the data into long format, include a couple of convenience variables, drop irrelevant columns
long_specimen_data <- raw_data %>%
  mutate("flo_age_in_wks"=as.numeric(date-dob)%/%7)%>%
  select(id, dob, date, flo_age_in_wks, mstatus, qPCRparsdens, ageinwks, SampleDate,BoxNumber1,PositionColumn1, PositionRow1, PBMC, Paxgene, Plasma, PlasmaPK, CellStabilizer, qPCR, visittype, withdrawaldate) %>%
  mutate("visit_id"=paste(id, date, sep="_"))%>%
  pivot_longer(cols = c(PBMC, Paxgene, Plasma, PlasmaPK, CellStabilizer, qPCR), names_to = "Specimen_Type", values_to = "Specimen_ID")%>%
  mutate(subject_id=id)%>%
  #Specimen_IDs are shared between specimen types, so let's create a unique code
  mutate(Specimen_ID_ID=paste(Specimen_Type, visit_id, sep="_"))%>%
  mutate("Timepoint_in_weeks"=if_else(
    flo_age_in_wks %in% sample_ages, flo_age_in_wks, ifelse(
      flo_age_in_wks %in% sample_ages_minus, flo_age_in_wks+1, if_else(
        flo_age_in_wks %in% sample_ages_plus, flo_age_in_wks-1, 999)
    )
  )
  )

# 461 individuals, 646 samples
pbmc_samples <- long_specimen_data %>%
  filter(BoxNumber1 %in% paste("MICD-", 3001:3008, sep=""))
  # filter(Specimen_Type =="PBMC") #%>%
  # select(id, date, Specimen_Type, ageinwks, BoxNumber1)

# 49 dropouts
dropout <- long_specimen_data %>%
  filter(!is.na(withdrawaldate))%>%
  mutate(age_at_withdrawal=withdrawaldate-dob)%>%
  select(id, withdrawaldate, age_at_withdrawal)%>%
  filter(!duplicated(id))

# 31 individuals; 43 tubes; 27 8 weeks; 15 24 weeks
pbmc_dropouts <- pbmc_samples%>%
  filter(id %in% dropout$id)

old_pbmc_dropouts <- pbmc_samples%>%
  filter(id %in% dropout$id, flo_age_in_wks>10)

kids_with_clinical <- long_specimen_data %>%
  filter(id %in% pbmc_samples$id, mstatus!=0 & ageinwks<24) %>%
  select(id, date, Specimen_Type, ageinwks, BoxNumber1, mstatus, withdrawaldate)

ex_vivo <- long_specimen_data %>%
  # make variable "days_since_malaria" that calculates the days a malaria episode and resets to 0 when mstatus!=0
  mutate(days_since_malaria=SampleDate)%>%
  filter(Specimen_Type == "PBMC")%>%
  filter()%>%
  select(id, date, , ageinwks, BoxNumber1, mstatus, withdrawaldate)
  

# 
# ex_vivo_pbmc <- long_specimen_data %>%
#   filter(Specimen_Type=="PBMC")%>%
#   group_by(id)%>%
#   mutate(
#     n_infection = cumsum(mstatus%in%c(1,2,3)),
#     date_of_malaria = if_else(mstatus%in%c(1,2,3), date, NA))%>%
#   group_by(id, n_infection)%>%
#   mutate(days_since_malaria=date-date_of_malaria)
# 
# n_inf_summary <- ex_vivo_pbmc %>%
#   group_by(id)%>%
#   summarise("n_inf"=max(n_infection))


mic_drop %>%
  filter(id %in% unique(pbmc_dropouts$id),
         !is.na(any_parsdens),
         # ageinwks < 25
         )%>%
  ggplot(., aes(x=date, y=as.numeric(any_parsdens)+0.001))+
  geom_point(aes(color=factor(mstatus, #levels=c("0",
                              #       "1",
                              #      "2",
                              #     "3")
  )))+
  geom_line(alpha=0.3, aes(group=id))+
  # ggrepel::geom_text_repel(data=label_df, aes(x=AGE, y=as.numeric(any_parsdens), label=age_in_days))+
  facet_wrap(~ id)+
  ylab("qPCR parasites / μl\n")+
  xlab("Date")+
  scale_y_log10(breaks=c(1/100, 1, 10^2, 10^4, 10^6))+
  scale_x_date(breaks="2 week"
               # labels = seq(0,24,by=8)
  )+
  theme_minimal()+
  geom_vline(xintercept = 24)+
  # scale_shape_manual(values=c(16,15))+
  # scale_color_manual(values=mstatus_pal)+
  guides(color=guide_legend(title=""))+
  theme(axis.text.x = element_text(size=5, angle=90, vjust=0.5))




mic_drop %>%
  filter(id %in% unique(old_pbmc_dropouts$id),
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
  # ggrepel::geom_text_repel(data=label_df, aes(x=AGE, y=as.numeric(any_parsdens), label=age_in_days))+
  facet_wrap(~ id)+
  ylab("qPCR parasites / μl\n")+
  xlab("Date")+
  scale_y_log10(breaks=c(1/100, 1, 10^2, 10^4, 10^6))+
  theme_minimal()+
  geom_vline(xintercept = 24)+
  geom_vline(xintercept = 52, linetype="dashed")+
  
  # scale_shape_manual(values=c(16,15))+
  # scale_color_manual(values=mstatus_pal)+
  guides(color=guide_legend(title=""))+
  theme(axis.text.x = element_text(size=5, angle=90, vjust=0.5))


ggsave("~/postdoc/stanford/clinical_data/MICDROP/visit_databases/2023_07/figures/21_days_symptoms.png", two_week_symptoms_plot, width = 24, height=8, bg="white")
