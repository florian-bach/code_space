library(dplyr)
library(tidyr)
library(xlsx)

`%notin%` <- Negate(`%in%`)
is.blank <- function(x){sapply(x, function(y) {ifelse(y=="", TRUE, FALSE)})}

raw_data <- haven::read_dta("~/postdoc/stanford/clinical_data/MICDROP/specimen_QC/2023_09/MICDSpecimenBoxSep23_withclinical.dta")


sample_ages <- c(8, 24, 52, 68, 84, 104, 120)
sample_ages_minus <- sample_ages-1
sample_ages_plus <- sample_ages+1
# 
sample_ranges <- sort(c(sample_ages, sample_ages_minus, sample_ages_plus))


long_specimen_data <- raw_data %>%
  mutate("flo_age_in_wks"=as.numeric(date-dob)%/%7)%>%
  select(id, dob, date, flo_age_in_wks, mstatus, qPCRparsdens, ageinwks, SampleDate, PBMC, Paxgene, Plasma, PlasmaPK, CellStabilizer, qPCR, visittype) %>%
  mutate("visit_id"=paste(id, date, sep="_"))%>%
  pivot_longer(cols = c(PBMC, Paxgene, Plasma, PlasmaPK, CellStabilizer, qPCR), names_to = "Specimen_Type", values_to = "Specimen_ID")%>%
  mutate(subject_id=id)%>%
  #Specimen_IDs are shared between specimen types, so let's create a unique code
  mutate(Specimen_ID_ID=paste(Specimen_Type, visit_id, sep="_"))


# add a timepoint label that allows for a week on either side
routine_ish_visits <- long_specimen_data %>%
  filter(flo_age_in_wks %in% sample_ranges) %>%
  mutate("Timepoint_in_weeks"=if_else(
    flo_age_in_wks %in% sample_ages, flo_age_in_wks, ifelse(
      flo_age_in_wks %in% sample_ages_minus, flo_age_in_wks+1, if_else(
        flo_age_in_wks %in% sample_ages_plus, flo_age_in_wks-1, 999)
    )
  )
  )



no_visit_type_recorded <- routine_ish_visits %>%
  filter(Specimen_ID != "" & visittype %notin% c(1,0) )

# subset routine visits with no sample code
samples_missing_visits <- routine_ish_visits %>%
  filter(Specimen_ID == "" & visittype==1)

# count the number of missing samples per visit, subset to only include visits where all are missing
all_samples_missing_visits <- samples_missing_visits %>%
  group_by(visit_id)%>%
  summarise("n_missing"=n(), Timepoint_in_weeks)%>%
  filter(n_missing==4, !duplicated(visit_id))

# identify visits where no samples where taken
blank_visits <- table(all_samples_missing_visits$Timepoint_in_weeks)


# how many individuals have attended each timepoint where samples were taken;
# sometimes samples are taken but visit-type is not recorded
# the 11 extra visits exist because the enrolment visit was like a day before the routine visit 
participant_summary <- routine_ish_visits %>%
  mutate("visittype_edit"=ifelse(!is.na(visittype), visittype, ifelse(visit_id %in% no_visit_type_recorded, 1, NA)))%>%
  filter(visittype_edit %notin% c(2, NA))%>%
  group_by(Timepoint_in_weeks) %>%
  summarise("Number_of_Individuals"=n_distinct(subject_id), "Number_of_Visits"=n_distinct(visit_id))



visits_with_samples_taken <- routine_ish_visits %>%
  filter(Specimen_ID != "") %>%
  group_by(Timepoint_in_weeks) %>%
  summarise("Number_of_Individuals"=n_distinct(subject_id),
            "Number_of_Visits"=n_distinct(visit_id))


# how many people had more than one visit around the time of a planned visit
routine_ish_visits %>% 
  group_by(id, Timepoint_in_weeks)%>%
  summarise("number_of_visits"=n_distinct(visit_id))%>%
  group_by(number_of_visits)%>%
  summarise("how_often"=n())


# one individual had PBMCs taken at the enrolment visit at 8 weeks
# routine_ish_visits %>%
#   filter(Timepoint_in_weeks %in% sample_ages)%>%
#   group_by(Specimen_Type, Timepoint_in_weeks) %>%
#   summarise("Number_of_Samples"=n_distinct(Specimen_ID_ID))

# number of samples of each type at each timepoint, including collection rate
sample_counts <- routine_visits %>%
  filter(Timepoint_in_weeks %in% sample_ages)%>%
  group_by(Specimen_Type, Timepoint_in_weeks) %>%
  summarise("Number_of_Samples"=n_distinct(Specimen_ID_ID))%>%
  ungroup()%>%
  mutate("Collection_Rate"=.$Number_of_Samples/participant_summary$Number_of_Individuals[match(.$Timepoint_in_weeks, participant_summary$Timepoint_in_weeks)],
         "Collection_Failure_Rate"=1-Collection_Rate,
         "Missed_Samples"=participant_summary$Number_of_Individuals[match(Timepoint_in_weeks, participant_summary$Timepoint_in_weeks)]-Number_of_Samples)


collection_rate_plot <- ggplot(sample_counts, aes(x=factor(Timepoint_in_weeks), y=Collection_Rate, fill=Specimen_Type))+
  geom_bar(stat="identity")+
  facet_wrap(~Specimen_Type, ncol=4)+
  scale_y_continuous(labels = scales::label_percent())+
  scale_fill_manual(values = colorspace::sequential_hcl(n=5, "RdPu")[1:4])+
  xlab("Timepoint (weeks)")+
  ylab("Collection Rate")+
  theme_minimal()+
  theme(legend.position = "none")

ggsave("~/postdoc/stanford/clinical_data/MICDROP/specimen_QC/2023_09/collection_rate_plot.png", collection_rate_plot, bg="white", dpi=444, width=8, height=4)


# subset visit database to be only "duplicate" visits and send IDRC


# how many individuals do we have n samples available
n_sample_summary <- routine_ish_visits %>%
  filter(Specimen_ID != "")%>%
  group_by(Specimen_Type, subject_id) %>%
  summarise("Number_of_Samples_of_Individual" = n()) %>%
  ungroup()%>%
  group_by(Specimen_Type,  Number_of_Samples_of_Individual) %>%
  summarise("Number_of_Individuals_with_n_Samples"=n())

n_sample_summary_plot  <- ggplot(n_sample_summary, aes(x=factor(Number_of_Samples_of_Individual), y=Number_of_Individuals_with_n_Samples, fill=Specimen_Type))+
  geom_bar(stat="identity")+
  geom_text(aes(label=Number_of_Individuals_with_n_Samples),vjust= -0.2, size=3)+
  facet_wrap(~Specimen_Type, ncol=4)+
  scale_fill_manual(values = colorspace::sequential_hcl(n=5, "RdPu")[1:4])+
  xlab("Number of Samples of Individual")+
  ylab("Number of Individuals with n Samples")+
  theme_minimal()+
  theme(legend.position = "none")

ggsave("~/postdoc/stanford/clinical_data/MICDROP/specimen_QC/2023_09/n_sample_summary_plot.png", n_sample_summary_plot, bg="white", dpi=444, width=8, height=4)


# visits without samples


samples_missing_visits <- routine_ish_visits %>%
  mutate("visittype_edit"=ifelse(!is.na(visittype), visittype, ifelse(visit_id %in% no_visit_type_recorded, 1, NA)))%>%
  filter(Specimen_ID == "" & visittype==1)

all_samples_missing_visits <- samples_missing_visits %>%
  group_by(visit_id)%>%
  summarise(id, "n_missing"=n(), Timepoint_in_weeks)%>%
  filter(n_missing==4, !duplicated(visit_id))

blank_visits <- table(all_samples_missing_visits$Timepoint_in_weeks)

# save summary-level info as exel sheet ####
wb = createWorkbook()

n_sample_summary_sheet = createSheet(wb, "N Sample Summary")
sample_summary_sheet = createSheet(wb, "Collection Summary")
participant_summary_sheet = createSheet(wb, "Participant Summary")

addDataFrame(as.data.frame(n_sample_summary), sheet=n_sample_summary_sheet, startColumn=1, row.names=FALSE)
addDataFrame(as.data.frame(sample_summary), sheet=sample_summary_sheet, startColumn=1, row.names=FALSE)
addDataFrame(as.data.frame(participant_summary), sheet=participant_summary_sheet, startColumn=1, row.names=FALSE)

saveWorkbook(wb, "/Users/fbach/Box Sync/MIC_DroP IPTc Study/Florian_Specimen_QC/2023_05/Specimen_Summary.xlsx")



# handling of QC errors / conflicts ####

# who has three cell stabilizer samples and why
the_threes <- long_specimen_data %>%
  group_by(Specimen_Type, id) %>%
  summarise("Number_of_Samples" = n()) %>%
  filter(Number_of_Samples==3,  Specimen_Type %in% c("PBMC", "Paxgene", "Cellstabilizer"))

# visits with erroneous? samples taken
err_visits <- long_specimen_data %>%
  filter(id %in% the_threes$id, Specimen_Type %in% c("PBMC", "Paxgene", "Cellstabilizer"), flo_age_in_wks %notin% sample_ranges)

write.csv(err_visits, "/Users/fbach/Box Sync/MIC_DroP IPTc Study/Florian_Specimen_QC/2023_05/third_sample_visits.csv")

# broken birthdays

broken_bb_visits <-long_specimen_data %>%
  filter(is.na(dob))
write.csv(broken_bb_visits, "/Users/fbach/Box Sync/MIC_DroP IPTc Study/Florian_Specimen_QC/2023_05/visits_with_missing_birthdays.csv")


broken_bb_individuals <- distinct(broken_bb_visits,id)
write.csv(broken_bb_individuals, "/Users/fbach/Box Sync/MIC_DroP IPTc Study/Florian_Specimen_QC/2023_05/individuals_with_missing_birthdays.csv")
# turn into pdf table ####



# wordlcoud ####


library(wordcloud)
library(tm)

make_cloudable <- function(corpus){
  
  corpus = corpus %>%
    tm_map(removeNumbers) %>%
    tm_map(removePunctuation) %>%
    tm_map(stripWhitespace)
  corpus = tm_map(corpus, content_transformer(tolower))
  corpus = tm_map(corpus, removeWords, stopwords("english"))
  
  matrix = as.matrix(TermDocumentMatrix(corpus))
  words = sort(rowSums(matrix),decreasing=TRUE) 
  df = data.frame(word = names(words),freq=words)
  return(df)
  
}


comments <- raw_data %>%
  dplyr::select(c(Comment1, Comment2, Comment5))%>%
  pivot_longer(cols=c(Comment1, Comment2, Comment5), values_to = "Comment", names_to = "Comment_Type")%>%
  filter(Comment != "")

all_corpus <-  Corpus(VectorSource(comments$Comment))
all_comments <- make_cloudable(all_corpus)


wordcloud(words = all_comments$word, freq = all_comments$freq, min.freq = 1, scale=c(3, .5),
          max.words=200, random.order=FALSE, rot.per=0.35,
          colors=brewer.pal(8, "Dark2"))


# wordlcoud split by comment type ####
# 
# 
# c1_corpus <- Corpus(VectorSource(
#   subset(raw_data, Comment1!= "", select = Comment1)
# ))
# 
# c2_corpus <- Corpus(VectorSource(
#   subset(raw_data, Comment2!= "", select = Comment2)
# ))
# 
# c5_corpus <- Corpus(VectorSource(
#   subset(raw_data, Comment5!= "", select = Comment5)
# ))
# 
# c1w <- make_cloudable(c1_corpus)
# c2w <- make_cloudable(c2_corpus)
# c5w <- make_cloudable(c5_corpus)
# 
# 
# 
# set.seed(1234) # for reproducibility
# 
# wordcloud(words = c1w$word, freq = c1w$freq, min.freq = 1, scale=c(3, .5),
#           max.words=200, random.order=FALSE, rot.per=0.35,
#           colors=brewer.pal(8, "Dark2"))
# 
# wordcloud(words = c2w$word, freq = c2w$freq, min.freq = 1, scale=c(3, .5),
#           max.words=200, random.order=FALSE, rot.per=0.35,
#           colors=brewer.pal(8, "Dark2"))
# 
# wordcloud(words = c5w$word, freq = c5w$freq, min.freq = 1, scale=c(3, .5),
#           max.words=200, random.order=FALSE, rot.per=0.35,
#           colors=brewer.pal(8, "Dark2"))

