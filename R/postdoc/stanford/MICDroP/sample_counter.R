library(dplyr)
library(tidyr)
`%notin%` <- Negate(`%in%`)


raw_data <- haven::read_dta("~/postdoc/stanford/clinical_data/MICDROP/MICDropSpecimensJan2023_withClinical.dta")


sample_ages <- seq(8, 88, by=16)
sample_ages_minus <- sample_ages-1
sample_ages_plus <- sample_ages+1
# 
sample_ranges <- sort(c(sample_ages, sample_ages_minus, sample_ages_plus))


long_specimen_data <- raw_data %>%
  mutate(dob=ifelse(id==10883, as.Date("2022-05-29"), dob))%>%
  mutate("flo_age_in_wks"=as.numeric(date-dob)%/%7)%>%
  select(id, dob, date, flo_age_in_wks, ageinwks, SampleDate, SpecimenID1, SpecimenID2, SpecimenID3, SpecimenID4, SpecimenID5, SpecimenID6, visittype) %>%
  mutate("visit_id"=paste(id, date, sep="_"))%>%
  pivot_longer(cols = c(SpecimenID1, SpecimenID2, SpecimenID3, SpecimenID4, SpecimenID5, SpecimenID6), names_to = "Specimen_Type", values_to = "Specimen_ID")%>%
  filter(Specimen_ID!="")%>%
  mutate(subject_id=id)%>%
  mutate(Specimen_Type=recode(Specimen_Type,
                            "SpecimenID1"="PBMC",
                            "SpecimenID2"="Paxgene",
                            "SpecimenID3"="Plasma",
                            "SpecimenID4"="PlasmaPK",
                            "SpecimenID5"="CellStabiliser",
                            "SpecimenID6"="qPCR"))%>%
  mutate(Specimen_ID_ID=paste(Specimen_Type, Specimen_ID, sep="_"))


#check for individuals who came at more than one visit / have more than one routine
routine_ish_visits <- long_specimen_data %>%
  filter(flo_age_in_wks %in% sample_ranges) %>%
  mutate("timepoint"=if_else(
    flo_age_in_wks %in% sample_ages, flo_age_in_wks, ifelse(
      flo_age_in_wks %in% sample_ages_minus, flo_age_in_wks+1, if_else(
        flo_age_in_wks %in% sample_ages_plus, flo_age_in_wks-1, 999)
      )
  )
  )

routine_visits <- routine_ish_visits %>%
  filter(visittype==1) 


# how many individuals have attended each planned visit
participant_summary <- routine_ish_visits %>%
  group_by(timepoint) %>%
  summarise("Number_of_Individuals"=n_distinct(id))


# how many people had more than one visit around the time of a planned visit
routine_ish_visits %>% 
  group_by(id, timepoint)%>%
  summarise("number_of_visits"=n_distinct(visit_id))%>%
  group_by(number_of_visits)%>%
  summarise("how_often"=n())


# one individual had PBMCs taken at the enrolment visit at 8 weeks
# routine_ish_visits %>%
#   filter(timepoint %in% sample_ages)%>%
#   group_by(Specimen_Type, timepoint) %>%
#   summarise("Number_of_Samples"=n_distinct(Specimen_ID_ID))

# number of samples of each type at each timepoint, including collection rate
sample_summary <- routine_visits %>%
  filter(timepoint %in% sample_ages)%>%
  group_by(Specimen_Type, timepoint) %>%
  summarise("Number_of_Samples"=n_distinct(Specimen_ID_ID))%>%
  ungroup()%>%
  mutate("Collection_Rate"=.$Number_of_Samples/participant_summary$Number_of_Individuals[match(.$timepoint, participant_summary$timepoint)],
         "Collection_Failure_Rate"=1-Collection_Rate,
         "Missed_Samples"=participant_summary$Number_of_Individuals[match(timepoint, participant_summary$timepoint)]-Number_of_Samples)



# who has three cell stabiliser samples and why
# subset visit database to be only "duplicate" visits and send IDRC


# how many individuals do we have n samples available
n_sample_summary <- long_specimen_data %>%
  group_by(Specimen_Type, id) %>%
  summarise("Number_of_Samples" = n()) %>%
  group_by(Specimen_Type,  Number_of_Samples) %>%
  summarise("Number_of_Individuals_with_n_Samples"=n())


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

