library(dplyr)
library(tidyr)

raw_data <- haven::read_dta("~/postdoc/stanford/clinical_data/MICDROP/MICDropSpecimensJan2023_withClinical.dta")

`%notin%` <- Negate(`%in%`)

sample_ages <- seq(8, 88, by=16)
sample_ages_minus <- sample_ages-1
sample_ages_plus <- sample_ages+1
# 
sample_ranges <- sort(c(sample_ages, sample_ages_minus, sample_ages_plus))


long_specimen_data <- raw_data %>%
  select(id, date, ageinwks, SampleDate, SpecimenID1, SpecimenID2, SpecimenID3, SpecimenID4, SpecimenID5, SpecimenID6, visittype) %>%
  mutate("visit_id"=paste(id, date))%>%
  pivot_longer(cols = c(SpecimenID1, SpecimenID2, SpecimenID3, SpecimenID4, SpecimenID5, SpecimenID6), names_to = "Specimen_Type", values_to = "Specimen_ID")%>%
  filter(Specimen_ID!="")%>%
  mutate(subject_id=id)%>%
  mutate(Specimen_Type=recode(Specimen_Type,
                            "SpecimenID1"="PBMC",
                            "SpecimenID2"="Paxgene",
                            "SpecimenID3"="Plasma",
                            "SpecimenID4"="PlasmaPK",
                            "SpecimenID5"="CellStabiliser",
                            "SpecimenID6"="qPCR"))


#check for individuals who came at more than one visit / have more than one routine
routine_visits <- long_specimen_data %>%
  #filter(ageinwks %in% sample_ranges & visittype==1) %>%
  filter(ageinwks %in% sample_ranges) %>%
  mutate("timepoint"=if_else(
    ageinwks %in% sample_ages, ageinwks, ifelse(
      ageinwks %in% sample_ages_minus, ageinwks+1, if_else(
        ageinwks %in% sample_ages_plus, ageinwks-1, 999)
      )
  )
  )



participant_summary <- routine_visits %>%
  group_by(timepoint, visittype) %>%
  summarise("Number_of_Individuals"=n_distinct(id))



routine_visits %>% 
  group_by(id, timepoint)%>%
  summarise("number_of_visits"=n_distinct(visit_id))%>%
  # ungroup()%>%
  # add_row(id=99999, timepoint=8, number_of_visits=4)%>%
  group_by(number_of_visits)%>%
  summarise("how_often"=n())

# one individual had PBMCs taken at the enrolment visit at 8 weeks
sample_summary <- routine_visits %>%
  filter(timepoint %in% sample_ages)%>%
  group_by(Specimen_Type, timepoint) %>%
  summarise("Number_of_Samples"=n_distinct(Specimen_ID))

  
sample_summary$Collection_Rate <- sample_summary$Number_of_Samples/participant_summary$Number_of_Individuals[match(sample_summary$timepoint, participant_summary$timepoint)]


sample_summary <- sample_summary %>%
  mutate("Collection_Failure_Rate"=1-Collection_Rate,
         "Missed_Samples"=participant_summary$Number_of_Individuals[match(timepoint, participant_summary$timepoint)]-Number_of_Samples)




n_sample_summary <- long_specimen_data %>%
  group_by(Specimen_Type, id) %>%
  summarise("Number_of_Samples" = n()) %>%
  group_by(Specimen_Type,  Number_of_Samples) %>%
  summarise("Number_of_Individuals_with_n_Samples"=n())


# make table ####



# wordlcoud ####


library(wordcloud)
library(tm)


c1_corpus <- Corpus(VectorSource(
  subset(raw_data, Comment1!= "", select = Comment1)
  ))

c2_corpus <- Corpus(VectorSource(
  subset(raw_data, Comment2!= "", select = Comment2)
))

c5_corpus <- Corpus(VectorSource(
  subset(raw_data, Comment5!= "", select = Comment5)
))


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


c1w <- make_cloudable(c1_corpus)
c2w <- make_cloudable(c2_corpus)
c5w <- make_cloudable(c5_corpus)



set.seed(1234) # for reproducibility
wordcloud(words = c5w$word, freq = df$freq, min.freq = 1, scale=c(3, .5),
          max.words=200, random.order=FALSE, rot.per=0.35,
          colors=brewer.pal(8, "Dark2"))

wordcloud(words = c1w$word, freq = df$freq, min.freq = 1, scale=c(3, .5),
          max.words=200, random.order=FALSE, rot.per=0.35,
          colors=brewer.pal(8, "Dark2"))
