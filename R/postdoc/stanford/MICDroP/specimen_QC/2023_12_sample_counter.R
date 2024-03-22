library(dplyr)
library(tidyr)
library(xlsx)

`%notin%` <- Negate(`%in%`)
is.blank <- function(x){sapply(x, function(y) {ifelse(y=="", TRUE, FALSE)})}

raw_data <- haven::read_dta("~/postdoc/stanford/clinical_data/MICDROP/specimen_QC/2023_12/MICDSpecimenBoxDec23_formerge.dta")
bdays <- haven::read_dta("~/postdoc/stanford/clinical_data/MICDROP/specimen_QC/2023_10/MICDSpecimenBoxOct23_withclinical.dta", col_select = c(id, dob))
bdays <- bdays[!duplicated(bdays),]
raw_data$dob <- bdays$id[match(raw_data$id, bdays$id)]


sample_ages <- c(8, 24, 52, 68, 84, 104, 120)
sample_ages_minus <- sample_ages-1
sample_ages_plus <- sample_ages+1

sample_ranges <- sort(c(sample_ages, sample_ages_minus, sample_ages_plus))

long_specimen_data <- raw_data %>%
  mutate("flo_age_in_wks"=as.numeric(date-dob)%/%7)%>%
  select(id, dob, date, flo_age_in_wks, SampleDate, PBMC, Paxgene, Plasma, PlasmaPK, CellStabilizer, qPCR) %>%
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


# add a timepoint label that allows for a week on either side
routine_ish_visits <- long_specimen_data %>%
  filter(flo_age_in_wks %in% sample_ranges)




# no_visit_type_recorded <- routine_ish_visits %>%
#   filter(Specimen_ID != "" & visittype %notin% c(0, 1, 2))

n_distinct(no_visit_type_recorded$visit_id)
n_distinct(no_visit_type_recorded$id)
#766; 739


# subset routine visits with no sample code
visits_missing_samples <- routine_ish_visits %>%
  filter(Specimen_ID == "" & visittype==1)

n_distinct(visits_missing_samples$visit_id)
n_distinct(visits_missing_samples$id)
#806 visits; 521 children


# count the number of missing samples per visit, subset to only include visits where all are missing
visits_missing_all_samples <- visits_missing_samples %>%
  group_by(visit_id)%>%
  mutate("n_missing"=n())%>%
  filter(n_missing==4, !duplicated(visit_id))
# 25 visits miss all samples 


# how many individuals have attended each timepoint where samples were taken;
# sometimes samples are taken but visit-type is not recorded
# the 11 extra visits exist because the enrolment visit was like a day before the routine visit 
participant_summary <- routine_ish_visits %>%
  # mutate("visittype_edit"=ifelse(!is.na(visittype), visittype, ifelse(visit_id %in% no_visit_type_recorded, 1, NA)))%>%
  filter(visittype %notin% c(2, NA))%>%
  group_by(Timepoint_in_weeks) %>%
  summarise("Number_of_Individuals"=n_distinct(subject_id), "Number_of_Visits"=n_distinct(visit_id))



participants_with_8 <- long_specimen_data %>%
  filter(Timepoint_in_weeks==8)

n_distinct(participants_with_8$id)

participants_with_no_8 <- raw_data$id[raw_data$id %notin% participants_with_8$id]

#correct for withdrawal
withdrawn <- raw_data %>%
  filter(!is.na(raw_data$withdrawaldate))%>%
  mutate("withdrawal_age" = withdrawaldate-dob)%>%
  dplyr::select(id, withdrawal_age)%>%
  distinct()

withdrawn_before_8 <- withdrawn%>%
  filter(withdrawal_age<65)
# ten people, exactly how many 8 week visits are missing

visits_with_samples_taken <- routine_ish_visits %>%
  filter(Specimen_ID != "") %>%
  group_by(Timepoint_in_weeks) %>%
  summarise("Number_of_Individuals"=n_distinct(subject_id),
            "Number_of_Visits"=n_distinct(visit_id))


# routine_ish_visits %>%
#   filter(Timepoint_in_weeks %in% sample_ages)%>%
#   group_by(Specimen_Type, Timepoint_in_weeks) %>%
#   summarise("Number_of_Samples"=n_distinct(Specimen_ID_ID))

# number of samples of each type at each timepoint, including collection rate
sample_counts <- routine_ish_visits %>%
  filter(Timepoint_in_weeks %in% sample_ages, Specimen_ID != "")%>%
  group_by(Specimen_Type, Timepoint_in_weeks) %>%
  summarise("Number_of_Samples"=n_distinct(Specimen_ID))%>%
  ungroup()%>%
  mutate("Collection_Rate"=.$Number_of_Samples/participant_summary$Number_of_Individuals[match(.$Timepoint_in_weeks, participant_summary$Timepoint_in_weeks)],
         "Collection_Failure_Rate"=1-Collection_Rate,
         "Missed_Samples"=participant_summary$Number_of_Individuals[match(Timepoint_in_weeks, participant_summary$Timepoint_in_weeks)]-Number_of_Samples)


collection_rate_plot <-   sample_counts%>%
  filter(Specimen_Type %in% c("CellStabilizer", "Paxgene", "PBMC", "Plasma"))%>%
  ggplot(aes(x=factor(Timepoint_in_weeks), y=Collection_Rate, fill=Specimen_Type))+
  geom_bar(stat="identity")+
  facet_grid(~Specimen_Type)+
  scale_y_continuous(labels = scales::label_percent())+
  scale_fill_manual(values = colorspace::sequential_hcl(n=7, "RdPu")[1:6])+
  xlab("Timepoint (weeks)")+
  ylab("Collection Rate")+
  theme_minimal()+
  theme(legend.position = "none")

ggsave("~/postdoc/stanford/clinical_data/MICDROP/specimen_QC/2023_10/collection_rate_plot.png", collection_rate_plot, bg="white", dpi=444, width=8, height=4)


indie_samples <- routine_ish_visits %>%
  filter(Timepoint_in_weeks %in% sample_ages, Specimen_ID!="")%>%
  group_by(id, Specimen_Type, Timepoint_in_weeks) %>%
  summarise("Number_of_Samples"=n_distinct(Specimen_ID))

table(indie_samples$Number_of_Samples, indie_samples$Specimen_Type)

routine_ish_visits %>%
  filter(id==11101, Timepoint_in_weeks==24)%>%
  select(id, date, Specimen_ID, Specimen_Type)%>%
  arrange(desc(Specimen_ID), date)

routine_ish_visits %>%
  filter(id==10680, Timepoint_in_weeks==52)%>%
  select(id, date, Specimen_ID, Specimen_Type)%>%
  arrange(desc(Specimen_ID), date)





# subset visit database to be only "duplicate" visits and send IDRC


# how many individuals do we have n samples available
n_sample_summary <- routine_ish_visits %>%
  filter(Specimen_ID != "")%>%
  group_by(Specimen_Type, subject_id) %>%
  summarise("Number_of_Samples_of_Individual" = n()) %>%
  ungroup()%>%
  group_by(Specimen_Type,  Number_of_Samples_of_Individual) %>%
  summarise("Number_of_Individuals_with_n_Samples"=n())



n_sample_summary_plot  <- n_sample_summary %>%
  filter(Specimen_Type %in% c("CellStabilizer", "Paxgene", "PBMC", "Plasma"))%>%
  ggplot(aes(x=factor(Number_of_Samples_of_Individual), y=Number_of_Individuals_with_n_Samples, fill=Specimen_Type))+
  geom_bar(stat="identity")+
  geom_text(aes(label=Number_of_Individuals_with_n_Samples),vjust= -0.2, size=3)+
  facet_grid(~Specimen_Type, scales = "free_x")+
  scale_fill_manual(values = colorspace::sequential_hcl(n=7, "RdPu")[1:6])+
  xlab("Number of Samples of Individual")+
  ylab("Number of Individuals with n Samples")+
  theme_minimal()+
  theme(legend.position = "none",
        # axis.text.x = element_text(angle=90)
  )

ggsave("~/postdoc/stanford/clinical_data/MICDROP/specimen_QC/2023_10/n_sample_summary_plot.png", n_sample_summary_plot, bg="white", dpi=444, width=8, height=4)


# visits without samples


samples_missing_visits <- routine_ish_visits %>%
  mutate("visittype_edit"=ifelse(!is.na(visittype), visittype, ifelse(visit_id %in% no_visit_type_recorded, 1, NA)))%>%
  filter(Specimen_ID == "" & visittype==1, !duplicated(visit_id))


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

saveWorkbook(wb, "/Users/fbach/Box Sync/MIC_DroP IPTc Study/Florian_Specimen_QC/2023_10/Specimen_Summary.xlsx")



# handling of QC errors / conflicts ####

# who has three cell stabilizer samples and why
# the_threes <- long_specimen_data %>%
#   group_by(Specimen_Type, id) %>%
#   summarise("Number_of_Samples" = n()) %>%
#   filter(Number_of_Samples==3,  Specimen_Type %in% c("PBMC", "Paxgene", "Cellstabilizer"))

# visits with erroneous? samples taken
err_visits <- long_specimen_data %>%
  filter(id %in% the_threes$id, Specimen_Type %in% c("PBMC", "Paxgene", "Cellstabilizer"), flo_age_in_wks %notin% sample_ranges)

write.csv(err_visits, "/Users/fbach/Box Sync/MIC_DroP IPTc Study/Florian_Specimen_QC/2023_10/third_sample_visits.csv")

# broken birthdays

broken_bb_visits <-long_specimen_data %>%
  filter(is.na(dob))
write.csv(broken_bb_visits, "/Users/fbach/Box Sync/MIC_DroP IPTc Study/Florian_Specimen_QC/2023_10/visits_with_missing_birthdays.csv")

broken_bb_individuals <- distinct(broken_bb_visits,id)
write.csv(broken_bb_individuals, "/Users/fbach/Box Sync/MIC_DroP IPTc Study/Florian_Specimen_QC/2023_10/individuals_with_missing_birthdays.csv")



# turn into pdf table ####
# wordlcoud ####
# 
# 
# library(wordcloud)
# library(tm)
# 
# make_cloudable <- function(corpus){
#   
#   corpus = corpus %>%
#     tm_map(removeNumbers) %>%
#     tm_map(removePunctuation) %>%
#     tm_map(stripWhitespace)
#   corpus = tm_map(corpus, content_transformer(tolower))
#   corpus = tm_map(corpus, removeWords, stopwords("english"))
#   
#   matrix = as.matrix(TermDocumentMatrix(corpus))
#   words = sort(rowSums(matrix),decreasing=TRUE) 
#   df = data.frame(word = names(words),freq=words)
#   return(df)
#   
# }
# 
# 
# comments <- raw_data %>%
#   dplyr::select(c(Comment1, Comment2, Comment5))%>%
#   pivot_longer(cols=c(Comment1, Comment2, Comment5), values_to = "Comment", names_to = "Comment_Type")%>%
#   filter(Comment != "")
# 
# all_corpus <-  Corpus(VectorSource(comments$Comment))
# all_comments <- make_cloudable(all_corpus)
# 
# 
# wordcloud(words = all_comments$word, freq = all_comments$freq, min.freq = 1, scale=c(3, .5),
#           max.words=200, random.order=FALSE, rot.per=0.35,
#           colors=brewer.pal(8, "Dark2"))


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

# sandbox
# 
# withdrawals <- data.frame(matrix(ncol=4))
# colnames(withdrawals) <- c("id", "withdrawal_age", "withdrawal_before_timepoint", "remaining_participants")
# 
# for(i in sample_ages){
#   
#   participants_with_8 <- long_specimen_data %>%
#     filter(Timepoint_in_weeks==i) %>%
#     distinct(id)
#   
#   participants_with_no_8 <- raw_data$id[raw_data$id %notin% participants_with_8$id]
#   
#   #correct for withdrawal
#   withdrawn <- raw_data %>%
#     filter(!is.na(raw_data$withdrawaldate))%>%
#     mutate("withdrawal_age" = withdrawaldate-dob,
#            "withdrawal_before_timepoint"=i)%>%
#     dplyr::select(id, withdrawal_age, withdrawal_before_timepoint)%>%
#     distinct()
#   
#   remaining_participants <-  nrow(participants_with_8)
#   
#   withdrawn_before_8 <- withdrawn %>%
#     filter(withdrawal_age<(i)*7)
#   
#   withdrawn_before_8$remaining_participants <- remaining_participants
#   withdrawals <- rbind(withdrawals, withdrawn_before_8)
#   
# }
# 
# withdrawals <- withdrawals %>%
#   distinct(id, .keep_all = TRUE)
# 
# left <- withdrawals %>%
#   group_by(withdrawal_before_timepoint)%>%
#   summarise(n())
# 
# right <- withdrawals %>%
#   group_by(withdrawal_before_timepoint)%>%
#   distinct(remaining_participants)
# 
# combo <- full_join(left, right)
# participant_summary
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# long_specimen_data %>%
#   filter(visittype==0)%>%
#   distinct(id)%>%
#   nrow()
# #924
# 
# participants_with_8 <- long_specimen_data %>%
#   filter(Timepoint_in_weeks==8) %>%
#   distinct(id)
# 
# nrow(participants_with_8)
# #919
# 
# 
# participants_with_no_8 <- raw_data$id[raw_data$id %notin% participants_with_8$id]
# 
# enrolment <- raw_data%>%
#   filter(visittype==0)
# 
eight <- long_specimen_data%>%
  filter(Timepoint_in_weeks==8)
# 
# n_distinct(enrolment$id)
# table(is.na(enrolment$dob), useNA = "ifany")
# table(enrolment$gender, useNA = "ifany")
# table(is.na(enrolment$dob), useNA = "ifany")
# table(is.na(enrolment$AGE), useNA = "ifany")
# 
# long_specimen_data <- enrolment %>%
#   mutate("flo_age_in_wks"=as.numeric(date-dob)%/%7)%>%
#   select(id, dob, date, flo_age_in_wks, mstatus, qPCRparsdens, ageinwks, SampleDate, PBMC, Paxgene, Plasma, PlasmaPK, CellStabilizer, qPCR, visittype) %>%
#   mutate("visit_id"=paste(id, date, sep="_"))%>%
#   pivot_longer(cols = c(PBMC, Paxgene, Plasma, PlasmaPK, CellStabilizer, qPCR), names_to = "Specimen_Type", values_to = "Specimen_ID")%>%
#   mutate(subject_id=id)%>%
#   #Specimen_IDs are shared between specimen types, so let's create a unique code
#   mutate(Specimen_ID_ID=paste(Specimen_Type, visit_id, sep="_"))%>%
#   mutate("Timepoint_in_weeks"=if_else(
#     flo_age_in_wks %in% sample_ages, flo_age_in_wks, ifelse(
#       flo_age_in_wks %in% sample_ages_minus, flo_age_in_wks+1, if_else(
#         flo_age_in_wks %in% sample_ages_plus, flo_age_in_wks-1, NA)
#     )
#   )
#   )
# 
# participant_summary <- long_specimen_data %>%
#   mutate("visittype_edit"=ifelse(!is.na(visittype), visittype, ifelse(visit_id %in% no_visit_type_recorded, 1, NA)))%>%
#   filter(visittype_edit %notin% c(2, NA))%>%
#   summarise("Number_of_Individuals"=n_distinct(subject_id), "Number_of_Visits"=n_distinct(visit_id))
# 
# 
# 
# sample_counts <- long_specimen_data %>%
#   group_by(Specimen_Type, Timepoint_in_weeks) %>%
#   summarise("Number_of_Samples"=n_distinct(Specimen_ID_ID))%>%
#   ungroup()%>%
#   mutate("Collection_Rate"=.$Number_of_Samples/participant_summary$Number_of_Individuals[match(.$Timepoint_in_weeks, participant_summary$Timepoint_in_weeks)],
#          "Collection_Failure_Rate"=1-Collection_Rate,
#          "Missed_Samples"=participant_summary$Number_of_Individuals[match(Timepoint_in_weeks, participant_summary$Timepoint_in_weeks)]-Number_of_Samples)
# 

# [1] "id"               "dob"              "gender"           "enroldate"        "ageinwksenroll"  
# [6] "withdrawaldate"   "withdrawalreason" "enddate"          "duration"         "date"            
# [11] "visittype"        "ageinwks"         "AGE"              "temp"             "fever"           
# [16] "febrile"          "paxgene"          "pcellstab"        "immunol"          "stool"           
# [21] "qPCRparsdens"     "qPCRdich"         "qPCRcat"          "hb"               "wbc"             
# [26] "mstatus"          "incidentmalaria"  "BoxNumber1"       "PositionRow1"     "PositionColumn1" 
# [31] "PBMC"             "RandomNumber1"    "Volume1"          "CellCount1"       "Comment1"        
# [36] "BoxNumber2"       "PositionRow2"     "PositionColumn2"  "Paxgene"          "RandomNumber2"   
# [41] "Volume2"          "Comment2"         "BoxNumber3"       "PositionRow3"     "PositionColumn3" 
# [46] "Plasma"           "Volume3"          "Comment3"         "BoxNumber4"       "PositionRow4"    
# [51] "PositionColumn4"  "PlasmaPK"         "RandomNumber4"    "Volume4"          "BoxNumber5"      
# [56] "PositionRow5"     "PositionColumn5"  "CellStabilizer"   "RandomNumber5"    "Volume5"         
# [61] "Comment5"         "BoxNumber6"       "PositionRow6"     "PositionColumn6"  "RandomNumber6"   
# [66] "Volume6"          "BoxNumber7"       "PositionRow7"     "PositionColumn7"  "qPCR"            
# [71] "RandomNumber7"    "Volume7"          "SampleDate"  


# visittype ==1 ####
routine_visits <- long_specimen_data %>%
  filter(visittype==1)


participant_summary2 <- routine_visits %>%
  group_by(Timepoint_in_weeks) %>%
  summarise("Number_of_Individuals"=n_distinct(subject_id), "Number_of_Visits"=n_distinct(visit_id))



sample_counts2 <- routine_visits %>%
  filter(Timepoint_in_weeks %in% sample_ages, Specimen_ID != "")%>%
  group_by(Specimen_Type, Timepoint_in_weeks) %>%
  summarise("Number_of_Samples"=n_distinct(Specimen_ID))%>%
  ungroup()%>%
  mutate("Collection_Rate"=.$Number_of_Samples/participant_summary$Number_of_Individuals[match(.$Timepoint_in_weeks, participant_summary$Timepoint_in_weeks)],
         "Collection_Failure_Rate"=1-Collection_Rate,
         "Missed_Samples"=participant_summary$Number_of_Individuals[match(Timepoint_in_weeks, participant_summary$Timepoint_in_weeks)]-Number_of_Samples)


collection_rate_plot2 <- ggplot(sample_counts, aes(x=factor(Timepoint_in_weeks), y=Collection_Rate, fill=Specimen_Type))+
  geom_bar(stat="identity")+
  facet_grid(~Specimen_Type)+
  scale_y_continuous(labels = scales::label_percent())+
  scale_fill_manual(values = colorspace::sequential_hcl(n=7, "RdPu")[1:6])+
  xlab("Timepoint (weeks)")+
  ylab("Collection Rate")+
  theme_minimal()+
  theme(legend.position = "none")

ggsave("~/postdoc/stanford/clinical_data/MICDROP/specimen_QC/2023_10/collection_rate_plot2.png", collection_rate_plot2, bg="white", dpi=444, width=8, height=4)


n_sample_summary2 <- routine_visits %>%
  filter(Specimen_ID != "")%>%
  group_by(Specimen_Type, subject_id) %>%
  summarise("Number_of_Samples_of_Individual" = n()) %>%
  ungroup()%>%
  group_by(Specimen_Type,  Number_of_Samples_of_Individual) %>%
  summarise("Number_of_Individuals_with_n_Samples"=n())



n_sample_summary_plot2  <- ggplot(n_sample_summary2, aes(x=factor(Number_of_Samples_of_Individual), y=Number_of_Individuals_with_n_Samples, fill=Specimen_Type))+
  geom_bar(stat="identity")+
  geom_text(aes(label=Number_of_Individuals_with_n_Samples),vjust= -0.2, size=3)+
  facet_wrap(~Specimen_Type, scales = "free_x")+
  scale_fill_manual(values = colorspace::sequential_hcl(n=7, "RdPu")[1:6])+
  xlab("Number of Samples of Individual")+
  ylab("Number of Individuals with n Samples")+
  theme_minimal()+
  theme(legend.position = "none",
        axis.text.x = element_text(angle=90))

ggsave("~/postdoc/stanford/clinical_data/MICDROP/specimen_QC/2023_10/n_sample_summary_plot2.png", n_sample_summary_plot2, bg="white", dpi=444, width=8, height=6)



putative_eight <- routine_ish_visits %>%
  # mutate("visittype_edit"=ifelse(!is.na(visittype), visittype, ifelse(visit_id %in% no_visit_type_recorded, 1, NA)))%>%
  filter(visittype %notin% c(2, NA))%>%
  filter(Timepoint_in_weeks==8) %>%
  filter(!duplicated(visit_id))%>%
  select(subject_id, date, visit_id)


definitive_eight <- routine_visits %>%
  filter(Timepoint_in_weeks==8) %>%
  # filter(!duplicated(visit_id))%>%
  select(subject_id, date, visit_id)



uncoded_visits <- putative_eight %>%
  filter(visit_id %notin% definitive_eight$visit_id)
"uncoded_visits.csv"

unknown_ids <- raw_data %>%
  filter(id %notin% putative_eight$subject_id)%>%
  distinct(id)
"unknown_ids.csv"

unknown_visits <- long_specimen_data %>%
  filter(id %notin% putative_eight$subject_id)%>%
  filter(!duplicated(visit_id))
"unknown_id_visits.csv"

missing_birthdays <- long_specimen_data %>%
  filter(is.na(dob))%>%
  filter(!duplicated(visit_id))%>%
  select(id, date)
"missing_birthdays.csv"

