`%notin%` <- Negate(`%in%`)


ltc_2018 <-read.csv("~/Downloads/ltc_2018_pos.csv")
ltc_2019 <- read.csv("~/Downloads/ltc_2019_pos.csv")
nhr_2018 <- read.csv("~/Downloads/nhr_2018_pos.csv")
nhr_2019 <- read.csv("~/Downloads/nhr_2019_pos.csv")
hcup_2018 <- readxl::read_excel("~/Downloads/combined_hcup_2018.xlsx")
hcup_2019 <- readxl::read_excel("~/Downloads/combined_hcup_2019.xlsx")


ltc_2018 <- ltc_2018%>%
  mutate(year=2018)

ltc_2019 <- ltc_2019%>%
  mutate(year=2019)

combo_ltc <- bind_rows(ltc_2018, ltc_2019)

wide_combo_ltc <- combo_ltc%>%
  group_by(State, County.Name, year)%>%
  reframe("number_of_facilities"=n())%>%
  pivot_wider(names_from = year, names_prefix = "year_", values_from = number_of_facilities, id_cols=c(County.Name, State))%>%
  mutate(year_change=year_2019-year_2018)




nhr_2018 <- nhr_2018%>%
  mutate(year=2018)%>%
  select(PROVNUM, STATE, COUNTY_NAME, OWNERSHIP, BEDCERT, SPRINKLER_STATUS, OVERALL_RATING, year)

nhr_2019 <- nhr_2019%>%
  mutate(year=2019)%>%
  mutate(COUNTY_NAME=County_name)
  select(PROVNUM, STATE, COUNTY_NAME, OWNERSHIP, BEDCERT, SPRINKLER_STATUS, OVERALL_RATING, year)

combo_nhr <- bind_rows(nhr_2018, nhr_2019)


wide_combo_nhr <- combo_nhr%>%
  group_by(STATE, COUNTY_NAME, year)%>%
  reframe("number_of_facilities"=n())%>%
  pivot_wider(names_from = year, names_prefix = "year_", values_from = number_of_facilities, id_cols=c(COUNTY_NAME, STATE))%>%
  mutate(year_change=year_2019-year_2018)
  




hcup_2018 <- hcup_2018%>%
  mutate(year=2018)

hcup_2019 <- hcup_2019%>%
  mutate(year=2019)

combo_hcup <- bind_rows(hcup_2018, hcup_2019)

wide_combo_hcup <- combo_hcup%>%
  mutate(
    state_short = case_when(
      state == "Alabama"               ~ "AL",
      state == "Alaska"                ~ "AK",
      state == "Arizona"               ~ "AZ",
      state == "Arkansas"              ~ "AR",
      state == "California"            ~ "CA",
      state == "Colorado"              ~ "CO",
      state == "Connecticut"           ~ "CT",
      state == "Delaware"              ~ "DE",
      state == "District of Columbia"  ~ "DC",
      state == "Florida"               ~ "FL",
      state == "Georgia"               ~ "GA",
      state == "Hawaii"                ~ "HI",
      state == "Idaho"                 ~ "ID",
      state == "Illinois"              ~ "IL",
      state == "Indiana"               ~ "IN",
      state == "Iowa"                  ~ "IA",
      state == "Kansas"                ~ "KS",
      state == "Kentucky"              ~ "KY",
      state == "Louisiana"             ~ "LA",
      state == "Maine"                 ~ "ME",
      state == "Maryland"              ~ "MD",
      state == "Massachusetts"         ~ "MA",
      state == "Michigan"              ~ "MI",
      state == "Minnesota"             ~ "MN",
      state == "Mississippi"           ~ "MS",
      state == "Missouri"              ~ "MO",
      state == "Montana"               ~ "MT",
      state == "Nebraska"              ~ "NE",
      state == "Nevada"                ~ "NV",
      state == "New Hampshire"         ~ "NH",
      state == "Newjersey"            ~ "NJ",
      state == "Newmexico"            ~ "NM",
      state == "New York"              ~ "NY",
      state == "Northcarolina"        ~ "NC",
      state == "Northdakota"          ~ "ND",
      state == "Ohio"                  ~ "OH",
      state == "Oklahoma"              ~ "OK",
      state == "Oregon"                ~ "OR",
      state == "Pennsylvania"          ~ "PA",
      state == "Rhode Island"          ~ "RI",
      state == "Southcarolina"        ~ "SC",
      state == "Southdakota"          ~ "SD",
      state == "Tennessee"             ~ "TN",
      state == "Texas"                 ~ "TX",
      state == "Utah"                  ~ "UT",
      state == "Vermont"               ~ "VT",
      state == "Virginia"              ~ "VA",
      state == "Washington"            ~ "WA",
      state == "Westvirginia"         ~ "WV",
      state == "Wisconsin"             ~ "WI",
      state == "Wyoming"               ~ "WY"))%>%
  mutate(`Average Length of Stay (in days)`=if_else(`Average Length of Stay (in days)`=="*", NA, `Average Length of Stay (in days)`))%>%
  mutate(`Rate of Discharges per 100,000 Population`=if_else(`Rate of Discharges per 100,000 Population`=="*", NA, `Rate of Discharges per 100,000 Population`))%>%
  mutate(`Average Length of Stay (in days)`=as.integer(`Average Length of Stay (in days)`))%>%
  mutate(`Rate of Discharges per 100,000 Population`=as.integer(`Rate of Discharges per 100,000 Population`))%>%
  
  pivot_wider(names_from = year, names_prefix = "year_",
              values_from = c(`Average Length of Stay (in days)`, `Rate of Discharges per 100,000 Population`), id_cols=c(state_short, County))%>%
  
  mutate("year_change_discharge"=`Rate of Discharges per 100,000 Population_year_2019`-`Rate of Discharges per 100,000 Population_year_2018`)%>%
  mutate("year_change_length_of_stay"=`Average Length of Stay (in days)_year_2019`-`Average Length of Stay (in days)_year_2018`)
 

# wide_combo_hcup$County[wide_combo_hcup$County%notin%wide_combo_nhr$COUNTY_NAME]
# 
# gsub(" Parish", "", wide_combo_hcup$County[wide_combo_hcup$County%notin%wide_combo_nhr$COUNTY_NAME])
#   
# mutate(County=gsub(" Parish", "", county))
# 

# filter(county %notin% c("State Total", "US Total"))

# inner_join()