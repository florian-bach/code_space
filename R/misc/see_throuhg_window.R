
table(combined_data$Stim, combined_data$Stim2, useNA = "ifany")







# test one func file in complete isolation
test_file <- readxl::read_excel(filenames[1], col_types = "text") |>
  filter(!is.na(ID))

# check Stim before pivoting
table(test_file$Stim, useNA = "ifany")

# then pivot and check alignment
test_long <- test_file |>
  pivot_longer(starts_with("Time subset"), names_to = "gate", values_to = "freq_or_count") |>
  mutate(Stim2 = ifelse("Stim" %in% colnames(test_file), Stim, "phenotypic"))

table(test_long$Stim, test_long$Stim2, useNA = "ifany")



