n_studies <- 32
n_coders <- 4

# Strategy: rotate coders so each study gets 2 coders,
# and every coder ends with exactly 16 assignments.

assignments <- data.frame(
  study = 1:n_studies,
  coder1 = ((1:n_studies - 1) %% n_coders) + 1,
  coder2 = ((1:n_studies - 1 + 1) %% n_coders) + 1
)

# Long form
library(tidyr)
assignments_long <- pivot_longer(assignments,
                                 cols = c(coder1, coder2),
                                 values_to = "coder",
                                 names_to = NULL)

# Check: each coder should have exactly 16
print(table(assignments_long$coder))  

# Now build the wide "coder × 16 studies" table
library(dplyr)

coder_schedule <- assignments_long %>%
  arrange(coder, study) %>%
  group_by(coder) %>%
  summarise(studies = list(sort(unique(study))), .groups = "drop") %>%
  tidyr::unnest_wider(studies, names_sep = "_")

names(coder_schedule) <- c("coder", paste0("study_", 1:(ncol(coder_schedule)-1)))

coder_schedule

write.csv(coder_schedule, "study_allocation.csv")
