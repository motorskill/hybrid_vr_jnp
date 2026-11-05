##### scr to run components for hybrid vr study #####
# created by cas 28 jul 2026 #
# modified by cas 17 oct 2026
# runs anovas 

library(openxlsx)
library(tidyverse)
library(afex)

base <- "/Users/motorskill/Desktop/R01s/vr/rep_vr/data_analyses/"

# Read data
cop_vr <- read.xlsx(file.path(base, "cop_vr.xlsx"))
cop_irl <- read.xlsx(file.path(base, "cop_irl.xlsx"))

lf_vr <- read.xlsx(file.path(base, "lf_vr.xlsx"))
lf_irl <- read.xlsx(file.path(base, "lf_irl.xlsx"))

gf_vr <- read.xlsx(file.path(base, "gf_vr.xlsx"))
gf_irl <- read.xlsx(file.path(base, "gf_irl.xlsx"))

rcop_vr <- read.xlsx(file.path(base, "rcop_vr.xlsx"))
rcop_irl <- read.xlsx(file.path(base, "rcop_irl2.xlsx"))

# Remove empty columns
remove_empty_columns <- function(data) {
  data[, colSums(!is.na(data)) > 0, drop = FALSE]
}

cop_vr <- remove_empty_columns(cop_vr)
cop_irl <- remove_empty_columns(cop_irl)

lf_vr <- remove_empty_columns(lf_vr)
lf_irl <- remove_empty_columns(lf_irl)

gf_vr <- remove_empty_columns(gf_vr)
gf_irl <- remove_empty_columns(gf_irl)

rcop_vr <- remove_empty_columns(rcop_vr)
rcop_irl <- remove_empty_columns(rcop_irl)

# Convert all data columns to numeric
convert_to_numeric <- function(data) {
  data %>%
    mutate(across(everything(), ~ as.numeric(as.character(.x))))
}

cop_vr <- convert_to_numeric(cop_vr)
cop_irl <- convert_to_numeric(cop_irl)

lf_vr <- convert_to_numeric(lf_vr)
lf_irl <- convert_to_numeric(lf_irl)

gf_vr <- convert_to_numeric(gf_vr)
gf_irl <- convert_to_numeric(gf_irl)

rcop_vr <- convert_to_numeric(rcop_vr)
rcop_irl <- convert_to_numeric(rcop_irl)

# Convert data to long format
make_long <- function(data, environment) {
  data %>%
    mutate(TrialNumber = row_number()) %>%
    pivot_longer(
      cols = -TrialNumber,
      names_to = "Subject",
      values_to = "Value"
    ) %>%
    mutate(
      Environment = environment,
      Block = ceiling(TrialNumber / 5),
      Trial = ((TrialNumber - 1) %% 5) + 1
    )
}

convert_to_numeric <- function(data) {
  data %>%
    mutate(across(everything(), ~ na_if(as.character(.x), ""))) %>%
    mutate(across(everything(), as.numeric))
}

# Combine environments and define factors
make_component_data <- function(vr_data, irl_data) {
  bind_rows(
    make_long(irl_data, "IRL"),
    make_long(vr_data, "VR")
  ) %>%
    mutate(
      Subject = factor(Subject),
      Environment = factor(Environment, levels = c("IRL", "VR")),
      Trial = factor(Trial, levels = 1:5),
      Block = factor(Block, levels = 1:5)
    )
}

cop_data <- make_component_data(cop_vr, cop_irl)
lf_data <- make_component_data(lf_vr, lf_irl)
gf_data <- make_component_data(gf_vr, gf_irl)
rcop_data <- make_component_data(rcop_vr, rcop_irl)

lf_data <- lf_data %>% filter(!is.na(Value))

# Run Trial x Environment ANOVA for each component
run_component_anova <- function(data, component_name) {
  cat("\n", component_name, "ANOVA\n")
  
  results <- aov_ez(
    id = "Subject",
    dv = "Value",
    data = data,
    within = c("Trial", "Environment"),
    type = 3
  )
  
  print(results)
  return(results)
}

cop_anova <- run_component_anova(cop_data, "COP")
lf_anova <- run_component_anova(lf_data, "LF")
gf_anova <- run_component_anova(gf_data, "GF")

# Run Trial x Environment x Block ANOVA for RCOP
cat("\nRCOP ANOVA\n")

rcop_anova <- aov_ez(
  id = "Subject",
  dv = "Value",
  data = rcop_data,
  within = c("Trial", "Environment", "Block"),
  type = 3,
  fun_aggregate = mean
)

print(rcop_anova)

warnings()

library(emmeans)

# COP: Compare trials
emmeans(cop_anova, pairwise ~ Trial, adjust = "tukey")

# LF: Compare trials
emmeans(lf_anova, pairwise ~ Trial, adjust = "tukey")

