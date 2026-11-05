###### POST-SWITCH MCom #####
# created by CAS 28 Jul 2026
# modified by CAS 13 Oct 2026
# Repeated measures ANOVA

library(openxlsx)
library(tidyverse)
library(afex)
library(emmeans)

base <- "/Users/motorskill/Desktop/R01s/vr/rep_vr/data_analyses/"

rep1_vr  <- read.xlsx(file.path(base, "post_switch_rep1_vr.xlsx"), colNames = FALSE)
rep5_vr  <- read.xlsx(file.path(base, "post_switch_rep5_vr.xlsx"), colNames = FALSE)
rep1_irl <- read.xlsx(file.path(base, "post_switch_rep1_irl.xlsx"), colNames = FALSE)
rep5_irl <- read.xlsx(file.path(base, "post_switch_rep5_irl.xlsx"), colNames = FALSE)

# Participants who completed VR first, then IRL
VR_first <- c("1", "2", "3", "4", "9", "11", "13", "15")

# Participants who completed IRL first, then VR
IRL_first <- c("5", "6", "7", "8", "10", "12", "14")

# Convert each file to long format
make_long <- function(data, repetition, environment) {
  # The first row contains participant IDs
  participant_ids <- as.character(unlist(data[1, ], use.names = FALSE))
  colnames(data) <- participant_ids

  # Remove the row containing participant IDs
  data <- data[-1, , drop = FALSE]

  # Convert observations to numeric
  data[] <- lapply(data, function(x) {
    as.numeric(as.character(x))
  })

  # Add block numbers and reshape to long format
  long_data <- data %>%
    mutate(Block = row_number()) %>%
    pivot_longer(
      cols = -Block,
      names_to = "Subject",
      values_to = "MCom"
    ) %>%
    mutate(
      Subject = as.character(Subject),
      Repetition = factor(repetition, levels = c(1, 5), labels = c("X1", "X5")),
      Environment = factor(environment, levels = c("IRL", "VR")),
      Order = case_when(
        Subject %in% VR_first ~ "VR_first",
        Subject %in% IRL_first ~ "IRL_first",
        TRUE ~ NA_character_
      ),
      Block = factor(Block, levels = 1:5, labels = paste0("X", 1:5)),
      Order = factor(Order, levels = c("IRL_first", "VR_first"))
    )

  return(long_data)
}

# Build the data set; each file is included once
post <- bind_rows(
  make_long(rep1_vr,  repetition = 1, environment = "VR"),
  make_long(rep5_vr,  repetition = 5, environment = "VR"),
  make_long(rep1_irl, repetition = 1, environment = "IRL"),
  make_long(rep5_irl, repetition = 5, environment = "IRL")
)

# Confirm number of observations and participants
cat("\nNumber of observations:", nrow(post), "\n")
cat("Number of participants:", n_distinct(post$Subject), "\n")

# Check for missing values
cat("\nMissing MCom observations:", sum(is.na(post$MCom)), "\n")
cat("Missing participant order assignments:", sum(is.na(post$Order)), "\n")

# Descriptive statistics
descriptive_stats <- post %>%
  group_by(Repetition, Block, Environment, Order) %>%
  summarise(
    Mean_MCom = mean(MCom, na.rm = TRUE),
    SD_MCom = sd(MCom, na.rm = TRUE),
    N = sum(!is.na(MCom)),
    .groups = "drop"
  )

print(descriptive_stats)

# Mixed repeated-measures ANOVA
# Within-subject factors: Repetition, Block, Environment
# Between-subject factor: Order
anova_model <- aov_ez(
  id = "Subject",
  dv = "MCom",
  data = post,
  within = c("Repetition", "Block", "Environment"),
  between = "Order",
  type = 3,
  return = "afex_aov"
)

# Display ANOVA results
print(anova_model)

anova_results <- nice(anova_model, correction = "GG", es = "pes")
print(anova_results)

# Main effect: Repetition
emm_repetition <- emmeans(anova_model, ~ Repetition)
print(emm_repetition)
pairs(emm_repetition, adjust = "tukey")

# Main effect: Block
emm_block <- emmeans(anova_model, ~ Block)
print(emm_block)
pairs(emm_block, adjust = "tukey")

# Main effect: Environment
emm_environment <- emmeans(anova_model, ~ Environment)
print(emm_environment)
pairs(emm_environment, adjust = "tukey")

# Main effect: Order
emm_order <- emmeans(anova_model, ~ Order)
print(emm_order)
pairs(emm_order, adjust = "tukey")

# Repetition × Environment
emm_rep_environment <- emmeans(
  anova_model,
  pairwise ~ Repetition | Environment,
  adjust = "tukey"
)
print(emm_rep_environment)

# Repetition × Block
emm_rep_block <- emmeans(
  anova_model,
  pairwise ~ Repetition | Block,
  adjust = "tukey"
)
print(emm_rep_block)

# Block × Environment
emm_block_environment <- emmeans(
  anova_model,
  pairwise ~ Block | Environment,
  adjust = "tukey"
)
print(emm_block_environment)

# Block × Order
emm_block_order <- emmeans(
  anova_model,
  pairwise ~ Block | Order,
  adjust = "tukey"
)
print(emm_block_order)

# Repetition × Block × Environment
emm_block_by_rep_environment <- emmeans(
  anova_model,
  pairwise ~ Block | Repetition * Environment,
  adjust = "tukey"
)
print(emm_block_by_rep_environment)
