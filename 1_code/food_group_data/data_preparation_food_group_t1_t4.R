no_funciton()
rm(list = ls())
library(tidyverse)

##load data
library(data.table)

setwd(r4projects::get_project_wd())

load("3-data_analysis/nutrition_t1_t4/data_preparation/sample_info")

sample_info_old <-
  sample_info

nutrition_data <-
  data.table::fread("2-data/Food_Cat_Weight_agg.csv")

colnames(nutrition_data)

##food group from Beverages
##impute NA with zero
sample_info <-
  nutrition_data %>%
  dplyr::select(Participant.ID:Time) %>% 
  dplyr::left_join(sample_info_old, by = c("Participant.ID", "Time"))

expression_data <-
  nutrition_data %>%
  dplyr::select(Beverages:Misc) %>% 
  as.data.frame()

expression_data[is.na(expression_data)] <- 0

variable_info =
  data.frame(variable_id = colnames(expression_data))

expression_data <-
  t(expression_data) %>%
  as.data.frame()

colnames(expression_data) <-
  sample_info$sample_id

colnames(expression_data) == sample_info$sample_id

rownames(expression_data) == variable_info$variable_id

dir.create("3-data_analysis/food_group/data_preparation", recursive = TRUE)

expression_data %>%
  apply(1, function(x) {
    sum(is.na(x)) * 100/ ncol(expression_data)
  }) %>%
  plot(xlab = "Food group", ylab = "Missing value percentage (%)")

expression_data %>%
  apply(2, function(x) {
    sum(is.na(x)) / nrow(expression_data)
  }) %>%
  plot(xlab = "Sample", ylab = "Missing value percentage (%)")

save(sample_info, file = "3-data_analysis/food_group/data_preparation/sample_info")
save(expression_data, file = "3-data_analysis/food_group/data_preparation/expression_data")
save(variable_info, file = "3-data_analysis/food_group/data_preparation/variable_info")
