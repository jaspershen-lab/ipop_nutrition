no_funciton()
rm(list = ls())
library(tidyverse)

##load data
library(data.table)
nutrition_data =
  data.table::fread("2-data/T1_Glucose_Markers_Diet_Demographics.csv")


colnames(nutrition_data)

sample_info = 
  nutrition_data %>% 
  dplyr::select(subjectid:DSD_Unix_Date, CollectionDate:baseline_ogtt_status)

expression_data = 
  nutrition_data %>% 
  dplyr::select(Energy:Chln)

variable_info = 
  data.frame(variable_id = colnames(expression_data))

expression_data = t(expression_data) %>%
  as.data.frame()

colnames(expression_data) = sample_info$subjectid

sample_info = 
sample_info %>% 
  dplyr::rename(subject_id = subjectid) %>% 
  dplyr::mutate(sample_id = subject_id) %>% 
  dplyr::select(sample_id, subject_id, everything())

colnames(expression_data) == sample_info$sample_id

rownames(expression_data) == variable_info$variable_id

save(sample_info, file = "3-data_analysis/nutrition/data_preparation/sample_info")
save(expression_data, file = "3-data_analysis/nutrition/data_preparation/expression_data")
save(variable_info, file = "3-data_analysis/nutrition/data_preparation/variable_info")

