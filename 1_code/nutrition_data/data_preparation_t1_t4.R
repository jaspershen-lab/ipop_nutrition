no_funciton()
rm(list = ls())
library(tidyverse)

##load data
library(data.table)

setwd(r4projects::get_project_wd())

nutrition_data =
  data.table::fread("2-data/Diet_Glucose_Marker_T1_through_T4_Healthy_Only.csv")

colnames(nutrition_data)

sample_info =
  nutrition_data %>%
  dplyr::select(subjectid:Energy, A1C_Date:sspg_status)

sample_info %>%
  dplyr::group_by(subjectid) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::filter(n == 4)

expression_data =
  nutrition_data %>%
  dplyr::select(Energy:Chln)

variable_info =
  data.frame(variable_id = colnames(expression_data))

expression_data = t(expression_data) %>%
  as.data.frame()

sample_info$sample_id = paste(sample_info$subjectid, sample_info$Time, sep = "_")

colnames(expression_data) = sample_info$sample_id

sample_info =
  sample_info %>%
  dplyr::rename(subject_id = subjectid) %>%
  dplyr::select(sample_id, subject_id, everything())

colnames(expression_data) == sample_info$sample_id

rownames(expression_data) == variable_info$variable_id

save(sample_info, file = "3-data_analysis/nutrition_t1_t4/data_preparation/sample_info")
save(expression_data, file = "3-data_analysis/nutrition_t1_t4/data_preparation/expression_data")
save(variable_info, file = "3-data_analysis/nutrition_t1_t4/data_preparation/variable_info")
