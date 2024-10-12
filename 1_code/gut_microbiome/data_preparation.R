no_funciton()
rm(list = ls())
library(tidyverse)
library(plyr)
tinyTools::setwd_project()

####load raw data
load("2-data/Revision_MultiOmes_0509.RData")

expression_data = st.df

sample_info = 
  expression_data %>% 
  dplyr::select(HostSampleID, collection.date:CollectionDate) %>% 
  dplyr::rename(sample_id = HostSampleID,
                subject_id = SubjectID) %>% 
  dplyr::left_join(sc, by = c("subject_id" = "SubjectID"))

expression_data = 
  expression_data %>% 
  dplyr::select(-c(HostSampleID, collection.date:CollectionDate)) %>% 
  t() %>% 
  as.data.frame()

colnames(expression_data) = sample_info$sample_id

variable_info = 
  data.frame(variable_id = rownames(expression_data)) %>%
  dplyr::mutate(level = unlist(lapply(stringr::str_split(variable_id, "_", n = 2), function(x)
    x[1])))

table(variable_info$level)

sum(expression_data$`69-001-1010`[which(variable_info$level == "genus")])

colSums(expression_data[which(variable_info$level == "genus"),]) %>% plot


dim(variable_info)
dim(expression_data)
dim(sample_info)

colnames(expression_data) = sample_info$sample_id

save(sample_info, file = "3_data_analysis/gut_microbiome/data_preparation/sample_info")
save(expression_data, file = "3_data_analysis/gut_microbiome/data_preparation/expression_data")
save(variable_info, file = "3_data_analysis/gut_microbiome/data_preparation/variable_info")

