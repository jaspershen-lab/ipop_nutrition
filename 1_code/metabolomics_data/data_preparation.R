no_funciton()
rm(list = ls())
library(tidyverse)
library(plyr)
tinyTools::setwd_project()

####load raw data
load("2-data/Revision_MultiOmes_0509.RData")

variable_info = metb.curated %>% 
  dplyr::select(Compounds_ID, dplyr::everything()) %>% 
  dplyr::rename(variable_id = Compounds_ID)

expression_data = metbcr.df

sample_info = 
  expression_data %>% 
  dplyr::select(SampleID, SubjectID:CL4) %>% 
  dplyr::rename(sample_id = SampleID,
                subject_id = SubjectID) %>% 
  dplyr::left_join(sc, by = c("subject_id" = "SubjectID"))

expression_data = 
  expression_data %>% 
  dplyr::select(-c(SampleID, SubjectID:CL4)) %>% 
  t() %>% 
  as.data.frame()

dim(variable_info)
dim(expression_data)
dim(sample_info)

colnames(expression_data) = sample_info$sample_id

setdiff(rownames(expression_data),
        variable_info$variable_id)

setdiff(variable_info$variable_id,
        rownames(expression_data))

length(variable_info$variable_id)
length(unique(variable_info$variable_id))

which(duplicated(variable_info$variable_id))
variable_info %>% 
  dplyr::filter(variable_id %in% variable_info$variable_id[c(383, 435)])

variable_info = 
  variable_info %>% 
  dplyr::distinct(variable_id, .keep_all = TRUE)

dim(variable_info)

variable_info$variable_id == rownames(expression_data)
sort(variable_info$variable_id) == sort(rownames(expression_data))
expression_data = 
  expression_data[variable_info$variable_id,]

rownames(expression_data) == variable_info$variable_id
colnames(expression_data) == sample_info$sample_id

variable_info = 
  variable_info %>% 
  dplyr::filter(!stringr::str_detect(Metabolite, "C[0-9]{1,2}H[0-9]{1,2}")) %>% 
  dplyr::filter(!stringr::str_detect(Metabolite, "Unknown")) %>% 
  dplyr::mutate(Metabolite = stringr::str_replace_all(Metabolite, '\"', "")) %>% 
  dplyr::mutate(true_name = 
                  stringr::str_replace_all(Metabolite, "\\([0-9]{1}\\)", ""))

mean_int = 
  expression_data[variable_info$variable_id,] %>% 
  apply(1, mean)

variable_info$mean_int = mean_int

variable_info = 
  variable_info %>% 
  plyr::dlply(.variables = .(true_name)) %>% 
  purrr::map(function(x){
    x %>% 
      dplyr::filter(MSMS == "YES") %>% 
      dplyr::filter(mean_int == max(mean_int))
  }) %>% 
  do.call(rbind, .) %>% 
  as.data.frame()

variable_info$Metabolite = variable_info$true_name

expression_data = expression_data[variable_info$variable_id,]

save(sample_info, file = "3_data_analysis/metabolomics/data_preparation/sample_info")
save(expression_data, file = "3_data_analysis/metabolomics/data_preparation/expression_data")
save(variable_info, file = "3_data_analysis/metabolomics/data_preparation/variable_info")

