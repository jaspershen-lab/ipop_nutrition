no_function()
setwd(r4projects::get_project_wd())
library(tidyverse)
library(plyr)
rm(list = ls())
source("1-code/tools.R")

###load data
{
  ##nutrition
  load("3-data_analysis/food_group/data_preparation/expression_data")
  load("3-data_analysis/food_group/data_preparation/sample_info")
  load("3-data_analysis/food_group/data_preparation/variable_info")
  
  nutrition_expression_data = expression_data
  nutrition_sample_info = sample_info
  nutrition_variable_info = variable_info
  
  filter_subject_id =
    nutrition_sample_info %>%
    dplyr::group_by(subject_id) %>%
    dplyr::summarise(n = n()) %>%
    dplyr::filter(n == 4) %>%
    pull(subject_id)
  
  nutrition_sample_info =
    nutrition_sample_info %>%
    dplyr::filter(subject_id %in% filter_subject_id)
  
  nutrition_expression_data =
    nutrition_expression_data[, nutrition_sample_info$sample_id]
  
  ##metabolomics
  load("3-data_analysis/metabolomics/data_preparation/expression_data")
  load("3-data_analysis/metabolomics/data_preparation/sample_info")
  load("3-data_analysis/metabolomics/data_preparation/variable_info")
  
  metabolomics_expression_data = expression_data
  metabolomics_sample_info = sample_info
  metabolomics_variable_info = variable_info
}

dir.create("3-data_analysis/food_group_vs_metabolomics_t1_t4/based_on_sspg/")
setwd("3-data_analysis/food_group_vs_metabolomics_t1_t4/based_on_sspg/")

###data preparation
nutrition_sample_info$Diet.Survey.Date
nutrition_sample_info$CollectionDate

####match metabolomics data from baseline
metabolomics_sample_info =
  metabolomics_sample_info %>%
  dplyr::mutate(CollectionDate = as.Date(CollectionDate, format = "%m/%d/%y"))

matched_sample <-
  nutrition_sample_info %>%
  t() %>%
  as.data.frame() %>%
  purrr::map(function(x) {
    subject_id1 = x[4]
    date1 = as.Date(x[5])
    temp =
      metabolomics_sample_info %>%
      dplyr::filter(subject_id == subject_id1)
    temp %>%
      dplyr::mutate(diff_days = as.numeric(CollectionDate - date1)) %>%
      dplyr::filter(abs(diff_days) < 7) %>%
      dplyr::filter(abs(diff_days) == min(abs(diff_days))) %>%
      head(1)
  })

matched_sample

matched_sample %>%
  lapply(nrow) %>%
  unlist() %>%
  plot

matched_sample =
  matched_sample %>%
  lapply(function(x) {
    temp = x$sample_id
    if (length(temp) == 0) {
      temp = NA
    }
    temp
  }) %>%
  unlist()

matched_idx =
  data.frame(sample_id1 = nutrition_sample_info$sample_id,
             sample_id2 = unname(matched_sample)) %>%
  dplyr::filter(!is.na(sample_id2))


matched_idx

cbind(nutrition_sample_info[match(matched_idx$sample_id1, nutrition_sample_info$sample_id),
                            c("sample_id", "Diet.Survey.Date")],
      metabolomics_sample_info[match(matched_idx$sample_id2, metabolomics_sample_info$sample_id),
                               c("sample_id", "CollectionDate")])

matched_idx =
  matched_idx %>%
  dplyr::left_join(nutrition_sample_info[, c("subject_id", "sample_id")],
                   by = c("sample_id1" = "sample_id"))

temp_subject_id =
  matched_idx %>%
  dplyr::group_by(subject_id) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::filter(n == 4) %>%
  dplyr::pull(subject_id)

matched_idx =
  matched_idx %>%
  dplyr::filter(subject_id %in% temp_subject_id) %>%
  dplyr::select(-subject_id)

nutrition_sample_info =
  nutrition_sample_info %>%
  dplyr::filter(sample_id %in% matched_idx$sample_id1)

nutrition_expression_data =
  nutrition_expression_data[, nutrition_sample_info$sample_id]

nutrition_expression_data = nutrition_expression_data[, matched_idx$sample_id1]
metabolomics_expression_data = metabolomics_expression_data[, matched_idx$sample_id2]

nutrition_sample_info =
  nutrition_sample_info[match(colnames(nutrition_expression_data),
                              nutrition_sample_info$sample_id), ]

metabolomics_sample_info =
  metabolomics_sample_info[match(colnames(metabolomics_expression_data),
                                 metabolomics_sample_info$sample_id), ]

dim(nutrition_expression_data)
dim(metabolomics_expression_data)

colnames(metabolomics_expression_data) =
  colnames(nutrition_expression_data)

# nutrition_expression_data[is.na(nutrition_expression_data)] <- 0

######calculate the correlation between nutrition and metabolites
metabolomics_expression_data =
  log(metabolomics_expression_data + 1, 2)

metabolomics_expression_data =
  metabolomics_expression_data %>%
  apply(1, function(x) {
    (x) / sd(x)
  }) %>%
  t() %>%
  as.data.frame()

###scale
nutrition_expression_data =
  nutrition_expression_data %>%
  apply(1, function(x) {
    (x) / sd(x)
  }) %>%
  t() %>%
  as.data.frame()

nutrition_expression_data[is.na(nutrition_expression_data)] <- 0

range(nutrition_expression_data)
range(metabolomics_expression_data)

###############################################################################
####based on the sspg status
####calculate the correlation between them
###cor_data for IS people
sample_info_IS =
  nutrition_sample_info %>%
  dplyr::filter(!is.na(sspg_status)) %>%
  dplyr::filter(sspg_status == "IS")

sample_info_IR =
  nutrition_sample_info %>%
  dplyr::filter(!is.na(sspg_status)) %>%
  dplyr::filter(sspg_status == "IR")

load("cor_data_IS")
load("cor_data_IR")

dim(cor_data_IS)
dim(cor_data_IR)

which(cor_data_IR$p_adjust < 0.05 & abs(cor_data_IR$cor) > 0.5)
which(cor_data_IS$p_adjust < 0.05& abs(cor_data_IS$cor) > 0.5)

cor_data_IS <- 
  cor_data_IS %>% 
  dplyr::filter(p_adjust < 0.05 & abs(cor_data_IS$cor) > 0.5)

cor_data_IR <- 
  cor_data_IR %>% 
  dplyr::filter(p_adjust < 0.05 & abs(cor_data_IR$cor) > 0.5)

kegg_id <-
  variable_info %>%
  dplyr::filter(variable_id %in% unique(cor_data_IS$data_set2)) %>%
  pull(KEGG)

kegg_id <-
  kegg_id[kegg_id != ""] %>%
  stringr::str_split(pattern = "\\|") %>%
  unlist() %>%
  unique()

data("kegg_hsa_pathway", package = "metpath")
kegg_hsa_pathway

# result_is <-
#   enrich_kegg(
#     query_id = kegg_id,
#     query_type = "compound",
#     id_type = "KEGG",
#     pathway_database = kegg_hsa_pathway,
#     p_cutoff = 0.05,
#     p_adjust_method = "BH",
#     threads = 3
#   )
# 
# save(result_is, file = "result_is")

load("result_is")

plot <-
  enrich_bar_plot(object = result_is,
                  x_axis = "p_value",
                  top = 10)
plot
ggsave(plot,
       filename = "enrichment_result_is.pdf",
       width = 7,
       height = 7)


write.csv(
  result_is@result %>% dplyr::filter(p_value < 0.05),
  file = "enrichment_result_is.csv",
  row.names = FALSE
)


kegg_id <-
  variable_info %>%
  dplyr::filter(variable_id %in% unique(cor_data_IR$data_set2)) %>%
  pull(KEGG)

kegg_id <-
  kegg_id[kegg_id != ""] %>%
  stringr::str_split(pattern = "\\|") %>%
  unlist() %>%
  unique()

# result_ir <-
#   enrich_kegg(
#     query_id = kegg_id,
#     query_type = "compound",
#     id_type = "KEGG",
#     pathway_database = kegg_hsa_pathway,
#     p_cutoff = 0.05,
#     p_adjust_method = "BH",
#     threads = 3
#   )
# 
# save(result_ir, file = "result_ir")

load("result_ir")

plot <-
  enrich_bar_plot(object = result_ir,
                  x_axis = "p_value",
                  top = 10)
plot
ggsave(plot,
       filename = "enrichment_result_ir.pdf",
       width = 7,
       height = 7)

write.csv(
  result_ir@result %>% dplyr::filter(p_value < 0.05),
  file = "enrichment_result_ir.csv",
  row.names = FALSE
)
