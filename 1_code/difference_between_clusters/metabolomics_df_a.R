no_function()
setwd(r4projects::get_project_wd())
library(tidyverse)
library(plyr)
rm(list = ls())
source("1-code/tools.R")

###load data
{
  ##metabolome
  load("3-data_analysis/metabolomics/data_preparation/expression_data")
  load("3-data_analysis/metabolomics/data_preparation/sample_info")
  load("3-data_analysis/metabolomics/data_preparation/variable_info")
  
  metabolomics_expression_data = expression_data
  metabolomics_sample_info = sample_info
  metabolomics_variable_info = variable_info
}

library(massdataset)

object = create_mass_dataset(
  expression_data = metabolomics_expression_data,
  sample_info = metabolomics_sample_info %>%
    dplyr::mutate(class = ""),
  variable_info = metabolomics_variable_info
)

object

df_a = readr::read_csv("2-data/DF_A_Energy.csv")

df_a2 = readr::read_csv("2-data/DF_A_Energy_SampleID.csv")

new_sample_info =
  df_a2 %>%
  dplyr::select(SampleID, clustering_tSNE) %>%
  dplyr::rename(sample_id = SampleID,
                cluster = clustering_tSNE)

object =
  object %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::left_join(new_sample_info, by = "sample_id") %>%
  dplyr::filter(!is.na(cluster))

object %>%
  extract_sample_info()

dir.create("3-data_analysis/metabolomics/different_metabolomics_dfa")
setwd("3-data_analysis/metabolomics/different_metabolomics_dfa")

###find markers
object %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::count(cluster)

library(massstat)

case_sample_id =
  object %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(cluster == 2) %>%
  pull(sample_id)

control_sample_id =
  object %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(cluster == 1) %>%
  pull(sample_id)

object =
  object %>%
  massstat::mutate_fc(
    control_sample_id = control_sample_id,
    case_sample_id = case_sample_id,
    mean_median = "median"
  ) %>%
  `+`(1) %>%
  log(2) %>%
  mutate_p_value(
    control_sample_id = control_sample_id,
    case_sample_id = case_sample_id,
    method = "wilcox",
    p_adjust_methods = "BH"
  )

plot =
  massstat::volcano_plot(
    object = object,
    fc_up_cutoff = 1,
    fc_down_cutoff = 1,
    p_value_column_name = "p_value_adjust",
    up_color = ggsci::pal_aaas()(n = 10)[2],
    down_color = ggsci::pal_aaas()(n = 10)[1],
    point_size_scale = "p_value_adjust",
    line_type = 2,
    add_text = TRUE,
    text_for = "marker",
    text_from = "Metabolite"
  )

plot

ggsave(plot,
       filename = "volocao_plot.pdf",
       width = 7,
       height = 7)

export_mass_dataset(object = object)
