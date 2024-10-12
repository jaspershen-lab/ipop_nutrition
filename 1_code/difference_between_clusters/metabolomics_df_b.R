no_function()
setwd(r4projects::get_project_wd())
library(tidyverse)
library(plyr)
rm(list = ls())
source("1_code/tools.R")

###load data
{
  ##metabolome
  load("3_data_analysis/metabolomics/data_preparation/expression_data")
  load("3_data_analysis/metabolomics/data_preparation/sample_info")
  load("3_data_analysis/metabolomics/data_preparation/variable_info")
  
  metabolomics_expression_data = expression_data
  metabolomics_sample_info = sample_info
  metabolomics_variable_info = variable_info
}

df_b = readr::read_csv("2-data/DF_B_Energy.csv")
df_a = readr::read_csv("2-data/DF_A_Energy_SampleID.csv")


dir.create("3_data_analysis/metabolomics/different_metabolomics_dfb")
setwd("3_data_analysis/metabolomics/different_metabolomics_dfb")

expression_data =
  df_b$subjectid %>%
  purrr::map(function(x) {
    cat(x, " ")
    idx = which(df_a$subjectid == x)
    sample_id = df_a$SampleID[idx]
    idx = match(sample_id, colnames(expression_data))
    idx = idx[!is.na(idx)]
    if (length(idx) == 0) {
      return(NA)
    } else{
      apply(metabolomics_expression_data[, idx, drop = FALSE], 1, mean)
    }
  }) %>%
  do.call(cbind, .) %>%
  as.data.frame()

colnames(expression_data) = df_b$subjectid
rownames(expression_data)

sample_info =
  df_b %>%
  dplyr::select(subjectid, clustering_tSNE, clustering_hierarchical) %>%
  dplyr::rename(sample_id = subjectid) %>%
  dplyr::mutate(subject_id = sample_id,
                class = "",
                group = "")

library(tidymass)
object = create_mass_dataset(
  expression_data = expression_data,
  sample_info = sample_info,
  variable_info = metabolomics_variable_info
)

object =
  object %>%
  mutate_sample_na_number() %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(na_number == 0)

###find markers based on clustering_tSNE
object %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::count(clustering_tSNE)

library(massstat)

case_sample_id =
  object %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(clustering_tSNE == 2) %>%
  pull(sample_id)

control_sample_id =
  object %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(clustering_tSNE == 1) %>%
  pull(sample_id)

object1 =
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

plot1 =
  massstat::volcano_plot(
    object = object1,
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
plot1

ggsave(plot1,
       filename = "volocao_plot_clustering_tSNE.pdf",
       width = 7,
       height = 7)

export_mass_dataset(object = object1, path = "clustering_tSNE")









###find markers based on clustering_tSNE
object %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::count(clustering_hierarchical)

library(massstat)

case_sample_id =
  object %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(clustering_hierarchical == 2) %>%
  pull(sample_id)

control_sample_id =
  object %>%
  activate_mass_dataset(what = "sample_info") %>%
  filter(clustering_hierarchical == 1) %>%
  pull(sample_id)

object2 =
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

plot2 =
  massstat::volcano_plot(
    object = object2,
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
plot2

ggsave(plot2,
       filename = "volocan_plot_clustering_hierarchical.pdf",
       width = 7,
       height = 7)

export_mass_dataset(object = object1, path = "clustering_hierarchical")
