no_function()
setwd(r4projects::get_project_wd())
library(tidyverse)
library(plyr)
rm(list = ls())
source("1-code/tools.R")

###load data
{
  ##nutrition
  load("3-data_analysis/nutrition/data_preparation/expression_data")
  load("3-data_analysis/nutrition/data_preparation/sample_info")
  load("3-data_analysis/nutrition/data_preparation/variable_info")
  
  nutrition_expression_data = expression_data
  nutrition_sample_info = sample_info
  nutrition_variable_info = variable_info
  
  ###remove some nutrients
  rownames(nutrition_expression_data)
  nutrition_variable_info <-
    nutrition_variable_info %>%
    dplyr::filter(
      !variable_id %in% c(
        "Fib",
        "Carb",
        "Pro",
        "Fat",
        "FatCals",
        "SatCals",
        "VitD_mcg",
        "Fol_DFE",
        "OCarb",
        "VitA_RAE"
      )
    )
  nutrition_expression_data <-
    nutrition_expression_data[nutrition_variable_info$variable_id, ]
  
  ##microbiome
  load("3-data_analysis/gut_microbiome/data_preparation/expression_data")
  load("3-data_analysis/gut_microbiome/data_preparation/sample_info")
  load("3-data_analysis/gut_microbiome/data_preparation/variable_info")
  
  microbiome_expression_data = expression_data
  microbiome_sample_info = sample_info
  microbiome_variable_info = variable_info
}

setwd("3-data_analysis/nutrition_vs_microbiome/based_on_ogtt/")

###data preparation
nutrition_sample_info$Diet.Survey.Date
nutrition_sample_info$CollectionDate

####match microbiome data from baseline
microbiome_sample_info =
  microbiome_sample_info %>%
  dplyr::mutate(CollectionDate = as.Date(CollectionDate, format = "%m/%d/%y"))

matched_sample =
  nutrition_sample_info %>%
  t() %>%
  as.data.frame() %>%
  purrr::map(function(x) {
    subject_id1 = x[1]
    date1 = as.Date(x[15])
    temp =
      microbiome_sample_info %>%
      dplyr::filter(subject_id == subject_id1)
    temp %>%
      dplyr::mutate(diff_days = as.numeric(CollectionDate - date1)) %>%
      dplyr::filter(abs(diff_days) < 7) %>%
      dplyr::filter(abs(diff_days) == min(abs(diff_days)))
  })

matched_sample

matched_sample %>% 
  lapply(nrow) %>% 
  unlist() %>% 
  plot

matched_sample = 
  matched_sample %>% 
  lapply(function(x){
    temp = x$sample_id
    if(length(temp) == 0){
      temp = NA
    }
    temp
  }) %>% 
  unlist()

matched_idx = 
  data.frame(sample_id1 = nutrition_sample_info$sample_id,
             sample_id2 = unname(matched_sample)) %>% 
  dplyr::filter(!is.na(sample_id2))

nutrition_expression_data = nutrition_expression_data[,matched_idx$sample_id1]
microbiome_expression_data = microbiome_expression_data[,matched_idx$sample_id2]

nutrition_sample_info = 
  nutrition_sample_info[match(colnames(nutrition_expression_data),
                              nutrition_sample_info$sample_id),]

microbiome_sample_info = 
  microbiome_sample_info[match(colnames(microbiome_expression_data), microbiome_sample_info$sample_id),]

dim(nutrition_expression_data)
dim(microbiome_expression_data)

colnames(microbiome_expression_data) =
  colnames(nutrition_expression_data)

######calculate the correlation between nutrition and microbiomes
####missing value
nutrition_expression_data %>% 
  apply(1, function(x){
    sum(is.na(x))
  }) %>% 
  plot()

library(impute)

nutrition_expression_data =
  impute::impute.knn(data = as.matrix(nutrition_expression_data))$data %>%
  as.data.frame()

nutrition_expression_data = 
  nutrition_expression_data %>% 
  apply(1, function(x){
    (x - mean(x)) / sd(x)
  }) %>% 
  t() %>% 
  as.data.frame()

###############################################################################
####based on the sspg status
####calculate the correlation between them
###cor_data for IS people
sample_info_Normal =
  nutrition_sample_info %>%
  dplyr::filter(!is.na(ogtt_t2d_status)) %>%
  dplyr::filter(ogtt_t2d_status == "Normal")

sample_info_preDM =
  nutrition_sample_info %>%
  dplyr::filter(!is.na(ogtt_t2d_status)) %>%
  dplyr::filter(ogtt_t2d_status == "preDM" | ogtt_t2d_status == "T2D")

# cor_data_Normal =
#   partial_cor(
#     data_set1 = nutrition_expression_data[, sample_info_Normal$sample_id],
#     data_set2 = microbiome_expression_data[, sample_info_Normal$sample_id],
#     sample_info = sample_info_Normal,
#     method = "spearman",
#     threads = 5
#   )
# 
# cor_data_Normal = cor_data_Normal[[1]]
# 
# cor_data_preDM =
#   partial_cor(
#     data_set1 = nutrition_expression_data[, sample_info_preDM$sample_id],
#     data_set2 = microbiome_expression_data[, sample_info_preDM$sample_id],
#     sample_info = sample_info_preDM,
#     method = "spearman",
#     threads = 5
#   )
# 
# cor_data_preDM = cor_data_preDM[[1]]
# 
# save(cor_data_Normal, file = "cor_data_Normal")
# save(cor_data_preDM, file = "cor_data_preDM")

load("cor_data_Normal")
load("cor_data_preDM")

dim(cor_data_Normal)
dim(cor_data_preDM)

library(openxlsx)
cor_data_Normal_output = 
  cor_data_Normal %>% 
  dplyr::left_join(microbiome_variable_info, by = c("data_set2" = "variable_id")) %>% 
  dplyr::select(-p_adjust2) %>% 
  dplyr::filter(p_adjust < 0.2)

cor_data_preDM_output = 
  cor_data_preDM %>% 
  dplyr::left_join(microbiome_variable_info, by = c("data_set2" = "variable_id")) %>% 
  dplyr::select(-p_adjust2) %>% 
  dplyr::filter(p_adjust < 0.2)

# openxlsx::write.xlsx(
#   cor_data_Normal_output,
#   "cor_data_Normal_output.xlsx",
#   asTable = TRUE,
#   overwrite = TRUE
# )
# 
# openxlsx::write.xlsx(
#   cor_data_preDM_output,
#   "cor_data_preDM_output.xlsx",
#   asTable = TRUE,
#   overwrite = TRUE
# )

which(cor_data_preDM$p_adjust < 0.05)
which(cor_data_Normal$p_adjust < 0.05)

plot(cor_data_preDM$cor)
plot(cor_data_Normal$cor)

####output results
library(openxlsx)

#####output the cor plot
idx = which(cor_data_preDM$p_adjust < 0.05)

cor_data_preDM %>%
  dplyr::arrange(desc(abs(cor))) %>% 
  head()

plot(as.numeric(nutrition_expression_data["Disacc", sample_info_preDM$sample_id]),
     as.numeric(microbiome_expression_data["genus_unclassified_Erysipelotrichaceae", sample_info_preDM$sample_id]))

cor(
  as.numeric(nutrition_expression_data["Disacc", sample_info_preDM$sample_id]),
  as.numeric(microbiome_expression_data["genus_unclassified_Erysipelotrichaceae", sample_info_preDM$sample_id]),
  method = "spearman"
)

temp_sample_info = 
  sample_info_preDM[,c("Sex", "Age")]

temp_sample_info$Sex[temp_sample_info$Sex == 'F'] = 0
temp_sample_info$Sex[temp_sample_info$Sex == 'M'] = 1
temp_sample_info$Sex = as.numeric(temp_sample_info$Sex)

ppcor::pcor.test(
  x = as.numeric(nutrition_expression_data["Disacc", sample_info_preDM$sample_id]),
  y = as.numeric(microbiome_expression_data["genus_unclassified_Erysipelotrichaceae", sample_info_preDM$sample_id]),
  z = temp_sample_info,
  method = "spearman"
)

cor_data_Normal %>%
  dplyr::arrange(desc(abs(cor))) %>% 
  head()

unique(cor_data_Normal$data_set2)
unique(cor_data_Normal$data_set1)

######all microbiome
normal_cor =
  cor_data_Normal %>%
  dplyr::select(-c(p:p_adjust2)) %>%
  tidyr::pivot_wider(names_from = data_set2, values_from = "cor") %>%
  tibble::column_to_rownames(var = "data_set1")

normal_p =
  cor_data_Normal %>%
  dplyr::select(-c(p, p_adjust2, cor)) %>%
  tidyr::pivot_wider(names_from = data_set2, values_from = "p_adjust") %>%
  tibble::column_to_rownames(var = "data_set1")

predm_cor =
  cor_data_preDM %>%
  dplyr::select(-c(p:p_adjust2)) %>%
  tidyr::pivot_wider(names_from = data_set2, values_from = "cor") %>%
  tibble::column_to_rownames(var = "data_set1")

predm_p =
  cor_data_preDM %>%
  dplyr::select(-c(p, p_adjust2, cor)) %>%
  tidyr::pivot_wider(names_from = data_set2, values_from = "p_adjust") %>%
  tibble::column_to_rownames(var = "data_set1")

colnames(predm_cor) = 
  colnames(normal_cor) = 
  microbiome_variable_info$variable_id[match(colnames(predm_cor), microbiome_variable_info$variable_id)]

library(ComplexHeatmap)
library(circlize)

col_fun = circlize::colorRamp2(
  breaks = seq(-1, 1, length.out = 11),
  colors = rev(RColorBrewer::brewer.pal(n = 11, name = "BrBG"))
)

library(wesanderson)

plot = 
  Heatmap(
    t(predm_cor),
    col = col_fun,
    border = TRUE,
    show_column_dend = TRUE,
    show_row_dend = TRUE,
    clustering_distance_rows = "euclidean",
    clustering_method_rows = "complete",
    clustering_distance_columns = "euclidean",
    clustering_method_columns = "complete",
    name = "Spearman correlation",
    cluster_columns = TRUE,
    cluster_rows = TRUE,
    column_names_gp = gpar(cex = 0.5, rot = 45),
    row_names_gp = gpar(cex = 0.3),
    column_names_rot = 45,
    cell_fun = function(j, i, x, y, w, h, fill) {
      gb1 = textGrob("*")
      gb_w1 = convertWidth(grobWidth(gb1), "mm")
      gb_h1 = convertHeight(grobHeight(gb1), "mm")
      gb2 = textGrob(".")
      gb_w2 = convertWidth(grobWidth(gb2), "mm")
      gb_h2 = convertHeight(grobHeight(gb2), "mm")
      if (t(predm_p)[i, j] < 0.05) {
        grid.text("*", x, y - gb_h1 * 0.5 + gb_w1 * 0.4, gp = gpar(col = "red"))
      }
      
      if (t(predm_p)[i, j] > 0.05 & t(predm_p)[i, j] < 0.2) {
        grid.points(pch = 20, x = x, y = y, size = unit(0.3, "char"), gp = gpar(col = "red"))
      }
    }
  )

plot = ggplotify::as.ggplot(plot)
plot
# ggsave(plot, filename = "IR_cor_all_microbiome.pdf", width = 10, height = 10)

plot = 
  Heatmap(
    t(normal_cor),
    col = col_fun,
    border = TRUE,
    show_column_dend = TRUE,
    show_row_dend = TRUE,
    clustering_distance_rows = "euclidean",
    clustering_method_rows = "complete",
    clustering_distance_columns = "euclidean",
    clustering_method_columns = "complete",
    name = "Spearman correlation",
    cluster_columns = TRUE,
    cluster_rows = TRUE,
    column_names_gp = gpar(cex = 0.5, rot = 45),
    row_names_gp = gpar(cex = 0.3),
    column_names_rot = 45,
    cell_fun = function(j, i, x, y, w, h, fill) {
      gb1 = textGrob("*")
      gb_w1 = convertWidth(grobWidth(gb1), "mm")
      gb_h1 = convertHeight(grobHeight(gb1), "mm")
      gb2 = textGrob(".")
      gb_w2 = convertWidth(grobWidth(gb2), "mm")
      gb_h2 = convertHeight(grobHeight(gb2), "mm")
      if (t(normal_p)[i, j] < 0.05) {
        grid.text("*", x, y - gb_h1 * 0.5 + gb_w1 * 0.4, gp = gpar(col = "red"))
      }
      
      if (t(normal_p)[i, j] > 0.05 & t(normal_p)[i, j] < 0.2) {
        grid.points(pch = 20, x = x, y = y, size = unit(0.3, "char"), gp = gpar(col = "red"))
        # grid.text(label = ".", x = x, y - gb_h2 * 0.5 + gb_w1 * 0.4)
      }
    }
  )

plot = ggplotify::as.ggplot(plot)
plot
# ggsave(plot, filename = "IS_cor_all_microbiome.pdf", width = 10, height = 10)


######only show the microbiomes with at least one significant correlation with nutrition
###remove the microbiomes which have no significant cor (p_adjust 0.2) with nutrition
remove_microbiome1 = 
  cor_data_Normal %>% 
  dplyr::group_by(data_set2) %>% 
  dplyr::summarise(n = sum(p_adjust < 0.2)) %>% 
  dplyr::filter(n == 0) %>% 
  dplyr::pull(data_set2)

remove_microbiome2 = 
  cor_data_preDM %>% 
  dplyr::group_by(data_set2) %>% 
  dplyr::summarise(n = sum(p_adjust < 0.2)) %>% 
  dplyr::filter(n == 0) %>% 
  dplyr::pull(data_set2)

remove_microbiome = 
  intersect(remove_microbiome1, remove_microbiome2)

length(remove_microbiome)

normal_cor =
  cor_data_Normal %>%
  dplyr::select(-c(p:p_adjust2)) %>%
  tidyr::pivot_wider(names_from = data_set2, values_from = "cor") %>%
  tibble::column_to_rownames(var = "data_set1")

normal_p =
  cor_data_Normal %>%
  dplyr::select(-c(p, p_adjust2, cor)) %>%
  tidyr::pivot_wider(names_from = data_set2, values_from = "p_adjust") %>%
  tibble::column_to_rownames(var = "data_set1")

predm_cor =
  cor_data_preDM %>%
  dplyr::select(-c(p:p_adjust2)) %>%
  tidyr::pivot_wider(names_from = data_set2, values_from = "cor") %>%
  tibble::column_to_rownames(var = "data_set1")

predm_p =
  cor_data_preDM %>%
  dplyr::select(-c(p, p_adjust2, cor)) %>%
  tidyr::pivot_wider(names_from = data_set2, values_from = "p_adjust") %>%
  tibble::column_to_rownames(var = "data_set1")

###remove the microbiomes which we need to remove
normal_cor =
  normal_cor %>%
  dplyr::select(-remove_microbiome)

normal_p =
  normal_p %>%
  dplyr::select(-remove_microbiome)

predm_cor =
  predm_cor %>%
  dplyr::select(-remove_microbiome)

predm_p =
  predm_p %>%
  dplyr::select(-remove_microbiome)

colnames(normal_cor)

colnames(predm_cor)

library(ComplexHeatmap)

colnames(predm_cor) = 
  colnames(normal_cor) = 
  microbiome_variable_info$variable_id[match(colnames(predm_cor), microbiome_variable_info$variable_id)]

library(circlize)

col_fun = circlize::colorRamp2(
  breaks = seq(-1, 1, length.out = 11),
  colors = rev(RColorBrewer::brewer.pal(n = 11, name = "BrBG"))
)

library(wesanderson)

plot = 
  Heatmap(
    t(predm_cor),
    col = col_fun,
    border = TRUE,
    show_column_dend = TRUE,
    show_row_dend = TRUE,
    clustering_distance_rows = "euclidean",
    clustering_method_rows = "complete",
    clustering_distance_columns = "euclidean",
    clustering_method_columns = "complete",
    name = "Spearman correlation",
    cluster_columns = TRUE,
    cluster_rows = TRUE,
    rect_gp = gpar(col = "white"),
    column_names_gp = gpar(cex = 0.5, rot = 45),
    row_names_gp = gpar(cex = 0.5),
    column_names_rot = 45,
    cell_fun = function(j, i, x, y, w, h, fill) {
      gb1 = textGrob("*")
      gb_w1 = convertWidth(grobWidth(gb1), "mm")
      gb_h1 = convertHeight(grobHeight(gb1), "mm")
      gb2 = textGrob(".")
      gb_w2 = convertWidth(grobWidth(gb2), "mm")
      gb_h2 = convertHeight(grobHeight(gb2), "mm")
      if (t(predm_p)[i, j] < 0.05) {
        grid.text("*", x, y - gb_h1 * 0.5 + gb_w1 * 0.4, gp = gpar(col = "red"))
      }
      
      if (t(predm_p)[i, j] > 0.05 & t(predm_p)[i, j] < 0.2) {
        grid.points(pch = 20, x = x, y = y, size = unit(0.3, "char"), gp = gpar(col = "red"))
        # grid.text(label = ".", x = x, y - gb_h2 * 0.5 + gb_w1 * 0.4)
      }
    }
  )

plot = ggplotify::as.ggplot(plot)
plot
# ggsave(plot, filename = "IR_cor.pdf", width = 10, height = 6)

plot = 
  Heatmap(
    t(normal_cor),
    col = col_fun,
    border = TRUE,
    show_column_dend = TRUE,
    show_row_dend = TRUE,
    rect_gp = gpar(col = "white"),
    clustering_distance_rows = "euclidean",
    clustering_method_rows = "complete",
    clustering_distance_columns = "euclidean",
    clustering_method_columns = "complete",
    name = "Spearman correlation",
    cluster_columns = TRUE,
    cluster_rows = TRUE,
    column_names_gp = gpar(cex = 0.5, rot = 45),
    row_names_gp = gpar(cex = 0.5),
    column_names_rot = 45,
    cell_fun = function(j, i, x, y, w, h, fill) {
      gb1 = textGrob("*")
      gb_w1 = convertWidth(grobWidth(gb1), "mm")
      gb_h1 = convertHeight(grobHeight(gb1), "mm")
      gb2 = textGrob(".")
      gb_w2 = convertWidth(grobWidth(gb2), "mm")
      gb_h2 = convertHeight(grobHeight(gb2), "mm")
      if (t(normal_p)[i, j] < 0.05) {
        grid.text("*", x, y - gb_h1 * 0.5 + gb_w1 * 0.4, gp = gpar(col = "red"))
      }
      
      if (t(normal_p)[i, j] > 0.05 & t(normal_p)[i, j] < 0.2) {
        grid.points(pch = 20, x = x, y = y, size = unit(0.3, "char"), gp = gpar(col = "red"))
        # grid.text(label = ".", x = x, y - gb_h2 * 0.5 + gb_w1 * 0.4)
      }
    }
  )

plot = ggplotify::as.ggplot(plot)
plot
# ggsave(plot, filename = "IS_cor.pdf", width = 10, height = 6)
