# no_function()
# setwd(r4projects::get_project_wd())
# library(tidyverse)
# library(plyr)
# rm(list = ls())
# source("1_code/tools.R")
# 
# ###load data
# {
#   ##nutrition
#   load("3_data_analysis/food_group/data_preparation/expression_data")
#   load("3_data_analysis/food_group/data_preparation/sample_info")
#   load("3_data_analysis/food_group/data_preparation/variable_info")
#   
#   nutrition_expression_data = expression_data
#   nutrition_sample_info = sample_info
#   nutrition_variable_info = variable_info
#   
#   load("3_data_analysis/nutrition/data_preparation/sample_info")
#   
#   nutrition_sample_info = 
#   nutrition_sample_info %>%
#     dplyr::left_join(sample_info[, c("subject_id", "a1c_t2d_status")], by = "subject_id")
#   
#   filter_subject_id = 
#   nutrition_sample_info %>% 
#     dplyr::group_by(subject_id) %>% 
#     dplyr::summarise(n = n()) %>% 
#     dplyr::filter(n == 4) %>% 
#     pull(subject_id)
#   
#   nutrition_sample_info =
#     nutrition_sample_info %>% 
#     dplyr::filter(subject_id %in% filter_subject_id)
#   
#   nutrition_expression_data = 
#     nutrition_expression_data[,nutrition_sample_info$sample_id]
#   
#   ##microbiome
#   load("3_data_analysis/gut_microbiome/data_preparation/expression_data")
#   load("3_data_analysis/gut_microbiome/data_preparation/sample_info")
#   load("3_data_analysis/gut_microbiome/data_preparation/variable_info")
#   
#   microbiome_expression_data = expression_data
#   microbiome_sample_info = sample_info
#   microbiome_variable_info = variable_info
# }
# 
# dir.create("3_data_analysis/food_group_vs_microbiome_t1_t4/based_on_a1c/")
# setwd("3_data_analysis/food_group_vs_microbiome_t1_t4/based_on_a1c/")
# 
# ###data preparation
# nutrition_sample_info$Diet.Survey.Date
# nutrition_sample_info$CollectionDate
# 
# ####match microbiome data from baseline
# microbiome_sample_info =
#   microbiome_sample_info %>%
#   dplyr::mutate(CollectionDate = as.Date(CollectionDate, format = "%m/%d/%y"))
# 
# matched_sample =
#   nutrition_sample_info %>%
#   t() %>%
#   as.data.frame() %>%
#   purrr::map(function(x) {
#     subject_id1 = x[4]
#     date1 = as.Date(x[5])
#     temp =
#       microbiome_sample_info %>%
#       dplyr::filter(subject_id == subject_id1)
#     temp %>%
#       dplyr::mutate(diff_days = as.numeric(CollectionDate - date1)) %>%
#       dplyr::filter(abs(diff_days) < 7) %>%
#       dplyr::filter(abs(diff_days) == min(abs(diff_days))) %>% 
#       head(1)
#   })
# 
# matched_sample
# 
# matched_sample %>% 
#   lapply(nrow) %>% 
#   unlist() %>% 
#   plot
# 
# matched_sample = 
#   matched_sample %>% 
#   lapply(function(x){
#     temp = x$sample_id
#     if(length(temp) == 0){
#       temp = NA
#     }
#     temp
#   }) %>% 
#   unlist()
# 
# matched_idx = 
#   data.frame(sample_id1 = nutrition_sample_info$sample_id,
#              sample_id2 = unname(matched_sample)) %>% 
#   dplyr::filter(!is.na(sample_id2))
# 
# 
# matched_idx
# 
# cbind(
#   nutrition_sample_info[match(matched_idx$sample_id1, nutrition_sample_info$sample_id),
#                         c("sample_id", "Diet.Survey.Date")],
#   microbiome_sample_info[match(matched_idx$sample_id2, microbiome_sample_info$sample_id),
#                            c("sample_id", "CollectionDate")]  
# )
# 
# matched_idx = 
# matched_idx %>% 
#   dplyr::left_join(nutrition_sample_info[,c("subject_id", "sample_id")],
#                    by = c("sample_id1" = "sample_id"))
#   
# temp_subject_id = 
# matched_idx %>% 
# dplyr::group_by(subject_id) %>% 
#   dplyr::summarise(n = n()) %>% 
#   dplyr::filter(n == 4) %>% 
#   dplyr::pull(subject_id)
#   
# matched_idx = 
#   matched_idx %>% 
#   dplyr::filter(subject_id %in% temp_subject_id) %>% 
#   dplyr::select(-subject_id)
# 
# nutrition_sample_info = 
#   nutrition_sample_info %>% 
#   dplyr::filter(sample_id %in% matched_idx$sample_id1)
# 
# nutrition_expression_data = 
#   nutrition_expression_data[,nutrition_sample_info$sample_id]
# 
# nutrition_expression_data = nutrition_expression_data[,matched_idx$sample_id1]
# microbiome_expression_data = microbiome_expression_data[,matched_idx$sample_id2]
# 
# nutrition_sample_info = 
#   nutrition_sample_info[match(colnames(nutrition_expression_data),
#                               nutrition_sample_info$sample_id),]
# 
# microbiome_sample_info = 
#   microbiome_sample_info[match(colnames(microbiome_expression_data), microbiome_sample_info$sample_id),]
# 
# dim(nutrition_expression_data)
# dim(microbiome_expression_data)
# 
# colnames(microbiome_expression_data) =
#   colnames(nutrition_expression_data)
# 
# ######calculate the correlation between nutrition and metabolites
# ####missing value
# nutrition_expression_data %>% 
#   apply(1, function(x){
#     sum(is.na(x))
#   }) %>% 
#   plot()
# 
# ###remove some food groups
# remove_idx <-
#   nutrition_expression_data %>% 
#   apply(1, sd) %>% 
#   `==`(0) %>% 
#   which()
# 
# if(length(remove_idx) > 0){
#   nutrition_expression_data <-
#     nutrition_expression_data[-remove_idx,]
#   nutrition_variable_info <-
#     nutrition_variable_info[-remove_idx,,drop = FALSE]
# }
# 
# library(impute)
# 
# nutrition_expression_data =
#   impute::impute.knn(data = as.matrix(nutrition_expression_data))$data %>%
#   as.data.frame()
# 
# ###scale
# nutrition_expression_data = 
#   nutrition_expression_data %>% 
#   apply(1, function(x){
#     (x) / sd(x)
#   }) %>% 
#   t() %>% 
#   as.data.frame()
# 
# range(nutrition_expression_data)
# range(microbiome_expression_data)
# 
# ###############################################################################
# ####based on the ogtt status
# ####calculate the correlation between them
# ###cor_data for IS people
# sample_info_normal =
#   nutrition_sample_info %>%
#   dplyr::filter(!is.na(a1c_t2d_status)) %>%
#   dplyr::filter(a1c_t2d_status == "Normal")
# 
# sample_info_predm =
#   nutrition_sample_info %>%
#   dplyr::filter(!is.na(a1c_t2d_status)) %>%
#   dplyr::filter(a1c_t2d_status == "preDM" | a1c_t2d_status == "T2D")
# 
# cor_data_IS =
#   partial_cor_pattern(
#     data_set1 = nutrition_expression_data[, sample_info_normal$sample_id],
#     data_set2 = microbiome_expression_data[, sample_info_normal$sample_id],
#     sample_info = sample_info_normal,
#     method = "spearman",
#     threads = 5
#   )
# 
# cor_data_IR =
#   partial_cor_pattern(
#     data_set1 = nutrition_expression_data[, sample_info_predm$sample_id],
#     data_set2 = microbiome_expression_data[, sample_info_predm$sample_id],
#     sample_info = sample_info_predm,
#     method = "spearman",
#     threads = 5
#   )
# 
# save(cor_data_IS, file = "cor_data_IS")
# save(cor_data_IR, file = "cor_data_IR")
# 
# load("cor_data_IS")
# load("cor_data_IR")
# 
# dim(cor_data_IS)
# dim(cor_data_IR)
# 
# which(cor_data_IR$p_adjust < 0.05)
# which(cor_data_IS$p_adjust < 0.05)
# 
# # which(cor_data_IR$p_adjust < 0.2)
# # which(cor_data_IS$p_adjust < 0.2)
# 
# ####output results
# library(openxlsx)
# cor_data_IS_output = 
#   cor_data_IS %>% 
#   dplyr::filter(p_adjust < 0.05) %>% 
#   dplyr::left_join(microbiome_variable_info, by = c("data_set2" = "variable_id")) %>% 
#   dplyr::select(-p_adjust2) %>% 
#   dplyr::arrange(p_adjust)
# 
# cor_data_IR_output = 
#   cor_data_IR %>% 
#   dplyr::filter(p_adjust < 0.05) %>% 
#   dplyr::left_join(microbiome_variable_info, by = c("data_set2" = "variable_id")) %>% 
#   dplyr::select(-p_adjust2) %>% 
#   dplyr::arrange(p_adjust)
# 
# openxlsx::write.xlsx(
#   cor_data_IS_output,
#   "cor_data_IS_output.xlsx",
#   asTable = TRUE,
#   overwrite = TRUE
# )
# 
# openxlsx::write.xlsx(
#   cor_data_IR_output,
#   "cor_data_IR_output.xlsx",
#   asTable = TRUE,
#   overwrite = TRUE
# )
# 
# #####output the cor plot
# idx = which(cor_data_IR$p_adjust < 0.05)
# 
# # for(i in idx) {
# # cat(i, " ")
# #   x = as.numeric(nutrition_expression_data[cor_data_IR$data_set1[i],])
# #   y = as.numeric(microbiome_expression_data[cor_data_IR$data_set2[i],])
# #   data.frame(x, y) %>% 
# #     ggplot(aes(x,y)) +
# #     geom_point()
# # }
# 
# cor_data_IR %>%
#   dplyr::arrange(desc(abs(cor))) %>% 
#   head()
# 
# plot(as.numeric(nutrition_expression_data["VitC", sample_info_predm$sample_id]),
#      as.numeric(microbiome_expression_data["nRPLC_441.3947_11.3", sample_info_predm$sample_id]))
# 
# abline(0,-1)
# 
# cor_data_IS %>%
#   dplyr::arrange(desc(abs(cor))) %>% 
#   head()
# 
# plot(as.numeric(nutrition_expression_data["VitC", sample_info_normal$sample_id]),
#      as.numeric(microbiome_expression_data["nRPLC_441.3947_11.3", sample_info_normal$sample_id]))
# 
# abline(0,-1)
# 
# cor(as.numeric(nutrition_expression_data["VitC", sample_info_normal$sample_id]),
#     as.numeric(microbiome_expression_data["nRPLC_441.3947_11.3", sample_info_normal$sample_id]),
#     method = "spearman")
# 
# temp_sample_info = 
#   sample_info_normal[,c("Sex", "Age")]
# 
# temp_sample_info$Sex[temp_sample_info$Sex == 'F'] = 0
# temp_sample_info$Sex[temp_sample_info$Sex == 'M'] = 1
# temp_sample_info$Sex = as.numeric(temp_sample_info$Sex)
# 
# ppcor::pcor.test(
#   x = as.numeric(nutrition_expression_data["VitE_a_Toco", sample_info_normal$sample_id]),
#   y = as.numeric(microbiome_expression_data["pHILIC_732.5525_5.1", sample_info_normal$sample_id]),
#   z = temp_sample_info,
#   method = "spearman"
# )
# 
# unique(cor_data_IS$data_set2)
# unique(cor_data_IS$data_set1)
# 
# ######all metabolites
# normal_cor =
#   cor_data_IS %>%
#   dplyr::select(-c(p:p_adjust2)) %>%
#   tidyr::pivot_wider(names_from = data_set2, values_from = "cor") %>%
#   tibble::column_to_rownames(var = "data_set1")
# 
# normal_p =
#   cor_data_IS %>%
#   dplyr::select(-c(p, p_adjust2, cor)) %>%
#   tidyr::pivot_wider(names_from = data_set2, values_from = "p_adjust") %>%
#   tibble::column_to_rownames(var = "data_set1")
# 
# predm_cor =
#   cor_data_IR %>%
#   dplyr::select(-c(p:p_adjust2)) %>%
#   tidyr::pivot_wider(names_from = data_set2, values_from = "cor") %>%
#   tibble::column_to_rownames(var = "data_set1")
# 
# predm_p =
#   cor_data_IR %>%
#   dplyr::select(-c(p, p_adjust2, cor)) %>%
#   tidyr::pivot_wider(names_from = data_set2, values_from = "p_adjust") %>%
#   tibble::column_to_rownames(var = "data_set1")
# 
# colnames(predm_cor) = 
#   colnames(normal_cor) = 
#   microbiome_variable_info$Metabolite[match(colnames(predm_cor), microbiome_variable_info$variable_id)]
# 
# library(ComplexHeatmap)
# library(circlize)
# 
# col_fun = circlize::colorRamp2(
#   breaks = seq(-1, 1, length.out = 11),
#   colors = rev(RColorBrewer::brewer.pal(n = 11, name = "BrBG"))
# )
# 
# library(wesanderson)
# 
# plot = 
#   Heatmap(
#     t(predm_cor),
#     col = col_fun,
#     border = TRUE,
#     show_column_dend = TRUE,
#     show_row_dend = TRUE,
#     clustering_distance_rows = "euclidean",
#     clustering_method_rows = "complete",
#     clustering_distance_columns = "euclidean",
#     clustering_method_columns = "complete",
#     name = "Spearman correlation",
#     cluster_columns = TRUE,
#     cluster_rows = TRUE,
#     column_names_gp = gpar(cex = 0.5, rot = 45),
#     row_names_gp = gpar(cex = 0.3),
#     column_names_rot = 45,
#     cell_fun = function(j, i, x, y, w, h, fill) {
#       gb1 = textGrob("*")
#       gb_w1 = convertWidth(grobWidth(gb1), "mm")
#       gb_h1 = convertHeight(grobHeight(gb1), "mm")
#       gb2 = textGrob(".")
#       gb_w2 = convertWidth(grobWidth(gb2), "mm")
#       gb_h2 = convertHeight(grobHeight(gb2), "mm")
#       if (t(predm_p)[i, j] < 0.05) {
#         grid.text("*", x, y - gb_h1 * 0.5 + gb_w1 * 0.4, gp = gpar(col = "red"))
#       }
#       
#       # if (t(predm_p)[i, j] > 0.05 & t(predm_p)[i, j] < 0.2) {
#       #   grid.points(pch = 20, x = x, y = y, size = unit(0.3, "char"), gp = gpar(col = "red"))
#       # }
#     }
#   )
# 
# plot = ggplotify::as.ggplot(plot)
# plot
# ggsave(plot, filename = "IR_cor_all_metabolite.pdf", width = 10, height = 10)
# 
# plot = 
#   Heatmap(
#     t(normal_cor),
#     col = col_fun,
#     border = TRUE,
#     show_column_dend = TRUE,
#     show_row_dend = TRUE,
#     clustering_distance_rows = "euclidean",
#     clustering_method_rows = "complete",
#     clustering_distance_columns = "euclidean",
#     clustering_method_columns = "complete",
#     name = "Spearman correlation",
#     cluster_columns = TRUE,
#     cluster_rows = TRUE,
#     column_names_gp = gpar(cex = 0.5, rot = 45),
#     row_names_gp = gpar(cex = 0.3),
#     column_names_rot = 45,
#     cell_fun = function(j, i, x, y, w, h, fill) {
#       gb1 = textGrob("*")
#       gb_w1 = convertWidth(grobWidth(gb1), "mm")
#       gb_h1 = convertHeight(grobHeight(gb1), "mm")
#       gb2 = textGrob(".")
#       gb_w2 = convertWidth(grobWidth(gb2), "mm")
#       gb_h2 = convertHeight(grobHeight(gb2), "mm")
#       if (t(normal_p)[i, j] < 0.05) {
#         grid.text("*", x, y - gb_h1 * 0.5 + gb_w1 * 0.4, gp = gpar(col = "red"))
#       }
#       
#       # if (t(normal_p)[i, j] > 0.05 & t(normal_p)[i, j] < 0.2) {
#       #   grid.points(pch = 20, x = x, y = y, size = unit(0.3, "char"), gp = gpar(col = "red"))
#       #   # grid.text(label = ".", x = x, y - gb_h2 * 0.5 + gb_w1 * 0.4)
#       # }
#     }
#   )
# 
# plot = ggplotify::as.ggplot(plot)
# plot
# ggsave(plot, filename = "IS_cor_all_metabolite.pdf", width = 10, height = 10)
# 
# ######only show the metabolites with at least one significant correlation with nutrition
# ###remove the metabolites which have no significant cor (p_adjust > 0.2) with nutrition
# remove_metabolite1 = 
#   cor_data_IS %>% 
#   dplyr::group_by(data_set2) %>% 
#   dplyr::summarise(n = sum(p_adjust < 0.05)) %>% 
#   dplyr::filter(n == 0) %>% 
#   dplyr::pull(data_set2)
# 
# remove_metabolite2 = 
#   cor_data_IR %>% 
#   dplyr::group_by(data_set2) %>% 
#   dplyr::summarise(n = sum(p_adjust < 0.05)) %>% 
#   dplyr::filter(n == 0) %>% 
#   dplyr::pull(data_set2)
# 
# remove_metabolite = 
#   intersect(remove_metabolite1, remove_metabolite2)
# 
# length(remove_metabolite)
# 
# normal_cor =
#   cor_data_IS %>%
#   dplyr::select(-c(p:p_adjust2)) %>%
#   tidyr::pivot_wider(names_from = data_set2, values_from = "cor") %>%
#   tibble::column_to_rownames(var = "data_set1")
# 
# normal_p =
#   cor_data_IS %>%
#   dplyr::select(-c(p, p_adjust2, cor)) %>%
#   tidyr::pivot_wider(names_from = data_set2, values_from = "p_adjust") %>%
#   tibble::column_to_rownames(var = "data_set1")
# 
# predm_cor =
#   cor_data_IR %>%
#   dplyr::select(-c(p:p_adjust2)) %>%
#   tidyr::pivot_wider(names_from = data_set2, values_from = "cor") %>%
#   tibble::column_to_rownames(var = "data_set1")
# 
# predm_p =
#   cor_data_IR %>%
#   dplyr::select(-c(p, p_adjust2, cor)) %>%
#   tidyr::pivot_wider(names_from = data_set2, values_from = "p_adjust") %>%
#   tibble::column_to_rownames(var = "data_set1")
# 
# ###remove the metabolites which we need to remove
# normal_cor =
#   normal_cor %>%
#   dplyr::select(-remove_metabolite)
# 
# normal_p =
#   normal_p %>%
#   dplyr::select(-remove_metabolite)
# 
# predm_cor =
#   predm_cor %>%
#   dplyr::select(-remove_metabolite)
# 
# predm_p =
#   predm_p %>%
#   dplyr::select(-remove_metabolite)
# 
# colnames(normal_cor)
# 
# colnames(predm_cor)
# 
# library(ComplexHeatmap)
# 
# colnames(predm_cor) = 
#   colnames(normal_cor) = 
#   microbiome_variable_info$Metabolite[match(colnames(predm_cor), microbiome_variable_info$variable_id)]
# 
# library(circlize)
# 
# col_fun = circlize::colorRamp2(
#   breaks = seq(-1, 1, length.out = 11),
#   colors = rev(RColorBrewer::brewer.pal(n = 11, name = "BrBG"))
# )
# 
# library(wesanderson)
# 
# plot = 
#   Heatmap(
#     t(predm_cor),
#     col = col_fun,
#     border = TRUE,
#     show_column_dend = TRUE,
#     show_row_dend = TRUE,
#     clustering_distance_rows = "euclidean",
#     clustering_method_rows = "complete",
#     clustering_distance_columns = "euclidean",
#     clustering_method_columns = "complete",
#     name = "Spearman correlation",
#     cluster_columns = TRUE,
#     cluster_rows = TRUE,
#     column_names_gp = gpar(cex = 0.5, rot = 45),
#     row_names_gp = gpar(cex = 0.5),
#     column_names_rot = 45,
#     cell_fun = function(j, i, x, y, w, h, fill) {
#       gb1 = textGrob("*")
#       gb_w1 = convertWidth(grobWidth(gb1), "mm")
#       gb_h1 = convertHeight(grobHeight(gb1), "mm")
#       gb2 = textGrob(".")
#       gb_w2 = convertWidth(grobWidth(gb2), "mm")
#       gb_h2 = convertHeight(grobHeight(gb2), "mm")
#       if (t(predm_p)[i, j] < 0.05) {
#         grid.text("*", x, y - gb_h1 * 0.5 + gb_w1 * 0.4, gp = gpar(col = "red"))
#       }
#       
#       if (t(predm_p)[i, j] > 0.05 & t(predm_p)[i, j] < 0.2) {
#         grid.points(pch = 20, x = x, y = y, size = unit(0.3, "char"), gp = gpar(col = "red"))
#         # grid.text(label = ".", x = x, y - gb_h2 * 0.5 + gb_w1 * 0.4)
#       }
#     }
#   )
# 
# plot = ggplotify::as.ggplot(plot)
# plot
# ggsave(plot, filename = "IR_cor.pdf", width = 10, height = 10)
# 
# plot = 
#   Heatmap(
#     t(normal_cor),
#     col = col_fun,
#     border = TRUE,
#     show_column_dend = TRUE,
#     show_row_dend = TRUE,
#     clustering_distance_rows = "euclidean",
#     clustering_method_rows = "complete",
#     clustering_distance_columns = "euclidean",
#     clustering_method_columns = "complete",
#     name = "Spearman correlation",
#     cluster_columns = TRUE,
#     cluster_rows = TRUE,
#     column_names_gp = gpar(cex = 0.5, rot = 45),
#     row_names_gp = gpar(cex = 0.5),
#     column_names_rot = 45,
#     cell_fun = function(j, i, x, y, w, h, fill) {
#       gb1 = textGrob("*")
#       gb_w1 = convertWidth(grobWidth(gb1), "mm")
#       gb_h1 = convertHeight(grobHeight(gb1), "mm")
#       gb2 = textGrob(".")
#       gb_w2 = convertWidth(grobWidth(gb2), "mm")
#       gb_h2 = convertHeight(grobHeight(gb2), "mm")
#       if (t(normal_p)[i, j] < 0.05) {
#         grid.text("*", x, y - gb_h1 * 0.5 + gb_w1 * 0.4, gp = gpar(col = "red"))
#       }
#       
#       if (t(normal_p)[i, j] > 0.05 & t(normal_p)[i, j] < 0.2) {
#         grid.points(pch = 20, x = x, y = y, size = unit(0.3, "char"), gp = gpar(col = "red"))
#         # grid.text(label = ".", x = x, y - gb_h2 * 0.5 + gb_w1 * 0.4)
#       }
#     }
#   )
# 
# plot = ggplotify::as.ggplot(plot)
# plot
# ggsave(plot, filename = "IS_cor.pdf", width = 10, height = 10)
