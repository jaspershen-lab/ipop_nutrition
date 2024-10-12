no_function()
setwd(r4projects::get_project_wd())
library(tidyverse)
library(plyr)
rm(list = ls())
source("1_code/tools.R")

###load data
{
  ##nutrition
  load("3_data_analysis/food_group/data_preparation/expression_data")
  load("3_data_analysis/food_group/data_preparation/sample_info")
  load("3_data_analysis/food_group/data_preparation/variable_info")
  
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
  
  ##microbiome
  load("3_data_analysis/gut_microbiome/data_preparation/expression_data")
  load("3_data_analysis/gut_microbiome/data_preparation/sample_info")
  load("3_data_analysis/gut_microbiome/data_preparation/variable_info")
  
  microbiome_expression_data = expression_data
  microbiome_sample_info = sample_info
  microbiome_variable_info = variable_info
}

dir.create("3_data_analysis/1_food_group_vs_microbiome_t1_t4/based_on_sspg/",
           recursive = TRUE)
setwd("3_data_analysis/1_food_group_vs_microbiome_t1_t4/based_on_sspg/")

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
    subject_id1 = x[4]
    date1 = as.Date(x[5])
    temp =
      microbiome_sample_info %>%
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
      microbiome_sample_info[match(matched_idx$sample_id2, microbiome_sample_info$sample_id),
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
microbiome_expression_data = microbiome_expression_data[, matched_idx$sample_id2]

nutrition_sample_info =
  nutrition_sample_info[match(colnames(nutrition_expression_data),
                              nutrition_sample_info$sample_id), ]

microbiome_sample_info =
  microbiome_sample_info[match(colnames(microbiome_expression_data),
                               microbiome_sample_info$sample_id), ]

dim(nutrition_expression_data)
dim(microbiome_expression_data)

colnames(microbiome_expression_data) =
  colnames(nutrition_expression_data)

######calculate the correlation between nutrition and metabolites
####missing value
nutrition_expression_data %>%
  apply(1, function(x) {
    sum(is.na(x))
  }) %>%
  plot()

library(impute)

remove_idx <-
  nutrition_expression_data %>%
  apply(1, sd) %>%
  `==`(0) %>%
  which()

if (length(remove_idx) > 0) {
  nutrition_expression_data <-
    nutrition_expression_data[-remove_idx,]
  nutrition_variable_info <-
    nutrition_variable_info[-remove_idx, , drop = FALSE]
}

nutrition_expression_data =
  impute::impute.knn(data = as.matrix(nutrition_expression_data))$data %>%
  as.data.frame()

###scale
nutrition_expression_data =
  nutrition_expression_data %>%
  apply(1, function(x) {
    (x) / sd(x)
  }) %>%
  t() %>%
  as.data.frame()

range(nutrition_expression_data)
range(microbiome_expression_data)

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

# cor_data =
#   partial_cor_pattern(
#     data_set1 = nutrition_expression_data,
#     data_set2 = microbiome_expression_data,
#     sample_info = nutrition_sample_info,
#     method = "spearman",
#     threads = 5
#   )
#
# cor_data_IS =
#   partial_cor_pattern(
#     data_set1 = nutrition_expression_data[, sample_info_IS$sample_id],
#     data_set2 = microbiome_expression_data[, sample_info_IS$sample_id],
#     sample_info = sample_info_IS,
#     method = "spearman",
#     threads = 5
#   )
#
# cor_data_IR =
#   partial_cor_pattern(
#     data_set1 = nutrition_expression_data[, sample_info_IR$sample_id],
#     data_set2 = microbiome_expression_data[, sample_info_IR$sample_id],
#     sample_info = sample_info_IR,
#     method = "spearman",
#     threads = 5
#   )
# save(cor_data, file = "cor_data")
# save(cor_data_IS, file = "cor_data_IS")
# save(cor_data_IR, file = "cor_data_IR")

load("cor_data")
load("cor_data_IS")
load("cor_data_IR")

head(cor_data_IR)
head(cor_data_IS)

plot(cor_data_IR$cor, cor_data_IS$cor)

dim(cor_data)
dim(cor_data_IS)
dim(cor_data_IR)

which(cor_data$p_adjust < 0.05)
which(cor_data_IR$p_adjust < 0.05)
which(cor_data_IS$p_adjust < 0.05)

# which(cor_data_IR$p_adjust < 0.2)
# which(cor_data_IS$p_adjust < 0.2)

####output results
library(openxlsx)
cor_data_output =
  cor_data %>%
  dplyr::filter(p_adjust < 0.05) %>%
  dplyr::left_join(microbiome_variable_info, by = c("data_set2" = "variable_id")) %>%
  dplyr::select(-p_adjust2) %>%
  dplyr::arrange(p_adjust)

cor_data_IS_output =
  cor_data_IS %>%
  dplyr::filter(p_adjust < 0.05) %>%
  dplyr::left_join(microbiome_variable_info, by = c("data_set2" = "variable_id")) %>%
  dplyr::select(-p_adjust2) %>%
  dplyr::arrange(p_adjust)

cor_data_IR_output =
  cor_data_IR %>%
  dplyr::filter(p_adjust < 0.05) %>%
  dplyr::left_join(microbiome_variable_info, by = c("data_set2" = "variable_id")) %>%
  dplyr::select(-p_adjust2) %>%
  dplyr::arrange(p_adjust)


# openxlsx::write.xlsx(
#   cor_data_output,
#   "cor_data_output.xlsx",
#   asTable = TRUE,
#   overwrite = TRUE
# )
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

data1 =
  cor_data_IS_output %>%
  dplyr::filter(p_adjust < 0.05) %>%
  dplyr::arrange(data_set1) %>%
  dplyr::mutate(class = "IS")

data2 =
  cor_data_IR_output %>%
  dplyr::filter(p_adjust < 0.05) %>%
  dplyr::arrange(data_set1) %>%
  dplyr::mutate(class = "IR")

data =
  rbind(data1, data2) %>%
  dplyr::arrange(data_set1)

openxlsx::write.xlsx(data,
                     file = "significant_cor_IR_IS.xlsx",
                     asTable = TRUE,
                     overwrite = TRUE)

#####output the cor plot
idx = which(cor_data_IR$p_adjust < 0.05)

# for (i in idx) {
#   cat(i, " ")
#   x = as.numeric(nutrition_expression_data[cor_data_IR$data_set1[i], ])
#   y = as.numeric(microbiome_expression_data[cor_data_IR$data_set2[i], ])
#   data.frame(x, y) %>%
#     ggplot(aes(x, y)) +
#     geom_point()
# }

cor_data_IR %>%
  dplyr::arrange(desc(abs(cor))) %>%
  head()

# plot(as.numeric(nutrition_expression_data["VitC", sample_info_IR$sample_id]),
#      as.numeric(microbiome_expression_data["nRPLC_441.3947_11.3", sample_info_IR$sample_id]))

# abline(0, -1)

cor_data_IS %>%
  dplyr::arrange(desc(abs(cor))) %>%
  head()

# plot(as.numeric(nutrition_expression_data["VitC", sample_info_IS$sample_id]),
#      as.numeric(microbiome_expression_data["nRPLC_441.3947_11.3", sample_info_IS$sample_id]))

# abline(0, -1)

cor(
  as.numeric(nutrition_expression_data["VitC", sample_info_IS$sample_id]),
  as.numeric(microbiome_expression_data["nRPLC_441.3947_11.3", sample_info_IS$sample_id]),
  method = "spearman"
)

# temp_sample_info =
#   sample_info_IS[, c("Sex", "Age")]

# temp_sample_info$Sex[temp_sample_info$Sex == 'F'] = 0
# temp_sample_info$Sex[temp_sample_info$Sex == 'M'] = 1
# temp_sample_info$Sex = as.numeric(temp_sample_info$Sex)

# ppcor::pcor.test(
#   x = as.numeric(nutrition_expression_data["VitE_a_Toco", sample_info_IS$sample_id]),
#   y = as.numeric(microbiome_expression_data["pHILIC_732.5525_5.1", sample_info_IS$sample_id]),
#   z = temp_sample_info,
#   method = "spearman"
# )

unique(cor_data_IS$data_set2)
unique(cor_data_IS$data_set1)

######all metabolites
###all participants
all_cor =
  cor_data %>%
  dplyr::select(-c(p:p_adjust2)) %>%
  tidyr::pivot_wider(names_from = data_set2, values_from = "cor") %>%
  tibble::column_to_rownames(var = "data_set1")

all_p =
  cor_data %>%
  dplyr::select(-c(p, p_adjust2, cor)) %>%
  tidyr::pivot_wider(names_from = data_set2, values_from = "p_adjust") %>%
  tibble::column_to_rownames(var = "data_set1")

###remove the diet all are zero
remove_idx <-
  apply(all_cor, 1, function(x) {
    sum(x == 0) / ncol(all_cor)
  }) %>%
  `>=`(0.9) %>%
  which()

if (length(remove_idx) > 0) {
  all_cor <- all_cor[-remove_idx, ]
  all_p <- all_p[-remove_idx, ]
}

####normal people, IS people
normal_cor =
  cor_data_IS %>%
  dplyr::select(-c(p:p_adjust2)) %>%
  tidyr::pivot_wider(names_from = data_set2, values_from = "cor") %>%
  tibble::column_to_rownames(var = "data_set1")

normal_p =
  cor_data_IS %>%
  dplyr::select(-c(p, p_adjust2, cor)) %>%
  tidyr::pivot_wider(names_from = data_set2, values_from = "p_adjust") %>%
  tibble::column_to_rownames(var = "data_set1")

###remove the diet all are zero
remove_idx <-
  apply(normal_cor, 1, function(x) {
    sum(x == 0) / ncol(normal_cor)
  }) %>%
  `>=`(0.9) %>%
  which()

if (length(remove_idx) > 0) {
  normal_cor <- normal_cor[-remove_idx, ]
  normal_p <- normal_p[-remove_idx, ]
}

####prediabetes, IR people
predm_cor =
  cor_data_IR %>%
  dplyr::select(-c(p:p_adjust2)) %>%
  tidyr::pivot_wider(names_from = data_set2, values_from = "cor") %>%
  tibble::column_to_rownames(var = "data_set1")

predm_p =
  cor_data_IR %>%
  dplyr::select(-c(p, p_adjust2, cor)) %>%
  tidyr::pivot_wider(names_from = data_set2, values_from = "p_adjust") %>%
  tibble::column_to_rownames(var = "data_set1")

###remove the diet all are zero
remove_idx <-
  apply(predm_cor, 1, function(x) {
    sum(x == 0) / ncol(predm_cor)
  }) %>%
  `>=`(0.9) %>%
  which()

if (length(remove_idx) > 0) {
  predm_cor <- predm_cor[-remove_idx, ]
  predm_p <- predm_p[-remove_idx, ]
}

colnames(all_cor) =
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
      
      # if (t(predm_p)[i, j] > 0.05 & t(predm_p)[i, j] < 0.2) {
      #   grid.points(pch = 20, x = x, y = y, size = unit(0.3, "char"), gp = gpar(col = "red"))
      # }
    }
  )

plot = ggplotify::as.ggplot(plot)
plot
# ggsave(plot,
#        filename = "IR_cor_all_metabolite.pdf",
#        width = 10,
#        height = 10)

plot <-
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
      
      # if (t(normal_p)[i, j] > 0.05 & t(normal_p)[i, j] < 0.2) {
      #   grid.points(pch = 20, x = x, y = y, size = unit(0.3, "char"), gp = gpar(col = "red"))
      #   # grid.text(label = ".", x = x, y - gb_h2 * 0.5 + gb_w1 * 0.4)
      # }
    }
  )

plot = ggplotify::as.ggplot(plot)
plot
# ggsave(plot,
#        filename = "IS_cor_all_metabolite.pdf",
#        width = 10,
#        height = 10)



plot <-
  Heatmap(
    t(all_cor),
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
      if (t(all_p)[i, j] < 0.05) {
        grid.text("*", x, y - gb_h1 * 0.5 + gb_w1 * 0.4, gp = gpar(col = "red"))
      }
      
      # if (t(all_p)[i, j] > 0.05 & t(all_p)[i, j] < 0.2) {
      #   grid.points(pch = 20, x = x, y = y, size = unit(0.3, "char"), gp = gpar(col = "red"))
      #   # grid.text(label = ".", x = x, y - gb_h2 * 0.5 + gb_w1 * 0.4)
      # }
    }
  )

plot = ggplotify::as.ggplot(plot)
plot
# ggsave(plot,
#        filename = "all_cor_all_metabolite.pdf",
#        width = 10,
#        height = 10)

######only show the metabolites with at least one significant correlation with nutrition
###remove the metabolites which have no significant cor (p_adjust > 0.2) with nutrition
remove_metabolite1 =
  cor_data_IS %>%
  dplyr::group_by(data_set2) %>%
  dplyr::summarise(n = sum(p_adjust < 0.05)) %>%
  dplyr::filter(n == 0) %>%
  dplyr::pull(data_set2)

remove_metabolite2 =
  cor_data_IR %>%
  dplyr::group_by(data_set2) %>%
  dplyr::summarise(n = sum(p_adjust < 0.05)) %>%
  dplyr::filter(n == 0) %>%
  dplyr::pull(data_set2)

remove_metabolite3 =
  cor_data %>%
  dplyr::group_by(data_set2) %>%
  dplyr::summarise(n = sum(p_adjust < 0.05)) %>%
  dplyr::filter(n == 0) %>%
  dplyr::pull(data_set2)

remove_metabolite =
  intersect(intersect(remove_metabolite1, remove_metabolite2),
            remove_metabolite3)

length(remove_metabolite)


####all participants, IS people
all_cor =
  cor_data %>%
  dplyr::select(-c(p:p_adjust2)) %>%
  tidyr::pivot_wider(names_from = data_set2, values_from = "cor") %>%
  tibble::column_to_rownames(var = "data_set1")

all_p =
  cor_data %>%
  dplyr::select(-c(p, p_adjust2, cor)) %>%
  tidyr::pivot_wider(names_from = data_set2, values_from = "p_adjust") %>%
  tibble::column_to_rownames(var = "data_set1")

###remove the diet all are zero
remove_idx <-
  apply(all_cor, 1, function(x) {
    sum(x == 0) / ncol(all_cor)
  }) %>%
  `>=`(0.9) %>%
  which()

if (length(remove_idx) > 0) {
  all_cor <- all_cor[-remove_idx, ]
  all_p <- all_p[-remove_idx, ]
}

####normal participants, IS people
normal_cor =
  cor_data_IS %>%
  dplyr::select(-c(p:p_adjust2)) %>%
  tidyr::pivot_wider(names_from = data_set2, values_from = "cor") %>%
  tibble::column_to_rownames(var = "data_set1")

normal_p =
  cor_data_IS %>%
  dplyr::select(-c(p, p_adjust2, cor)) %>%
  tidyr::pivot_wider(names_from = data_set2, values_from = "p_adjust") %>%
  tibble::column_to_rownames(var = "data_set1")

###remove the diet all are zero
remove_idx <-
  apply(normal_cor, 1, function(x) {
    sum(x == 0) / ncol(normal_cor)
  }) %>%
  `>=`(0.9) %>%
  which()

if (length(remove_idx) > 0) {
  normal_cor <- normal_cor[-remove_idx, ]
  normal_p <- normal_p[-remove_idx, ]
}


predm_cor =
  cor_data_IR %>%
  dplyr::select(-c(p:p_adjust2)) %>%
  tidyr::pivot_wider(names_from = data_set2, values_from = "cor") %>%
  tibble::column_to_rownames(var = "data_set1")

predm_p =
  cor_data_IR %>%
  dplyr::select(-c(p, p_adjust2, cor)) %>%
  tidyr::pivot_wider(names_from = data_set2, values_from = "p_adjust") %>%
  tibble::column_to_rownames(var = "data_set1")



###remove the diet all are zero
remove_idx <-
  apply(predm_cor, 1, function(x) {
    sum(x == 0) / ncol(predm_cor)
  }) %>%
  `>=`(0.9) %>%
  which()

if (length(remove_idx) > 0) {
  predm_cor <- predm_cor[-remove_idx, ]
  predm_p <- predm_p[-remove_idx, ]
}


###remove the metabolites which we need to remove
all_cor =
  all_cor %>%
  dplyr::select(-remove_metabolite)

all_p =
  all_p %>%
  dplyr::select(-remove_metabolite)

normal_cor =
  normal_cor %>%
  dplyr::select(-remove_metabolite)

normal_p =
  normal_p %>%
  dplyr::select(-remove_metabolite)

predm_cor =
  predm_cor %>%
  dplyr::select(-remove_metabolite)

predm_p =
  predm_p %>%
  dplyr::select(-remove_metabolite)

colnames(all_cor)

colnames(normal_cor)

colnames(predm_cor)

library(ComplexHeatmap)

colnames(all_cor) =
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
    t(all_cor),
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
    row_names_gp = gpar(cex = 0.5),
    column_names_rot = 45,
    cell_fun = function(j, i, x, y, w, h, fill) {
      gb1 = textGrob("*")
      gb_w1 = convertWidth(grobWidth(gb1), "mm")
      gb_h1 = convertHeight(grobHeight(gb1), "mm")
      gb2 = textGrob(".")
      gb_w2 = convertWidth(grobWidth(gb2), "mm")
      gb_h2 = convertHeight(grobHeight(gb2), "mm")
      if (t(all_p)[i, j] < 0.05) {
        grid.text("*", x, y - gb_h1 * 0.5 + gb_w1 * 0.4, gp = gpar(col = "red"))
      }
      
      if (t(all_p)[i, j] > 0.05 & t(all_p)[i, j] < 0.2) {
        grid.points(
          pch = 20,
          x = x,
          y = y,
          size = unit(0.3, "char"),
          gp = gpar(col = "red")
        )
        # grid.text(label = ".", x = x, y - gb_h2 * 0.5 + gb_w1 * 0.4)
      }
    }
  )

plot = ggplotify::as.ggplot(plot)
plot
# ggsave(plot,
#        filename = "all_cor.pdf",
#        width = 10,
#        height = 10)


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
        grid.points(
          pch = 20,
          x = x,
          y = y,
          size = unit(0.3, "char"),
          gp = gpar(col = "red")
        )
        # grid.text(label = ".", x = x, y - gb_h2 * 0.5 + gb_w1 * 0.4)
      }
    }
  )

plot = ggplotify::as.ggplot(plot)
plot
# ggsave(plot,
#        filename = "IR_cor.pdf",
#        width = 10,
#        height = 10)

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
        grid.points(
          pch = 20,
          x = x,
          y = y,
          size = unit(0.3, "char"),
          gp = gpar(col = "red")
        )
        # grid.text(label = ".", x = x, y - gb_h2 * 0.5 + gb_w1 * 0.4)
      }
    }
  )

plot = ggplotify::as.ggplot(plot)
plot
# ggsave(plot,
#        filename = "IS_cor.pdf",
#        width = 10,
#        height = 10)

cor_data_IR


temp_cor_is <-
  normal_cor  %>%
  tibble::rownames_to_column(var = "food_group") %>%
  tidyr::pivot_longer(cols = -food_group,
                      names_to = "microbiome",
                      values_to = "cor")

temp_p_is <-
  normal_p  %>%
  tibble::rownames_to_column(var = "food_group") %>%
  tidyr::pivot_longer(cols = -food_group,
                      names_to = "microbiome",
                      values_to = "p")


temp_is <-
  temp_cor_is  %>%
  dplyr::left_join(temp_p_is, by = c("food_group", "microbiome"))


temp_cor_ir <-
  predm_cor  %>%
  tibble::rownames_to_column(var = "food_group") %>%
  tidyr::pivot_longer(cols = -food_group,
                      names_to = "microbiome",
                      values_to = "cor")

temp_p_ir <-
  predm_cor  %>%
  tibble::rownames_to_column(var = "food_group") %>%
  tidyr::pivot_longer(cols = -food_group,
                      names_to = "microbiome",
                      values_to = "p")

temp_ir <-
  temp_cor_ir  %>%
  dplyr::left_join(temp_p_ir, by = c("food_group", "microbiome"))

dim(temp_is)
dim(temp_ir)

###use histgram to show the distribution of the cor
library(ggplot2)
library(ggpubr)
temp_is  %>%
  ggplot(aes(cor)) +
  geom_histogram(binwidth = 0.1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


###use the density plot to show the distribution of the cor

temp <-
  rbind(data.frame(temp_is, class = "IS"),
        data.frame(temp_ir, class = "IR"))

###the density plot below use the color to fill the density plot
plot <-
  temp %>%
  ggplot(aes(cor)) +
  geom_density(aes(color = class)) +
  theme_bw() +
  scale_color_manual(values = ir_is_color) +
  labs(x = "Correlation", y = "Density")

plot

ggsave(plot,
       filename = "cor_density_comparison.pdf",
       width = 8,
       height = 6)
