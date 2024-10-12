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

df_a <-
  readr::read_csv("2-data/DF_A_Energy_SampleID (1)_with_DC.csv")
df_a <-
  df_a %>%
  dplyr::select(SubjectID, Participant.ID, cluster_801) %>%
  dplyr::distinct(SubjectID, .keep_all = TRUE)


dir.create(
  "3_data_analysis/1_food_group_vs_microbiome_t1_t4/based_on_cluster801/",
  recursive = TRUE
)
setwd("3_data_analysis/1_food_group_vs_microbiome_t1_t4/based_on_cluster801/")

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
####based on the cluster801 status
####calculate the correlation between them
###cor_data for cluster0 people
sample_info_cluster0 =
  nutrition_sample_info %>%
  dplyr::left_join(df_a, by = c("subject_id" = "SubjectID")) %>%
  dplyr::filter(!is.na(cluster_801)) %>%
  dplyr::filter(cluster_801 == "0")

sample_info_cluster1 =
  nutrition_sample_info %>%
  dplyr::left_join(df_a, by = c("subject_id" = "SubjectID")) %>%
  dplyr::filter(!is.na(cluster_801)) %>%
  dplyr::filter(cluster_801 == "1")

# cor_data =
#   partial_cor_pattern(
#     data_set1 = nutrition_expression_data,
#     data_set2 = microbiome_expression_data,
#     sample_info = nutrition_sample_info,
#     method = "spearman",
#     threads = 5
#   )
#
# cor_data_cluster0 =
#   partial_cor_pattern(
#     data_set1 = nutrition_expression_data[, sample_info_cluster0$sample_id],
#     data_set2 = microbiome_expression_data[, sample_info_cluster0$sample_id],
#     sample_info = sample_info_cluster0,
#     method = "spearman",
#     threads = 5
#   )
#
# cor_data_cluster1 =
#   partial_cor_pattern(
#     data_set1 = nutrition_expression_data[, sample_info_cluster1$sample_id],
#     data_set2 = microbiome_expression_data[, sample_info_cluster1$sample_id],
#     sample_info = sample_info_cluster1,
#     method = "spearman",
#     threads = 5
#   )
# save(cor_data, file = "cor_data")
# save(cor_data_cluster0, file = "cor_data_cluster0")
# save(cor_data_cluster1, file = "cor_data_cluster1")

load("cor_data")
load("cor_data_cluster0")
load("cor_data_cluster1")

head(cor_data_cluster1)
head(cor_data_cluster0)

plot(cor_data_cluster1$cor, cor_data_cluster0$cor)

dim(cor_data)
dim(cor_data_cluster0)
dim(cor_data_cluster1)

which(cor_data$p_adjust < 0.05)
which(cor_data_cluster1$p_adjust < 0.05)
which(cor_data_cluster0$p_adjust < 0.05)

# which(cor_data_cluster1$p_adjust < 0.2)
# which(cor_data_cluster0$p_adjust < 0.2)

####output results
library(openxlsx)
cor_data_output =
  cor_data %>%
  dplyr::filter(p_adjust < 0.05) %>%
  dplyr::left_join(microbiome_variable_info, by = c("data_set2" = "variable_id")) %>%
  dplyr::select(-p_adjust2) %>%
  dplyr::arrange(p_adjust)

cor_data_cluster0_output =
  cor_data_cluster0 %>%
  dplyr::filter(p_adjust < 0.05) %>%
  dplyr::left_join(microbiome_variable_info, by = c("data_set2" = "variable_id")) %>%
  dplyr::select(-p_adjust2) %>%
  dplyr::arrange(p_adjust)

cor_data_cluster1_output =
  cor_data_cluster1 %>%
  dplyr::filter(p_adjust < 0.05) %>%
  dplyr::left_join(microbiome_variable_info, by = c("data_set2" = "variable_id")) %>%
  dplyr::select(-p_adjust2) %>%
  dplyr::arrange(p_adjust)


openxlsx::write.xlsx(cor_data_output,
                     "cor_data_output.xlsx",
                     asTable = TRUE,
                     overwrite = TRUE)
openxlsx::write.xlsx(
  cor_data_cluster0_output,
  "cor_data_cluster0_output.xlsx",
  asTable = TRUE,
  overwrite = TRUE
)

openxlsx::write.xlsx(
  cor_data_cluster1_output,
  "cor_data_cluster1_output.xlsx",
  asTable = TRUE,
  overwrite = TRUE
)

data1 =
  cor_data_cluster0_output %>%
  dplyr::filter(p_adjust < 0.05) %>%
  dplyr::arrange(data_set1) %>%
  dplyr::mutate(class = "Cluster 0")

data2 =
  cor_data_cluster1_output %>%
  dplyr::filter(p_adjust < 0.05) %>%
  dplyr::arrange(data_set1) %>%
  dplyr::mutate(class = "Cluster 1")

data =
  rbind(data1, data2) %>%
  dplyr::arrange(data_set1)

openxlsx::write.xlsx(data,
                     file = "significant_cor_cluster0_1.xlsx",
                     asTable = TRUE,
                     overwrite = TRUE)

#####output the cor plot
idx = which(cor_data_cluster1$p_adjust < 0.05)

# for (i in idx) {
#   cat(i, " ")
#   x = as.numeric(nutrition_expression_data[cor_data_cluster1$data_set1[i], ])
#   y = as.numeric(microbiome_expression_data[cor_data_cluster1$data_set2[i], ])
#   data.frame(x, y) %>%
#     ggplot(aes(x, y)) +
#     geom_point()
# }

cor_data_cluster1 %>%
  dplyr::arrange(desc(abs(cor))) %>%
  head()

# plot(as.numeric(nutrition_expression_data["VitC", sample_info_cluster1$sample_id]),
#      as.numeric(microbiome_expression_data["nRPLC_441.3947_11.3", sample_info_cluster1$sample_id]))

# abline(0, -1)

cor_data_cluster0 %>%
  dplyr::arrange(desc(abs(cor))) %>%
  head()

# plot(as.numeric(nutrition_expression_data["VitC", sample_info_cluster0$sample_id]),
#      as.numeric(microbiome_expression_data["nRPLC_441.3947_11.3", sample_info_cluster0$sample_id]))

# abline(0, -1)

# cor(
#   as.numeric(nutrition_expression_data["VitC", sample_info_cluster0$sample_id]),
#   as.numeric(microbiome_expression_data["nRPLC_441.3947_11.3", sample_info_cluster0$sample_id]),
#   method = "spearman"
# )

# temp_sample_info =
#   sample_info_cluster0[, c("Sex", "Age")]

# temp_sample_info$Sex[temp_sample_info$Sex == 'F'] = 0
# temp_sample_info$Sex[temp_sample_info$Sex == 'M'] = 1
# temp_sample_info$Sex = as.numeric(temp_sample_info$Sex)

# ppcor::pcor.test(
#   x = as.numeric(nutrition_expression_data["VitE_a_Toco", sample_info_cluster0$sample_id]),
#   y = as.numeric(microbiome_expression_data["pHILIC_732.5525_5.1", sample_info_cluster0$sample_id]),
#   z = temp_sample_info,
#   method = "spearman"
# )

unique(cor_data_cluster0$data_set2)
unique(cor_data_cluster0$data_set1)

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

####normal people, cluster0 people
cluster0_cor =
  cor_data_cluster0 %>%
  dplyr::select(-c(p:p_adjust2)) %>%
  tidyr::pivot_wider(names_from = data_set2, values_from = "cor") %>%
  tibble::column_to_rownames(var = "data_set1")

cluster0_p =
  cor_data_cluster0 %>%
  dplyr::select(-c(p, p_adjust2, cor)) %>%
  tidyr::pivot_wider(names_from = data_set2, values_from = "p_adjust") %>%
  tibble::column_to_rownames(var = "data_set1")

###remove the diet all are zero
remove_idx <-
  apply(cluster0_cor, 1, function(x) {
    sum(x == 0) / ncol(cluster0_cor)
  }) %>%
  `>=`(0.9) %>%
  which()

if (length(remove_idx) > 0) {
  cluster0_cor <- cluster0_cor[-remove_idx, ]
  cluster0_p <- cluster0_p[-remove_idx, ]
}

####prediabetes, cluster1 people
cluster1_cor =
  cor_data_cluster1 %>%
  dplyr::select(-c(p:p_adjust2)) %>%
  tidyr::pivot_wider(names_from = data_set2, values_from = "cor") %>%
  tibble::column_to_rownames(var = "data_set1")

cluster1_p =
  cor_data_cluster1 %>%
  dplyr::select(-c(p, p_adjust2, cor)) %>%
  tidyr::pivot_wider(names_from = data_set2, values_from = "p_adjust") %>%
  tibble::column_to_rownames(var = "data_set1")

###remove the diet all are zero
remove_idx <-
  apply(cluster1_cor, 1, function(x) {
    sum(x == 0) / ncol(cluster1_cor)
  }) %>%
  `>=`(0.9) %>%
  which()

if (length(remove_idx) > 0) {
  cluster1_cor <- cluster1_cor[-remove_idx, ]
  cluster1_p <- cluster1_p[-remove_idx, ]
}

colnames(all_cor) =
  colnames(cluster1_cor) =
  colnames(cluster0_cor) =
  microbiome_variable_info$variable_id[match(colnames(cluster1_cor),
                                             microbiome_variable_info$variable_id)]

library(ComplexHeatmap)
library(circlize)

col_fun = circlize::colorRamp2(
  breaks = seq(-1, 1, length.out = 11),
  colors = rev(RColorBrewer::brewer.pal(n = 11, name = "BrBG"))
)

library(wesanderson)

plot =
  Heatmap(
    t(cluster1_cor),
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
      if (t(cluster1_p)[i, j] < 0.05) {
        grid.text("*", x, y - gb_h1 * 0.5 + gb_w1 * 0.4, gp = gpar(col = "red"))
      }
      
      # if (t(cluster1_p)[i, j] > 0.05 & t(cluster1_p)[i, j] < 0.2) {
      #   grid.points(pch = 20, x = x, y = y, size = unit(0.3, "char"), gp = gpar(col = "red"))
      # }
    }
  )

plot = ggplotify::as.ggplot(plot)
plot

ggsave(plot,
       filename = "cluster1_cor_all_metabolite.pdf",
       width = 10,
       height = 10)

plot <-
  Heatmap(
    t(cluster0_cor),
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
      if (t(cluster0_p)[i, j] < 0.05) {
        grid.text("*", x, y - gb_h1 * 0.5 + gb_w1 * 0.4, gp = gpar(col = "red"))
      }
      
      # if (t(cluster0_p)[i, j] > 0.05 & t(cluster0_p)[i, j] < 0.2) {
      #   grid.points(pch = 20, x = x, y = y, size = unit(0.3, "char"), gp = gpar(col = "red"))
      #   # grid.text(label = ".", x = x, y - gb_h2 * 0.5 + gb_w1 * 0.4)
      # }
    }
  )

plot = ggplotify::as.ggplot(plot)
plot
ggsave(plot,
       filename = "cluster0_cor_all_metabolite.pdf",
       width = 10,
       height = 10)



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
ggsave(plot,
       filename = "all_cor_all_metabolite.pdf",
       width = 10,
       height = 10)

######only show the metabolites with at least one significant correlation with nutrition
###remove the metabolites which have no significant cor (p_adjust > 0.2) with nutrition
remove_metabolite1 =
  cor_data_cluster0 %>%
  dplyr::group_by(data_set2) %>%
  dplyr::summarise(n = sum(p_adjust < 0.05)) %>%
  dplyr::filter(n == 0) %>%
  dplyr::pull(data_set2)

remove_metabolite2 =
  cor_data_cluster1 %>%
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


####all participants, cluster0 people
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

####normal participants, cluster0 people
cluster0_cor =
  cor_data_cluster0 %>%
  dplyr::select(-c(p:p_adjust2)) %>%
  tidyr::pivot_wider(names_from = data_set2, values_from = "cor") %>%
  tibble::column_to_rownames(var = "data_set1")

cluster0_p =
  cor_data_cluster0 %>%
  dplyr::select(-c(p, p_adjust2, cor)) %>%
  tidyr::pivot_wider(names_from = data_set2, values_from = "p_adjust") %>%
  tibble::column_to_rownames(var = "data_set1")

###remove the diet all are zero
remove_idx <-
  apply(cluster0_cor, 1, function(x) {
    sum(x == 0) / ncol(cluster0_cor)
  }) %>%
  `>=`(0.9) %>%
  which()

if (length(remove_idx) > 0) {
  cluster0_cor <- cluster0_cor[-remove_idx, ]
  cluster0_p <- cluster0_p[-remove_idx, ]
}

cluster1_cor =
  cor_data_cluster1 %>%
  dplyr::select(-c(p:p_adjust2)) %>%
  tidyr::pivot_wider(names_from = data_set2, values_from = "cor") %>%
  tibble::column_to_rownames(var = "data_set1")

cluster1_p =
  cor_data_cluster1 %>%
  dplyr::select(-c(p, p_adjust2, cor)) %>%
  tidyr::pivot_wider(names_from = data_set2, values_from = "p_adjust") %>%
  tibble::column_to_rownames(var = "data_set1")

###remove the diet all are zero
remove_idx <-
  apply(cluster1_cor, 1, function(x) {
    sum(x == 0) / ncol(cluster1_cor)
  }) %>%
  `>=`(0.9) %>%
  which()

if (length(remove_idx) > 0) {
  cluster1_cor <- cluster1_cor[-remove_idx, ]
  cluster1_p <- cluster1_p[-remove_idx, ]
}

###remove the metabolites which we need to remove
all_cor =
  all_cor %>%
  dplyr::select(-remove_metabolite)

all_p =
  all_p %>%
  dplyr::select(-remove_metabolite)

cluster0_cor =
  cluster0_cor %>%
  dplyr::select(-remove_metabolite)

cluster0_p =
  cluster0_p %>%
  dplyr::select(-remove_metabolite)

cluster1_cor =
  cluster1_cor %>%
  dplyr::select(-remove_metabolite)

cluster1_p =
  cluster1_p %>%
  dplyr::select(-remove_metabolite)

colnames(all_cor)

colnames(cluster0_cor)

colnames(cluster1_cor)

library(ComplexHeatmap)

colnames(all_cor) =
  colnames(cluster1_cor) =
  colnames(cluster0_cor) =
  microbiome_variable_info$variable_id[match(colnames(cluster1_cor),
                                             microbiome_variable_info$variable_id)]

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
ggsave(plot,
       filename = "all_cor.pdf",
       width = 10,
       height = 10)

plot =
  Heatmap(
    t(cluster1_cor),
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
      if (t(cluster1_p)[i, j] < 0.05) {
        grid.text("*", x, y - gb_h1 * 0.5 + gb_w1 * 0.4, gp = gpar(col = "red"))
      }
      
      if (t(cluster1_p)[i, j] > 0.05 & t(cluster1_p)[i, j] < 0.2) {
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
ggsave(plot,
       filename = "cluster1_cor.pdf",
       width = 10,
       height = 10)

plot =
  Heatmap(
    t(cluster0_cor),
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
      if (t(cluster0_p)[i, j] < 0.05) {
        grid.text("*", x, y - gb_h1 * 0.5 + gb_w1 * 0.4, gp = gpar(col = "red"))
      }
      
      if (t(cluster0_p)[i, j] > 0.05 & t(cluster0_p)[i, j] < 0.2) {
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
ggsave(plot,
       filename = "cluster0_cor.pdf",
       width = 10,
       height = 10)

cor_data_cluster1

temp_cor_cluster0 <-
  cluster0_cor  %>%
  tibble::rownames_to_column(var = "food_group") %>%
  tidyr::pivot_longer(cols = -food_group,
                      names_to = "microbiome",
                      values_to = "cor")

temp_p_cluster0 <-
  cluster0_p  %>%
  tibble::rownames_to_column(var = "food_group") %>%
  tidyr::pivot_longer(cols = -food_group,
                      names_to = "microbiome",
                      values_to = "p")

temp_cluster0 <-
  temp_cor_cluster0  %>%
  dplyr::left_join(temp_p_cluster0, by = c("food_group", "microbiome"))

temp_cor_cluster1 <-
  cluster1_cor  %>%
  tibble::rownames_to_column(var = "food_group") %>%
  tidyr::pivot_longer(cols = -food_group,
                      names_to = "microbiome",
                      values_to = "cor")

temp_p_cluster1 <-
  cluster1_cor  %>%
  tibble::rownames_to_column(var = "food_group") %>%
  tidyr::pivot_longer(cols = -food_group,
                      names_to = "microbiome",
                      values_to = "p")

temp_cluster1 <-
  temp_cor_cluster1  %>%
  dplyr::left_join(temp_p_cluster1, by = c("food_group", "microbiome"))

dim(temp_cluster0)
dim(temp_cluster1)

###use histgram to show the distribution of the cor
library(ggplot2)
library(ggpubr)
temp_cluster0  %>%
  ggplot(aes(cor)) +
  geom_histogram(binwidth = 0.1) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


###use the density plot to show the distribution of the cor

temp <-
  rbind(
    data.frame(temp_cluster0, class = "Cluster 0"),
    data.frame(temp_cluster1, class = "Cluster 1")
  )

###the density plot below use the color to fill the density plot
plot <-
  temp %>%
  ggplot(aes(cor)) +
  geom_density(aes(color = class)) +
  theme_bw() +
  scale_color_manual(values = cluster0_1_color) +
  labs(x = "Correlation", y = "Density")

plot

ggsave(plot,
       filename = "cor_density_comparison.pdf",
       width = 8,
       height = 6)






