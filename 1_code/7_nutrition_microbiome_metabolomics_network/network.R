no_function()
setwd(r4projects::get_project_wd())
library(tidyverse)
library(plyr)
rm(list = ls())
source("1_code/tools.R")

###load variable_info
load("3_data_analysis/nutrition_group/data_preparation/variable_info")
load("3_data_analysis/nutrition/data_preparation/variable_info")
nutrition_variable_info <-
  variable_info

load("3_data_analysis/gut_microbiome/data_preparation/variable_info")
microbiome_variable_info <-
  variable_info

load("3_data_analysis/metabolomics/data_preparation/variable_info")
metabolomics_variable_info <-
  variable_info

###load correlation data
load("3_data_analysis/nutrition_vs_microbiome_t1_t4/based_on_sspg/cor_data_IS")
load("3_data_analysis/nutrition_vs_microbiome_t1_t4/based_on_sspg/cor_data_IR")

nutrition_microbiome_cor_is <-
  cor_data_IS %>%
  dplyr::filter(p_adjust < 0.05) %>%
  dplyr::arrange(data_set1) %>%
  dplyr::mutate(class = "IS") %>%
  dplyr::mutate(class1 = "food",
                class2 = "microbiome") %>%
  dplyr::mutate(node_name1 = data_set1,
                node_name2 = data_set2) %>%
  dplyr::rename(node1 = data_set1,
                node2 = data_set2)

nutrition_microbiome_cor_ir <-
  cor_data_IR %>%
  dplyr::filter(p_adjust < 0.05) %>%
  dplyr::arrange(data_set1) %>%
  dplyr::mutate(class = "IR") %>%
  dplyr::mutate(class1 = "food",
                class2 = "microbiome") %>%
  dplyr::mutate(node_name1 = data_set1,
                node_name2 = data_set2) %>%
  dplyr::rename(node1 = data_set1,
                node2 = data_set2)

###load correlation data
load("3_data_analysis/nutrition_vs_metabolomics_t1_t4/based_on_sspg/cor_data_IS")
load("3_data_analysis/nutrition_vs_metabolomics_t1_t4/based_on_sspg/cor_data_IR")

nutrition_metabolomics_cor_is <-
  cor_data_IS %>%
  dplyr::filter(p_adjust < 0.05) %>%
  dplyr::arrange(data_set1) %>%
  dplyr::mutate(class = "IS") %>%
  dplyr::mutate(class1 = "food",
                class2 = "metabolomics") %>%
  dplyr::left_join(metabolomics_variable_info[, c("variable_id", "true_name")],
                   by = c("data_set2" = "variable_id")) %>%
  dplyr::rename(node_name1 = true_name) %>%
  dplyr::mutate(node_name2 = data_set2) %>%
  dplyr::rename(node1 = data_set1,
                node2 = data_set2)

nutrition_metabolomics_cor_ir <-
  cor_data_IR %>%
  dplyr::filter(p_adjust < 0.05) %>%
  dplyr::arrange(data_set1) %>%
  dplyr::mutate(class = "IR") %>%
  dplyr::mutate(class1 = "food",
                class2 = "metabolomics") %>%
  dplyr::left_join(metabolomics_variable_info[, c("variable_id", "true_name")],
                   by = c("data_set2" = "variable_id")) %>%
  dplyr::rename(node_name1 = true_name) %>%
  dplyr::mutate(node_name2 = data_set2) %>%
  dplyr::rename(node1 = data_set1,
                node2 = data_set2)

dir.create("3_data_analysis/7_nutrition_microbiome_metabolomics_network",
           recursive = TRUE)
setwd("3_data_analysis/7_nutrition_microbiome_metabolomics_network")

edge_data_is <-
  rbind(nutrition_microbiome_cor_is,
        nutrition_metabolomics_cor_is)

edge_data_ir <-
  rbind(nutrition_microbiome_cor_ir,
        nutrition_metabolomics_cor_ir)

edge_data_is$p_adjust[edge_data_is$p_adjust == 0] <-
  min(edge_data_is$p_adjust[edge_data_is$p_adjust != 0])

edge_data_ir$p_adjust[edge_data_ir$p_adjust == 0] <-
  min(edge_data_ir$p_adjust[edge_data_ir$p_adjust != 0])

node_data_is <-
  rbind(
    edge_data_is[, c("node1", "node_name1", "class1")] %>%
      dplyr::rename(
        node = node1,
        node_name = node_name1,
        class = class1
      ),
    edge_data_is[, c("node2", "node_name2", "class2")] %>%
      dplyr::rename(
        node = node2,
        node_name = node_name2,
        class = class2
      )
  ) %>%
  dplyr::distinct(node, .keep_all = TRUE)

node_data_ir <-
  rbind(
    edge_data_ir[, c("node1", "node_name1", "class1")] %>%
      dplyr::rename(
        node = node1,
        node_name = node_name1,
        class = class1
      ),
    edge_data_ir[, c("node2", "node_name2", "class2")] %>%
      dplyr::rename(
        node = node2,
        node_name = node_name2,
        class = class2
      )
  ) %>%
  dplyr::distinct(node, .keep_all = TRUE)

dim(node_data_ir)
dim(node_data_is)

intersect(node_data_ir$node,
          node_data_is$node)

node_data <-
  rbind(node_data_ir,
        node_data_is) %>%
  dplyr::distinct(node, .keep_all = TRUE)

library(ggraph)
library(igraph)
library(tidygraph)

edge_data_ir <-
  edge_data_ir %>%
  dplyr::rowwise() %>%
  dplyr::mutate(edge_name = paste(node1, node2, sep = "-")) %>%
  dplyr::mutate(class = case_when(cor > 0 ~ "pos",
                                  cor < 0 ~ "neg"))
edge_data_is <-
  edge_data_is %>%
  dplyr::rowwise() %>%
  dplyr::mutate(edge_name = paste(node1, node2, sep = "-")) %>%
  dplyr::mutate(class = case_when(cor > 0 ~ "pos",
                                  cor < 0 ~ "neg"))
edge_data <-
  rbind(edge_data_ir,
        edge_data_is) %>%
  dplyr::distinct(edge_name, .keep_all = TRUE) %>%
  dplyr::mutate(class = case_when(cor > 0 ~ "pos",
                                  cor < 0 ~ "neg"))

graph_all <-
  tidygraph::tbl_graph(nodes = node_data,
                       edges = edge_data,
                       directed = TRUE) %>%
  dplyr::mutate(Degree = centrality_degree(mode = 'all'))

graph_ir <-
  tidygraph::tbl_graph(nodes = node_data,
                       edges = edge_data_ir,
                       directed = TRUE) %>%
  dplyr::mutate(Degree = centrality_degree(mode = 'all'))

graph_is <-
  tidygraph::tbl_graph(nodes = node_data,
                       edges = edge_data_is,
                       directed = TRUE) %>%
  dplyr::mutate(Degree = centrality_degree(mode = 'all'))

pal <-
  wesanderson::wes_palette(name = "Zissou1", n = 100, type = "continuous")

g <- graph_all

V(g)$type <- igraph::bipartite_mapping(g)$type

coords <-
  create_layout(g, layout = "bipartite")

coords$y[coords$class == "food"] = 0
coords$y[coords$class == "microbiome"] = 0.6
coords$y[coords$class == "metabolomics"] = 1

coords <-
  coords %>%
  dplyr::mutate(
    theta = x / (max(x) + 1) * 2 * pi,
    r = y + 1,
    x = r * cos(theta),
    y = r * sin(theta)
  )

cor_color <-
  rev(RColorBrewer::brewer.pal(n = 11, name = "BrBG"))

network_is <-
  ggraph(graph_is,
         layout = 'manual',
         x = coords$x,
         y = coords$y) +
  geom_edge_diagonal(show.legend = TRUE,
                     aes(color = cor,
                         width = -log(p_adjust, 10))) +
  geom_node_point(aes(fill = class,
                      size = Degree),
                  shape = 21,
                  show.legend = TRUE) +
  geom_node_text(
    aes(
      x = coords$x,
      y = coords$y,
      label = ifelse(class == "food", node_name, NA)
    ),
    size = 2,
    show.legend = FALSE
  ) +
  scale_edge_color_gradientn(colors = cor_color) +
  scale_fill_manual(values = data_class_color) +
  scale_size_continuous(range = c(1, 5)) +
  scale_edge_width_continuous(range = c(0.1, 2)) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "bottom"
  )
network_is
extrafont::loadfonts()
ggsave(network_is,
       filename = "network_is.pdf",
       width = 6,
       height = 6.21)

network_ir <-
  ggraph(graph_ir,
         layout = 'manual',
         x = coords$x,
         y = coords$y) +
  geom_edge_diagonal(show.legend = TRUE,
                     aes(color = cor,
                         width = -log(p_adjust, 10))) +
  geom_node_point(aes(fill = class,
                      size = Degree),
                  shape = 21,
                  show.legend = TRUE) +
  scale_edge_color_gradientn(colors = cor_color) +
  scale_fill_manual(values = data_class_color) +
  scale_size_continuous(range = c(1, 5)) +
  scale_edge_width_continuous(range = c(0.1, 2)) +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "bottom"
  )
network_ir
ggsave(network_ir,
       filename = "network_ir.pdf",
       width = 6,
       height = 6.21)
