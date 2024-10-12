library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

##read data
data <- readr::read_csv("3_data_analysis/coefficients.csv")

dir.create("3_data_analysis/cytokine_coefficients", recursive = TRUE)
setwd("3_data_analysis/cytokine_coefficients")

data

data <-
data %>% 
  dplyr::select(Y_Col, "Diet Pattern", SSPG) %>% 
  tibble::column_to_rownames("Y_Col")

# Load required libraries
library(ggplot2)

library(ComplexHeatmap)
library(circlize)

col_fun1 = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

circos.clear()

circos.par(start.degree = 90, gap.degree = 20)

circos.heatmap(data,  col = col_fun1,
               rownames.side = "inside",
               dend.side = "outside",
               cell.border = "black")

