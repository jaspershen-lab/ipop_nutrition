####calculate the partial spearman correlation for data_set1 and metabolome
partial_cor <-
  function(data_set1,
           data_set2,
           sample_info,
           method = c("spearman", "pearson", "all"),
           threads = 5) {
    # browser()
    method = match.arg(method)
    library(future)
    library(furrr)
    tryCatch(
      expr = plan(strategy = multisession(workers = threads)),
      error = function(x) {
        cat("something is wrong.\n")
      }
    )
    
    cor =
      rownames(data_set1) %>%
      furrr::future_map(function(name) {
        x = as.numeric(data_set1[name,])
        temp_partial_cor =
          purrr::map(
            as.data.frame(t(data_set2)),
            .f = function(y) {
              temp_data =
                data.frame(x = x,
                           y = y,
                           sample_info)
              # dplyr::filter(!is.na(sspg_status))
              temp_data$Sex[temp_data$Sex == 'F'] = 0
              temp_data$Sex[temp_data$Sex == 'M'] = 1
              temp_data$Sex = as.numeric(temp_data$Sex)
              
              # temp_data$Ethnicity[temp_data$Ethnicity == 'C'] = 0
              # temp_data$Ethnicity[temp_data$Ethnicity == 'H'] = 1
              # temp_data$Ethnicity[temp_data$Ethnicity == 'B'] = 2
              # temp_data$Ethnicity[temp_data$Ethnicity == 'A'] = 3
              # temp_data$Ethnicity = as.numeric(temp_data$Ethnicity)
              
              ##partial correlation
              if (method == "all") {
                cor_value1 =
                  ppcor::pcor.test(
                    x = temp_data$x,
                    y = temp_data$y,
                    z = temp_data[, c("Sex", "Age")],
                    method = "spearman"
                  )
                
                result1 =
                  c(cor = cor_value1$estimate,
                    p = cor_value1$p.value)
                
                if (is.na(result1[1])) {
                  result1[1] = 0
                }
                
                if (is.na(result1[2])) {
                  result1[2] = 1
                }
                
                cor_value2 =
                  ppcor::pcor.test(
                    x = temp_data$x,
                    y = temp_data$y,
                    z = temp_data[, c("Sex", "Age")],
                    method = "pearson"
                  )
                
                result2 =
                  c(cor = cor_value2$estimate,
                    p = cor_value2$p.value)
                
                if (is.na(result2[1])) {
                  result2[1] = 0
                }
                
                if (is.na(result2[2])) {
                  result2[2] = 1
                }
                
                list(result1 = result1,
                     result2 = result2)
                
              } else{
                cor_value =
                  ppcor::pcor.test(
                    x = temp_data$x,
                    y = temp_data$y,
                    z = temp_data[, c("Sex", "Age")],
                    method = method
                  )
                
                result =
                  c(cor = cor_value$estimate,
                    p = cor_value$p.value)
                
                if (is.na(result[1])) {
                  result[1] = 0
                }
                
                if (is.na(result[2])) {
                  result[2] = 1
                }
                list(result)
              }
            }
          )
        # do.call(rbind, .) %>%
        # as.data.frame()
        
        if (method == "all") {
          temp_partial_cor1 =
            temp_partial_cor %>%
            purrr::map(function(x) {
              x[[1]]
            }) %>%
            do.call(rbind, .) %>%
            as.data.frame() %>%
            tibble::rownames_to_column(var = "data_set2") %>%
            dplyr::mutate(data_set1 = name) %>%
            dplyr::select(data_set1, data_set2, dplyr::everything())
          
          temp_partial_cor1$p_adjust = p.adjust(temp_partial_cor1$p, method = "BH")
          
          temp_partial_cor2 =
            temp_partial_cor %>%
            purrr::map(function(x) {
              x[[2]]
            }) %>%
            do.call(rbind, .) %>%
            as.data.frame() %>%
            tibble::rownames_to_column(var = "data_set2") %>%
            dplyr::mutate(data_set1 = name) %>%
            dplyr::select(data_set1, data_set2, dplyr::everything())
          
          temp_partial_cor2$p_adjust = p.adjust(temp_partial_cor2$p, method = "BH")
          list(temp_partial_cor1 = temp_partial_cor1,
               temp_partial_cor2 = temp_partial_cor2)
        } else{
          temp_partial_cor =
            temp_partial_cor %>%
            purrr::map(function(x) {
              x[[1]]
            }) %>%
            do.call(rbind, .) %>%
            as.data.frame() %>%
            tibble::rownames_to_column(var = "data_set2") %>%
            dplyr::mutate(data_set1 = name) %>%
            dplyr::select(data_set1, data_set2, dplyr::everything())
          temp_partial_cor$p_adjust = p.adjust(temp_partial_cor$p, method = "BH")
          list(temp_partial_cor)
        }
      }, .progress = TRUE)
    # do.call(rbind, .) %>%
    # as.data.frame()
    
    if (method == "all") {
      cor1 =
        cor %>%
        purrr::map(function(x) {
          x[[1]]
        }) %>%
        do.call(rbind, .) %>%
        as.data.frame()
      cor1$p_adjust2 =
        p.adjust(cor1$p, method = "BH")
      
      cor2 =
        cor %>%
        purrr::map(function(x) {
          x[[2]]
        }) %>%
        do.call(rbind, .) %>%
        as.data.frame()
      
      cor2$p_adjust2 =
        p.adjust(cor2$p, method = "BH")
      list(cor1 = cor1, cor2 = cor2)
    } else{
      cor =
        cor %>%
        lapply(function(x) {
          x[[1]]
        }) %>%
        do.call(rbind, .) %>%
        as.data.frame()
      cor$p_adjust2 =
        p.adjust(cor$p, method = "BH")
      list(cor)
    }
    
    
  }





####calculate the partial spearman correlation for data_set1 and metabolome
partial_cor_pattern = function(data_set1,
                               data_set2,
                               sample_info,
                               method = c("spearman", "pearson", "all"),
                               threads = 5) {
  # browser()
  method = match.arg(method)
  library(future)
  library(furrr)
  
  tryCatch(
    expr = plan(strategy = multisession(workers = threads)),
    error = function(x) {
      cat("something is wrong.\n")
    }
  )
  
  cor =
    rownames(data_set1) %>%
    furrr::future_map(function(name) {
      x = as.numeric(data_set1[name,])
      temp_partial_cor =
        purrr::map(
          as.data.frame(t(data_set2)),
          .f = function(y) {
            # cat(y[1], " ")
            temp_data =
              data.frame(x = x,
                         y = y,
                         sample_info)
            library(plyr)
            
            temp_data =
              temp_data %>%
              plyr::dlply(.variables = .(subject_id)) %>%
              purrr::map(function(x) {
                # cat(x$subject_id[1])
                x =
                  x %>%
                  dplyr::arrange(Time)
                
                x$x = x$x - x$x[1]
                x$y = x$y - x$y[1]
                
                x1 = get_slope_data(x$x)
                x1_auc = get_auc(value = x$x)
                
                y1 = get_slope_data(x$y)
                y1_auc = get_auc(value = x$y)
                
                value1 = c(x$x, x1, x1_auc)
                value2 = c(x$y, y1, y1_auc)
                
                data.frame(x = value1,
                           y = value2,
                           subject_id = x$subject_id[1])
              }) %>%
              dplyr::bind_rows()
            
            result =
              cor.test(temp_data$x, temp_data$y, method = method)
            
            data.frame(
              data_set1 = name,
              cor = ifelse(is.na(unname(
                result$estimate
              )), 0, unname(result$estimate)),
              p = ifelse(is.na(unname(
                result$p.value
              )), 1, unname(result$p.value))
            )
          }
        ) %>%
        dplyr::bind_rows()
      
      temp_partial_cor$data_set2 =
        rownames(data_set2)
      
      temp_partial_cor$p_adjust =
        p.adjust(temp_partial_cor$p, method = "BH")
      
      temp_partial_cor =
        temp_partial_cor %>%
        dplyr::select(data_set1, data_set2, everything())
      
    }, .progress = TRUE) %>%
    dplyr::bind_rows()
  
  cor$p_adjust2 = p.adjust(cor$p, method = "BH")
  cor
}



get_slope_data = function(value = c(1, 2, 3, 4, 5)) {
  value = value - value[1]
  purrr::map(
    1:(length(value) - 1),
    .f = function(i) {
      value[(i + 1):length(value)] - value[i]
    }
  ) %>%
    unlist()
}


get_auc = function(value = c(1, 2, 3, 4, 5)) {
  value = value - value[1]
  pracma::trapz(x = 1:length(value), y = value)
}


data_class_color <-
  c(
    "food" = ggsci::pal_lancet()(n = 9)[1],
    "microbiome" = ggsci::pal_lancet()(n = 9)[4],
    "metabolomics" = ggsci::pal_lancet()(n = 9)[5]
  )


ir_is_color <-
  c("IR" = "#EE0000FF",
    "IS" = "#008280FF")

cluster0_1_color <-
  c(
    "Cluster 0" = "#EE0000FF",
    "Cluster 1" = "#008280FF"
  )
