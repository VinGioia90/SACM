###############################################################################
# Needed function for visualising the performance on the test set and the     #
# stepwise-related subset of selected models for the Application of Section 7 #
###############################################################################



###############################################################
# LogScores:                                                  #
# VALIDATION: Currently no rolling origin forecasting         #
# This function takes in input                                #
# -) obj a list of fitted model with decreasing number        #
#    of effects modelling the entries of the Theta matrix     #
# -) dgrid: step length considered                            #
# -) d: dimensionsion of the outcome                          #
# -) data_train: train set (to obtain in-sample logScore)     #
# -) data_test: test set (to obtain out-of-sample logScore)   #
#                                                             #
# The function returns a list including in- and out-of-sample #
# logScore                                                    #
###############################################################
#
# res_perf <- function(obj, dgrid, d, data_train, data_test, param){
#   grid <- seq(d + d * (d + 1)/2, d, -dgrid)
#   lpi_pred_in <- lapply(1 : length(grid), function(x) predict(obj[[x]]))
#   lpi_pred_out <- lapply(1: length(grid), function(x) predict(obj[[x]], newdata = data_test ))
#
#   if(param == "mcd"){
#     logSin <- unlist(lapply(1: length(lpi_pred_in), function(x) -SCM::internal()$ll_mcd(lpi_pred_in[[x]],as.matrix(GEF14_data[1 : nrow(data_train), 1:d]))))
#     logSout <- unlist(lapply(1: length(lpi_pred_out), function(x) -SCM::internal()$ll_mcd(lpi_pred_out[[x]],as.matrix(GEF14_data[1 : nrow(data_test), 1:d]))))
#   } else {
#     logSin <- unlist(lapply(1: length(lpi_pred_in), function(x) -SCM::internal()$ll_mcd(lpi_pred_in[[x]],as.matrix(GEF14_data[1 : nrow(data_train), 1:d]))))
#     logSout <- unlist(lapply(1: length(lpi_pred_out), function(x) -SCM::internal()$ll_mcd(lpi_pred_out[[x]],as.matrix(GEF14_data[1 : nrow(data_test), 1:d]))))
#   }
#   return(list(logSin = logSin, logSout = logSout))
# }

res_perf <- function(obj, d, data, param, sets_eval){
  ncv_rof <- length(obj)
  nmodel <- length(obj[[1]]$gen)

  eta_hat <- lapply(1 : nmodel,
                    function(ii){ do.call("rbind", lapply(1 : ncv_rof, function(x) obj[[x]]$gen[[ii]]$lpi_pred_out)) } )

  y <- as.matrix(data[(sets_eval[1]+1): sets_eval[length(sets_eval)], 1:d])
  # Get log-likelihood at a function of number of effects
  if(param == "mcd"){
    nloglik <- sapply(1 : nmodel, function(ii) return( -SCM::internal()$ll_mcd(eta_hat[[ii]], y) ))
  } else{
    nloglik <- sapply(1 : nmodel, function(ii) return( -SCM::internal()$ll_logm(eta_hat[[ii]], y) ))
  }

  return(nloglik)
}


res_perf_MSE_mar <- function(obj, d, data, param, sets_eval){
  ncv_rof <- length(obj)
  nmodel <- length(obj[[1]]$gen)

  eta_hat <- lapply(1 : nmodel,
                    function(ii){ do.call("rbind", lapply(1 : ncv_rof, function(x) obj[[x]]$gen[[ii]]$lpi_pred_out[,1:d])) } )

  y <- as.matrix(data[(sets_eval[1]+1): sets_eval[length(sets_eval)], 1:d])

  MSE_mar <- sapply(1 : nmodel, function(ii) return( apply((y-eta_hat[[ii]])^2, 2, mean )))

  return(MSE_mar)
}


############################################################################################
# This function takes in input:                                                            #
# -) the list of fitted models                                                             #
# -) a vector of the name of effects included in the covariance matrix model formula       #
#   (actually does not work since it is set to a single effect)                            #
# -) dimension of the outcome                                                              #
# and return a matrix including the name of effects and their position in Theta matrix     #
############################################################################################
# get_eff_idx <- function(obj, name_eff, d){
#   out <- data.frame(matrix(NA, length(name_eff) * d * (d + 1)/2, 6))        # NA for the element of Theta modelled via intercepts
#   out[, 1] <- rep(name_eff, each = d * (d + 1)/2) # name of the effects
#
#   count <- 1
#   for(j in 1 : length(name_eff)){                               # for loop over the effects
#     for(k in (d + 1) : length(obj$foo_print)){                  # for loop over the elements of the covariance matrix model formula
#       if(str_detect(deparse(obj$foo_print[[k]]), name_eff[j])){ # if the effect is in the model formula save the results in the out matrix
#         out[count, 2] <- lhs.vars(obj$foo_print[[k]])           # lhs of the model formula Th_<something>
#         length_string <- nchar(out[count, 2])
#         out[count, 3] <- as.numeric(substr(out[count, 2],       # get the row of Theta
#                                            start = 4 ,
#                                            stop =   gregexpr("\\.",  out[count, 2])[[1]][1]-1))
#         out[count, 4] <- as.numeric(substr(out[count, 2],       # get the column of Theta
#                                            start = gregexpr("\\.",  out[count, 2])[[1]][1]+1,
#                                            stop = length_string))
#         out[count, 5] <- out[count, 3] + 0.32 * cos(2 * pi * j/length(name_eff)) - 1   # jittering of the rows and columns for graphical purposes
#         out[count, 6] <- out[count, 4] - 0.32 * cos(2 * pi * j/length(name_eff)) - 1
#         count <- count + 1
#       }
#     }
#   }
#   colnames(out) <- c("name_eff", "Th_el", "row", "col", "jit_row", "jit_col")
#   return(out)
# }

get_eff_idx <- function(foo_vcov, name_eff, d){
  out <- data.frame(matrix(NA, length(name_eff) * (length(foo_vcov) - d), 7))        # NA for the element of Theta modelled via intercepts
  out[, 1] <- rep(name_eff, each = (length(foo_vcov) - d)) # name of the effects

  count <- 1
  for(j in 1 : length(name_eff)){                               # for loop over the effects
    for(k in (d + 1) : length(foo_vcov)){                  # for loop over the elements of the covariance matrix model formula
      if(str_detect(deparse(foo_vcov[[k]]), name_eff[j])){ # if the effect is in the model formula save the results in the out matrix
        out[count, 2] <- lhs.vars(foo_vcov[[k]])           # lhs of the model formula Th_<something>
        length_string <- nchar(out[count, 2])
        out[count, 3] <- as.numeric(substr(out[count, 2],       # get the row of Theta
                                           start = 4 ,
                                           stop =   gregexpr("\\.",  out[count, 2])[[1]][1]-1))
        out[count, 4] <- as.numeric(substr(out[count, 2],       # get the column of Theta
                                           start = gregexpr("\\.",  out[count, 2])[[1]][1]+1,
                                           stop = length_string))
        out[count, 5] <- out[count, 3] + 0.32 * cos(2 * pi * j/length(name_eff)) - 1   # jittering of the rows and columns for graphical purposes
        out[count, 6] <- out[count, 4] - 0.32 * cos(2 * pi * j/length(name_eff)) - 1
        count <- count + 1
      }
    }
  }
  colnames(out) <- c("name_eff", "Th_el", "row", "col", "jit_row", "jit_col")
  return(out)
}


get_eff_idx2 <- function(foo_vcov, name_eff, d, Param){
  out <- data.frame(matrix(NA, length(name_eff) * (length(foo_vcov) - d), 7))        # NA for the element of Theta modelled via intercepts
  out[, 1] <- rep(name_eff, each = (length(foo_vcov) - d)) # name of the effects

  count <- 1
  for(j in 1 : length(name_eff)){                               # for loop over the effects
    for(k in (d + 1) : length(foo_vcov)){                  # for loop over the elements of the covariance matrix model formula
      if(str_detect(deparse(foo_vcov[[k]]), name_eff[j])){ # if the effect is in the model formula save the results in the out matrix
        out[count, 2] <- lhs.vars(foo_vcov[[k]])           # lhs of the model formula Th_<something>
        length_string <- nchar(out[count, 2])
        out[count, 3] <- as.numeric(substr(out[count, 2],       # get the row of Theta
                                           start = 4 ,
                                           stop =   gregexpr("\\.",  out[count, 2])[[1]][1]-1))
        out[count, 4] <- as.numeric(substr(out[count, 2],       # get the column of Theta
                                           start = gregexpr("\\.",  out[count, 2])[[1]][1]+1,
                                           stop = length_string))
        if(Param == "MCD"){
          out[count, 5] <- out[count, 3] + 0.32/2 * cos(2 * pi * j/length(name_eff)) - 1   # jittering of the rows and columns for graphical purposes
          out[count, 6] <- out[count, 4] - 0.32/2 * cos(2 * pi * j/length(name_eff)) - 1
        } else {
          out[count, 5] <- out[count, 3] + 0.32/2 * cos(2 * pi * j/length(name_eff)) - 0.5 - 1   # jittering of the rows and columns for graphical purposes
          out[count, 6] <- out[count, 4] - 0.32/2 * cos(2 * pi * j/length(name_eff)) + 0.5 - 1
        }
        out[count, 7] <- Param
        count <- count + 1
      }
    }
  }
  colnames(out) <- c("name_eff", "Th_el", "row", "col", "jit_row", "jit_col", "Param")
  return(out)
}

############################################################################################
# This function produces a list of plots related to selection process                      #
# on the effects acting on Theta matrix                                                    #
# The arguments are                                                                        #
# -) the list of fitted models                                                             #
# -) a vector of the name of effects included in the covariance matrix model formula       #
#   (actually does not work since it is set to a single effect)                            #
# -) dimension of the outcome                                                              #
# and return a matrix including the name of effects and their position in Theta matrix     #
############################################################################################

# amounf of effects included in covariance matrix modelling              #
##########################################################################
get_plots <- function(obj, name_eff, d, grid_length, param){
  if(param == "mcd"){
    param2 <- "MCD"
  }  else {
    param2 <- "logM"
  }

  #removed first (full) and last (static)
  res_plot <- lapply(2 : (length(obj$foo) - 1), function(x) get_eff_idx(foo_vcov = obj$foo[[x]], name_eff = name_eff, d = d))

  pl <- list() # save the plots in a list

  segment_data = data.frame(  # for drawing the segments
    x = c(rep(-0.5, d + 1), (seq(0.5, d + 0.5, by = 1))),
    xend = c(seq(1.5, d + 1.5, by = 1),  (seq(0.5, d + 0.5, by = 1))),
    y =  c(seq(0.5, d + 0.5, by = 1), rep(d - 0.5, d + 1)),
    yend =  c(seq(0.5, d + 0.5, by = 1), seq(-0.5, d - 0.5, by = 1))
  )



  for(j in 1 : length(res_plot)){
    eff_grid <- nrow(res_plot[[j]])
    pl[[j]] <-ggplot(res_plot[[j]], aes(x = jit_col, y = jit_row, shape = name_eff))+
      geom_point(size = 3) +
      scale_y_continuous(breaks = 0 : (d - 1), expand = c(0, 0), limits = c(d , -0.5), name = "Hours", trans = reverse_trans()) +
      scale_x_continuous(breaks = 0 : (d - 1), expand = c(0, 0), limits = c(-1, d), name = "Hours") +
      geom_segment(data = segment_data, mapping = aes(x = x, y = y, xend = xend, yend = yend),
                   inherit.aes = FALSE, colour = "gray") +
      theme_bw() +
      labs(shape = "Effect") +
      labs(title=paste(param2, " - Covariance matrix model with", eff_grid, "effects"))+
      theme(legend.title = element_text(size = 14),
            legend.text = element_text(size = 14),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.title.x = element_text(size = 16, colour = "black"),
            axis.title.y = element_text(size = 16, colour = "black"),
            axis.text.x = element_text(size = 14, colour = "black"),
            axis.text.y = element_text(size = 14, colour = "black"),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_blank(),
            plot.title = element_text(hjust = 0.5))
    #ggsave(paste0("Covmod_with", eff_grid[j], "Effects.pdf"),  plot=pl[[j]], width = 20, height = 20, units = "cm")
  }
  return(invisible(pl))
}


# This function allows for plotting the covariance matrix model selected for the MCD and logM
get_plots2 <- function(obj1_mcd, obj2_logm, name_eff, d, grid_length, neff1_mcd, neff2_logm){
  #removed first (full) and last (static)
  res_plot_mcd <- lapply(2 : (length(obj1_mcd$foo) - 1), function(x) get_eff_idx2(foo_vcov = obj1_mcd$foo[[x]],
                                                                                  name_eff = name_eff, d = d, Param = "MCD"))
  res_plot_logm <- lapply(2 : (length(obj2_logm$foo) - 1), function(x) get_eff_idx2(foo_vcov = obj2_logm$foo[[x]],
                                                                                    name_eff = name_eff, d = d, Param = "logM"))

  segment_data = data.frame(  # for drawing the segments
    x = c(rep(-0.5, d + 1), (seq(0.5, d + 0.5, by = 1))),
    xend = c(seq(1.5, d + 1.5, by = 1),  (seq(0.5, d + 0.5, by = 1))),
    y =  c(seq(0.5, d + 0.5, by = 1), rep(d - 0.5, d + 1)),
    yend =  c(seq(0.5, d + 0.5, by = 1), seq(-0.5, d - 0.5, by = 1))
  )


  seq_mcd <- seq(nrow(res_plot_mcd[[1]]), nrow(res_plot_mcd[[length(res_plot_mcd)]]), by = -grid_length)
  seq_logm <- seq(nrow(res_plot_logm[[1]]), nrow(res_plot_logm[[length(res_plot_logm)]]), by = -grid_length)


  jElem_mcd <-  which(neff1_mcd == seq_mcd)
  jElem_logm <-  which(neff2_logm == seq_logm)


  pl <- ggplot(res_plot_mcd[[jElem_mcd]], aes(x = jit_col, y = jit_row, shape = Param, col = Param))+
    geom_point(size = 3) +
    geom_point(data = res_plot_logm[[jElem_logm]], aes(x = jit_col, y = jit_row, shape = Param), size = 3) +
    scale_y_continuous(breaks = 0 : (d - 1), expand = c(0, 0), limits = c(d , -0.5), name = "Hours", trans = reverse_trans()) +
    scale_x_continuous(breaks = 0 : (d - 1), expand = c(0, 0), limits = c(-1, d), name = "Hours") +
    geom_segment(data = segment_data, mapping = aes(x = x, y = y, xend = xend, yend = yend),
                 inherit.aes = FALSE, colour = "gray") +
    scale_color_manual(name = "Parametrisation", values = c("MCD" = "#F8766D", "logM" = "#619CFF")) +
    scale_shape_manual(name = "Parametrisation", values = c("MCD" = 17, "logM" = 16)) +
    theme_bw() +
    labs(shape = "Parametrisation") +
    #labs(title=paste("MCD:", neff1_mcd, "effects; logM:", neff2_logm, "effects"))+
    theme(legend.title = element_text(size = 16),
          legend.text = element_text(size = 16),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.x = element_text(size = 16, colour = "black"),
          axis.title.y = element_text(size = 16, colour = "black"),
          axis.text.x = element_text(size = 14, colour = "black", angle = 90, vjust = 0.5),
          axis.text.y = element_text(size = 14, colour = "black"),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "bottom",
          plot.title = element_text(hjust = 0.5))
  return(invisible(pl))
}

heatmap_FitCov <- function(PredCov, PredCov2 = NULL, d, range_var = NULL, range_corr = c(0, 1), label = "h",
                           label_xaxis, label_yaxis, col_var, col_cor){

  diagDf <- data.frame(
    var1 = c(paste0("h0", 1:9), paste0("h", 10:d)),
    var2 = c(paste0("h", (d):10), paste0("h0", 9:1)),
    if(is.null(PredCov2)){
      Stdev =  sqrt(diag(PredCov))
    } else {
      Stdev =  sqrt(diag(PredCov )) - sqrt(diag(PredCov2))
    }

  )

  sDf <- data.frame(X1 = rep(NA, d ^ 2), X2 = rep(NA, d ^ 2), Corr = rep(NA, d ^ 2))

  count <- 1
  for(i in 1 : d){
    for(j in 1 : d){
      if(i >= 10 & ((d - j) + 2>= 10)){
        sDf$X1[count] <- paste0(label, i)
        sDf$X2[count] <- paste0(label, (d - j) + 1)
      }
      if(i >= 10 & ((d - j) + 2 <= 10)){
        sDf$X1[count] <- paste0(label, i)
        sDf$X2[count] <- paste0(label, "0", (d - j) + 1)
      }

      if(i < 10 & ((d - j) + 2  >= 10)){
        sDf$X1[count] <- paste0(label, "0", i)
        sDf$X2[count] <- paste0(label, (d - j) + 1)
      }

      if(i < 10 & ((d - j) + 2  <= 10)){
        sDf$X1[count] <- paste0(label, "0", i)
        sDf$X2[count] <- paste0(label, "0", (d - j) + 1)
      }

      if(i == j){
        sDf$Corr[count] <- 0
      }
      if(j < i){
        sDf$Corr[count] <- NA
      }
      if(j > i){
        if(is.null(PredCov2)){
          sDf$Corr[count] <- PredCov[j, i]
        } else {
          sDf$Corr[count] <- PredCov[j, i] - PredCov2[j, i]
        }
      }
      count <- count + 1
    }
  }


  gg1 <- ggplot(sDf, aes(X1, X2)) +
    geom_tile(aes(fill = Corr)) +
    scale_fill_gradientn(colors =  col_cor, limits = range_corr, na.value = "white") +
    new_scale_fill() +
    geom_tile(data = diagDf, aes(var1, var2, fill = Stdev)) +
    scale_fill_gradientn(colors = col_var, limits = range_var) +
    theme(aspect.ratio = 1,
          axis.title.x = element_text(size=23),
          axis.title.y = element_text(size=23),
          axis.text.x=element_text(size=23, angle = 90, vjust = 0.5, hjust=1),
          axis.text.y=element_text(size=23),
          legend.key.width  = unit(1, "lines"),
          legend.key.height = unit(1, "lines"),
          legend.text=element_text(size=16),
          legend.title=element_text(size=16))+
    scale_x_discrete(labels = label_xaxis,name = "Hours")+
    scale_y_discrete(labels = label_yaxis,name = "Hours")

  return(invisible(gg1))

}

# At first we fit the model and save the results
fit_model_AllData <- function(data, flag_res = FALSE, param = c("mcd", "logm"),
                              model_type = c("reduced", "full"), neff_reduced = NULL,
                              grid_length = 5, grid_d = grid_d, d = 24){
  param <- match.arg(param)
  model_type <- match.arg(model_type)

  setwd(root_dir)
  setwd("content/Section6/Results")
  if(flag_res == TRUE){
    load( paste0("res_stepwise_param", param, "_d_", d, "_lstep_", grid_length, "_residuals.RData"))
  } else {
    load( paste0("res_stepwise_param", param, "_d_", d, "_lstep_", grid_length, "_response.RData"))
  }
  if(param == "mcd"){
    res_fit <- res_mcd
    rm(res_mcd)
  } else {
    res_fit <- res_logm
    rm(res_logm)
  }



  # Static ACM model
  res_fixed <- gam_scm(res_fit$foo[[length(grid_d)]], family = mvn_scm(d = d, param = param), data = data)

  if(model_type == "full"){
    res_final <- gam_scm(res_fit$foo[[1]], family = mvn_scm(d = d, param = param), data = data)
  }
  if(model_type == "reduced"){
    jElem <-  which(neff_reduced == grid_d)
    res_final <- gam_scm(res_fit$foo[[jElem]], family = mvn_scm(d = d, param = param), data = data)
  }

  if(flag_res == TRUE){
    save(res_fixed, file = paste0("fit_AllData_param", param, "_d_", d, "_model_Static_residuals.RData"))
    save(res_final, file = paste0("fit_AllData_param", param, "_d_", d, "_model_", model_type, "_residuals.RData"))
  } else {
    save(res_fixed, file = paste0("fit_AllData_param", param, "_d_", d, "_model_Static_response.RData"))
    save(res_final, file = paste0("fit_AllData_param", param, "_d_", d, "_model_", model_type, "_response.RData"))
  }
  return(NULL)
}




plot_Heatmap <- function(data, flag_res = TRUE, model_type = c("reduced", "full"),
                         param = c("mcd", "logm"), neff_reduced = NULL,
                         grid_length = 5, grid_d = grid_d, d = 24){
  param <- match.arg(param)
  model_type <- match.arg(model_type)

  setwd(root_dir)
  setwd("content/Section6/Results")

  if(flag_res == TRUE){
    load(paste0("fit_AllData_param", param, "_d_", d, "_model_Static_residuals.RData"))
    load(paste0("fit_AllData_param", param, "_d_", d, "_model_", model_type, "_residuals.RData"))
  } else {
    load(paste0("fit_AllData_param", param, "_d_", d, "_model_Static_response.RData"))
    load(paste0("fit_AllData_param", param, "_d_", d, "_model_", model_type, "_response.RData"))
  }


  Sigma_pred_fixed <- predict(res_fixed, type = "response")
  Sigma_predMat_fixed <-  Sigma_mat(Sigma_pred_fixed[,-c(1 : d)])

  Sigma_pred_model <- predict(res_final, type = "response")
  Sigma_predMat_model <-  Sigma_mat(Sigma_pred_model[,-c(1 : d)])

  idx_min_Corr <- which.min(unlist(lapply(1 : length(Sigma_predMat_model), function(x) unlist(mean(Sigma_predMat_model[[x]][lower.tri(Sigma_predMat_model[[x]])])))))
  idx_max_Corr <- which.max(unlist(lapply(1 : length(Sigma_predMat_model), function(x) unlist(mean(Sigma_predMat_model[[x]][lower.tri(Sigma_predMat_model[[x]])])))))

  col_cor <- rev(colorspace::sequential_hcl(palette = "Blues 3", n = 100))
  col_var <- rev(colorspace::sequential_hcl(palette = "Red", n = 10)[1 : 5])

  label_xaxis <- c(paste0("h0", 0:9), paste0("h", 10:(d-1)))
  label_yaxis <- c(paste0("h", (d-1):10), paste0("h0", 9:0))
  days <- c("Sunday", "Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday")

  setwd(root_dir)
  setwd("content/Section6/Plots")

  # Static
  pl1_fixed <- heatmap_FitCov(Sigma_predMat_fixed[[1]], d = d, range_var = NULL, range_corr = c(0, 1), label = "h",
                              label_xaxis = label_xaxis, label_yaxis = label_yaxis,  col_var = col_var, col_cor = col_cor)
  print(pl1_fixed)

  # Min corr using the model predictions
  pl1 <- heatmap_FitCov(Sigma_predMat_model[[idx_min_Corr]], d = d, range_var = NULL, range_corr = c(0, 1), label = "h",
                        label_xaxis = label_xaxis, label_yaxis = label_yaxis,  col_var = col_var, col_cor = col_cor)

  date <- as.Date(paste(data[idx_min_Corr, "year"], data[idx_min_Corr, "doy"]), format = "%Y %j")
  dow <- days[wday(date)]
  pl1 <- pl1 + annotate("text",  x = (d-1) + 0.25, y = (d-1) + 0.25, label = paste(dow, date), vjust=1, hjust=1, cex = 15)
  print(pl1)

  # Max corr using the model predictions
  pl2 <- heatmap_FitCov(Sigma_predMat_model[[idx_max_Corr]], d = d, range_var = NULL, range_corr = c(0, 1), label = "h",
                        label_xaxis = label_xaxis, label_yaxis = label_yaxis,  col_var = col_var, col_cor = col_cor)
  date <- as.Date(paste(data[idx_max_Corr,"year"], data[idx_max_Corr,"doy"]), format = "%Y %j")
  dow <- days[wday(date)]
  pl2 <- pl2 +  annotate("text",  x = (d-1) + 0.25, y = (d-1) + 0.25, label = paste(dow, date), vjust = 1, hjust = 1, cex = 15)
  print(pl2)

  ############################################################################################################################
  col_var2 <- rev(colorspace::diverging_hcl(palette = "Blue - Red 2", n = 100))
  col_cor2 <- rev(colorspace::divergingx_hcl(palette = "RdBu", n = 100))[100:1]

  minV_min_Corr <- min((sqrt(diag(Sigma_predMat_model[[idx_min_Corr]])) - sqrt(diag(Sigma_predMat_fixed[[1]]))))
  maxV_min_Corr <- max((sqrt(diag(Sigma_predMat_model[[idx_min_Corr]])) - sqrt(diag(Sigma_predMat_fixed[[1]]))))

  lowTri <- lower.tri(Sigma_predMat_fixed[[1]])
  minC_min_Corr <- min((Sigma_predMat_model[[idx_min_Corr]] - Sigma_predMat_fixed[[1]])[lowTri])
  maxC_min_Corr <- max((Sigma_predMat_model[[idx_min_Corr]] - Sigma_predMat_fixed[[1]])[lowTri])


  minV_max_Corr <- min((sqrt(diag(Sigma_predMat_model[[idx_max_Corr]])) - sqrt(diag(Sigma_predMat_fixed[[1]]))))
  maxV_max_Corr <- max((sqrt(diag(Sigma_predMat_model[[idx_max_Corr]])) - sqrt(diag(Sigma_predMat_fixed[[1]]))))


  lowTri <- lower.tri(Sigma_predMat_fixed[[1]])
  minC_max_Corr <- min((Sigma_predMat_model[[idx_max_Corr]] - Sigma_predMat_fixed[[1]])[lowTri])
  maxC_max_Corr <- max((Sigma_predMat_model[[idx_max_Corr]] - Sigma_predMat_fixed[[1]])[lowTri])

  pl1_diff <- heatmap_FitCov(PredCov = Sigma_predMat_model[[idx_min_Corr]], PredCov2 = Sigma_predMat_fixed[[1]],
                             d = d, range_var = c(min(minV_min_Corr, minV_max_Corr) , max(maxV_min_Corr ,   maxV_max_Corr )),
                             range_corr = c(min(minC_min_Corr, minC_max_Corr) , max(maxC_min_Corr ,   maxC_max_Corr )), label = "h",
                             label_xaxis = label_xaxis, label_yaxis = label_yaxis,  col_var = col_var2, col_cor = col_cor2)

  date_min_Corr <- as.Date(paste(data[idx_min_Corr, "year"], data[idx_min_Corr, "doy"]), format = "%Y %j")
  dow_min_Corr <- days[wday(date_min_Corr)]
  pl1_diff  <- pl1_diff  + annotate("text",  x = (d-1) + 0.25, y = (d-1) + 0.25, label = paste(dow_min_Corr, date_min_Corr), vjust=1, hjust=1, cex = 7)
  print(pl1_diff)

  pl2_diff  <- heatmap_FitCov(Sigma_predMat_model[[idx_max_Corr]], PredCov2 = Sigma_predMat_fixed[[1]],
                              d = d,range_var = c(min(minV_min_Corr, minV_max_Corr) , max(maxV_min_Corr ,   maxV_max_Corr )),
                              range_corr = c(min(minC_min_Corr, minC_max_Corr) , max(maxC_min_Corr ,   maxC_max_Corr )), label = "h",
                              label_xaxis = label_xaxis, label_yaxis = label_yaxis,  col_var = col_var2, col_cor = col_cor2)

  date_max_Corr <- as.Date(paste(data[idx_max_Corr,"year"], data[idx_max_Corr,"doy"]), format = "%Y %j")
  dow_max_Corr <- days[wday(date_max_Corr)]
  pl2_diff <- pl2_diff +  annotate("text",  x = (d-1) + 0.25, y = (d-1) + 0.25, label = paste(dow_max_Corr, date_max_Corr), vjust = 1, hjust = 1, cex = 7)
  print(pl2_diff)

  if(flag_res == TRUE){
    ggsave(paste0("Residuals_Vcormat_param", param, "_model_fixed.pdf"),  plot =  pl1_fixed, width = 20, height = 20, units = "cm")
    ggsave(paste0("Residuals_Vcormat_lowCorr_param", param, "_model_", model_type, ".pdf"),  plot =  pl1, width = 20, height = 20, units = "cm")
    ggsave(paste0("Residuals_Vcormat_HighCorr_param", param, "_model_", model_type, ".pdf"),  plot =  pl2, width = 20, height = 20, units = "cm")
    ggsave(paste0("Residuals_Vcormat_lowCorr_param", param, "_model_", model_type, "_diffStatic.pdf"),  plot =  pl1_diff, width = 20, height = 20, units = "cm")
    ggsave(paste0("Residuals_Vcormat_HighCorr_param", param, "_model_", model_type, "_diffStatic.pdf"),  plot =  pl2_diff, width = 20, height = 20, units = "cm")
  } else {
    ggsave(paste0("Response_Vcormat_param", param, "_model_fixed.pdf"),  plot =  pl1_fixed, width = 20, height = 20, units = "cm")
    ggsave(paste0("Response_Vcormat_lowCorr_param", param, "_model_", model_type, ".pdf"),  plot =  pl1, width = 20, height = 20, units = "cm")
    ggsave(paste0("Response_Vcormat_HighCorr_param", param, "_model_", model_type, ".pdf"),  plot =  pl2, width = 20, height = 20, units = "cm")
    ggsave(paste0("Response_Vcormat_lowCorr_param", param, "_model_", model_type, "_diffStatic.pdf"),  plot =  pl1_diff, width = 20, height = 20, units = "cm")
    ggsave(paste0("Response_Vcormat_HighCorr_param", param, "_model_", model_type, "_diffStatic.pdf"),  plot =  pl2_diff, width = 20, height = 20, units = "cm")
  }
  return(NULL)
}




plot_Heatmap2 <- function(data, flag_res = TRUE, model_type = c("reduced", "full"),
                          neff_reduced = NULL, grid_length = 5, grid_d = grid_d, d = 24){
  model_type <- match.arg(model_type)

  setwd(root_dir)
  setwd("content/Section6/Results")

  if(flag_res == TRUE){
    load(paste0("fit_AllData_parammcd_d_", d, "_model_Static_residuals.RData"))
    load(paste0("fit_AllData_parammcd_d_", d, "_model_", model_type, "_residuals.RData"))
  } else {
    load(paste0("fit_AllData_parammcd_d_", d, "_model_Static_response.RData"))
    load(paste0("fit_AllData_parammcd_d_", d, "_model_", model_type, "_response.RData"))
  }


  Sigma_pred_fixed_MCD <- predict(res_fixed, type = "response")
  Sigma_predMat_fixed_MCD <-  Sigma_mat(Sigma_pred_fixed_MCD[,-c(1 : d)])

  Sigma_pred_model_MCD <- predict(res_final, type = "response")
  Sigma_predMat_model_MCD <-  Sigma_mat(Sigma_pred_model_MCD[,-c(1 : d)])

  idx_min_Corr_MCD <- which.min(unlist(lapply(1 : length(Sigma_predMat_model_MCD), function(x) unlist(mean(Sigma_predMat_model_MCD[[x]][lower.tri(Sigma_predMat_model_MCD[[x]])])))))
  idx_max_Corr_MCD <- which.max(unlist(lapply(1 : length(Sigma_predMat_model_MCD), function(x) unlist(mean(Sigma_predMat_model_MCD[[x]][lower.tri(Sigma_predMat_model_MCD[[x]])])))))

  if(flag_res == TRUE){
    load(paste0("fit_AllData_paramlogm_d_", d, "_model_Static_residuals.RData"))
    load(paste0("fit_AllData_paramlogm_d_", d, "_model_", model_type, "_residuals.RData"))
  } else {
    load(paste0("fit_AllData_paramlogm_d_", d, "_model_Static_response.RData"))
    load(paste0("fit_AllData_paramlogm_d_", d, "_model_", model_type, "_response.RData"))
  }

  Sigma_pred_fixed_logM <- predict(res_fixed, type = "response")
  Sigma_predMat_fixed_logM <-  Sigma_mat(Sigma_pred_fixed_logM[,-c(1 : d)])

  Sigma_pred_model_logM <- predict(res_final, type = "response")
  Sigma_predMat_model_logM <-  Sigma_mat(Sigma_pred_model_logM[,-c(1 : d)])

  col_cor <- rev(colorspace::sequential_hcl(palette = "Blues 3", n = 100))
  col_var <- rev(colorspace::sequential_hcl(palette = "Red", n = 10)[1 : 5])

  label_xaxis <- c(paste0("h0", 0:9), paste0("h", 10:(d-1)))
  label_yaxis <- c(paste0("h", (d-1):10), paste0("h0", 9:0))
  days <- c("Sunday", "Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday")

  setwd(root_dir)
  setwd("content/Section6/Plots")


  # Min corr using the model predictions
  pl1_MCD <- heatmap_FitCov(Sigma_predMat_model_MCD[[idx_min_Corr_MCD]], d = d, range_var = NULL, range_corr = c(0, 1), label = "h",
                            label_xaxis = label_xaxis, label_yaxis = label_yaxis,  col_var = col_var, col_cor = col_cor)

  date <- as.Date(paste(data[idx_min_Corr_MCD, "year"], data[idx_min_Corr_MCD, "doy"]), format = "%Y %j")
  dow <- days[wday(date)]
  pl1_MCD <- pl1_MCD + annotate("text",  x = (d-1) + 0.25, y = (d-1) + 0.25, label = paste(dow, date), vjust=1, hjust=1, cex = 7)
  print(pl1_MCD)

  pl1_logM <- heatmap_FitCov(Sigma_predMat_model_logM[[idx_min_Corr_MCD]], d = d, range_var = NULL, range_corr = c(0, 1), label = "h",
                             label_xaxis = label_xaxis, label_yaxis = label_yaxis,  col_var = col_var, col_cor = col_cor)

  date <- as.Date(paste(data[idx_min_Corr_MCD, "year"], data[idx_min_Corr_MCD, "doy"]), format = "%Y %j")
  dow <- days[wday(date)]
  pl1_logM <- pl1_logM + annotate("text",  x = (d-1) + 0.25, y = (d-1) + 0.25, label = paste(dow, date), vjust=1, hjust=1, cex = 7)
  print(pl1_logM)


  # Max corr using the model predictions
  pl2_MCD <- heatmap_FitCov(Sigma_predMat_model_MCD[[idx_max_Corr_MCD]], d = d, range_var = NULL, range_corr = c(0, 1), label = "h",
                            label_xaxis = label_xaxis, label_yaxis = label_yaxis,  col_var = col_var, col_cor = col_cor)
  date <- as.Date(paste(data[idx_max_Corr_MCD,"year"], data[idx_max_Corr_MCD,"doy"]), format = "%Y %j")
  dow <- days[wday(date)]
  pl2_MCD <- pl2_MCD +  annotate("text",  x = (d-1) + 0.25, y = (d-1) + 0.25, label = paste(dow, date), vjust = 1, hjust = 1, cex = 7)
  print(pl2_MCD)

  pl2_logM <- heatmap_FitCov(Sigma_predMat_model_logM[[idx_max_Corr_MCD]], d = d, range_var = NULL, range_corr = c(0, 1), label = "h",
                             label_xaxis = label_xaxis, label_yaxis = label_yaxis,  col_var = col_var, col_cor = col_cor)
  date <- as.Date(paste(data[idx_max_Corr_MCD,"year"], data[idx_max_Corr_MCD,"doy"]), format = "%Y %j")
  dow <- days[wday(date)]
  pl2_logM <- pl2_logM +  annotate("text",  x = (d-1) + 0.25, y = (d-1) + 0.25, label = paste(dow, date), vjust = 1, hjust = 1, cex = 7)
  print(pl2_logM)


  ############################################################################################################################
  col_var2 <- rev(colorspace::diverging_hcl(palette = "Blue - Red 2", n = 100))
  col_cor2 <- rev(colorspace::divergingx_hcl(palette = "RdBu", n = 100))[100:1]

  minV <- min((sqrt(diag(Sigma_predMat_model_MCD[[idx_min_Corr_MCD]])) - sqrt(diag(Sigma_predMat_fixed_MCD[[1]]))))
  maxV <- max((sqrt(diag(Sigma_predMat_model_MCD[[idx_min_Corr_MCD]])) - sqrt(diag(Sigma_predMat_fixed_MCD[[1]]))))

  lowTri <- lower.tri(Sigma_predMat_fixed_MCD[[1]])
  minC <- min((Sigma_predMat_model_MCD[[idx_min_Corr_MCD]] - Sigma_predMat_fixed_MCD[[1]])[lowTri])
  maxC <- max((Sigma_predMat_model_MCD[[idx_min_Corr_MCD]] - Sigma_predMat_fixed_MCD[[1]])[lowTri])

  pl1_diff_MCD <- heatmap_FitCov(Sigma_predMat_model_MCD[[idx_min_Corr_MCD]], PredCov2 = Sigma_predMat_fixed_MCD[[1]],
                                 d = d, range_var = c(minV,maxV), range_corr = c(minC, maxC), label = "h",
                                 label_xaxis = label_xaxis, label_yaxis = label_yaxis,  col_var = col_var2, col_cor = col_cor2)

  date <- as.Date(paste(data[idx_min_Corr_MCD, "year"], data[idx_min_Corr_MCD, "doy"]), format = "%Y %j")
  dow <- days[wday(date)]
  pl1_diff_MCD  <- pl1_diff_MCD  + annotate("text",  x = (d-1) + 0.25, y = (d-1) + 0.25, label = paste(dow, date), vjust=1, hjust=1, cex = 7)
  print(pl1_diff_MCD )


  minV <- min((sqrt(diag(Sigma_predMat_model_logM[[idx_min_Corr_MCD]])) - sqrt(diag(Sigma_predMat_fixed_logM[[1]]))))
  maxV <- max((sqrt(diag(Sigma_predMat_model_logM[[idx_min_Corr_MCD]])) - sqrt(diag(Sigma_predMat_fixed_logM[[1]]))))

  lowTri <- lower.tri(Sigma_predMat_fixed_logM[[1]])
  minC <- min((Sigma_predMat_model_logM[[idx_min_Corr_MCD]] - Sigma_predMat_fixed_logM[[1]])[lowTri])
  maxC <- max((Sigma_predMat_model_logM[[idx_min_Corr_MCD]] - Sigma_predMat_fixed_logM[[1]])[lowTri])

  pl1_diff_logM <- heatmap_FitCov(Sigma_predMat_model_logM[[idx_min_Corr_MCD]], PredCov2 = Sigma_predMat_fixed_logM[[1]],
                                  d = d, range_var = c(minV,maxV), range_corr = c(minC, maxC), label = "h",
                                  label_xaxis = label_xaxis, label_yaxis = label_yaxis,  col_var = col_var2, col_cor = col_cor2)

  date <- as.Date(paste(data[idx_min_Corr_MCD, "year"], data[idx_min_Corr_MCD, "doy"]), format = "%Y %j")
  dow <- days[wday(date)]
  pl1_diff_logM  <- pl1_diff_logM  + annotate("text",  x = (d-1) + 0.25, y = (d-1) + 0.25, label = paste(dow, date), vjust=1, hjust=1, cex = 7)
  print(pl1_diff_logM )

  ####################################
  minV <- min((sqrt(diag(Sigma_predMat_model_MCD[[idx_max_Corr_MCD]])) - sqrt(diag(Sigma_predMat_fixed_MCD[[1]]))))
  maxV <- max((sqrt(diag(Sigma_predMat_model_MCD[[idx_max_Corr_MCD]])) - sqrt(diag(Sigma_predMat_fixed_MCD[[1]]))))


  lowTri <- lower.tri(Sigma_predMat_fixed_MCD[[1]])
  minC <- min((Sigma_predMat_model_MCD[[idx_max_Corr_MCD]] - Sigma_predMat_fixed_MCD[[1]])[lowTri])
  maxC <- max((Sigma_predMat_model_MCD[[idx_max_Corr_MCD]] - Sigma_predMat_fixed_MCD[[1]])[lowTri])

  pl2_diff_MCD  <- heatmap_FitCov(Sigma_predMat_model_MCD[[idx_max_Corr_MCD]], PredCov2 = Sigma_predMat_fixed_MCD[[1]],
                                  d = d, range_var = c(minV,maxV), range_corr = c(minC, maxC), label = "h",
                                  label_xaxis = label_xaxis, label_yaxis = label_yaxis,  col_var = col_var2, col_cor = col_cor2)

  date <- as.Date(paste(data[idx_max_Corr_MCD,"year"], data[idx_max_Corr_MCD,"doy"]), format = "%Y %j")
  dow <- days[wday(date)]
  pl2_diff_MCD <- pl2_diff_MCD +  annotate("text",  x = (d-1) + 0.25, y = (d-1) + 0.25, label = paste(dow, date), vjust = 1, hjust = 1, cex = 7)
  print(pl2_diff_MCD)

  lowTri <- lower.tri(Sigma_predMat_fixed_logM[[1]])
  minC <- min((Sigma_predMat_model_logM[[idx_max_Corr_MCD]] - Sigma_predMat_fixed_logM[[1]])[lowTri])
  maxC <- max((Sigma_predMat_model_logM[[idx_max_Corr_MCD]] - Sigma_predMat_fixed_logM[[1]])[lowTri])

  pl2_diff_logM  <- heatmap_FitCov(Sigma_predMat_model_logM[[idx_max_Corr_MCD]], PredCov2 = Sigma_predMat_fixed_logM[[1]],
                                   d = d, range_var = c(minV,maxV), range_corr = c(minC, maxC), label = "h",
                                   label_xaxis = label_xaxis, label_yaxis = label_yaxis,  col_var = col_var2, col_cor = col_cor2)

  date <- as.Date(paste(data[idx_max_Corr_MCD,"year"], data[idx_max_Corr_MCD,"doy"]), format = "%Y %j")
  dow <- days[wday(date)]
  pl2_diff_logM <- pl2_diff_logM +  annotate("text",  x = (d-1) + 0.25, y = (d-1) + 0.25, label = paste(dow, date), vjust = 1, hjust = 1, cex = 7)
  print(pl2_diff_logM)

  plot1 <- ggarrange(pl1_diff_MCD + ggtitle("MCD"), pl1_diff_logM + ggtitle("logM"), nrow = 1)
  plot2 <- ggarrange(pl2_diff_MCD + ggtitle("MCD"), pl2_diff_logM + ggtitle("logM"), nrow = 1)


  if(flag_res == TRUE){
    # ggsave(paste0("Residuals_Vcormat_param", param, "_model_fixed.pdf"),  plot =  pl1_fixed, width = 20, height = 20, units = "cm")
    # ggsave(paste0("Residuals_Vcormat_lowCorr_param", param, "_model_", model_type, ".pdf"),  plot =  pl1, width = 20, height = 20, units = "cm")
    # ggsave(paste0("Residuals_Vcormat_HighCorr_param", param, "_model_", model_type, ".pdf"),  plot =  pl2, width = 20, height = 20, units = "cm")
    ggsave(paste0("Residuals_Vcormat_lowCorr_param_MCDvslogM_model_", model_type, "_diffStatic_FIXDAY.pdf"),  plot = plot1 , width = 40, height = 20, units = "cm")
    ggsave(paste0("Residuals_Vcormat_lowCorr_param_MCDvslogM_model_", model_type, "_diffStatic_FIXDAY.pdf"),  plot = plot2, width = 40, height = 20, units = "cm")

  } else {
    # ggsave(paste0("Response_Vcormat_param", param, "_model_fixed.pdf"),  plot =  pl1_fixed, width = 20, height = 20, units = "cm")
    # ggsave(paste0("Response_Vcormat_lowCorr_param", param, "_model_", model_type, ".pdf"),  plot =  pl1, width = 20, height = 20, units = "cm")
    # ggsave(paste0("Response_Vcormat_HighCorr_param", param, "_model_", model_type, ".pdf"),  plot =  pl2, width = 20, height = 20, units = "cm")
    ggsave(paste0("Response_Vcormat_lowCorr_param_MCDvslogM_model_", model_type, "_diffStatic_FIXDAY.pdf"),  plot =  plot1, width =40, height = 20, units = "cm")
    ggsave(paste0("Response_Vcormat_lowCorr_param_MCDvslogM_model_", model_type, "_diffStatic_FIXDAY.pdf"),  plot =  plot2, width = 40, height = 20, units = "cm")

  }
  return(NULL)
}

plot_heat_and_traj <- function(mu_sigma, yobs, idx1, idx2, dat, nsim){

  # Extract dates
  date1 <- format(as.Date(paste(dat[idx1, "doy"], dat[idx1, "year"]), format = "%j %Y"), format="%d-%m-%Y")
  date2 <- format(as.Date(paste(dat[idx2, "doy"], dat[idx2, "year"]), format = "%j %Y"), format="%d-%m-%Y")

  # Convert to MW
  d <- ncol(yobs)
  conv_fact <- 0.1
  mu_sigma[ , 1:d] <- mu_sigma[ , 1:d] * conv_fact
  mu_sigma[ , (d+1):(2*d)] <- mu_sigma[ , (d+1):(2*d)] * conv_fact^2
  yobs <- yobs * conv_fact

  # Get covariance matrices: note SIGMA decomposed in SD and correlation hence these matrices
  # are not positive definite but that's OK here!!
  Sigma_list <-  Sigma_mat(mu_sigma[,-c(1 : d)])

  # Simulate trajectories for each date
  set.seed(5151)
  COV <- Sigma_list[[idx1]]
  sdev <- sqrt(diag(COV))
  COV <- diag(sdev) %*% COV %*% diag(sdev)
  diag(COV) <- sdev^2
  X <- rmvn(nsim, rep(0, 24), COV)
  my_dat <- data.frame("hour" = rep(1:24, nsim), "y" = as.vector(t(X)))

  set.seed(5151)
  COV <- Sigma_list[[idx2]]
  sdev <- sqrt(diag(COV))
  COV <- diag(sdev) %*% COV %*% diag(sdev)
  diag(COV) <- sdev^2
  X2 <- rmvn(nsim, rep(0, 24), COV)

  # Created data frame for ggplot
  my_dat <- rbind(my_dat,
                  data.frame("hour" = rep(1:24, nsim), "y" = as.vector(t(X2))))
  my_dat$day <- rep(c(1, 2), each = 24 * nsim)
  my_dat$ID <- rep(1:(2*nsim), each = 24)
  my_dat$Day <- as.character(my_dat$day)
  my_dat$Day[my_dat$day == "1"] <- date1
  my_dat$Day[my_dat$day == "2"] <- date2
  my_dat$Day <- as.factor(my_dat$Day)

  # Heatmap for first date
  plt <- list()
  col_cor <- rev(colorspace::sequential_hcl(palette = "Blues 3", n = 100))
  col_var <- rev(colorspace::sequential_hcl(palette = "Red", n = 100)[1:90])
  label_xaxis <- c(paste0(0:9), paste0(10:(d-1))) #c(paste0("h0", 0:9), paste0("h", 10:(d-1)))
  label_yaxis <- c(paste0((d-1):10), paste0(9:0))
  days <- c("Sunday", "Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday")
  range_var <- NULL #range(c(sqrt(diag(Sigma_list[[idx1]])), sqrt(diag(Sigma_list[[idx2]]))))
  pl1_MCD <- heatmap_FitCov(Sigma_list[[idx1]], d = d, range_var = range_var, range_corr = c(0, 1), label = "h",
                            label_xaxis = label_xaxis, label_yaxis = label_yaxis,  col_var = col_var, col_cor = col_cor)
  pl1_MCD <- pl1_MCD +
    annotate("text",  x = (d-1) + 0.25, y = (d-1) + 0.25, label = date1, vjust=1, hjust=1, cex = 7) +
    theme(legend.position=c(.9,.55))
  plt[[1]] <- pl1_MCD

  # Heatmap for 2nd date
  pl1_MCD <- heatmap_FitCov(Sigma_list[[idx2]], d = d, range_var = range_var, range_corr = c(0, 1), label = "h",
                            label_xaxis = label_xaxis, label_yaxis = label_yaxis,  col_var = col_var, col_cor = col_cor)
  pl1_MCD <- pl1_MCD + annotate("text",  x = (d-1) + 0.25, y = (d-1) + 0.25, label = date2,
                                vjust=1, hjust=1, cex = 7) +
    theme(legend.position=c(.9,.55))
  plt[[2]] <- pl1_MCD

  # Trajectories plot
  my_dat$hour <- my_dat$hour - 1
  cols <- c("blue", "red")
  names(cols) <- as.character(c(date1, date2))
  plt[[3]] <-  ggplot(my_dat, mapping = aes(y = y, x = hour, group = ID, colour=Day, linetype = Day)) +
    geom_line(alpha = 0.3) +
    scale_color_manual(values=cols) +
    theme_bw() +
    geom_line(data = data.frame(hour = 0:23, y = unlist(yobs[idx1, ] - mu_sigma[idx1, 1:24])),
              aes(y = y, x = hour), inherit.aes = FALSE, colour = "blue", linewidth = 1.5) +
    geom_line(data = data.frame(hour = 0:23, y = unlist(yobs[idx2, ] - mu_sigma[idx2, 1:24])),
              aes(y = y, x = hour), inherit.aes = FALSE, colour = "red", linewidth = 1.5) +
    ylab("Residuals (GW)") + xlab("Hour") +
    theme(legend.position = c(.2,.1)) +
    guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.5, linetype = 1))) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.text.x = element_text(size = 20, vjust = 0.5),
          axis.text.y = element_text(size = 20),
          text = element_text(size = 23),
          legend.text=element_text(size=20), strip.text.x = element_text(size = 23)) #+

  return( plt )

}


get_data_4_heat <- function(cv_results, neff, tot_eff, grid_step, param, d = 24){

  index <- (tot_eff - neff) / grid_step + 1

  if( sum(apply(cv_results[[1]][[1]][[index]]$lpi_pred_out, 2, sd) > 1e-6) != (neff + d) ){
    stop("Wrong number of modelled effects")
  }

  eta <- do.call("rbind", sapply(cv_results, function(o) o[[1]][[index]]$lpi_pred_out))
  mu <- eta * NA

  if(param == "mcd") pred_fun <- internal()$pred_mcd
  if(param == "logm") pred_fun <- internal()$pred_logm

  pred_fun(eta = eta, pred = mu, d = d, cor_flag = TRUE)

  return(mu)
}
