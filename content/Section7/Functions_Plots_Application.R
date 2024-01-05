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
  out <- data.frame(matrix(NA, length(name_eff) * (length(foo_vcov) - d), 6))        # NA for the element of Theta modelled via intercepts
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



