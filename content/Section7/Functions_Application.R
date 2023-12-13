################################################################################
# Needed functions for fitting the models of Section  7: GEFCom14 Application  #
################################################################################


#######################################################################
# Function for extracting the string of Theta elements to be selected #
#######################################################################
get_substr_rhs_lpi_Vcov <- function(name_Vcov_eff){
  string_eff <- rep(NA, length(name_Vcov_eff))
  for(i in 1 : length(name_Vcov_eff)){
    det_Th <- gregexpr("Th_", name_Vcov_eff[i])[[1]][1]
    if(det_Th != -1){
      string_eff[i] <- substring(name_Vcov_eff[i], first = det_Th)
    }
  }
  return(string_eff)
}

################################################################
# Function for extracting the string of effects to be selected #
################################################################
get_substr_lhs_lpi_Vcov <- function(name_Vcov_eff = NULL){
  string_eff <- rep(NA, length(name_Vcov_eff))
  for(i in 1 : length(name_Vcov_eff)){
    det_Th <- gregexpr("Th_", name_Vcov_eff[i])[[1]][1]
    if(det_Th != -1){
      string_eff[i] <- substring(name_Vcov_eff[i], first =   1, last = det_Th - 2)
    }
  }
  return(string_eff)
}

#################################################################################
# Function for extracting the position indices of Theta elements to be selected #
#################################################################################
get_idx_lpi_Vcov <- function(string_eff){
  idx <- which(grepl("Th_", string_eff))
  return(idx)
}

################################################################################
# Functions for checking if an object is not an integer(0) or not a logical(0) #
################################################################################
is.integer0 <- function(x){
  is.integer(x) && length(x) == 0L
}

is.logical0 <- function(x){
  is.logical(x) && length(x) == 0L
}

#################################################################################################################################################################
# Function for building the model formula:
# It starts from a fitted model, get the summary and order the  the smooth (and linear) effects
# involved in covariance matrix modelling according to the p-value, build the model formula for theta and append such formula to the mean model
#################################################################################################################################################################
mod_foo_building <- function(summary_Peff, summary_Seff, mean_model_formula,  K_covmod, grid_length, metric){

  #!!! Problem: if set a different number of basis you cannot extract the number of basis from the summary print (this is related to the mgcv summary function)
  # Thus, maybe it should be better to consider the default k=10

  name_Peff <- row.names(summary_Peff)
  name_Seff <- row.names(summary_Seff)

  string_Peff <- get_substr_rhs_lpi_Vcov(name_Peff)
  string_Seff <- get_substr_rhs_lpi_Vcov(name_Seff)

  idx_Peff <- get_idx_lpi_Vcov(string_Peff)
  idx_Seff <- get_idx_lpi_Vcov(string_Seff)

  mat_eff <- data.frame(matrix(0, nrow = length(idx_Seff) + length(idx_Peff), 7))
  names(mat_eff) <- c("Eff_Lpi", "StTest", "pvalue", "lpvalue", "rank_ST", "rank_p", "rank_lp")
  if(!(is.integer0(idx_Peff) | is.logical0(idx_Peff)) & (!is.integer0(idx_Seff) | !is.logical0(idx_Seff))){
    mat_eff[, 1] <- c(name_Peff[idx_Peff], name_Seff[idx_Seff])
    mat_eff[, 2] <- c(summary_Peff[idx_Peff, 3], summary_Seff[idx_Seff, 3])
    mat_eff[, 3] <- c(summary_Peff[idx_Peff, 4], summary_Seff[idx_Seff, 4])
    mat_eff[, 4] <- c(log(abs(summary_Peff[idx_Peff, 4])),log(abs(summary_Seff[idx_Seff, 4])))
  } else if((is.integer0(idx_Peff) | is.logical0(idx_Peff)) & (!is.integer0(idx_Seff) | !is.logical0(idx_Seff))){
    mat_eff[, 1] <- name_Seff[idx_Seff]
    mat_eff[, 2] <- summary_Seff[idx_Seff, 3]
    mat_eff[, 3] <- summary_Seff[idx_Seff, 4]
    mat_eff[, 4] <- log(abs(summary_Seff[idx_Seff, 4]))
  } else if((!is.integer0(idx_Peff) | !is.logical0(idx_Peff)) & (is.integer0(idx_Seff) | is.logical0(idx_Seff))){
    mat_eff[, 1] <- name_Peff[idx_Peff]
    mat_eff[, 2] <- summary_Peff[idx_Peff, 3]
    mat_eff[, 3] <- summary_Peff[idx_Peff, 4]
    mat_eff[, 4] <- log(abs(summary_Peff[idx_Peff, 4]))
  } else {
    stop("Nothing to do")
  }

  if(metric == "p"){
    order.pvalue <- order(mat_eff$pvalue, mat_eff$Eff_Lpi)
    mat_eff$rank_p[order.pvalue] <- 1 : nrow(mat_eff)
    get_eff <- mat_eff[!(mat_eff$rank_p %in% ((K_covmod - grid_length + 1) : K_covmod)),]
  } else if(metric == "logp"){
    order.logpvalue <- order(mat_eff$lpvalue, mat_eff$Eff_Lpi)
    mat_eff$rank_lp[order.logpvalue] <- 1 : nrow(mat_eff)
    get_eff <- mat_eff[!(mat_eff$rank_lp %in% ((K_covmod - grid_length + 1) : K_covmod)),]
  } else if(metric == "ST"){
    order.ST <- order(mat_eff$StTest, mat_eff$Eff_Lpi, decreasing = TRUE)
    mat_eff$rank_ST[order.ST] <- 1 : nrow(mat_eff)
    get_eff <- mat_eff[!(mat_eff$rank_ST %in% ((K_covmod - grid_length + 1) : K_covmod)),]
  }

  get_eff$Effect <- get_substr_lhs_lpi_Vcov(get_eff[, 1])
  get_eff$Lpi <- get_substr_rhs_lpi_Vcov(get_eff[, 1])
  bool_duplicated <- duplicated(get_eff$Lpi)


  theta_foo <- list()
  count <- 1
  for(j in 1 : nrow(get_eff)){
    if(bool_duplicated[j] == FALSE){
      theta_foo[[count]] <- as.formula(paste0(get_eff$Lpi[j], "~",  get_eff$Effect[j]))
      count <- count + 1
    } else {
      string_dupl_lpi <- get_eff$Lpi[j]
      idx_dup <- unlist(lapply(1 : (count - 1), function(.x)  grepl(string_dupl_lpi, theta_foo[[.x]])[2]))
      theta_foo[[which(idx_dup)]] <- as.formula(paste0(deparse(theta_foo[[idx_dup]]), "+", get_eff$Effect[j]))
    }
  }

  global_formula <- c(mean_model_formula,  theta_foo)
  return(global_formula)
}





######################################################################################################################################
# This function                                                                                                                      #
# -) At first, fit the full models                                                                                                   #
# -) Get the summary and order the  the smooth (and linear) effects involved in covariance matrix modelling according to the p-value #
# -) Remove some not useful objects (to save memory) from the fitted model                                                           #
# -) Then return the model formula that must be sequentially used for model fitting                                                  #
# -) Fit the models by removing grid_length effects for each step                                                                    #
######################################################################################################################################
# param: parametrisation
# d: dimension of the outcome
# grid_length: number of effects sequentially removed
# mean_model_formula: unchanged during the fitting workflow
# data: dataset

stepw_res <- function(param, d, grid_length, mean_model_formula, data, sets_eval, eff_vcov = NULL,
                      metric = c("p", "logp", "ST"), ncores = 1, save.gam = NULL){

  metric <- match.arg(metric)
  eff_vcov <- eff_vcov

  # Set the covariance matrix model formula for the full model
  theta_foo <- list()
  for(j in 1 : d){
    theta_foo[[j]] <- as.formula(paste0("Th_", j, ".", j, "~", eff_vcov))
  }
  count <- d + 1
  for(j in 2 : d){
    for(k in 1 : (j - 1)){
      theta_foo[[count]] <- as.formula(paste0("Th_", j, ".", k, "~", eff_vcov))
      count <- count + 1
    }
  }

  global_formula_full <- c(mean_formula,  theta_foo)

  eff_grid <- seq(d * (d + 1)/2, 0, by = -grid_length)



  res_sim <- function(param, d,  grid_length, eff_grid, sets_eval,
                      mean_model_formula, global_formula_full,
                      data, data_test, eff_vcov,  metric, save.gam){
    out <- list()
    #dss <- lapply(1 : length(eff_grid), function(jj){
    for(jj in 1: length(eff_grid)){
      out[[jj]] <- list()
      if(save.gam){
        if(jj == 1){
          #############################################################################################################################
          # Full model
          time <- microbenchmark(
            out[[jj]]$stepw_model <- gam_scm(global_formula_full,
                                       family = mvn_scm(d = d), optimizer = "efs",
                                       data = data,
                                       aGam = list(control = list(trace = FALSE))), times = 1L)
            out[[jj]]$time_fit <- time$time
            out[[jj]]$sum_Vcov_Peff <- summary(out[[jj]]$stepw_model, print = FALSE)$p.table # summary table for linear effects (if there are any)
            out[[jj]]$sum_Vcov_Seff <- summary(out[[jj]]$stepw_model, print = FALSE)$s.table # table of smooth effects
            # Reduce memory footprint by removing gam.object not useful for the procedure
            out[[jj]]$stepw_model$R <- out[[jj]]$stepw_model$H <- out[[jj]]$stepw_model$Vc <- out[[jj]]$stepw_model$Vp <- out[[jj]]$stepw_model$Ve <- out[[jj]]$stepw_model$lbb <- out[[jj]]$stepw_model$L <- out[[jj]]$stepw_model$St <- NULL
            #############################################################################################################################
        } else if (jj > 1 & jj < length(eff_grid)){
           K_covmod <- d * (d + 1)/2
           global_formula <- mod_foo_building(summary_Peff = out[[jj - 1]]$sum_Vcov_Peff,
                                              summary_Seff = out[[jj - 1]]$sum_Vcov_Seff,
                                              mean_model_formula, K_covmod, grid_length, metric)
           time <- microbenchmark(
             out[[jj]]$stepw_model <- gam_scm(global_formula,
                                        family = mvn_scm(d = d, param = param),
                                        optimizer = "efs",
                                        data = data,
                                        aGam = list(control = list(trace = FALSE))), times = 1L)
             out[[jj]]$time_fit <- time$time

             out[[jj]]$sum_Vcov_Peff <- summary(out[[jj]]$stepw_model, print = FALSE)$p.table # summary table for linear effects (if there are any)
             out[[jj]]$sum_Vcov_Seff <- summary(out[[jj]]$stepw_model, print = FALSE)$s.table

             out[[jj]]$stepw_model$R <- out[[jj]]$stepw_model$H <- out[[jj]]$stepw_model$Vc <- out[[jj]]$stepw_model$Vp <- out[[jj]]$stepw_model$Ve <- out[[jj]]$stepw_model$lbb <- out[[jj]]$stepw_model$L <- out[[jj]]$stepw_model$St <- NULL
            #############################################################################################################################
        } else if(jj == length(eff_grid)){
             global_formula <- mean_model_formula
             time <- microbenchmark(out[[jj]]$stepw_model <- gam_scm(global_formula,
                                                               family = mvn_scm(d = d, param = param),
                                                               optimizer = "efs",
                                                               data = data,
                                                               aGam = list(control = list(trace = FALSE))), times = 1L)
             out[[jj]]$time_fit <- time$time
             out[[jj]]$sum_Vcov_Peff <- summary(out[[jj]]$stepw_model, print = FALSE)$p.table # summary table for linear effects (if there are any)
             out[[jj]]$sum_Vcov_Seff <- summary(out[[jj]]$stepw_model, print = FALSE)$s.table
             # Reduce memory footprint by removing gam.object not useful for the procedure
             out[[jj]]$stepw_model$R <- out[[jj]]$stepw_model$H <-out[[jj]]$stepw_model$Vc <- out[[jj]]$stepw_model$Vp <- out[[jj]]$stepw_model$Ve <- out[[jj]]$stepw_model$lbb <- out[[jj]]$stepw_model$L <- out[[jj]]$stepw_model$St <- NULL
        }
      } else {
        if(jj == 1){
          #############################################################################################################################
          # Full model
          time <- microbenchmark(
            stepw_model <- gam_scm(global_formula_full,
                                   family = mvn_scm(d = d), optimizer = "efs",
                                   data = data,
                                   aGam = list(control = list(trace = FALSE))), times = 1L)
          out[[jj]]$time_fit <- time$time
          out[[jj]]$lpi_pred_in <- predict(stepw_model)
          out[[jj]]$lpi_pred_out <- predict(stepw_model, new.data = data_test)
          out[[jj]]$mformula <- stepw_model$formula
          out[[jj]]$sum_Vcov_Peff <- summary(stepw_model, print = FALSE)$p.table # summary table for linear effects (if there are any)
          out[[jj]]$sum_Vcov_Seff <- summary(stepw_model, print = FALSE)$s.table # table of smooth effects
          # Reduce memory footprint by removing gam.object not useful for the procedure
          stepw_model$R <- stepw_model$H <- stepw_model$Vc <- stepw_model$Vp <- stepw_model$Ve <- stepw_model$lbb <- stepw_model$L <- stepw_model$St <- NULL
          #############################################################################################################################
        } else if (jj > 1 & jj < length(eff_grid)){
          K_covmod <- d * (d + 1)/2
          global_formula <- mod_foo_building(summary_Peff = out[[jj - 1]]$sum_Vcov_Peff,
                                             summary_Seff = out[[jj - 1]]$sum_Vcov_Seff,
                                             mean_model_formula, K_covmod, grid_length, metric)

          time <- microbenchmark(
          stepw_model <- gam_scm(global_formula,
                                 family = mvn_scm(d = d, param = param),
                                 optimizer = "efs",
                                 data = data,
                                 aGam = list(control = list(trace = FALSE))), times = 1L)
          out[[jj]]$time_fit <- time$time
          out[[jj]]$lpi_pred_in <- predict(stepw_model)
          out[[jj]]$lpi_pred_out <- predict(stepw_model, new.data = data_test)
          out[[jj]]$mformula <- stepw_model$formula
          out[[jj]]$sum_Vcov_Peff <- summary(stepw_model, print = FALSE)$p.table # summary table for linear effects (if there are any)
          out[[jj]]$sum_Vcov_Seff <- summary(stepw_model, print = FALSE)$s.table

          stepw_model$R <- stepw_model$H <- stepw_model$Vc <- stepw_model$Vp <- stepw_model$Ve <- stepw_model$lbb <- stepw_model$L <- stepw_model$St <- NULL
          #############################################################################################################################
        } else if(jj == length(eff_grid)){
          global_formula <- mean_model_formula
          time <- microbenchmark(stepw_model <- gam_scm(global_formula,
                                                        family=mvn_scm(d = d, param = param),
                                                        optimizer = "efs",
                                                        data = data,
                                                        aGam = list(control = list(trace = FALSE))), times = 1L)
          out[[jj]]$time_fit <- time$time
          out[[jj]]$lpi_pred <- predict(stepw_model)
          out[[jj]]$lpi_pred_out <- predict(stepw_model, new.data = data_test)
          out[[jj]]$mformula <- stepw_model$formula
          out[[jj]]$sum_Vcov_Peff <- summary(stepw_model, print = FALSE)$p.table # summary table for linear effects (if there are any)
          out[[jj]]$sum_Vcov_Seff <- summary(stepw_model, print = FALSE)$s.table
          # Reduce memory footprint by removing gam.object not useful for the procedure
          stepw_model$R <- stepw_model$H <- stepw_model$Vc <- stepw_model$Vp <- stepw_model$Ve <- stepw_model$lbb <- stepw_model$L <- stepw_model$St <- NULL
        }
      }
    }
    return(out)
  }
  data_test <- data

  cl <- makePSOCKcluster(ncores)
  setDefaultCluster(cl)
  clusterExport(NULL, c("get_substr_rhs_lpi_Vcov", "get_substr_lhs_lpi_Vcov", "get_idx_lpi_Vcov",
                        "is.integer0", "is.logical0", "mod_foo_building", "save.gam",
                        "param", "d", "grid_length", "eff_grid", "mean_model_formula", "data", "data_test",
                        "sets_eval", "eff_vcov",
                        "global_formula_full", "metric", "res_sim"), envir = environment())

  clusterEvalQ(NULL, {
    library(SCM)
    library(microbenchmark)
    library(stringr)
    library(BMisc)
  })

  out_res_sim <- function(.x){
    out2 <- res_sim(param = param, d = d,  grid_length = grid_length, eff_grid = eff_grid, sets_eval = sets,
                    mean_model_formula = mean_model_formula, global_formula_full = global_formula_full,
                    data = data[1: sets_eval[.x], ], data_test =  data_test[(sets_eval[.x] + 1): sets_eval[.x + 1], ],
                    eff_vcov = eff_vcov,  metric = metric, save.gam = save.gam)
    return(list(gen = out2))
  }


  environment(out_res_sim) <- .GlobalEnv

  res <- list()
  res <- parLapply(NULL, 1 : (length(sets_eval) - 1), out_res_sim)

  stopCluster(cl)
  rm(cl)
  return(res)
}




