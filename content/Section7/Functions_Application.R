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
      string_eff[i] <-substring(name_Vcov_eff[i],
                                first =   det_Th)
    }
  }
  return(string_eff)
}

sourceCpp(code = "
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
void over_writediag(arma::mat& E, const arma::vec& rho, const arma::vec& ind, const int k_sp) {
  int kk = k_sp - 1;
  int dd;
  for (int i = 0; i < ind.n_elem; ++i) {
    dd = ind(i) - 1;
    E(dd, dd) = exp(rho(kk) * 0.5);
  }
}
")


################################################################
# Function for extracting the string of effects to be selected #
################################################################
get_substr_lhs_lpi_Vcov <- function(name_Vcov_eff = NULL){
  string_eff <- rep(NA, length(name_Vcov_eff))
  for(i in 1 : length(name_Vcov_eff)){
    det_Th <- gregexpr("Th_", name_Vcov_eff[i])[[1]][1]
    if(det_Th != -1){
      string_eff[i] <-substring(name_Vcov_eff[i],
                                first =   1, last = det_Th - 2)
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
mod_foo_building <- function(summary_Peff, summary_Seff, mean_model_formula,  K_covmod, grid_length, metric = metric){

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
    mat_eff[,1] <- c(name_Peff[idx_Peff], name_Seff[idx_Seff])
    mat_eff[,2] <- c(summary_Peff[idx_Peff ,3], summary_Seff[idx_Seff, 3])
    mat_eff[,3] <- c(summary_Peff[idx_Peff ,4], summary_Seff[idx_Seff, 4])
    mat_eff[,4] <- c(log(abs(summary_Peff[idx_Peff ,4])),log(abs(summary_Seff[idx_Seff, 4])))
  } else if((is.integer0(idx_Peff) | is.logical0(idx_Peff)) & (!is.integer0(idx_Seff) | !is.logical0(idx_Seff))){
    mat_eff[,1] <- name_Seff[idx_Seff]
    mat_eff[,2] <- summary_Seff[idx_Seff, 3]
    mat_eff[,3] <- summary_Seff[idx_Seff, 4]
    mat_eff[,4] <- log(abs(summary_Seff[idx_Seff, 4]))
  } else if((!is.integer0(idx_Peff) | !is.logical0(idx_Peff)) & (is.integer0(idx_Seff) | is.logical0(idx_Seff))){
    mat_eff[,1] <- name_Peff[idx_Peff]
    mat_eff[,2] <- summary_Peff[idx_Peff, 3]
    mat_eff[,3] <- summary_Peff[idx_Peff, 4]
    mat_eff[,4] <- log(abs(summary_Peff[idx_Peff, 4]))
  } else {
    stop("Nothing to do")
  }


  if(metric == "p"){
    order.pvalue <- order(mat_eff$pvalue, mat_eff$Eff_Lpi)
    numb_p0 <- sum(mat_eff$pvalue <= 0)
    mat_eff$rank_p[order.pvalue] <- 1 : nrow(mat_eff)
    get_eff <- mat_eff[!(mat_eff$rank_p %in% ((K_covmod - grid_length + 1) : K_covmod)),]
  } else if(metric == "logp"){
    order.logpvalue <- order(mat_eff$lpvalue, mat_eff$Eff_Lpi)
    numb_p0 <- sum(mat_eff$pvalue <= 0)
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
    if( bool_duplicated[j] == FALSE){
      theta_foo[[count]] <- as.formula(paste0(get_eff$Lpi[j], "~",  get_eff$Effect[j]))
      count <- count + 1
    } else {
      string_dupl_lpi <- get_eff$Lpi[j]
      idx_dup <- unlist(lapply(1 : (count - 1), function(.x)  grepl(string_dupl_lpi, theta_foo[[.x]])[2]))
      theta_foo[[which(idx_dup)]] <- as.formula(paste0(deparse(theta_foo[[idx_dup]]), "+", get_eff$Effect[j]))
    }
  }

  global_formula <- c(mean_model_formula,  theta_foo)
  return(list(global_formula = global_formula, numb_p0 = numb_p0))
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
# data_train: dataset (only train set where we perform model selection)
# eff_vcov: name of the effects to be used for "covariance matrix" model formula
# metric: one among p (p-values), logp ("log p-values"), ST ("test staistic")
stepw_res <- function(param, d, grid_length, mean_model_formula, data_train, eff_vcov = NULL,
                      metric = c("p", "logp", "ST"), save.gam = NULL){

  metric <- match.arg(metric)
  if(is.null(save.gam)){
    save.gam <- FALSE
  } else{
    save.gam <- save.gam
  }

  if(save.gam){
    stepw_model <- list()      # list of fitted model
  } else {
    model_formula <- list()
  }

  time_fit <- list()         #list of model fitting times
  num_iter <- list()
  sum_Vcov_Peff <- list()
  sum_Vcov_Seff <- list()
  betas <- list()
  sps <- list()

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

  global_formula_full <- c(mean_model_formula,  theta_foo)
  K_covmod <- d * (d + 1)/2

  if(K_covmod %% grid_length != 0){
    stop("Consider for grid_length a divisor of d(d+1)/2")
  }

  #First model is the full model
  if(save.gam){
    tmp <- gam_scm(global_formula_full,
                   family = mvn_scm(d = d, param = param), optimizer = "efs",
                   data = data_train,
                   aGam = list(fit = FALSE))
    time <- microbenchmark(
      stepw_model[[1]] <- gam_scm(optimizer = "efs",
                                  aGam = list(G = tmp, control = list(trace = FALSE))),
      times = 1L)
    num_iter[[1]] <-  stepw_model[[1]]$family$getNC() - 1
    sum_Vcov_Peff[[1]] <- summary(stepw_model[[1]], print = FALSE)$p.table # summary table for linear effects (if there are any)
    sum_Vcov_Seff[[1]] <- summary(stepw_model[[1]], print = FALSE)$s.table # table of smooth effects
    stepw_model[[1]]$R <- stepw_model[[1]]$H <- stepw_model[[1]]$Vc <- stepw_model[[1]]$Vp <- stepw_model[[1]]$Ve <- stepw_model[[1]]$lbb <- stepw_model[[1]]$L <- stepw_model[[1]]$St <- NULL
  } else {
    model_formula[[1]] <- global_formula_full
    tmp <- gam_scm(global_formula_full,
                   family = mvn_scm(d = d, param = param), optimizer = "efs",
                   data = data_train,
                   aGam = list(fit = FALSE))
    time <- microbenchmark(
      stepw_model <- gam_scm(optimizer = "efs",
                             aGam = list(G = tmp, control = list(trace = FALSE))), times = 1L)
    num_iter[[1]] <-  stepw_model$family$getNC() - 1
    sum_Vcov_Peff[[1]] <- summary(stepw_model, print = FALSE)$p.table # summary table for linear effects (if there are any)
    sum_Vcov_Seff[[1]] <- summary(stepw_model, print = FALSE)$s.table # table of smooth effects
    betas[[1]] <- stepw_model$coefficients
    sps[[1]] <- stepw_model$sp
    stepw_model$R <- stepw_model$H <- stepw_model$Vc <- stepw_model$Vp <- stepw_model$Ve <- stepw_model$lbb <- stepw_model$L <- stepw_model$St <- NULL
  }
  time_fit[[1]] <- time$time
  np0 <- c(-99, rep(0, K_covmod/grid_length - 1), -99) # Note that we are excluding the full and the null covariance matrix model


  # if(save.gam){
  #   tmp_save <- list(fit = stepw_model, time_fit = time_fit, summary_Peff = sum_Vcov_Peff, summary_Seff = sum_Vcov_Seff, np0 = np0, betas = betas, sps = sps)
  #   save(file = "tmp_save.RData", tmp_save)
  # } else {
  #   tmp_save <- list(foo = model_formula, time_fit = time_fit, summary_Peff = sum_Vcov_Peff, summary_Seff = sum_Vcov_Seff, np0 = np0, betas = betas, sps = sps)
  #   save(file = "tmp_save.RData", tmp_save)
  # }



  # Build the model formula of the first nested model

  formula_func <- mod_foo_building(summary_Peff = sum_Vcov_Peff[[1]], summary_Seff = sum_Vcov_Seff[[1]],
                                   mean_model_formula, K_covmod, grid_length, metric)
  global_formula_old <- global_formula_full
  global_formula <- formula_func$global_formula


  old_nams <- sapply(global_formula_old, function(x) as.character(x[2]))
  full_nams <- old_nams
  new_nams <- sapply(formula_func$global_formula, function(x) as.character(x[2]))
  rem_ind_lp <- which(!(old_nams %in% new_nams))
  rem_ind_beta <- unlist(sapply(rem_ind_lp, function(.kk) attr(stepw_model$formula, "lpi")[[.kk]][-1]))
  beta_init <- betas[[1]][-rem_ind_beta]
  rem_smooth <- paste0("s.", rem_ind_lp-1, "(")
  rem_sp_ind <- sapply(rem_smooth, function(nam) which(grepl(nam, names(sps[[1]]), fixed = TRUE)))
  sp_init <- sps[[1]][-rem_sp_ind]

  dgrid <- seq((d * (d + 1)/2 - grid_length), grid_length, -grid_length)

  #K_covmod <- d * (d + 1)/2 - grid_length

  counter <- 2
  for(j in 1:length(dgrid)){
    #for(j in 1:3){

    np0[counter] <- formula_func$numb_p0
    if(save.gam){
      tmp <- gam_scm(global_formula,
                     family=mvn_scm(d = d, param = param),
                     optimizer = "efs",
                     data = data_train,
                     aGam = list(fit = FALSE))
      time <- microbenchmark(
        stepw_model[[counter]] <- gam_scm(optimizer = "efs",
                                          aGam = list(G = tmp, control = list(trace = FALSE)))
        , times=1L)
      num_iter[[counter]] <-  stepw_model[[counter]]$family$getNC() - 1
      sum_Vcov_Peff[[counter]] <- summary(stepw_model[[counter]], print = FALSE)$p.table # summary table for linear effects (if there are any)
      sum_Vcov_Seff[[counter]] <- summary(stepw_model[[counter]], print = FALSE)$s.table
      stepw_model[[counter]]$R <- stepw_model[[counter]]$H <- stepw_model[[counter]]$Vc <- stepw_model[[counter]]$Vp <- stepw_model[[counter]]$Ve <- stepw_model[[counter]]$lbb <- stepw_model[[counter]]$L <- stepw_model[[counter]]$St <- NULL
    } else {
      if(np0[counter] == (length(model_formula[[counter-1]])-d)){
        print("All the remaining smooth effects have p-value equal to zero")
        break
      } else if (np0[counter] > (length(model_formula[[counter-1]])-d -grid_length)){
        formula_func <- mod_foo_building(summary_Peff= sum_Vcov_Peff[[counter-1]], summary_Seff = sum_Vcov_Seff[[counter-1]],
                                         mean_model_formula, length(model_formula[[counter-1]]) - d, length(model_formula[[counter-1]]) - d - np0[counter], metric)
        global_formula <- formula_func$global_formula
        model_formula[[counter]] <- global_formula
        tmp <- gam_scm(global_formula,
                       family=mvn_scm(d = d, param = param),
                       optimizer = "efs",
                       data = data_train,
                       aGam = list(fit = FALSE))
        time <- microbenchmark(
          stepw_model <- gam_scm(optimizer = "efs",
                                 aGam = list(G = tmp, control = list(trace = FALSE)), in.out = list(sp = sp_init, scale = 1), start = beta_init), times=1L)
        num_iter[[counter]] <-  stepw_model$family$getNC() - 1
        sum_Vcov_Peff[[counter]] <- summary(stepw_model, print = FALSE)$p.table # summary table for linear effects (if there are any)
        sum_Vcov_Seff[[counter]] <- summary(stepw_model, print = FALSE)$s.table
        betas[[counter]] <- stepw_model$coefficients
        sps[[counter]] <- stepw_model$sp
        stepw_model$R <- stepw_model$H <- stepw_model$Vc <- stepw_model$Vp <- stepw_model$Ve <- stepw_model$lbb <- stepw_model$L <- stepw_model$St <- NULL
      } else {
        model_formula[[counter]] <- global_formula
        tmp <- gam_scm(global_formula,
                       family=mvn_scm(d = d, param = param),
                       optimizer = "efs",
                       data = data_train,
                       aGam = list(fit = FALSE))
        time <- microbenchmark(
          stepw_model <- gam_scm(optimizer = "efs",
                                 aGam = list(G = tmp, control = list(trace = FALSE), in.out = list(sp = sp_init, scale = 1), start = beta_init)), times=1L)
        par(mfrow = c(1, 2))
        plot(stepw_model$coef, beta_init)
        abline(0, 1, col = 2)
        plot(stepw_model$sp, sp_init)
        abline(0, 1, col = 2)
        num_iter[[counter]] <-  stepw_model$family$getNC() - 1
        sum_Vcov_Peff[[counter]] <- summary(stepw_model, print = FALSE)$p.table # summary table for linear effects (if there are any)
        sum_Vcov_Seff[[counter]] <- summary(stepw_model, print = FALSE)$s.table
        betas[[counter]] <- stepw_model$coefficients
        sps[[counter]] <- stepw_model$sp
        stepw_model$R <- stepw_model$H <- stepw_model$Vc <- stepw_model$Vp <- stepw_model$Ve <- stepw_model$lbb <- stepw_model$L <- stepw_model$St <- NULL
      }
    }
    time_fit[[counter]] <- time$time

    if(j == length(dgrid)){
      global_formula_old <- global_formula
      global_formula <- mean_model_formula
    } else {
      formula_func <- mod_foo_building(summary_Peff= sum_Vcov_Peff[[counter]], summary_Seff = sum_Vcov_Seff[[counter]],
                                       mean_model_formula, dgrid[j], grid_length, metric)  # Build the model formula of subsequent steps
      global_formula_old <- global_formula
      global_formula <- formula_func$global_formula
    }

    old_nams <- sapply(global_formula_old, function(x) as.character(x[2]))
    new_nams <- sapply(global_formula, function(x) as.character(x[2]))
    rem_ind_th <- which(!(old_nams %in% new_nams))
    rem_ind_lp <- which(full_nams %in% old_nams[rem_ind_th])
    rem_ind_beta <- unlist(sapply(rem_ind_lp, function(.kk) attr(stepw_model$formula, "lpi")[[.kk]][-1]))
    beta_init <- betas[[counter]][-rem_ind_beta]
    rem_smooth <- paste0("s.", rem_ind_lp-1, "(")
    rem_sp_ind <- sapply(rem_smooth, function(nam) which(grepl(nam, names(sps[[counter]]), fixed = TRUE)))
    sp_init <- sps[[counter]][-rem_sp_ind]
    counter <- counter + 1

    # if(save.gam){
    #   tmp_save <- list(fit = stepw_model, time_fit = time_fit, summary_Peff = sum_Vcov_Peff, summary_Seff = sum_Vcov_Seff, np0 = np0, betas = betas, sps = sps)
    #   save(file = "tmp_save.RData", tmp_save)
    # } else {
    #   tmp_save <- list(foo = model_formula, time_fit = time_fit, summary_Peff = sum_Vcov_Peff, summary_Seff = sum_Vcov_Seff, np0 = np0, betas = betas, sps = sps)
    #   save(file = "tmp_save.RData", tmp_save)
    # }
    par(mfrow = c(1, 2))
    plot(rev(unlist(time_fit)))
    plot(rev(unlist(time_fit)/unlist(num_iter)))
  }

  if(save.gam){
    time <- microbenchmark(
      stepw_model[[counter]] <- gam_scm(global_formula,
                                        family=mvn_scm(d = d, param = param),
                                        optimizer = "efs",
                                        data = data_train,
                                        aGam = list(control = list(trace = FALSE))), times=1L)
    num_iter[[counter]] <-  stepw_model[[counter]]$family$getNC() - 1
    sum_Vcov_Peff[[counter]] <- summary(stepw_model[[counter]], print = FALSE)$p.table # summary table for linear effects (if there are any)
    sum_Vcov_Seff[[counter]] <- summary(stepw_model[[counter]], print = FALSE)$s.table
    stepw_model[[counter]]$R <- stepw_model[[counter]]$H <- stepw_model[[counter]]$Vc <- stepw_model[[counter]]$Vp <- stepw_model[[counter]]$Ve <- stepw_model[[counter]]$lbb <- stepw_model[[counter]]$L <- stepw_model[[counter]]$St <- NULL
  } else {
    model_formula[[counter]] <- global_formula
    tmp <-  gam_scm(global_formula,
                    family=mvn_scm(d = d, param = param),
                    optimizer = "efs",
                    data = data_train,
                    aGam = list(fit = FALSE, control = list(trace = FALSE)))
    time <- microbenchmark(
      stepw_model <- gam_scm(optimizer = "efs",
                             aGam = list(G = tmp, control = list(trace = FALSE), in.out = list(sp = sp_init, scale = 1), start = beta_init)), times=1L)
    num_iter[[counter]] <-  stepw_model$family$getNC() - 1
    sum_Vcov_Peff[[counter]] <- summary(stepw_model, print = FALSE)$p.table # summary table for linear effects (if there are any)
    sum_Vcov_Seff[[counter]] <- summary(stepw_model, print = FALSE)$s.table
    betas[[counter]] <- stepw_model$coefficients
    sps[[counter]] <- stepw_model$sp
    stepw_model$R <- stepw_model$H <- stepw_model$Vc <- stepw_model$Vp <- stepw_model$Ve <- stepw_model$lbb <- stepw_model$L <- stepw_model$St <- NULL
  }

  time_fit[[counter]] <- time$time

  if(save.gam){
    #tmp_save <- list(fit = stepw_model, time_fit = time_fit, summary_Peff = sum_Vcov_Peff, summary_Seff = sum_Vcov_Seff, np0 = np0, betas = betas, sps = sps)
    #save(file = "tmp_save.RData", tmp_save)
    return(list(fit = stepw_model, time_fit = time_fit, num_iter = num_iter, summary_Peff = sum_Vcov_Peff, summary_Seff = sum_Vcov_Seff, np0 = np0, betas = betas, sps = sps))
  } else {
    #tmp_save <- list(foo = model_formula, time_fit = time_fit, summary_Peff = sum_Vcov_Peff, summary_Seff = sum_Vcov_Seff, np0 = np0, betas = betas, sps = sps)
    #save(file = "tmp_save.RData", tmp_save)
    return(list(foo =    model_formula, time_fit = time_fit, num_iter = num_iter, summary_Peff = sum_Vcov_Peff, summary_Seff = sum_Vcov_Seff, np0 = np0, betas = betas, sps = sps))
  }
}



###################################################################
cross_val <- function(obj, param, d, data, sets_eval, idx_vcov, ncores = 1, save.gam = NULL, root_dir){

  res_sim <- function(obj, param, d, sets_eval, idx_vcov, data, save.gam){
    dss <- lapply(idx_vcov, function(jj){
      out <- list()
      if(save.gam){
        print("not yet implemented")
      } else {
        mod_foo <- obj$foo[[jj]]
        time <- microbenchmark(
          stepw_model <- gam_scm(mod_foo,
                                 family = mvn_scm(d = d, param = param), optimizer = "efs",
                                 data = data[1 : sets_eval[1],],
                                 aGam = list(start = obj$betas[[jj]], in.out = list(sp = obj$sps[[jj]], scale = 1))), times = 1L)
        out$lpi_pred_in <- predict(stepw_model)
        out$lpi_pred_out <- predict(stepw_model, newdata = data[(sets_eval[1] + 1): sets_eval[2], ])
      }
      return(out)
    })
    return(dss)
  }

  cl <- makePSOCKcluster(ncores)
  setDefaultCluster(cl)

  clusterExport(NULL, c("root_dir"), envir = environment())

  clusterEvalQ(NULL, {
    library("mgcv", lib.loc=paste0(root_dir, "/my_library"))
  })

  clusterExport(NULL, c("obj", "save.gam", "param", "d", "idx_vcov", "data", "res_sim", "sets_eval"), envir = environment())

  clusterEvalQ(NULL, {
    library(SCM)
    library(microbenchmark)
    library(stringr)
    library(BMisc)
    if(packageVersion("mgcv") != "9.0"){
      stop("Wrong version of mgcv!!")
    }
    library(Rcpp)
    sourceCpp(code = "
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
void over_writediag(arma::mat& E, const arma::vec& rho, const arma::vec& ind, const int k_sp) {
  int kk = k_sp - 1;
  int dd;
  for (int i = 0; i < ind.n_elem; ++i) {
    dd = ind(i) - 1;
    E(dd, dd) = exp(rho(kk) * 0.5);
  }
}
")
  })


  out_res_sim <- function(.x){
    out2 <- res_sim(obj = obj, param = param, d = d, data = data, idx_vcov = idx_vcov, sets_eval = sets_eval[.x : (.x + 1)], save.gam = save.gam)
    return(list(gen = out2))
  }


  environment(out_res_sim) <- .GlobalEnv

  res <- list()
  res <- parLapply(NULL, 1 : (length(sets_eval) - 1), out_res_sim)

  stopCluster(cl)
  rm(cl)
  return(res)
}





