################################################################################
# Needed functions for fitting the models of Section  7: GEFCom14 Application  #
################################################################################


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


stepw_res <- function(param, d, grid_length, mean_model_formula, data, eff_vcov = NULL, metric = c("p", "logp", "ST")){

  metric <- match.arg(metric)
  stepw_model <- list()      # list of fitted model
  time_fit <- list()         #list of model fitting times
  sum_Vcov_Peff <- list()
  sum_Vcov_Seff <- list()

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


  #First model is the full model
  time <- microbenchmark(
    stepw_model[[1]] <- gam_scm(global_formula_full,
                                family = mvn_scm(d = d), optimizer = "efs",
                                data = data,
                                aGam = list(control = list(trace = FALSE))),
    times=1L)
  time_fit[[1]] <- time$time

  # This function extract the string of Theta elements to be selected
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


  # This function extract the string of effects to be selected
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

  # This function extract the position indices of Theta elements to be selected
  get_idx_lpi_Vcov <- function(string_eff){
    idx <- which(grepl("Th_", string_eff))
    return(idx)
  }

  # Functions for checking if an object is not an integer(0) or not a logical(0)
  is.integer0 <- function(x){
    is.integer(x) && length(x) == 0L
  }

  is.logical0 <- function(x){
    is.logical(x) && length(x) == 0L
  }



  # Function for building the model formula:
  # It starts from a fitted model, get the summary and order the  the smooth (and linear) effects
  # involved in covariance matrix modelling according to the p-value, build the model formula for theta and append such formula to the mean model
  mod_foo_building <- function(summary_Peff, summary_Seff, mean_model_formula,  K_covmod, grid_length, metric = metric){

    #!!! Problem: if set a different number of basis you cannot extract the number of basis from the summary print (this is related to the mgcv summary function)
    # Thus, maybe it should be better to consider the default k=10
    #sum_Vcov_Peff <- summary(model, print = FALSE)$p.table # summary table for linear effects (if there are any)
    #sum_Vcov_Seff <- summary(model, print = FALSE)$s.table # table of smooth effects

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
    return(global_formula)
  }


  sum_Vcov_Peff[[1]] <- summary(stepw_model[[1]], print = FALSE)$p.table # summary table for linear effects (if there are any)
  sum_Vcov_Seff[[1]] <- summary(stepw_model[[1]], print = FALSE)$s.table # table of smooth effects

  tmp_save <- list(fit = stepw_model, time_fit = time_fit, summary_Peff = sum_Vcov_Peff, summary_Seff = sum_Vcov_Seff)
  save(file = "tmp_save.RData", tmp_save)


  K_covmod <- d * (d + 1)/2
  # Build the model formula of the first nested model
  global_formula <- mod_foo_building(summary_Peff = sum_Vcov_Peff[[1]], summary_Seff = sum_Vcov_Seff[[1]],
                                     mean_model_formula, K_covmod, grid_length, metric)


  # Reduce memory footprint by removing gam.object not useful for the procedure
  stepw_model[[1]]$R <- stepw_model[[1]]$H <- stepw_model[[1]]$Vc <- stepw_model[[1]]$Vp <- stepw_model[[1]]$Ve <- stepw_model[[1]]$lbb <- stepw_model[[1]]$L <- stepw_model[[1]]$St <- NULL

  K_covmod <- d * (d + 1)/2 - grid_length
  counter <- 2
  while(K_covmod > 0){
    time <- microbenchmark(
      stepw_model[[counter]] <- gam_scm(global_formula,
                                        family=mvn_scm(d = d, param = param),
                                        optimizer = "efs",
                                        data = data,
                                        aGam = list(control = list(trace = FALSE)))
      , times=1L)
    time_fit[[counter]] <- time$time

    sum_Vcov_Peff[[counter]] <- summary(stepw_model[[counter]], print = FALSE)$p.table # summary table for linear effects (if there are any)
    sum_Vcov_Seff[[counter]] <- summary(stepw_model[[counter]], print = FALSE)$s.table

    K_covmod <- K_covmod -  grid_length
    global_formula <- mod_foo_building( sum_Vcov_Peff[[counter]], sum_Vcov_Seff[[counter]] , mean_model_formula, K_covmod, grid_length, metric)  # Build the model formula of subsequent steps
    stepw_model[[counter]]$R <- stepw_model[[counter]]$H <- stepw_model[[counter]]$Vc <- stepw_model[[counter]]$Vp <- stepw_model[[counter]]$Ve <- stepw_model[[counter]]$lbb <- stepw_model[[counter]]$L <- stepw_model[[counter]]$St <- NULL
    tmp_save <- list(fit = stepw_model, time_fit = time_fit, summary_Peff = sum_Vcov_Peff, summary_Seff = sum_Vcov_Seff)
    save(file = "tmp_save.RData", tmp_save)
    counter <- counter + 1
  }

  global_formula <- mean_model_formula
  time <- microbenchmark(
    stepw_model[[counter]] <- gam_scm(global_formula,
                                      family=mvn_scm(d = d, param = param),
                                      optimizer = "efs",
                                      data = data,
                                      aGam = list(control = list(trace = FALSE)))
    , times=1L)
  time_fit[[counter]] <- time$time
  sum_Vcov_Peff[[counter]] <- summary(stepw_model[[counter]], print = FALSE)$p.table # summary table for linear effects (if there are any)
  sum_Vcov_Seff[[counter]] <- summary(stepw_model[[counter]], print = FALSE)$s.table
  stepw_model[[counter]]$R <- stepw_model[[counter]]$H <- stepw_model[[counter]]$Vc <- stepw_model[[counter]]$Vp <- stepw_model[[counter]]$Ve <- stepw_model[[counter]]$lbb <- stepw_model[[counter]]$L <- stepw_model[[counter]]$St <- NULL
  tmp_save <- list(fit = stepw_model, time_fit = time_fit, summary_Peff = sum_Vcov_Peff, summary_Seff = sum_Vcov_Seff)
  save(file = "tmp_save.RData", tmp_save)

  return(list(fit = stepw_model, time_fit = time_fit, summary_Peff = sum_Vcov_Peff, summary_Seff = sum_Vcov_Seff))
}




