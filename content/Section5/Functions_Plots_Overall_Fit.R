##############################################################################
# Needed function for visualising the computational times and the            #
# iterations of the model fit for comparisons under the MCD parametrisation  #
##############################################################################


############################################################################
# This function extracts and summarise the computational times. Arguments: #
#  obj: an object including the simulation results                         #
#  param: parametrisation (her only "mcd" is set)                          #
#  dgrid: grid of outcome vector dimensions                                #
#  nrun: number of evaluations                                             #
############################################################################

fit_time <- function(obj, param = "mcd", dgrid = NULL, nrun){
  if(param == "mcd") param2 <- "MCD"


  ###################
  # Time in minutes #
  ###################
  data_time_efs <-  data.frame(unlist(lapply(1 : length(dgrid),
                                             function(z) unlist(lapply(1 : nrun, function(x) obj[[x]]$gen[[z]]$time_efs / (1e9 * 60))))),
                               rep(dgrid, each = nrun))
  data_time_exact_efs <-  data.frame(unlist(lapply(1 : length(dgrid),
                                             function(z) unlist(lapply(1 : nrun, function(x) obj[[x]]$gen[[z]]$time_exact_efs / (1e9 * 60))))),
                               rep(dgrid, each = nrun))

  data_time_exact_efs_initialised <-  data.frame(unlist(lapply(1 : length(dgrid),
                                                        function(z) unlist(lapply(1 : nrun, function(x) obj[[x]]$gen[[z]]$time_exact_efs_initialised / (1e9 * 60))))),
                                                 rep(dgrid, each = nrun))

  data_time_bfgs <-  data.frame(unlist(lapply(1 : length(dgrid),
                                             function(z) unlist(lapply(1 : nrun, function(x) obj[[x]]$gen[[z]]$time_bfgs / (1e9 * 60))))),
                               rep(dgrid, each = nrun))
  data_time_bfgsinit <-  data.frame(unlist(lapply(1 : length(dgrid),
                                              function(z) unlist(lapply(1 : nrun, function(x) obj[[x]]$gen[[z]]$time_bfgs_efs / (1e9 * 60))))),
                                rep(dgrid, each = nrun))
  data_time_bamlss <-  data.frame(unlist(lapply(1 : length(dgrid),
                                                  function(z) unlist(lapply(1 : nrun, function(x) obj[[x]]$gen[[z]]$time_bamlss / (1e9 * 60))))),
                                    rep(dgrid, each = nrun))

  res <- cbind(data_time_efs[,1], data_time_exact_efs[,1], data_time_exact_efs_initialised[,1], data_time_bfgs[,1], data_time_bfgsinit[,1], data_time_bamlss[,1], data_time_efs[,2])
  colnames(res) <- c("time_efs", "time_exact_efs", "time_exact_efs_initialised", "time_bfgs", "time_bfgsinit", "time_bamlss", "d")
  res <- as.data.frame(res)

  ###################
  # Summary indices #
  ###################
  mean_efs_time <- aggregate(res$time_efs, list(res$d), FUN = mean)[,2]
  mean_exact_efs_time <- aggregate(res$time_exact_efs, list(res$d), FUN = mean)[,2]
  mean_exact_efs_time_initialised <- aggregate(res$time_exact_efs_initialised, list(res$d), FUN = mean)[,2]
  mean_bfgs_time <- aggregate(res$time_bfgs, list(res$d), FUN = mean)[,2]
  mean_bfgsinit_time <- aggregate(res$time_bfgsinit, list(res$d), FUN = mean)[,2]
  mean_bamlss_time <- aggregate(res$time_bamlss, list(res$d), FUN = mean)[,2]

  out <- cbind(mean_efs_time, mean_exact_efs_time, mean_exact_efs_time_initialised, mean_bfgs_time, mean_bfgsinit_time, mean_bamlss_time, dgrid)
  colnames(out) <- c("efs", "efsExact", "efsExact_initalised", "bfgs", "bfgsinit", "bamlss", "d")
  sum_res <- as.data.frame(out)

  return(list(sum_res = data.frame(sum_res),
              res = data.frame(res)))
}



fit_time_red <- function(obj, param = "mcd", dgrid = NULL, nrun){
  if(param == "mcd") param2 <- "MCD"


  ###################
  # Time in minutes #
  ###################
  data_time_efs <-  data.frame(unlist(lapply(1 : length(dgrid),
                                             function(z) unlist(lapply(1 : nrun, function(x) obj[[x]]$gen[[z]]$time_efs / (1e9 * 60))))),
                               rep(dgrid, each = nrun))
  data_time_exact_efs <-  data.frame(unlist(lapply(1 : length(dgrid),
                                                   function(z) unlist(lapply(1 : nrun, function(x) obj[[x]]$gen[[z]]$time_exact_efs / (1e9 * 60))))),
                                     rep(dgrid, each = nrun))


  res <- cbind(data_time_efs[,1], data_time_exact_efs[,1], data_time_efs[,2])
  colnames(res) <- c("time_efs", "time_exact_efs",  "d")
  res <- as.data.frame(res)

  ###################
  # Summary indices #
  ###################
  mean_efs_time <- aggregate(res$time_efs, list(res$d), FUN = mean)[,2]
  mean_exact_efs_time <- aggregate(res$time_exact_efs, list(res$d), FUN = mean)[,2]

  out <- cbind(mean_efs_time, mean_exact_efs_time,  dgrid)
  colnames(out) <- c("efs", "efsExact",  "d")
  sum_res <- as.data.frame(out)

  return(list(sum_res = data.frame(sum_res),
              res = data.frame(res)))
}
########################################################
# The idx_bamlss() function allows to get the indices  #
# of bamlss in the ordering used in our implementation #
########################################################
idx_bamlss <- function(d){
  mat_idx <- matrix(0, d, d)
  count <- 2 * d + 1
  for(i in 1 : d){
    mat_idx[i, i] <- i + d
    if(i < d){
      for(j in (i + 1) : d){
        mat_idx[i, j] <- count
        count <- count + 1
      }
    }
  }
  idx1 <-   c(diag(mat_idx), mat_idx[upper.tri(mat_idx, diag = FALSE)])
  idx2 <- (d + 1) : (d + d * (d + 1)/2)
  return(list(idx1 = idx1, idx2 = idx2) )
}

####################################################################
# This function computes the logScore on the train set. Arguments: #
#  obj: an object including the simulation results                 #
#  nrun: number of evaluations                                     #
#  dgrid: grid of outcome vector dimensions                        #
#  nobs: number of observations (of the train set)                 #
#  param: parametrisation (only "mcd")                             #
####################################################################

log_Score_train <- function(obj, nrun, dgrid, nobs, param = "mcd"){

  res_efs <- lapply(1 : length(dgrid),
                    function(z) lapply(1 : nrun,
                                   function(x) SCM::internal()$ll_mcd(obj[[x]]$gen[[z]]$lpi_pred_efs, as.matrix(obj[[x]]$gen[[z]]$sim$data_train[, c(1 : dgrid[z])]))))
  res_exact_efs <- lapply(1 : length(dgrid),
                    function(z) lapply(1 : nrun,
                                       function(x) SCM::internal()$ll_mcd(obj[[x]]$gen[[z]]$lpi_pred_exact_efs, as.matrix(obj[[x]]$gen[[z]]$sim$data_train[, c(1 : dgrid[z])]))))
  res_exact_efs_initialised <- lapply(1 : length(dgrid),
                                      function(z) lapply(1 : nrun,
                                             function(x) SCM::internal()$ll_mcd(obj[[x]]$gen[[z]]$lpi_pred_exact_efs_initialised, as.matrix(obj[[x]]$gen[[z]]$sim$data_train[, c(1 : dgrid[z])]))))

  res_bfgs <- lapply(1 : length(dgrid),
                     function(z) lapply(1 : nrun,
                                       function(x) SCM::internal()$ll_mcd(obj[[x]]$gen[[z]]$lpi_pred_bfgs, as.matrix(obj[[x]]$gen[[z]]$sim$data_train[, c(1 : dgrid[z])]))))
  res_bfgs_efs <- lapply(1 : length(dgrid),
                         function(z) lapply(1 : nrun,
                                        function(x) SCM::internal()$ll_mcd(obj[[x]]$gen[[z]]$lpi_pred_bfgs_efs, as.matrix(obj[[x]]$gen[[z]]$sim$data_train[, c(1 : dgrid[z])]))))

  lpi_pred_bamlss <- list()
  for(j in 1:length(dgrid)){
    lpi_pred_bamlss[[j]]<-list()
    for(k in 1:nrun){
      lpi_pred_bamlss[[j]][[k]] <- matrix(0,nobs,dgrid[j]+dgrid[j]*(dgrid[j]+1)/2)
      idx <- idx_bamlss(dgrid[j])
      for(l in 1:dgrid[j]){
        lpi_pred_bamlss[[j]][[k]][,l] <-  obj[[k]]$gen[[j]]$lpi_pred_bamlss[[l]]
      }
      for(l in 1:(dgrid[j]*(dgrid[j]+1)/2)){
        if(idx$idx2[l] > 2*dgrid[j]){
          lpi_pred_bamlss[[j]][[k]][,idx$idx2[l]] <-  -obj[[k]]$gen[[j]]$lpi_pred_bamlss[[idx$idx1[l]]]
        } else {
          lpi_pred_bamlss[[j]][[k]][,idx$idx2[l]] <-  obj[[k]]$gen[[j]]$lpi_pred_bamlss[[idx$idx1[l]]]
        }
      }
    }
  }

  res_bamlss <- lapply(1 : length(dgrid),
                       function(z) lapply(1 : nrun,
                                          function(x) SCM::internal()$ll_mcd(lpi_pred_bamlss[[z]][[x]], as.matrix(obj[[x]]$gen[[z]]$sim$data[, c(1 : dgrid[z])]))))

  lpi_pred_gen <- list()
  for(i in 1:length(dgrid)){
    lpi_pred_gen[[i]] <- list()
    for(j in 1:nrun){
      lpi_pred_gen[[i]][[j]] <- matrix(NA, nobs, dgrid[i] + dgrid[i]*(dgrid[i] + 1)/2)
      for(k in 1: nobs){
        lpi_pred_gen[[i]][[j]][k, ] <- c(obj[[j]]$gen[[i]]$sim$mu_train[[k]],   diag(obj[[j]]$gen[[i]]$sim$Theta_train[[k]]),
                                         obj[[j]]$gen[[i]]$sim$Theta_train[[k]][upper.tri(obj[[j]]$gen[[i]]$sim$Theta_train[[k]], diag=FALSE)])
      }
    }
  }


  res_gen <- lapply(1 : length(dgrid),
                function(z) lapply(1 : nrun,
                                   function(x) SCM::internal()$ll_mcd(lpi_pred_gen[[z]][[x]], as.matrix(obj[[x]]$gen[[z]]$sim$data_train[, c(1 : dgrid[z])]))))
  res <- list(res_gen, res_efs, res_exact_efs, res_exact_efs_initialised, res_bfgs, res_bfgs_efs, res_bamlss)

  out <- list()
  for(i in 1 : 7){
    out[[i]] <- matrix(0, nrun, length(dgrid))
    for(j in 1:length(dgrid)){
      for(k in 1:nrun){
        out[[i]][k,j] <- -res[[i]][[j]][[k]]
      }
    }
  }

  return(out)
}

####################################################################
# This function computes the logScore on the test set. Arguments:  #
#  obj: an object including the simulation results                 #
#  nrun: number of evaluations                                     #
#  dgrid: grid of outcome vector dimensions                        #
#  nobs: number of observations (of the test set)                  #
#  param: parametrisation (only "mcd")                             #
####################################################################
log_Score_test <- function(obj, nrun, dgrid, nobs, param = "mcd"){

  res_efs <- lapply(1 : length(dgrid),
                    function(z) lapply(1 : nrun,
                                       function(x) SCM::internal()$ll_mcd(obj[[x]]$gen[[z]]$lpi_pred_test_efs, as.matrix(obj[[x]]$gen[[z]]$sim$data_test[, c(1 : dgrid[z])]))))
  res_exact_efs <- lapply(1 : length(dgrid),
                    function(z) lapply(1 : nrun,
                                       function(x) SCM::internal()$ll_mcd(obj[[x]]$gen[[z]]$lpi_pred_test_exact_efs, as.matrix(obj[[x]]$gen[[z]]$sim$data_test[, c(1 : dgrid[z])]))))
  res_exact_efs_initialised <- lapply(1 : length(dgrid),
                                      function(z) lapply(1 : nrun,
                                             function(x) SCM::internal()$ll_mcd(obj[[x]]$gen[[z]]$lpi_pred_test_exact_efs_initialised, as.matrix(obj[[x]]$gen[[z]]$sim$data_test[, c(1 : dgrid[z])]))))

  res_bfgs <- lapply(1 : length(dgrid),
                     function(z) lapply(1 : nrun,
                                        function(x) SCM::internal()$ll_mcd(obj[[x]]$gen[[z]]$lpi_pred_test_bfgs, as.matrix(obj[[x]]$gen[[z]]$sim$data_test[, c(1 : dgrid[z])]))))
  res_bfgs_efs <- lapply(1 : length(dgrid),
                         function(z) lapply(1 : nrun,
                                            function(x) SCM::internal()$ll_mcd(obj[[x]]$gen[[z]]$lpi_pred_test_bfgs_efs, as.matrix(obj[[x]]$gen[[z]]$sim$data_test[, c(1 : dgrid[z])]))))

  lpi_pred_test_bamlss <- list()
  for(j in 1:length(dgrid)){
    lpi_pred_test_bamlss[[j]]<-list()
    for(k in 1:nrun){
      lpi_pred_test_bamlss[[j]][[k]] <- matrix(0,nobs,dgrid[j]+dgrid[j]*(dgrid[j]+1)/2)
      idx <- idx_bamlss(dgrid[j])
      for(l in 1:dgrid[j]){
        lpi_pred_test_bamlss[[j]][[k]][,l] <-  obj[[k]]$gen[[j]]$lpi_pred_test_bamlss[[l]]
      }
      for(l in 1:(dgrid[j]*(dgrid[j]+1)/2)){
        if(idx$idx2[l] > 2*dgrid[j]){
          lpi_pred_test_bamlss[[j]][[k]][,idx$idx2[l]] <-  -obj[[k]]$gen[[j]]$lpi_pred_test_bamlss[[idx$idx1[l]]]
        } else {
          lpi_pred_test_bamlss[[j]][[k]][,idx$idx2[l]] <-  obj[[k]]$gen[[j]]$lpi_pred_test_bamlss[[idx$idx1[l]]]
        }
      }
    }
  }

  res_bamlss <- lapply(1 : length(dgrid),
                       function(z) lapply(1 : nrun,
                                          function(x) SCM::internal()$ll_mcd(lpi_pred_test_bamlss[[z]][[x]], as.matrix(obj[[x]]$gen[[z]]$sim$data_test[, c(1 : dgrid[z])]))))

  lpi_pred_test_gen <- list()
  for(i in 1:length(dgrid)){
    lpi_pred_test_gen[[i]] <- list()
    for(j in 1:nrun){
      lpi_pred_test_gen[[i]][[j]] <- matrix(NA, nobs, dgrid[i] + dgrid[i]*(dgrid[i] + 1)/2)
      for(k in 1: nobs){
        lpi_pred_test_gen[[i]][[j]][k, ] <- c(obj[[j]]$gen[[i]]$sim$mu_test[[k]],   diag(obj[[j]]$gen[[i]]$sim$Theta_test[[k]]),
                                              obj[[j]]$gen[[i]]$sim$Theta_test[[k]][upper.tri(obj[[j]]$gen[[i]]$sim$Theta_test[[k]], diag=FALSE)])
      }
    }
  }

  res_gen <- lapply(1 : length(dgrid),
                    function(z) lapply(1 : nrun,
                                       function(x) SCM::internal()$ll_mcd(lpi_pred_test_gen[[z]][[x]], as.matrix(obj[[x]]$gen[[z]]$sim$data_test[, c(1 : dgrid[z])]))))
  res <- list(res_gen, res_efs, res_exact_efs, res_exact_efs_initialised, res_bfgs, res_bfgs_efs, res_bamlss)

  out <- list()
  for(i in 1 : 7){
    out[[i]] <- matrix(0, nrun, length(dgrid))
    for(j in 1:length(dgrid)){
      for(k in 1:nrun){
        out[[i]][k,j] <- -res[[i]][[j]][[k]]
      }
    }
  }

  return(out)
}


##################################################################################
# This function extract the LAML (on the train set) for efs and bfgs. Arguments: #
#  obj: an object including the simulation results                               #
#  nrun: number of evaluations                                                   #
#  dgrid: grid of outcome vector dimensions                                      #
#  nobs: number of observations (of the train set)                               #
#  param: parametrisation (only "mcd")                                           #
#################################################################################

LAML_extraction <- function(obj, nrun, dgrid, nobs, param = "mcd"){

  res_efs <- lapply(1 : length(dgrid),
                    function(z) lapply(1 : nrun,
                                       function(x) sim_mcd_fit[[x]]$gen[[z]]$LAML_efs))
  res_exact_efs <- lapply(1 : length(dgrid),
                     function(z) lapply(1 : nrun,
                                        function(x) sim_mcd_fit[[x]]$gen[[z]]$LAML_exact_efs))

  res_exact_efs_initialisation <- lapply(1 : length(dgrid),
                          function(z) lapply(1 : nrun,
                                             function(x) sim_mcd_fit[[x]]$gen[[z]]$LAML_exact_efs_initialised))

  res_bfgs <- lapply(1 : length(dgrid),
                     function(z) lapply(1 : nrun,
                                        function(x) sim_mcd_fit[[x]]$gen[[z]]$LAML_bfgs))
  res_bfgs_efs <- lapply(1 : length(dgrid),
                         function(z) lapply(1 : nrun,
                                            function(x) sim_mcd_fit[[x]]$gen[[z]]$LAML_bfgs_efs))

  res <- list(res_efs, res_exact_efs, res_exact_efs_initialisation, res_bfgs, res_bfgs_efs)

  out <- list()
  for(i in 1 : 5){
    out[[i]] <- matrix(0, nrun, length(dgrid))
    for(j in 1:length(dgrid)){
      for(k in 1:nrun){
        out[[i]][k,j] <- res[[i]][[j]][[k]]
      }
    }
  }

  return(out)
}
