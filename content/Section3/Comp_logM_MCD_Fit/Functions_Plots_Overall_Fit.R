###########################################################
# Needed function for visualising the computational times #
# and the iterations of the model fit                     #
###########################################################

# This function extracts and summarise the computational times, as well as the number of inner iterations. Arguments:
#  obj: an object including the simulation results
#  param: parametrisation
#  dgrid: grid of outcome vector dimensions
#  nrun: number of evaluations

fit_time <- function(obj, param = c("mcd", "logm"), dgrid = NULL, nrun){
  if(param == "mcd") param2 <- "MCD" else param2 <- "logM"


  # Number of inner iterations
  data_iter <- data.frame(unlist(lapply(1 : length(dgrid),
                                        function(z) unlist(lapply(1 : nrun,
                                                  function(x) obj[[x]]$gen[[z]]$inner)))),
                          rep(dgrid, each = nrun))

  data_iter <- cbind(data_iter, rep(param2, each = length(dgrid)))
  colnames(data_iter) <- c("iter", "d", "Type")

  data_sum_iter <- aggregate(data_iter$iter, list(data_iter$d), FUN = mean)
  data_sum_iter_min <- aggregate(data_iter$iter, list(data_iter$d), FUN = min)
  data_sum_iter_max <- aggregate(data_iter$iter, list(data_iter$d), FUN = max)

  out_iter <- cbind(data_sum_iter, data_sum_iter_min[,2], data_sum_iter_max[,2])

  sum_res_iter <- as.data.frame(out_iter)
  sum_res_iter <- cbind(sum_res_iter, rep(param2, each = length(dgrid)))
  colnames(sum_res_iter) <- c("d", "mean_iter", "min_iter", "max_iter", "Type")

  # Time in minutes
  data_time <-  data.frame(unlist(lapply(1 : length(dgrid),
                                         function(z) unlist(lapply(1 : nrun, function(x) obj[[x]]$gen[[z]]$time / (1e9 * 60))))),
                           rep(dgrid, each = nrun))
  res <- data_time

  res <- cbind(res, rep(param2, each = length(dgrid)))
  colnames(res) <- c("time", "d", "Type")

  data <- aggregate(res$time, list(res$d), FUN = mean)
  data_min <- aggregate(res$time, list(res$d), FUN = min)
  data_max <- aggregate(res$time, list(res$d), FUN = max)

  out <- cbind(data,data_min[,2], data_max[,2])

  sum_res <- as.data.frame(out)
  sum_res <- cbind(sum_res, rep(param2, each = length(dgrid)))
  colnames(sum_res) <- c("d", "mean_time", "min_time", "max_time", "Type")

  return(list(sum_res = data.frame(sum_res),
              res = data.frame(res),
              iter = data.frame(data_iter),
              sum_iter = data.frame(sum_res_iter)))
}



# This function computes the logScore. Arguments:
#  obj: an object including the simulation results
#  nrun: number of evaluations
#  dgrid: grid of outcome vector dimensions
#  nobs: number of observations
#  param: parametrisation

log_Score <- function(obj, nrun, dgrid, nobs, param = c("logm","mcd"), ncov = NULL){
  param <- match.arg(param)
  if(param == "mcd"){
    res <- lapply(1 : length(dgrid),
                  function(z) lapply(1 : nrun,
                                     function(x) SCM::internal()$ll_mcd(obj[[x]]$gen[[z]]$lpi_pred, as.matrix(obj[[x]]$gen[[z]]$sim$data[, c(1 : dgrid[z])]))))
  } else {
    res <- lapply(1 : length(dgrid),
                  function(z) lapply(1 : nrun,
                                     function(x) SCM::internal()$ll_logm(obj[[x]]$gen[[z]]$lpi_pred, as.matrix(obj[[x]]$gen[[z]]$sim$data[, c(1 : dgrid[z])]))))
  }

  out <- matrix(0, nrun, length(dgrid))
  for(j in 1:length(dgrid)){
    for(k in 1:nrun){
      out[k,j] <- -res[[j]][[k]]
    }
  }
  return(out)
}

