##########################################################
# Needed function for visualisng the computational times #
# of evaluating the Hessian w.r.t. eta and beta          # 
##########################################################

# Log Scale
myscale_trans <- function()
{
  trans_new("myscale", log, exp) 
}

# Size
scaleFUN <- function(x) sprintf("%.2f", x)


# Function for extracting and summarising the computational times. Arguments: 
#  obj: an object includin the computational times
#  param: parametrisation ("logm" or "mcd") 
#  dgrid: grid of outcome vector dimensions
#  nrun: number of evaluations
#  type: vector including the strategise to compare 
#        (so for the hessian w.r.t . eta type = c("eff","TMB"), 
#         while for the hessian w.r.t . beta type = c("",""))
#  type1: vector labelling the strategies to compare 
#         (so for the hessian w.r.t . eta type = c("eff","TMB"), 
#         while for the hessian w.r.t . beta type = c("",""))
#  beta: flag for indicating if we are comparing the derivatives w.r.t. eta or w.r.t. beta

time_hessian <- function(obj, param = c("logm", "mcd"),  
                         dgrid = NULL, nrun, type, type1, beta = FALSE){
  param <- match.arg(param)
  
  # Summary of the times 
  out <- lapply(1 : length(type), function(jj){
    if(beta = = FALSE){
      str_type <- paste0("obj[[x]]$time$tD2_", type[jj])
    } else {
      str_type <- paste0("obj[[x]]$time$thessian_", type[jj])
    }
    
    # Time in minutes
    data_time <- data.frame(cbind(unlist(lapply(1 : nrun, 
                                                function(x) eval(parse(text = str_type)) / (1e9 * 60))),
                                  rep(dgrid, times = nrun)))
    colnames(data_time) <- c("time", "d")
    data <- aggregate(data_time$time, list(data_time$d), FUN = mean)
    data_min <- aggregate(data_time$time, list(data_time$d), FUN = min) #not used
    data_max <- aggregate(data_time$time, list(data_time$d), FUN = max) #not used
    data <- cbind(data, data_min[, 2], data_max[, 2]) 
    return(data)
  })
  
  sum_res <- matrix(0, length(dgrid) * length(out), 4)
  for(j in 1 : length(out)){
    sum_res[c((length(dgrid) * (j - 1) + 1) : (length(dgrid) * j)),] <- as.matrix(out[[j]])
  }  
  sum_res <- as.data.frame(sum_res)
  sum_res <- cbind(sum_res, 
                   rep(type1, each = length(dgrid)),
                   rep(param, each = length(dgrid)))
  colnames(sum_res) <- c("d", "mean_time", "min_time", "max_time", "Type", "param")
  
  # Overall times: repetition of (part of) the code  above (improbvable)
  out <- lapply(1 : length(type), function(jj){
    if(beta ==FALSE){
      str_type <- paste0("obj[[x]]$time$tD2_", type[jj])
    } else {
      str_type <- paste0("obj[[x]]$time$thessian_", type[jj])
    }
    # Times in minutes
    data_time <- data.frame(cbind(unlist(lapply(1 : nrun, 
                                                function(x) eval(parse(text = str_type)) / (1e9 * 60))),
                                  rep(dgrid, times = nrun)))
    colnames(data_time) <- c("time", "d")
    return(data_time)
  })
  
  res <- out[[1]]
  for(j in 2 : length(out)){
    res <- rbind(res, out[[j]]) # not beautiful
  }
  
  res <- cbind(res, 
               rep(type1, each = length(dgrid) * nrun), 
               rep(param, each = length(dgrid) * nrun))
  
  colnames(res) <- c("time", "d", "Type", "param") 
  
  return(list(sum_res = data.frame(sum_res),
              res = data.frame(res)))
}

