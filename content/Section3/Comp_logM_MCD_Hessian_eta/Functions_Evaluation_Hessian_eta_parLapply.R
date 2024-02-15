#########################################################################################
# Needed functions for reproducing the results of Section 3.3 - 2nd derivatives wrt eta #
#########################################################################################
source("idxHess_no0.R") # List of indices to deal with the sparsity of 2nd derivatives

#########################################################################
# Compile TMB functions for evaluating the 2nd order derivatives via AD #
#########################################################################
compile("nll_MCD_TMB.cpp")
compile("nll_logM_TMB.cpp")

#######################################################
# Main function for obtaining the computational times #
#######################################################

# The time_Deta() function evaluate the 2nd derivatives w.r.t. eta for Efficient and TMB approaches. Arguments:
# nobs: number of rows
# dgrid: grid of outcome vector dimension
# nrun: number of runs
# ncores: number of cores
# param: type of parametrisation ("mcd" or "logm")

time_Deta <- function(nobs, dgrid,  nrun, ncores, param = c("mcd", "logm")){

  if(is.null(param)){
    param <- "mcd"
  } else {
    param <- match.arg(param)
  }

  sim_time <- function(dgrid, nobs, param){
    sourceCpp("idx_zwGt.cpp") # Function for obtaining the indices z, w, G (which are functions of d)

    # out is a list containing the computational times for Efficient (tD2_eff) and TMB (tD2_TMB)
    out <- list()
    out$tD2_eff <- out$tD2_TMB <- rep(0, length(dgrid))
    for(jj in 1 : length(dgrid)){
      d <- dgrid[jj]
      q <- d + d * (d + 1)/2   # Number of linear predictors

      # eta and y could be filled with 0 (but this initialisation is not considered in the computational times)
      eta <- matrix(rnorm(nobs * q), nobs, q)
      y <- matrix(rnorm(nobs * d), nobs, d)

      # Indices z, w,
      z <- w <- t <- rep(0, (d * (d - 1)/2))
      Gm <- matrix(0, d - 1, d - 1)
      mode(Gm) <- mode(z) <- mode(w) <- mode(t) <- "integer"
      idx_zwGt(d, z, w, Gm, t)

      # Hessian MCD parametrisation
      if(param == "mcd"){
        sourceCpp("d2_mcd_eta_row.cpp")   # Efficient implementation
        dyn.load(dynlib("nll_MCD_TMB"))   # TMB implementation

        # Number of Hessian elements different from zero
        nHel <- d * (d ^ 2 + 15 * d + 2)/6
        if(d > 2) nHel <- nHel + d * (d - 1) * (d - 2)/3

        # Ger the indices of Hessian elements different from zero
        idx_jk <- idxHess_no0(q, z, w, param = 1)

        # Initialisation of the 2nd derivatives vector
        res <- matrix(0, 1, nHel)

        # Computational times for evaluating the Hessian derivatives (Efficient)
        time <- microbenchmark({
          for(i in 1 : nobs){
            d2_mcd_eta_row(eta[i,], y[i,], res, z, w, Gm, t, idx_jk)
          }
        }, times = 1L)
        out$tD2_eff[jj] <- time$time


        # Computational times for evaluating the Hessian derivatives (AD)
        obj <- list()
        for(i in 1 : nobs){
          obj[[i]] <- MakeADFun(list(Y = matrix(y[i,], ncol = d)),
                                list(eta = matrix(eta[i,], ncol = q)),
                                DLL = "nll_MCD_TMB", silent = TRUE)
        }

        time <- microbenchmark({
          for(i in 1 : nobs){
            HTMBi <- obj[[i]]$he()
          }
        }, times = 1L)
        out$tD2_TMB[jj] <- time$time

      } else {        # Hessian logM parametrisation
        sourceCpp("d2_logm_eta_row.cpp") # Efficient
        dyn.load(dynlib("nll_logM_TMB")) # TMB

        # Number of elements of the hessian matrix (no sparsity)
        nHel <- q * (q + 1)/2

        # Initialisation of the 2nd derivatives vector
        res <- matrix(0, 1, nHel)

        # Computational times for Efficient
        time <- microbenchmark({
          for(i in 1 : nobs){
            d2_logm_eta_row(eta[i,], y[i,], res, z, w)
          }
        }, times = 1L)
        out$tD2_eff[jj] <- time$time

        # Computational times for TMB
        obj <- list()
        for(i in 1 : nobs){
          obj[[i]] <- MakeADFun(list(Y = matrix(y[i,], ncol = d)),
                                list(eta = matrix(eta[i,], ncol = q)),
                                DLL = "nll_logM_TMB", silent = TRUE)
        }

        time <- microbenchmark({
          for(i in 1 : nobs){
            HTMBi <- obj[[i]]$he()
          }
        }, times = 1L)
        out$tD2_TMB[jj] <- time$time
      }
    }
    return(out)

  }

  # Setting parallel computation
  cl <- makePSOCKcluster(ncores)
  setDefaultCluster(cl)
  clusterExport(NULL, c("nobs", "sim_time", "dgrid", "param","idxHess_no0"), envir = environment())
  clusterEvalQ(NULL, {
    library(TMB)
    library(microbenchmark)
    library(Rcpp)
  })

  out_time <- function(.x){
    out2 <- sim_time(dgrid = dgrid, nobs = nobs, param = param)
    return(list(time = out2))
  }

  res <- list()
  res <- parLapply(NULL, 1 : nrun, out_time)

  stopCluster(cl)
  rm(cl)
  return(res)
}
