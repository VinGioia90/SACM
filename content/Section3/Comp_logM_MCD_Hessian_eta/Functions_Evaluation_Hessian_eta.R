#########################################################################################
# useful functions for reproducing the results of Section 3.3 - 2nd derivatives wrt eta #
#########################################################################################
source("idxHess_no0.R") # List of indices to deal with the sparsity of 2nd derivatives

# Removed for loop over i (due to the TMB comparison) and final for loop for saving
sourceCpp("d2_mcd_eta_row.cpp") # MCD
sourceCpp("d2_logm_eta_row.cpp") # logM

sourceCpp("idx_zwGt.cpp") # Indices z,w,G

#######
# TMB #
#######
compile("nll_MCD_TMB.cpp")
dyn.load(dynlib("nll_MCD_TMB"))

compile("nll_logM_TMB.cpp")
dyn.load(dynlib("nll_logM_TMB"))

########################################################
# Main functions for obtaining the computational times #
########################################################

# This function evaluate the 2nd derivatives w.r.t. eta for Efficient and TMB. Arguments
# nobs: number of rows
# dgrid: grid of outcome vector dimension
# nrun: number of runs
# ncores: number of cores
# param: "mcd" or "logm" parametrisation

time_Deta <- function(nobs, dgrid,  nrun, ncores, param = c("mcd", "logm")){
  param <- match.arg(param)

  time <- mclapply(1 : nrun, function(ii){ #parLapply??? maybe it's not a big problem using mclapply here
    out <- list()
    out$tD2_eff <- out$tD2_TMB <- rep(0, length(dgrid))

    for(jj in 1 : length(dgrid)){
      d <- dgrid[jj]
      q <- d + d * (d + 1)/2

      # eta and y could be filled with 0
      eta <- matrix(rnorm(nobs * q), nobs, q)
      y <- matrix(rnorm(nobs * d), nobs, d)

      # Indices
      z <- w <- t <- rep(0, (d * (d - 1)/2))
      Gm <- matrix(0, d - 1, d - 1)
      mode(Gm) <- mode(z) <- mode(w) <- mode(t) <- "integer"
      idx_zwGt(d, z, w, Gm, t)

      # Hessian
      if(param == "mcd"){
        nHel <- d * (d ^ 2 + 15 * d + 2)/6
        if(d > 2) nHel <- nHel + d * (d - 1) * (d - 2)/3
        idx_jk <- idxHess_no0(q, z, w, param = 1)
        res <- matrix(0, 1, nHel) # Consider that we can avoid passing res (as in the no allocation version we implemented last time)
        time <- microbenchmark({
          for(i in 1 : nobs){
          d2_mcd_eta_row(eta[i,], y[i,], res, z, w, Gm, t, idx_jk)
          }
        }, times = 1L)
        out$tD2_eff[jj] <- time$time


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

      } else {
        nHel <- q * (q + 1)/2
        res<- matrix(0, 1, nHel)

        time <- microbenchmark({
          for(i in 1 : nobs){
          d2_logm_eta_row(eta[i,], y[i,], res, z, w)
          }
        }, times = 1L)
        out$tD2_eff[jj] <- time$time

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

    return(list("time" = out))
  },  mc.cores = ncores)
  return(time)
}




