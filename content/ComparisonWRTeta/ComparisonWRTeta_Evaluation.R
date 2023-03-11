#####################################################################################################################
# Code for reproducibility of the results of Section 2.2.3 - Computational comparison of model-specific quantities: #
#####################################################################################################################
library(rstudioapi)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("loadPackages.R")



#rm(list = ls())
#setwd("C:/Users/Gioia/Desktop/ScalableAdditiveCovarianceMatrixModels")
#setwd("C:/Users/Gioia/Desktop/ScalableAdditiveCovarianceMatrixModels/ComparisonWRTeta")

source("idxHess_no0.R") # List of indices to deal the sparsity of the Hessian

sourceCpp("d2_mcd_eta.cpp") # 2nd derivatives MCD (Software version)
sourceCpp("d2_mcd_eta_row.cpp") # 2nd derivatives MCD (removed for loop over i and final for loop for saving)
sourceCpp("d2_mcd_eta_row_noalloc.cpp") # 2nd derivatives MCD (removed for loop over i and final for loop for saving)

sourceCpp("d2_logm_eta.cpp") # 2nd derivatives logM  (Software version)
sourceCpp("d2_logm_eta_row.cpp") # 2nd derivatives logM (removed for loop over i and final for loop for saving)
sourceCpp("d2_logm_eta_row_noalloc.cpp") # 2nd derivatives MCD (removed for loop over i and final for loop for saving)

sourceCpp("idx_zwGt.cpp") # Indices z,w,G

#######
# TMB #
#######
compile("nll_MCD_TMB.cpp")
dyn.load(dynlib("nll_MCD_TMB"))

##############
# TMB - logM #
##############
compile("nll_logM_TMB.cpp")
dyn.load(dynlib("nll_logM_TMB"))



###############################################
time_Deta <- function(nobs, dgrid,  nrun, ncores, param = c("mcd", "logm"), TMB=TRUE){
  param <- match.arg( param )

  time <- mclapply(1:nrun, function(ii){
    out <- list()
    out$tD2_eff<-rep(0,length(dgrid))
    out$tD2_eff_noalloc<-rep(0,length(dgrid))
    out$tD2_TMB<-rep(0,length(dgrid))

    for(jj in 1: length(dgrid)){
      d <- dgrid[jj]
      q<-d+d*(d+1)/2

      eta<- matrix(rnorm(nobs * q), nobs, q) # eventually set to 0
      y<- matrix(rnorm(nobs * d), nobs, d) # eventually set to 0

      # Indices
      z <- w <- t <- rep(0, (d * (d - 1)/2))
      Gm <- matrix(0, d - 1, d - 1)
      mode(Gm) <- mode(z) <- mode(w) <- mode(t) <- "integer"
      idx_zwGt(d, z, w, Gm, t)

      # Hessian
      if(param == "mcd"){
        nHel <- d * (d^2 + 15 * d + 2)/6
        if ( d > 2 ) nHel <- nHel + d * (d - 1) * (d - 2)/3
        idx_jk <- idxHess_no0(q, z, w, 1)
        res<- matrix(0,1,nHel)
        time<- microbenchmark({
          for(i in 1: nobs){
          d2_mcd_eta_row(eta[i,],y[i,],res,z,w,Gm,t,idx_jk)
          }
        } ,times=1L)
        out$tD2_eff[jj] <- time$time

        time<- microbenchmark({
          for(i in 1: nobs){
            d2_mcd_eta_row_noalloc(eta[i,],y[i,],z,w,Gm,t,idx_jk)
          }
        } ,times=1L)
        out$tD2_eff_noalloc[jj] <- time$time

        obj <- list()
        for(i in 1:nobs){
          obj[[i]] <- MakeADFun(list(Y=matrix(y[i,],ncol=d)), list(eta = matrix(eta[i,],ncol=q)), DLL = "nll_MCD_TMB", silent=TRUE)
        }

        #HTMB <- matrix(0, nobs, nHel)
        time<-microbenchmark({
          for(i in 1:nobs){
            HTMBi <- obj[[i]]$he()
          }
        }, times=1L)
        out$tD2_TMB[jj] <- time$time

      } else {
        res<- matrix(0,1,q*(q+1)/2)

        time<- microbenchmark({
          for(i in 1: nobs){
          d2_logm_eta_row(eta[i,],y[i,],res,z,w)
          }
        } ,times=1L)
        out$tD2_eff[jj] <- time$time


        time<- microbenchmark({
          for(i in 1: nobs){
            d2_logm_eta_row_noalloc(eta[i,],y[i,],z,w)
          }
        } ,times=1L)
        out$tD2_eff_noalloc[jj] <- time$time


        if(TMB == TRUE){
        obj <- list()
        for(i in 1:nobs){
          obj[[i]] <- MakeADFun(list(Y=matrix(y[i,],ncol=d)), list(eta = matrix(eta[i,],ncol=q)), DLL = "nll_logM_TMB", silent=TRUE)
        }

        time<-microbenchmark({
          for(i in 1:nobs){
            HTMBi <- obj[[i]]$he()
          }
        }, times=1L)
        out$tD2_TMB[jj] <- time$time
        }
      }
    }

    return(list("time" = out))
  },  mc.cores = ncores)
  return(time)
}




############################################################
# Simulation
############################################################
nobs<-1000
dgrid <- seq(5,50,ny=5)

nrun <- 1  # Set the number of runs
ncores <- 1 # Set the number of cores
TIME_MCD_D2eta <- time_Deta(nobs, dgrid,  nrun, ncores, param="mcd", TMB = TRUE)
save(TIME_MCD_D2eta, file=paste0("TIME_mcd_D2eta_dgrid_min_", min(dgrid),"_max_",max(dgrid),".RData"))

TIME_logM_D2eta <- time_Deta(nobs, dgrid,  nrun, ncores, param ="logm", TMB = TRUE)
save(TIME_logM_D2eta, file=paste0("TIME_logm_D2eta_dgrid_min_", min(dgrid),"_max_",max(dgrid),".RData"))


