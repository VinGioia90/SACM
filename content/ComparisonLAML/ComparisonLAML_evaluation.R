###########################################################################################################################
# Code for reproducibility of the results of Section 2.4:                                                                 #
# Computational comparison between the old and the new formulation of the term involved in LAML-derivative based approach #
# Time for fomputing the third derivatives is also included                                                               #
###########################################################################################################################
library(parallel)
library(Rcpp)
library(microbenchmark)

library(rstudioapi)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
sourceCpp("d3_mcd_eta.cpp")
sourceCpp("idx_zwGt.cpp")

source("aux_idx_l3.R")
source("il3.R")
source("il3_no0_mcd.R")

time_d3_eta_dh_drho <- function(nobs, dgrid,  nrun, ncores, ncoef = 10){

  time <- mclapply(1:nrun, function(ii){
    out <- list()
    out$td3_mcd_eta<-rep(0,length(dgrid)) #time for evaluating 3rd derivatives
    out$td3_complete<-rep(0,length(dgrid)) #time for evaluating OLD version of tr(H^{-1}l^{bb\rho_u})
    out$td3_new<-rep(0,length(dgrid)) #time for evaluating NEW version of tr(H^{-1}l^{bb\rho_u})

    for(kk in 1: length(dgrid)){
      d <- dgrid[kk]

      #Set of indices
      z <- w <- t <- rep(0, (d * (d - 1)/2))
      Gm <- matrix(0, d - 1, d - 1)
      mode(Gm) <- mode(z) <- mode(w) <- mode(t) <- "integer"
      idx_zwGt(d, z, w, Gm, t)

      # Dimension of lpi
      no_eta <- d + d * (d + 1)/2

      # Build lpi indices
      jj <-list()
      count <- 1
      for(j in 1 : no_eta){
        jj[[j]] <- count : (count + ncoef - 1)
        count <- count + ncoef
      }

      idx_l3 <- il3(d)
      idx_neq0 <- aux_idx_l3(d,z,w,Gm)
      idx_jkq <- il3_no0_mcd(d, z, w)
      l3 <- matrix(0, nobs, d * (4 * d^2 + 3 * d + 2)/3)

      eta <- matrix(0, nobs, no_eta)
      y <- matrix(0, nobs, d)
      X <- matrix(0, nobs, no_eta*ncoef)
      p <- ncol(X)
      m <- ncol(eta) ## number of smoothing parameters: in this case there is a smooth effect for each eta
      d1eta <- matrix(0,nobs*no_eta,m)

      time <- microbenchmark({
        d3_mcd_eta(eta, y, l3, z, w,  Gm,  t,  idx_l3$h, idx_l3$h2, idx_l3$h3, idx_l3$idx3_1, idx_l3$idx3_2, idx_l3$idx3_3, idx_l3$idx3_4, idx_l3$idx3_5, idx_l3$idx3_6, idx_neq0)
      }, times=1L)
      out$td3_mcd_eta[kk] <- time$time

      d2b <- matrix(0,p, p) # this replace the call to fh in the original function

      #New R version

      time <- microbenchmark({
        d1H <- rep(0,m)
        for(t in 1:length(idx_jkq$idxl3)){
          a <- rowSums((X[,jj[[idx_jkq$idxj[t]]],drop=FALSE] %*% d2b[jj[[idx_jkq$idxj[t]]],jj[[idx_jkq$idxk[t]]]]) * X[,jj[[idx_jkq$idxk[t]]],drop=FALSE])
          mult <- if (idx_jkq$idxk[t]==idx_jkq$idxj[t]) 1 else 2
          for (l in 1:m) { ## sp loop
            v <- rep(0,nobs)
            for (q in 1:length(idx_jkq$idxl3[[t]])) {
              v <- v + l3[,idx_jkq$idxl3[[t]][q]] * d1eta[((idx_jkq$idxq[[t]][q]-1)*nobs + 1):(idx_jkq$idxq[[t]][q]*nobs),l]
            }
            d1H[l] <- d1H[l] + mult * (sum(a*v)) ## accumulate tr(Hp^{-1}dH/drho_l) ??part of??
          }
        }
      }, times=1L)
      out$td3_new[kk] <- time$time


      #old R -- complete version
      d2b <- matrix(0,p, p) # this was managed internally by mgcv
      time <- microbenchmark({
        d1H <- list()
        d1H_v <- rep(0,m)
        for (l in 1:m) {
          d1H[[l]] <- matrix(0,p,p)

          for(j in 1:length(idx_jkq$idxl3)){
            v <- rep(0,nobs)

            for (q in 1:length(idx_jkq$idxl3[[j]])) {
              v <- v + l3[,idx_jkq$idxl3[[j]][q]] * d1eta[((idx_jkq$idxq[[j]][q]-1)*nobs + 1):(idx_jkq$idxq[[j]][q]*nobs),l]
            }
            A <- crossprod(X[,jj[[idx_jkq$idxj[j]]],drop=FALSE],v*X[,jj[[idx_jkq$idxk[j]]],drop=FALSE])
            d1H[[l]][jj[[idx_jkq$idxj[j]]],jj[[idx_jkq$idxk[j]]]] <- d1H[[l]][jj[[idx_jkq$idxj[j]]],jj[[idx_jkq$idxk[j]]]] + A
            if (idx_jkq$idxk[j]>idx_jkq$idxj[j]) d1H[[l]][jj[[idx_jkq$idxk[j]]],jj[[idx_jkq$idxj[j]]]] <- d1H[[l]][jj[[idx_jkq$idxk[j]]],jj[[idx_jkq$idxj[j]]]] + t(A)
          }
          d1H_v[l] <- sum(diag(crossprod(d2b,d1H[[l]])))
        }
      }, times=1L)

      out$td3_complete[kk] <- time$time
    }
    return(list("time" = out))
  },  mc.cores = ncores)
  return(time)
}

nobs <- 1000
dgrid <- c(2,3,5)#,10,15,20)
nrun <- 1 # to set
ncores <- 1 # to set

TIME_d3 <- time_d3_eta_dh_drho(nobs, dgrid,  nrun,ncores)
save(TIME_d3, file="TIME_d3.RData")

## !!! Memory problems with large dimensions albeit with 128 GB and 7 cores
#nobs <- 1000
#dgrid <- 25
#nrun <- 7
#ncores <- 7

#TIME_d3 <- time_d3_eta_dh_drho (nobs, dgrid,  nrun,ncores)
#save(TIME_d3, file="TIME_d3_d25.RData")
