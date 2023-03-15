##################################################################################################
# Code for reproducibility of Section 2.5.2 - Comparison between different estimation approaches #
##################################################################################################
rm(list=ls())
library(SCM)
library(parallel)
library(microbenchmark)
library(bamlss)
library(mvnchol)

library(rstudioapi)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#######################################################################################
# Function for data generation:                                                       #
# -) input: d - dimension of the outcome; n - sample size; pint - % intercepts; seed; #
#           param - type of parametrisation, one between "mcd" and "logm"             #
# -) output: list of                                                                  #
#           -) data - (y,x) (matrix)                                                  #
#           -) Sigma - simulated vcov (list of matrix)                                #
#           -) Theta - simulated Theta (list of matrix)                               #
#           -) mu - simulated mu (list of vectors)
#           -) idxint - indices intercepts (vector)                                   #
#######################################################################################
datagen <- function(d = 2, n = 1000, pint = 0, seed = 13, param = NULL){

  set.seed(seed)
  data <- matrix(0, n, d)
  x0 <- runif(n, 0, 1)
  sx0 <- sort(x0)

  intercepts <- sample(0 : 1, d * (d + 1)/2, prob= c( pint, 1 - pint), replace = TRUE)
  fmean <- function(x, m1, m2, m3, m4) m1 * x^11 * (m2 * (1 - x))^6 + m3 * (m4 * x)^3 * (1 - x)^10
  fvcov <- function(x, c1, c2, c3, beta0c, ind) beta0c + ind * (c1 * sin(2 * pi * (x + c2)) + c3 * cos(2 * pi * (x+c2)))

  c1 <- runif(d * (d + 1)/2, -0.5, 0.5)
  c2 <- runif(d * (d + 1)/2, -0.5, 0.5)
  c3 <- runif(d * (d + 1)/2, -0.25, 0.25)
  beta0c <- runif(d * (d + 1)/2, -0.25, 0.25)

  m1 <- runif(d, 0, 0.5)
  m2 <- runif(d, 9, 11)
  m3 <- runif(d, 9, 11)
  m4 <- runif(d, 9, 11)

  mu <- list()
  Sigma<- list()
  Theta <- list()

  for(i in 1:n){
    mu[[i]] <- rep(0,d)
    for(j in 1:d) mu[[i]][j] <- fmean(sx0[i], m1[j], m2[j], m3[j], m4[j])

    Theta[[i]] <- matrix(0,d,d)
    Theta[[i]][1,1] <- fvcov(sx0[i], c1[1], c2[1], c3[1], beta0c[1],  intercepts[1])
    count <- d + 1
    for(j in 2:d){
      Theta[[i]][j,j] <- fvcov(sx0[i], c1[j], c2[j], c3[j], beta0c[j],  intercepts[j])
      for(k in 1 : (j - 1)){
        Theta[[i]][j,k] <- Theta[[i]][k,j] <- fvcov(sx0[i], c1[count], c2[count], c3[count], beta0c[count],  intercepts[count])
        count <- count + 1
      }
    }
    if ( param ==  "mcd" ) {
      D <- exp(0.5*diag(Theta[[i]]))
      Dm2 <-diag(D^(-2))

      T <- matrix(0, d, d)
      for(j in 2:d){
        for(k in 1:(j-1)){
          T[j,k] <- Theta[[i]][j,k]
        }
      }
      diag(T) <- rep(1,d)
      iSigma <- t(T)%*%Dm2%*%T
      Sigma[[i]] <- solve(iSigma)
      L <- solve(T)
      u <- mvnfast::rmvn(1, rep(0, d), diag(rep(1, d)))
      data[i,] <- t(mu[[i]] + L %*% t(u*D))
    }

    if(param ==  "logm" ) {
      lpi <- c(mu[[i]],diag(Theta[[i]]), Theta[[i]][upper.tri(Theta[[i]], diag=FALSE)])
      Sigma[[i]] <- SCM::internal()$logM_Sigma(lpi, d)
      C <- chol(Sigma[[i]])
      u <- mvnfast::rmvn(1, rep(0, d), diag(rep(1, d)))
      data[i,] <- t(mu[[i]] + t(C) %*% t(u))
    }
  }

  colnames(data) <-  paste0("y_", 1 : d)

  return(list(data = cbind(data, sx0), Sigma = Sigma, Theta=Theta, mu=mu, idxint = intercepts))
}

#########################################################################################
# Function for model formula  generation:                                               #
# -) input: d - dimension of the outcome; expl - string of the predictor (only one now) #
#           idxcov_int - indices for covariance modelling via intercepts                #
# -) output: list containing the model formula                                          #
#########################################################################################

mformula <- function(d=2, expl, idxcov_int=rep(1,d)){

  meanv_f <- paste0("y_", 1:(d-1), sep = "|", collapse="")
  meanv_l <- paste0("y_",d , "~ s(", expl,")")
  mean_foo <- formula(paste0(meanv_f, meanv_l))
  labelTh <- rep(0, d*(d+1)/2)

  for(j in (d+1):(d+d*(d+1)/2)){
    labelTh[j-d] <- SCM:::labTh(d,j)
  }

  labelTh_i <- labelTh[idxcov_int!=1]
  labelTh_x <- labelTh[idxcov_int==1]
  lint <- sum(idxcov_int!=1)
  lnoint <- sum(idxcov_int==1)
  covi_foo_f <- covx_foo_f <- covx_foo_l <- covi_foo_l <- c()
  if( lnoint > 0){
    if(lnoint > 1) {covx_foo_f <- paste0( labelTh_x[1:(lnoint - 1)], sep = "|", collapse = "" )}
    covx_foo_l <- paste0( labelTh_x[lnoint], "~ s(", expl,")")
    covx_foo <- formula(paste0(covx_foo_f , covx_foo_l))
  }
  if(lint > 0){
    if(lint > 1){ covi_foo_f <- paste0( labelTh_i[1:(lint - 1)], sep = "|", collapse = "" )}
    covi_foo_l <- paste0( labelTh_i[lint], "~ 1" )
    covi_foo <- formula(paste0(covi_foo_f , covi_foo_l))
  }

  if( lnoint > 0 & lint == 0) cov_foo <- unlist(covx_foo)
  if(lnoint == 0 & lint > 0)    cov_foo <- unlist(covi_foo)
  if(lnoint >0 & lint > 0)    cov_foo <- c(unlist(covx_foo), unlist(covi_foo))

  return(c(mean_foo,cov_foo))

}

form_mvnchol <- function(d){
  if(d==2) foo <-  make_formula(y_1 | y_2 ~ s(sx0) | s(sx0) | s(sx0))
  if(d==3) foo <-  make_formula(y_1 | y_2 | y_3 ~ s(sx0) | s(sx0) | s(sx0))
  if(d==5) foo <-  make_formula(y_1 | y_2 | y_3 | y_4 | y_5 ~ s(sx0) | s(sx0) | s(sx0))
  if(d==10) foo <- make_formula(y_1 | y_2 | y_3 | y_4 | y_5 | y_6 | y_7 | y_8 | y_9 | y_10 ~ s(sx0) | s(sx0) | s(sx0))
  if(d==15) foo <- make_formula(y_1 | y_2 | y_3 | y_4 | y_5 | y_6 | y_7 | y_8 | y_9 | y_10 | y_11 | y_12 | y_13 | y_14 | y_15~ s(sx0) | s(sx0) | s(sx0))
  return(foo)
}

#######################################################################################################################
# Function for fitting the models with efs and bfgs for the mcd parametrisation:                                      #
# -) input: nobs - sample size; dgrid - vector of dimension of the outcome; nrun - number of replications for each d; #
#           save.gam: if "TRUE" a gam object is saved; otherwise only the lpi and the vcov predicted are saved        #
#           blockN: a vector of the same dimension of dgrid, allowing to specify the number of obesrvations' blocks   #
# -) output: a list containing:                                                                                       #
#            -) the objects resulting from the data generation                                                        #
#            -) the model formula                                                                                     #
#            -) the fitting time                                                                                      #
#            -) the gam objects or the lpi and covpred predicted acccording to the specification of save.gam in input #
#               (here we return both the efs and the bfgs fitting)                                                    #
#######################################################################################################################

sim_est_mcd <- function(nobs, dgrid, nrun, ncores, save.gam = "TRUE", blockN = rep(1,length(dgrid)), pureBFGS="TRUE"){
  sim1 <- mclapply(1:nrun, function(ii){

    dss <- lapply(1:length(dgrid), function(jj){
      out <- list()
      out$sim <- datagen(d = dgrid[jj], n = nobs, pint = 0, seed= 23*dgrid[jj]*(ii^2), param = "mcd")
      out$foo <-  mformula(d = dgrid[jj], expl="sx0", idxcov_int = rep(1,dgrid[jj]*(dgrid[jj]+1)/2))
      foo_mvnchol <-  form_mvnchol(d = dgrid[jj])
      if(save.gam == "TRUE"){
        time_bamlss <- microbenchmark(out$bamlss_fit <- bamlss(foo_mvnchol, family = mvnchol_bamlss(k = dgrid[jj], type="modified"), data = as.data.frame(out$sim$data), sampler=FALSE), times=1L)
        out$time_bamlss <- time_bamlss$time

        time_efs <- microbenchmark(out$fit_efs <- gam_scm(out$foo, family=mvn_scm(d = dgrid[jj], param = "mcd", nb = blockN[jj]), optimizer= "efs",  data=as.data.frame(out$sim$data)), times=1L)
        out$time_efs <- time_efs$time
        out$outer_efs <- out$fit_efs$iter
        out$inner_efs <- out$fit_efs$family$getNC()-1
        if(pureBFGS == "TRUE"){
          time_bfgs <- microbenchmark(out$fit_bfgs <- gam_scm(out$foo, family=mvn_scm(d = dgrid[jj], param = "mcd", nb = blockN[jj]), optimizer= "bfgs",  data=as.data.frame(out$sim$data)), times=1L)
          out$inner_bfgs <- out$fit_bfgs$family$getNC()-1
          out$time_bfgs <- time_bfgs$time
          out$outer_bfgs <- out$fit_bfgs$iter

        }
        time_bfgs_efs <- microbenchmark(out$fit_bfgs_efs <- gam_scm(out$foo, family=mvn_scm(d = dgrid[jj], param = "mcd", nb = blockN[jj]), optimizer= "bfgs",  data=as.data.frame(out$sim$data),
                                                                    aGam=list(start=out$fit_efs$coef, in.out=list(sp=out$fit_efs$sp, scale=1))), times=1L)
        out$time_bfgs_efs <- time_bfgs_efs$time
        out$outer_bfgs_efs <- out$fit_bfgs_efs$iter
        out$inner_bfgs_efs <- out$fit_bfgs_efs$family$getNC()-1
      } else {
        time_bamlss <- microbenchmark(bamlss_fit <- bamlss(foo_mvnchol, family = mvnchol_bamlss(k = dgrid[jj], type="modified"), data = as.data.frame(out$sim$data), sampler=FALSE), times=1L)
        out$time_bamlss <- time_bamlss$time

        time_efs <- microbenchmark(fit_efs <- gam_scm(out$foo, family=mvn_scm(d = dgrid[jj], param = "mcd", nb = blockN[jj]), optimizer= "efs",  data=as.data.frame(out$sim$data)), times=1L)
        out$time_efs <- time_efs$time
        out$outer_efs  <- fit_efs$iter
        out$inner_efs  <- fit_efs$family$getNC()-1
        if(pureBFGS == "TRUE"){
          time_bfgs <- microbenchmark(fit_bfgs <- gam_scm(out$foo, family=mvn_scm(d = dgrid[jj], param = "mcd", nb = blockN[jj]), optimizer= "bfgs",  data=as.data.frame(out$sim$data)), times=1L)
          out$time_bfgs <- time_bfgs$time
          out$outer_bfgs <- fit_bfgs$iter
          out$inner_bfgs <- fit_bfgs$family$getNC()-1
          out$lpi_pred_bfgs <- predict(fit_bfgs)
          fit_bfgs$family$put_cflag(FALSE)
          out$vcov_pred_bfgs <- predict(fit_bfgs, type = "response")

        }
        # For bfgs estimation we use the starting values of the efs
        time_bfgs_efs <- microbenchmark(fit_bfgs_efs <- gam_scm(out$foo, family=mvn_scm(d = dgrid[jj], param = "mcd", nb = blockN[jj]), optimizer= "bfgs",  data=as.data.frame(out$sim$data),
                                                                aGam=list(start=fit_efs$coef, in.out=list(sp=fit_efs$sp, scale=1))), times=1L)

        out$time_bfgs_efs <- time_bfgs_efs$time
        out$outer_bfgs_efs <- fit_bfgs_efs$iter
        out$inner_bfgs_efs <- fit_bfgs_efs$family$getNC()-1

        out$lpi_pred_bamlss <- predict(bamlss_fit)
        out$lpi_pred_efs <- predict(fit_efs)
        #out$lpi_pred_bfgs <- predict(fit_bfgs)
        out$lpi_pred_bfgs_efs <- predict(fit_bfgs_efs) #it was wrong, so use vcov_pred_bfgs_efs

        fit_efs$family$put_cflag(FALSE)
        #fit_bfgs$family$put_cflag(FALSE)
        fit_bfgs_efs$family$put_cflag(FALSE)
        out$vcov_pred_efs <- predict(fit_efs, type = "response")
        #out$vcov_pred_bfgs <- predict(fit_bfgs, type = "response")
        out$vcov_pred_bfgs_efs <- predict(fit_bfgs_efs, type = "response")
      }
      return(out)
    })

    return(list("gen" = dss))
  },  mc.cores = ncores)
  return(sim1)
}



##################################
# small dimension setting
#
dgrid <- c(2,3,5)
nrun <- 1 # vary according to the number of the cores
ncores <- 1 # vary according to the number of the cores
nobs <- 1000

sg <- "FALSE"

old <- Sys.time()
sim_mcd_bfgs <- sim_est_mcd(nobs, dgrid, nrun, ncores, save.gam=sg)
new <- Sys.time() -old

save(sim_mcd_bfgs,
    file=paste0("MCD_bfgs_nrun_",nrun,"_n_",nobs,"_d_",paste0(dgrid, collapse="_"), "savegam",sg ,".RData"))

#  dgrid <- c(2,3,5)
#  nrun <- 7 # vary according to the number of the cores
#  ncores <- 7 # vary according to the number of the cores
#  nobs <- 5000
#  sg <- "FALSE"
#
#  old <- Sys.time()
#  sim_mcd_bfgs <- sim_est_mcd(nobs, dgrid, nrun, ncores, save.gam=sg)
#  new <- Sys.time() -old
#
#
#  save(sim_mcd_bfgs,
#       file=paste0("MCD_bfgs_nrun_",nrun,"_n_",nobs,"_d_",paste0(dgrid, collapse="_"), "savegam",sg ,".RData"))
#
#
# dgrid <- 10
# nrun <-  7
# ncores <- 7 # vary according to the number of the cores
# nobs <- 5000
# sg <- "FALSE"
# #
# old <- Sys.time()
# sim_mcd_bfgs <- sim_est_mcd(nobs, dgrid, nrun, ncores, save.gam=sg, blockN = rep(1,length(dgrid)))
# new <- Sys.time() -old
#
#
#  save(sim_mcd_bfgs,
#        file=paste0("MCD_bfgs_nrun_",nrun,"_n_",nobs,"_d_",paste0(dgrid, collapse="_"), "savegam",sg ,".RData"))
#
# dgrid <- 10
# nrun <-  7
# ncores <- 7 # vary according to the number of the cores
# nobs <- 25000
# sg <- "FALSE"
#
# old <- Sys.time()
# sim_mcd_bfgs <- sim_est_mcd(nobs, dgrid, nrun, ncores, save.gam=sg, blockN = rep(5,length(dgrid)))
# new <- Sys.time() -old
# save(sim_mcd_bfgs,
#      file=paste0("MCD_bfgs_nrun_",nrun,"_n_",nobs,"_d_",paste0(dgrid, collapse="_"), "savegam",sg ,".RData"))



