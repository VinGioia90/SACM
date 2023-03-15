###############################################
# Code for reproducibility of Section 2.4.1   #
###############################################
rm(list=ls())

library(SCM)
library(parallel)
library(microbenchmark)
library(matrixcalc)

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
#           -) idxint - indices intercepts (vector)                                   #
#######################################################################################
datagen_Sigma <- function(d = 2, n = 1000, pint = 0, seed = 13, range.var=c(1,2), range.cov=c(-0.5,0.5)){

  set.seed(seed)
  data <- matrix(0, n, d)
  x0 <- runif(n, 0, 1)
  sx0 <- sort(x0)

  intercepts <- sample(0 : 1, d * (d + 1)/2, prob= c( pint, 1 - pint), replace = TRUE)
  fmean <- function(x, m1, m2, m3, m4) m1 * x^11 * (m2 * (1 - x))^6 + m3 * (m4 * x)^3 * (1 - x)^10
  fvar <- function(x, c1, c2, c3, beta0var, ind) beta0var + ind * (c1 * sin(2 * pi * (x + c2)) + c3 * cos(2 * pi * (x+c2)))
  fcov <- function(x, c1, c2, c3, beta0cov, ind) beta0cov + ind * (c1 * sin(2 * pi * (x + c2)) + c3 * cos(2 * pi * (x+c2)))

  c1 <- runif(d * (d + 1)/2, -0.5, 0.5)
  c2 <- runif(d * (d + 1)/2, -0.5, 0.5)
  c3 <- runif(d * (d + 1)/2, -0.25, 0.25)
  beta0vcov <- c(runif(d, range.var[1], range.var[2]),runif(d * (d - 1)/2, range.cov[1], range.cov[2]))

  m1 <- runif(d, 0, 0.5)
  m2 <- runif(d, 9, 11)
  m3 <- runif(d, 9, 11)
  m4 <- runif(d, 9, 11)

  mu <- list()
  Sigma<- list()
  Theta_mcd <- list()
  Theta_logm <- list()
  pos_def <- rep(0, n)



  for(i in 1:n){
    mu[[i]] <- rep(0,d)
    for(j in 1:d) mu[[i]][j] <- fmean(sx0[i], m1[j], m2[j], m3[j], m4[j])

    Theta_mcd[[i]] <- matrix(0,d,d)
    Theta_logm[[i]] <- matrix(0,d,d)

    Sigma[[i]] <- matrix(0,d,d)
    Sigma[[i]][1,1] <- fvar(sx0[i], c1[1], c2[1], c3[1], beta0vcov[1],  intercepts[1])
    count <- d + 1
    for(j in 2:d){
      Sigma[[i]][j,j] <- fvar(sx0[i], c1[j], c2[j], c3[j],beta0vcov[j],  intercepts[j])
      for(k in 1 : (j - 1)){
        Sigma[[i]][j,k] <- Sigma[[i]][k,j] <- fcov(sx0[i], c1[count], c2[count], c3[count], beta0vcov[count],  intercepts[count])
        count <- count + 1
      }
    }

    pos_def[i] <-  is.positive.definite( Sigma[[i]], tol=1e-8)

    if(pos_def[i] == 1){
      C <- chol(Sigma[[i]])
      u <- mvnfast::rmvn(1, rep(0, d), diag(rep(1, d)))
      data[i,] <- t(mu[[i]] + t(C) %*% t(u))
      eigd <- eigen(Sigma[[i]])
      Theta_logm[[i]] <- eigd$vectors %*% diag(log(eigd$values)) %*%  t(eigd$vectors)

      D <- diag(C)
      L <- t(C)%*%diag(1/D)
      T <- solve(L)
      diag(Theta_mcd[[i]]) <- log(D^2)
      for(j in 2:d){
        for(k in 1:(j-1)){
          Theta_mcd[[i]][j,k] <- Theta_mcd[[i]][k,j] <- T[j,k]
        }
      }
    }

  }

  colnames(data) <-  paste0("y_", 1 : d)

  return(list(data = cbind(data, sx0), Sigma = Sigma,  Theta_mcd=Theta_mcd, Theta_logm = Theta_logm, mu=mu, idxint = intercepts, pos_def=pos_def))
}


#PROVA<-datagen_Sigma(d = 20, n = 50000, seed = 1, range.var=c(3,5), range.cov=c(-0.5,0.5)
#sum(PROVA$pos_def)


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



#######################################################################################################################
# Function for fitting the models with efs:                                                                           #
# -) input: nobs - sample size; dgrid - vector of dimension of the outcome;                                           #
#           nrun - number of replications for each d; ncores = number of cores                                        #
#           param1: type of parametrization ("mcd" or "logm") for data generation                                     #
#           param2: type of parametrization ("mcd" or "logm") for model fitting                                       #
#           save.gam: if "TRUE" a gam object is saved; otherwise the predicted variances and covariances  are saved   #
#           blockN: a vector of the same dimension of dgrid, allowing to specify the number of obesrvations' blocks   #
#           pint = percentage of intercepts                                                                           #
# -) output: a list containing:                                                                                       #
#            -) the object resulting from the data generation                                                         #
#            -) the model formula                                                                                     #
#            -) the fitting time                                                                                      #
#            -) the gam object or the predicted variances and covariances according to the flag save.gam              #
#######################################################################################################################
sim_est_efs <- function(nobs, dgrid,  nrun, ncores, param2=NULL, save.gam = "TRUE", blockN = rep(1,length(dgrid)), pint=0){
  sim1 <- mclapply(1:nrun, function(ii){

    dss <- lapply(1:length(dgrid), function(jj){
      out <- list()
      out$sim <- datagen_Sigma(d = dgrid[jj], n = nobs, pint = pint, seed= 13*dgrid[jj]*(ii^2),range.var=c(3,5), range.cov=c(-0.5,0.5))
      out$foo <-  mformula(d = dgrid[jj], expl="sx0", idxcov_int = out$sim$idxint )# rep(1,dgrid[jj]*(dgrid[jj]+1)/2))
      if(save.gam == "TRUE"){
        time <- microbenchmark(out$fit <- gam_scm(out$foo, family=mvn_scm(d = dgrid[jj], param = param2, nb = blockN[jj]), optimizer= "efs",  data=as.data.frame(out$sim$data)), times=1L)
        out$time <- time$time
        out$outer<- out$fit$iter
        out$inner <- out$fit$family$getNC()-1
      } else {
        time <- microbenchmark(fit <- gam_scm(out$foo, family=mvn_scm(d = dgrid[jj], param = param2, nb = blockN[jj]), optimizer= "efs",  data=as.data.frame(out$sim$data)), times=1L)
        out$time <- time$time
        out$outer<- fit$iter
        out$inner <- fit$family$getNC()-1
        fit$family$put_cflag(FALSE)
        out$vcov_pred <- predict(fit, type = "response")
      }
      return(out)
    })
    return(list("gen" = dss))
  },  mc.cores = ncores)
  return(sim1)
}

#############################################
# Generation from MCD
#############################################

# Fit with MCD
##################################
# small dimension setting
##################################
 dgrid <- c(2,3,5)
 nrun <-  1
 ncores <- 1  # vary according to the number of the cores
 nobs <- 5000
 sg <- FALSE
#
# # Fit with MCD
old <- Sys.time()
sim_mcd_small <- sim_est_efs(nobs, dgrid, nrun, ncores,  param2="mcd", save.gam=sg)
new <- Sys.time() - old
print(new)

save(sim_mcd_small,
     file=paste0("sim_Sigmamcd_nrun_",nrun,"_n_",nobs,"_d_",paste0(dgrid, collapse="_"),".RData"))

# Fit with logM
  old <- Sys.time()
  sim_logm_small <- sim_est_efs(nobs, dgrid, nrun, ncores, param2="logm", save.gam=sg)
  new <- Sys.time() - old
print(new)

save(sim_logm_small,
     file=paste0("sim_Sigmalogm_nrun_",nrun,"_n_",nobs,"_d_",paste0(dgrid, collapse="_"),".RData"))




##################################
# moderate dimension setting
##################################

#dgrid <- c(10,15)
#nrun <-  7
#ncores <- 7  # vary according to the number of the cores
#nobs <- 25000

#sg <- FALSE

 # Fit with MCD
 # old <- Sys.time()
 # sim_mcd_moderate <- sim_est_efs(nobs, dgrid, nrun, ncores, param2="mcd", save.gam=sg, blockN = rep(5,length(dgrid)))
 # new <- Sys.time() - old
 #
 # save(sim_mcd_moderate,
 #      file=paste0("sim_Sigmamcd_nrun_",nrun,"_n_",nobs,"_d_",paste0(dgrid, collapse="_"),".RData"))

#
# # Fit with logM
#old <- Sys.time()
#sim_logm_moderate <- sim_est_efs(nobs, dgrid, nrun, ncores,  param2="logm", save.gam=sg, blockN = rep(5,length(dgrid)))
#new <- Sys.time() - old

#save(sim_logm_moderate,
#      file=paste0("sim_Sigmalogm_nrun_",nrun,"_n_",nobs,"_d_",paste0(dgrid, collapse="_"),".RData"))

# dgrid <- 20
# nrun <-  7
# ncores <- 7  # vary according to the number of the cores
# nobs <- 25000
#
# sg <- FALSE

# # # Fit with MCD
# old <- Sys.time()
# sim_mcd_moderate <- sim_est_efs(nobs, dgrid, nrun, ncores, param2="mcd", save.gam=sg, blockN = rep(5,length(dgrid)))
# new <- Sys.time() - old

# save(sim_mcd_moderate,
#      file=paste0("sim_Sigmamcd_nrun_",nrun,"_n_",nobs,"_d_",paste0(dgrid, collapse="_"),".RData"))

#
# # Fit with logM
# old <- Sys.time()
# sim_logm_moderate <- sim_est_efs(nobs, dgrid, nrun, ncores,  param2="logm", save.gam=sg, blockN = rep(5,length(dgrid)))
# new <- Sys.time() - old
#
# save(sim_logm_moderate,
#     file=paste0("sim_Sigmalogm_nrun_",nrun,"_n_",nobs,"_d_",paste0(dgrid, collapse="_"),".RData"))
