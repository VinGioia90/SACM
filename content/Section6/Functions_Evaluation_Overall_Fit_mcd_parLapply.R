#############################################################################
# Needed functions for generating the data and fit the models of Section  6 #
#############################################################################

########################################################################################
# This function is meant for data generation under the MCD parametrisation. Arguments: #
#  d: dimension of the outcome                                                         #
#  nobs_train: number of observations for the training set                             #
#  nobs_test: number of observations for the test set                                  #
#  pint: percentage of intercepts in covariance matrix modelling                       #
#  seed: seed for generations                                                          #
#  param - type of parametrisation, only "mcd" here                                    #
#                                                                                      #
# Output: for the training and  test set list of                                       #
#  data_train: (y,x) (matrix)                                                          #
#  mu_train: simulated mean vectors                                                    #
#  Sigma_train: simulated vcov (list of matrix)                                        #
#  Theta_train: simulated Theta (list of matrix)                                       #
#  data_test: (y,x) (matrix)                                                           #
#  mu_test: simulated mean vectors                                                     #
#  Sigma_test: simulated vcov (list of matrix)                                         #
#  Theta_test: simulated Theta (list of matrix)                                        #
#  idxint: indices intercepts (vector)                                                 #
########################################################################################
datagen <- function(d = 2, nobs_train = 1000, nobs_test = 1000, pint = 0, seed = 13, param = "mcd"){

  set.seed(seed)
  y <- matrix(0, nobs_train + nobs_test, d)
  x <- list()
  x$x1 <- runif(nobs_train + nobs_test)
  x$x2 <- runif(nobs_train + nobs_test)
  x$x3 <- runif(nobs_train + nobs_test)

  intercepts <- sample(0 : 1, d * (d + 1)/2, prob= c( pint, 1 - pint), replace = TRUE)

  ## Generating functions for the mean vectors and the Theta matrix
  f1m <- function(x1, m1) m1 * sin(pi * x1)
  f2m <- function(x2, m2) exp(m2 * x2)
  f3m <- function(x3, m3, m4, m5, m6) m3 * x3 ^ 11 * (m4 * (1 - x3)) ^ 6 + m5 * (m6 * x3) ^ 3 * (1 - x3) ^ 10
  fvc <- function(x1,x2,  vc0, vc1, vc2, vc3, vc4, vc5, vc6, ind){
    vc0 + ind * (vc1 * sin(2 * pi * (x1 + vc2)) + vc3 * cos(2 * pi * (x1 + vc2))
                 + vc4 * sin(2 * pi * (x2 + vc5)) + vc6 * cos(2 * pi * (x2 + vc5)))
  }

  # Setting the generating parameter values
  # For the meanvectors
  m1 <- runif(d, 1, 3)
  m2 <- runif(d, 1, 3)
  m3 <- runif(d, 0, 0.5)
  m4 <- runif(d, 9, 11)
  m5 <- runif(d, 9, 11)
  m6 <- runif(d, 9, 11)

  # For the entries of the (reparametrised) covariance matrix, that is for genereting the entries of Theta
  vc0 <- runif(d * (d + 1)/2, -0.25, 0.25)
  vc1 <- runif(d * (d + 1)/2, -0.5, 0.5)
  vc2 <- runif(d * (d + 1)/2, -0.5, 0.5)
  vc3 <- runif(d * (d + 1)/2, -0.25, 0.25)
  vc4 <- runif(d * (d + 1)/2, -1, 1)
  vc5 <- runif(d * (d + 1)/2, -1, 1)
  vc6 <- runif(d * (d + 1)/2, -0.5, 0.5)

  # Generation
  mu <- list()
  Sigma <- list()
  Theta <- list()

  for(i in 1 : (nobs_train + nobs_test)){
    mu[[i]] <- rep(0, d)
    for(j in 1 : d) mu[[i]][j] <- f1m(x$x1[i], m1[j]) + f2m(x$x2[i], m2[j]) +  f3m(x$x3[i], m3[j], m4[j], m5[j], m6[j])

    Theta[[i]] <- matrix(0, d, d)

    Theta[[i]] <- matrix(0, d, d)
    count <- d + 1
    for(j in 1 : d){
      Theta[[i]][j, j] <- fvc(x$x1[i], x$x2[i], vc0[j], vc1[j], vc2[j], vc3[j], vc4[j], vc5[j], vc6[j], intercepts[j])
      if(j > 1){
        for(k in 1 : (j - 1)){
          Theta[[i]][j, k] <- Theta[[i]][k, j] <- fvc(x$x1[i], x$x2[i], vc0[count], vc1[count], vc2[count], vc3[count], vc4[count], vc5[count], vc6[count], intercepts[count])
          count <- count + 1
        }
      }
    }
    # Converting the generated Theta entries into the sensible lpi vector ordering
    lpi <- c(mu[[i]],diag(Theta[[i]]), Theta[[i]][upper.tri(Theta[[i]], diag=FALSE)])

    # MCD parametrisation: data generation
    LD <- SCM::internal()$mcd_LD(lpi, d)
    L <- LD
    diag(L) <- rep(1,d)
    D <- diag(LD)
    Sigma[[i]] <- crossprod(t(L) * D) # L%*%diag(D^2)%*%t(L)
    u <- mvnfast::rmvn(1, rep(0, d), diag(rep(1, d)))
    y[i,] <- t(mu[[i]] + L %*% t(u * D))
  }

  colnames(y) <-  paste0("y_", 1 : d)
  X <- matrix(0, nobs_train + nobs_test, length(x))
  for(j in 1 : length(x)) X[,j] <-  x[[j]]
  colnames(X) <- names(x)
  data <- cbind(y, X)

  data_train <- data[1 : nobs_train,]
  data_test <- data[(nobs_train + 1) : (nobs_train + nobs_test), ]

  mu_train <- lapply(1 : nobs_train, function(x) mu[[x]])
  mu_test <- lapply((nobs_train + 1) : (nobs_train + nobs_test), function(x) mu[[x]])

  Theta_train <- lapply(1 : nobs_train, function(x) Theta[[x]])
  Theta_test <- lapply((nobs_train + 1) : (nobs_train + nobs_test), function(x) Theta[[x]])

  Sigma_train <- lapply(1 : nobs_train, function(x) Sigma[[x]])
  Sigma_test <- lapply((nobs_train + 1) : (nobs_train + nobs_test), function(x) Sigma[[x]])

  return(list(data_train = data.frame(data_train),
              mu_train = mu_train,
              Theta_train = Theta_train,
              Sigma_train = Sigma_train,
              data_test = data.frame(data_test),
              mu_test = mu_test,
              Theta_test = Theta_test,
              Sigma_test = Sigma_test,
              idxint = intercepts))
}


###################################################################
# Function for setting model formulas. Arguments;                 #
# d: dimension of the outcome                                     #
# expl_mean - vector of strings for the mean vector predictors    #
# expl_Theta - vector of strings for the mean vector predictors   #
# idxcov_int - indices for covariance modelling using intercepts  #
#                                                                 #
#  Output: list containing the model formulas                     #
###################################################################

mformula <- function(d = 2, expl_mean, expl_Theta, idxcov_int = rep(1, d * (d + 1)/2)){

  meanv_f <- paste0("y_", 1 : (d - 1), sep = "|", collapse = "") # First d - 1 outcome included in the model formula
  meanv_l <- paste0("y_", d, "~ s(", expl_mean[1], ") + s(", expl_mean[2],") + s(", expl_mean[3],")") # last part of the model formula (not beautiful:
  # only limited to this specific case)
  mean_foo <- formula(paste0(meanv_f, meanv_l)) # Complete mean model formula specification

  # Create the Theta labels (discriminating the intercepts from the covariate-dependent case)
  labelTh <- rep(0, d * (d + 1)/2)
  for(j in (d + 1) : (d + d * (d + 1)/2)) labelTh[j - d] <- SCM:::labTh(d, j)
  labelTh_i <- labelTh[idxcov_int != 1]
  labelTh_x <- labelTh[idxcov_int == 1]

  lint <- sum(idxcov_int != 1)
  lnoint <- sum(idxcov_int == 1)

  # f and l are meant for the first and lst part of the model formula, respectively
  # i and x are designed for capturing the possibility to consider intercepts and covariate-dependent formulas, respectively
  covi_foo_f <- covx_foo_f <- covx_foo_l <- covi_foo_l <- c()

  # Covariate-dependent case
  if( lnoint > 0){
    if(lnoint > 1) {
      covx_foo_f <- paste0(labelTh_x[1 : (lnoint - 1)], sep = "|", collapse = "")
    }
    covx_foo_l <- paste0(labelTh_x[lnoint], "~ s(", expl_Theta[1], ") + s(", expl_Theta[2], ")") #not beautiful: only limited to this specific case
    covx_foo <- formula(paste0(covx_foo_f , covx_foo_l))
  }

  # Intercept case
  if(lint > 0){
    if(lint > 1){
      covi_foo_f <- paste0(labelTh_i[1 : (lint - 1)], sep = "|", collapse = "")
    }
    covi_foo_l <- paste0(labelTh_i[lint], "~ 1" )
    covi_foo <- formula(paste0(covi_foo_f , covi_foo_l))
  }

  # Carrying out the final Theta model formula
  if(lnoint > 0 & lint == 0) cov_foo <- unlist(covx_foo)
  if(lnoint == 0 & lint > 0) cov_foo <- unlist(covi_foo)
  if(lnoint > 0 & lint > 0) cov_foo <- c(unlist(covx_foo), unlist(covi_foo))

  return(c(mean_foo, cov_foo))
}

################################################################
# Function for building the model formula for the bamlss model #
################################################################
# New version: to be improved
mformula_mvnchol <- function(d, expl_mean, expl_Theta){
  foo <- make_formula(as.formula(paste(paste(paste0("y_", 1 : d),  collapse = " | "),
                                       "~ s(", expl_mean[1], ") + s(", expl_mean[2],") + s(", expl_mean[3],") | s(", expl_Theta[1], ") + s(", expl_Theta[2], ") | s(", expl_Theta[1], ") + s(", expl_Theta[2], ")")))
  return(foo)
}


############################################################################################
# Function for fitting the models using efs/bfgs/bfgs-efs methods and bamlss. Arguments:   #
#   nobs_train: number of observations (training set)                                      #
#   nobs_test: number of observations (test set)                                           #
#   dgrid: grid of the outcome vector dimensions                                           #
#   nrun: number of replications for each d                                                #
#   ncores: number of cores                                                                #
#   param: type of parametrization, only "mcd" here                                        #
#   expl_mean: vector including the explanatory variables for mean modelling               #
#   expl_Theta: vector including the explanatory variables for covariance matrix modelling #
#   save.gam: if "TRUE" a gam object is saved; otherwise the predicted lpi are saved       #
#   blockN: a vector of the same dimension of dgrid,                                       #
#           allowing to specify the number of obesrvations' blocks                         #
#   pint: perc. of intercepts in covariance matrix modelling (not used in the simulations) #
#   pureBFGS: set to FALSE to consider only BFGS with the initialised version of the efs   #
#                                                                                          #
# Output: a list containing:                                                               #
#  -) the object resulting from the data generation                                        #
#  -) the model formulas                                                                   #
#  -) the fitting time                                                                     #
#  -) the number of iterations                                                             #
#  -) the gam object or the predicted variances and                                        #
#     covariances according to the flag save.gam                                           #
#     (both for the test and training set)                                                 #
#  -) the LAML for efs/bfgs/bfgs-efs (training)                                            #
############################################################################################
sim_est_efs_bfgs_bamlss <- function(nobs_train, nobs_test, dgrid,
                                    nrun, ncores, param = "mcd", expl_mean, expl_Theta,
                                    save.gam = "TRUE", blockN = rep(1, length(dgrid)),
                                    pint = 0, pureBFGS = "TRUE", root_dir){

  res_sim <- function(dgrid, nobs_train, nobs_test, param,  expl_mean, expl_Theta, save.gam, blockN, pint, seq_seed, pureBFGS){
    dss <- lapply(1 : length(dgrid), function(jj){
      out <- list()
      out$sim <- datagen(d = dgrid[jj], nobs_train = nobs_train, nobs_test = nobs_test,
                         pint = pint, seed = 13 * dgrid[jj] * seq_seed, param = param)
      out$foo <-  mformula(d = dgrid[jj], expl_mean = expl_mean, expl_Theta = expl_Theta, idxcov_int = out$sim$idxint)
      foo_mvnchol <-  mformula_mvnchol(d = dgrid[jj], expl_mean = expl_mean, expl_Theta = expl_Theta)
      if(save.gam == "TRUE"){
        time_bamlss <- microbenchmark(out$bamlss_fit <- bamlss(foo_mvnchol, family = mvnchol_bamlss(k = dgrid[jj], type = "modified"),
                                                               data = as.data.frame(out$sim$data_train), sampler = FALSE), times = 1L)
        out$time_bamlss <- time_bamlss$time


        time_efs <- microbenchmark(out$fit_efs <- gam_scm(out$foo,
                                                          family = mvn_scm(d = dgrid[jj], param = param, nb = blockN[jj]),
                                                          optimizer= "efs",  data = as.data.frame(out$sim$data_train)), times = 1L)
        out$time_efs <- time_efs$time
        out$outer_efs <- out$fit_efs$iter
        out$inner_efs <- out$fit_efs$family$getNC() - 1

        if(pureBFGS){
          time_bfgs <- microbenchmark(out$fit_bfgs <- gam_scm(out$foo,
                                                              family = mvn_scm(d = dgrid[jj], param = param, nb = blockN[jj]),
                                                              optimizer = "bfgs",  data = as.data.frame(out$sim$data_train)), times = 1L)
          out$time_bfgs <- time_bfgs$time
          out$outer_bfgs <- out$fit_bfgs$iter
          out$inner_bfgs <- out$fit_bfgs$family$getNC() - 1
        }

        time_bfgs_efs <- microbenchmark(out$fit_bfgs_efs <- gam_scm(out$foo, family = mvn_scm(d = dgrid[jj], param = param, nb = blockN[jj]),
                                                                    optimizer = "bfgs",  data = as.data.frame(out$sim$data_train),
                                                                    aGam = list(start = out$fit_efs$coef, in.out = list(sp = out$fit_efs$sp, scale = 1))), times = 1L)
        out$time_bfgs_efs <- time_bfgs_efs$time
        out$outer_bfgs_efs <- out$fit_bfgs_efs$iter
        out$inner_bfgs_efs <- out$fit_bfgs_efs$family$getNC() - 1

      } else {
        time_bamlss <- microbenchmark(bamlss_fit <- bamlss(foo_mvnchol, family = mvnchol_bamlss(k = dgrid[jj], type = "modified"),
                                                           data = as.data.frame(out$sim$data_train), sampler = FALSE), times = 1L)
        out$time_bamlss <- time_bamlss$time
        out$lpi_pred_bamlss <- predict(bamlss_fit)
        out$lpi_pred_test_bamlss <- predict(bamlss_fit, newdata = as.data.frame(out$sim$data_test))

        time_efs <- microbenchmark(fit_efs <- gam_scm(out$foo,
                                                      family = mvn_scm(d = dgrid[jj], param = param, nb = blockN[jj]),
                                                      optimizer= "efs",  data = as.data.frame(out$sim$data_train)), times = 1L)
        out$time_efs <- time_efs$time
        out$outer_efs<- fit_efs$iter
        out$inner_efs <- fit_efs$family$getNC() - 1
        out$LAML_efs <- fit_efs$gcv.ubre
        out$lpi_pred_efs <- predict(fit_efs)
        out$lpi_pred_test_efs <- predict(fit_efs, newdata = as.data.frame(out$sim$data_test))

        if(pureBFGS){
          time_bfgs <- microbenchmark(fit_bfgs <- gam_scm(out$foo,
                                                          family = mvn_scm(d = dgrid[jj], param = param, nb = blockN[jj]),
                                                          optimizer = "bfgs",  data = as.data.frame(out$sim$data_train)), times = 1L)

          out$time_bfgs <- time_bfgs$time
          out$outer_bfgs<- fit_bfgs$iter
          out$inner_bfgs <- fit_bfgs$family$getNC() - 1
          out$LAML_bfgs <- fit_bfgs$gcv.ubre

          out$lpi_pred_bfgs <- predict(fit_bfgs)
          out$lpi_pred_test_bfgs <- predict(fit_bfgs, newdata = as.data.frame(out$sim$data_test))
        }

        time_bfgs_efs <- microbenchmark(fit_bfgs_efs <- gam_scm(out$foo,
                                                                family = mvn_scm(d = dgrid[jj], param = param, nb = blockN[jj]),
                                                                optimizer = "bfgs",  data = as.data.frame(out$sim$data_train),
                                                                aGam = list(start = fit_efs$coef, in.out = list(sp = fit_efs$sp, scale = 1))), times = 1L)
        out$time_bfgs_efs <- time_bfgs_efs$time
        out$outer_bfgs_efs <- fit_bfgs_efs$iter
        out$inner_bfgs_efs <- fit_bfgs_efs$family$getNC() - 1
        out$LAML_bfgs_efs <- fit_bfgs_efs$gcv.ubre

        out$lpi_pred_bfgs_efs  <- predict(fit_bfgs_efs)
        out$lpi_pred_test_bfgs_efs  <- predict(fit_bfgs_efs, newdata = as.data.frame(out$sim$data_test))
      }
      return(out)
    })
    return(dss)
  }

  seq_seed <- (1:nrun)^2

  cl <- makePSOCKcluster(ncores)
  setDefaultCluster(cl)
  clusterExport(NULL, c("nobs_train", "nobs_test", "dgrid", "param", "expl_mean", "expl_Theta",
                        "save.gam", "blockN", "pint", "datagen", "mformula",
                        "mformula_mvnchol", "seq_seed", "pureBFGS", "res_sim", "root_dir"), envir = environment())
  clusterEvalQ(NULL, {
    library("mgcv", lib.loc=paste0(root_dir, "/my_library"))
    library(bamlss)
    library(mvnchol)
    library(SCM)
    library(microbenchmark)
    if(packageVersion("mgcv") != "9.0"){
      stop("Wrong version of mgcv!!")
    }
    library(Rcpp)
    sourceCpp(code = "
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
void over_writediag(arma::mat& E, const arma::vec& rho, const arma::vec& ind, const int k_sp) {
  int kk = k_sp - 1;
  int dd;
  for (int i = 0; i < ind.n_elem; ++i) {
    dd = ind(i) - 1;
    E(dd, dd) = exp(rho(kk) * 0.5);
  }
}
")
  })

  out_res_sim <- function(.x){
    out2 <- res_sim(dgrid = dgrid, nobs_train = nobs_train, nobs_test = nobs_test, param = param,
                    expl_mean = expl_mean, expl_Theta = expl_Theta,
                    save.gam = save.gam, blockN = blockN, pint = pint, seq_seed = seq_seed[.x], pureBFGS = pureBFGS)
    return(list(gen = out2))
  }

  environment(out_res_sim) <- .GlobalEnv

  res <- list()
  res <- parLapply(NULL, 1 : nrun, out_res_sim)

  stopCluster(cl)
  rm(cl)
  return(res)
}


