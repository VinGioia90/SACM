##############################################################################
# Needed functions for generating the data and fit the models of Section 3.3 #
##############################################################################

###################################################################
# Function for data generation. Arguments:                        #
#  d: dimension of the outcome                                    #
#  n: number of observations                                      #
#  pint: percentage of intercepts in covariance matrix modelling  #
#  seed: seed for generations                                     #
#  param - type of parametrisation, one between "mcd" and "logm"  #
#                                                                 #
# Output: list of                                                 #
#  data: (y,x) (matrix)                                           #
#  mu: simulated mu vectors (list of vectors)                     #
#  Sigma: simulated vcov (list of matrix)                         #
#  Theta: simulated Theta (list of matrix)                        #
#  idxint: indices intercepts (vector)                            #
###################################################################
datagen <- function(d = 2, n = 1000, pint = 0, seed = 13, param = c("mcd", "logm")){
  param <- match.arg(param)

  set.seed(seed)
  y <- matrix(0, n, d)
  x <- list()
  x$x1 <- runif(n)
  x$x2 <- runif(n)
  x$x3 <- runif(n)

  intercepts <- sample(0 : 1, d * (d + 1)/2, prob = c( pint, 1 - pint), replace = TRUE)

  ## Generating functions for the mean vectors and the Theta matrix
  f1m <- function(x1, m1) m1 * sin(pi * x1)
  f2m <- function(x2, m2) exp(m2 * x2)
  f3m <- function(x3, m3, m4, m5, m6) m3 * x3 ^ 11 * (m4 * (1 - x3)) ^ 6 + m5 * (m6 * x3) ^ 3 * (1 - x3) ^ 10
  fvc <- function(x1,x2,  vc0, vc1, vc2, vc3, vc4, vc5, vc6, ind){
    vc0 + ind * (vc1 * sin(2 * pi * (x1 + vc2)) + vc3 * cos(2 * pi * (x1 + vc2)) + vc4 * sin(2 * pi * (x2 + vc5)) + vc6 * cos(2 * pi * (x2 + vc5)))
  }

  # Setting the generating parameter values
  # For the mean vectors
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

  for(i in 1 : n){
    mu[[i]] <- rep(0, d)
    for(j in 1 : d) mu[[i]][j] <- f1m(x$x1[i], m1[j]) + f2m(x$x2[i], m2[j]) +  f3m(x$x3[i], m3[j], m4[j], m5[j], m6[j])

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
    lpi <- c(mu[[i]], diag(Theta[[i]]), Theta[[i]][upper.tri(Theta[[i]], diag = FALSE)])

    # MCD parametrisation: data generation
    if (param ==  "mcd" ) {
      LD <- SCM::internal()$mcd_LD(lpi, d)
      L <- LD
      diag(L) <- rep(1, d)
      D <- diag(LD)
      Sigma[[i]] <- crossprod(t(L) * D) # L%*%diag(D^2)%*%t(L)
      u <- mvnfast::rmvn(1, rep(0, d), diag(rep(1, d)))
      y[i,] <- t(mu[[i]] + L %*% t(u * D))
    }

    # logM parametrisation: data generation
    if(param ==  "logm" ) {
      Sigma[[i]] <- SCM::internal()$logM_Sigma(lpi, d)
      C <- chol(Sigma[[i]])
      u <- mvnfast::rmvn(1, rep(0, d), diag(rep(1, d)))
      y[i,] <- t(mu[[i]] + t(C) %*% t(u))
    }
  }

  # Build the data data.frame
  colnames(y) <-  paste0("y_", 1 : d)
  X <- matrix(0, n, length(x))
  for(j in 1 : length(x)) X[,j] <-  x[[j]]
  colnames(X) <- names(x)
  data <- cbind(y, X)

  return(list(data = data.frame(data),
              mu = mu,
              Theta = Theta,
              Sigma = Sigma,
              idxint = intercepts))
}


#################################################################
# Function for setting model formulas. Arguments:               #
# d: dimension of the outcome                                   #
# expl_mean - vector of strings for the mean vector predictors  #
# expl_Theta - vector of strings for the mean vector predictors #
# idxcov_int - indices for covariance modelling via intercepts  #
#                                                               #
# Output: list containing the model formulas                    #
#################################################################

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


#################################################################################
# Function for (simulating data and) fitting the models. Arguments:  a          #
# nobs: number of observations                                                  #
# dgrid: grid of orucome dimension values                                       #
# param1: parametrisation for data generation ("mcd" or "logm")                 #
# param2: parametrisation for model fitting ("mcd" or "logm")                   #
# expl_mean: vector of strings for the mean vector predictors                   #
# expl_Theta: vector of strings for the mean vector predictors                  #
# save.gam: flag for saving the gam objects or return only the predicted values #
# blockN: number of chunking blocks                                             #
# pint: percentage of intercepts                                                #
#                                                                               #
# Output: list containing                                                       #
# time: fitting times                                                           #
# outer: number of outer iterations                                             #
# inner: number of inner iterations                                             #
# fit: a gam object if save.gam = "TRUE"                                        #
# lpi_pred: the predicted linear predictor values                               #
#################################################################################

sim_est_efs <- function(nobs, dgrid,  nrun, ncores, param1, param2 = NULL, expl_mean, expl_Theta,
                         save.gam = "TRUE", blockN = rep(1, length(dgrid)), pint = 0, root_dir){
  if(is.null(param2)){
    param2 <- param1
  }

  res_sim <- function(dgrid, nobs, param1, param2,  expl_mean, expl_Theta, save.gam, blockN, pint, seq_seed){
    dss <- lapply(1 : length(dgrid), function(jj){
      out <- list()
      # Generating data
      out$sim <- datagen(d = dgrid[jj], n = nobs, pint = pint, seed = 13 * dgrid[jj] * seq_seed, param = param1)
      # Setting the model formula
      out$foo <-  mformula(d = dgrid[jj], expl_mean = expl_mean, expl_Theta = expl_Theta, idxcov_int = out$sim$idxint)
      if(save.gam == "TRUE"){ # Return a gam.object
        time <- microbenchmark(out$fit <- gam_scm(out$foo,
                                                  family = mvn_scm(d = dgrid[jj], param = param2, nb = blockN[jj]),
                                                  optimizer= "efs",  data = as.data.frame(out$sim$data)), times = 1L) #here we also saved gam fit
        out$time <- time$time                     # Computational times
        out$outer <- out$fit$iter                 # Outer iterations
        out$inner <- out$fit$family$getNC() - 1   # Inner iterations
      } else { ## Return only the predicted values
        time <- microbenchmark(fit <- gam_scm(out$foo,
                                              family=mvn_scm(d = dgrid[jj], param = param2, nb = blockN[jj]),
                                              optimizer= "efs",  data = as.data.frame(out$sim$data)), times = 1L)
        out$time <- time$time                     # Computational times
        out$outer<- fit$iter                      # Outer iterations
        out$inner <- fit$family$getNC() - 1       # Inner iterations
        out$lpi_pred <- predict(fit)              # Predicted linear predictors
      }
      return(out)
    })
    return(dss)
  }

  # Setting parallel computing
  seq_seed <- (1:nrun)^2

  cl <- makePSOCKcluster(ncores)
  setDefaultCluster(cl)
  clusterExport(NULL, c("nobs", "dgrid", "param1", "param2", "expl_mean", "expl_Theta",
                        "save.gam", "blockN", "pint", "datagen","mformula","seq_seed", "root_dir"), envir = environment())
  clusterEvalQ(NULL, {
    library("mgcv", lib.loc=paste0(root_dir, "/my_library"))
    library(SCM)
    library(microbenchmark)
    if(packageVersion("mgcv") != "9.0"){
      stop("Wrong version of mgcv!!")
    }
  })

  out_res_sim <- function(.x){
    out2 <- res_sim(dgrid = dgrid, nobs = nobs, param1 = param1, param2 = param2,
                    expl_mean = expl_mean, expl_Theta = expl_Theta,
                    save.gam = save.gam, blockN = blockN, pint = pint, seq_seed = seq_seed[.x])
    return(list(gen = out2))
  }

  res <- list()
  res <- parLapply(NULL, 1 : nrun, out_res_sim)

  stopCluster(cl)
  rm(cl)
  return(res)
}

