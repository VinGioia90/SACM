#########################################################
# Code for reproducibility of the results of Section 4: #
# Hessian w.r.t. beta                                   #
#########################################################
sourceCpp("d2_beta_noeta.cpp")
sourceCpp("idx_zwGt.cpp")
source("idxHess_no0.R")
#new version: an ifelse statement allows to consider the intercepts' blocks strategy (as implemented in the family) or not
source("aux_idx_new.R")

# This function generates the quantities used in the function time_hessian_beta. Arguments
# d: dimension of the outcome vector
# pint: percentage of intercepts involved in covariance matrix modelling
# nobs: number of rows
# ncoef: number of coefficients for each linear predictor (excluding the intercepts' case)
# nb: number of observations' block
# param: mcd (1) and logm (2)
# block: intercepts block strategy (TRUE) and no block strategy (FALSE)
# seed: is used to get fixed some quantities during the comparison
get_quantity <- function(d, pint, nobs, ncoef, nb, param, block, seed = 123){

  #Set of indices
  z <- w <- t <- rep(0, (d * (d - 1)/2))
  Gm <- matrix(0, d - 1, d - 1)
  mode(Gm) <- mode(z) <- mode(w) <- mode(t) <- "integer"
  idx_zwGt(d, z, w, Gm, t)

  # Dimension of lpi
  no_eta <- d + d * (d + 1)/2

  # generate eta, y: it can be filled with zeros
  set.seed(seed)
  eta <- matrix(rnorm(nobs * no_eta), nobs, no_eta)
  y <- matrix(rnorm(nobs * d), nobs, d)

  # random generation of the intercepts in the covariance matrix model
  intercepts <- sample(0 : 1, d * (d + 1)/2, prob= c(pint, 1 - pint), replace = TRUE)
  intercepts <- c(rep(1, d), intercepts)

  # X generation: it can be filled with zeros (or considering ones in the interpcept indices)
  dimX <- sum(intercepts == 0) + ncoef * sum(intercepts == 1)
  X <- matrix(rnorm(nobs * dimX), nobs, dimX)

  # Build lpi indices
  jj <- list()
  count <- 1
  for(j in 1 : no_eta){
    if(intercepts[j] == 0){
      jj[[j]] <- count
      #X[,count]<-rep(1,nobs) # We avoid to place intercepts and we retain the generated X
      count <- count + 1
    } else {
      jj[[j]] <- count : (count + ncoef - 1)
      count <- count + ncoef
    }
  }


  #number of (non-redundant) elements different from zero in the hessian matrix
  if(param == 1){ # See the computational paper
    nHel <- d * (d ^ 2 + 15 * d + 2)/6
    if ( d > 2 ) nHel <- nHel + d * (d - 1) * (d - 2)/3
  }
  if(param == 2) nHel <- no_eta * (no_eta + 1)/2 #logm: no sparsity

  # Get the indices of the derivatives different from zero
  idx_jk <- idxHess_no0(no_eta, z, w, param)
  idx_aux <- aux_idx_new(jj, idx_jk, no_eta, block)

  # Set the number of blocks if the allocation of the 2nd derivatives matrix is more than 1GB
  check.size <- ((nobs %/% nb) * (nHel - idx_aux$llls) * 8/1e9) < 1
  while(check.size == FALSE){
    nb <- nb + 1
    check.size <- ((nobs %/% nb) * (nHel - idx_aux$llls) * 8/1e9) < 1
  }
  if(check.size == FALSE) print(nb)

  nlast <- nobs %% nb
  nset <- nobs %/% nb

  # Derivatives w.r.t. eta initialization
  if(nb > 1){
    l2 <- matrix(0, nset, nHel - idx_aux$llls)
  } else { # Such trick allows to pass the matrix l2 in the case nb = 1, avoiding cpp issues
    l2 <- matrix(0, 1, 1)
  }
  l2_v <- rep(0, idx_aux$llls)
  l2_v_l <- rep(0, idx_aux$llls)

  # Create the last block quantities
  if ( nlast == 0 ) {
    nobs_b <- nset
    idx_b_seq <- cumsum(rep(nset, nb - 1)) - 1
  } else {
    nobs_b <- nlast
    idx_b_seq <- cumsum(rep(nset, nb)) - 1
  }
  l2_l <- matrix(0, nobs_b, nHel - idx_aux$llls)
  idx_b <- c(-1, idx_b_seq, nobs - 1)

  lpi <- lapply(jj, function(x) x - 1) # trick for using the c++ functions

  return(list(X = X, eta = eta, y = y, lpi = lpi, K = no_eta,
              l2 = l2, l2_l = l2_l, l2_v = l2_v, l2_v_l = l2_v_l,
              idx_b = idx_b, z = z, w = w, Gm = Gm, t = t,
              idx_aux = idx_aux,  param = param, intercepts = intercepts, nb = nb))

}

############################################################################################
# Computational time comparison between the intercepts block strategy and the standard one #
############################################################################################
# This function calls the generated quantities function and then evaluate the derivatives w.r.t. beta for each scenario. Arguments:
# nobs: number of rows
# dgrid: grid of outcome vector dimensions
# nrun: number of runs
# ncores: number of cores
# pint: scenario of covariance matrix modelling sparsity
# ncoef: number of parameters for each linear predictor
# nb: number of blocks
# param: parametrisation (MCD, 1; logM, 2)
# pint_value: percentage of intercepts for the const scenario

time_hessian_beta <- function(nobs, dgrid,  nrun, ncores,
                              pint = c("dm05", "dm1", "dm2", "const"),
                              ncoef = 10, nb = 1, param = 1, pint_value = 1){
  pint <- match.arg(pint)
  pint_function <- switch(pint,
                          const = function(d) pint_value,
                          dm05 = function(d) 1 - 1/sqrt(d),
                          dm1 = function(d) 1 - 1/d,
                          dm2 = function(d) 1 - 1/(d ^ 2) )

  time <- mclapply(1 : nrun, function(ii){ #parLapply??? maybe it's not a big problem using mclapply here
    out <- list()

    out$thessian_block <- rep(0,length(dgrid))
    out$thessian_noblock <- rep(0,length(dgrid))

    for(jj in 1 : length(dgrid)){
       d <- dgrid[jj]
       seed <- sample(1 : nobs, 1)
       getq1 <- get_quantity(d = d, pint = pint_function(d), nobs = nobs, ncoef = ncoef,
                             nb = nb, param = param, block = TRUE, seed = seed) # via intercept blocks
       getq2 <- get_quantity(d = d, pint = pint_function(d), nobs = nobs, ncoef = ncoef,
                             nb = nb, param = param, block = FALSE, seed = seed) #
       attr(out, "nb1") <- getq1$nb
       p <- ncol(getq1$X)
       lbb <- matrix(0, p, p)

       time<- microbenchmark({
         d2_beta_noeta(getq1$X, getq1$eta, getq1$y, getq1$lpi, getq1$K,
                       lbb, getq1$l2, getq1$l2_v, getq1$l2_l, getq1$l2_v_l,
                       getq1$idx_b, getq1$z, getq1$w, getq1$Gm, getq1$t,
                       getq1$idx_aux$b1_eta, getq1$idx_aux$b1, getq1$idx_aux$b2,
                       getq1$idx_aux$b3, getq1$idx_aux$idx_b1_eta, getq1$idx_aux$idx_b3,
                       getq1$idx_aux$l2_el, getq1$idx_aux$l2_el2 , getq1$param)

      }, times = 1L)
      out$thessian_block[jj] <- time$time

      lbb <- matrix(0, p, p)
      time<- microbenchmark({
        d2_beta_noeta(getq2$X, getq2$eta, getq2$y, getq2$lpi, getq2$K,
                      lbb, getq2$l2, getq2$l2_v, getq2$l2_l, getq2$l2_v_l,
                      getq2$idx_b, getq2$z, getq2$w, getq2$Gm, getq2$t,
                      getq2$idx_aux$b1_eta, getq2$idx_aux$b1, getq2$idx_aux$b2,
                      getq2$idx_aux$b3, getq2$idx_aux$idx_b1_eta, getq2$idx_aux$idx_b3,
                      getq2$idx_aux$l2_el, getq2$idx_aux$l2_el2, getq2$param)

      } ,times = 1L)
      out$thessian_noblock[jj] <- time$time
    }
    return(list("time" = out))
  },  mc.cores = ncores)
  return(time)
}

# This function evaluate the derivatives w.r.t. beta for all the scenarios scenario. Arguments:
# nobs: number of rows
# dgrid: grid of outcome vector dimensions
# nrun: number of runs
# ncores: number of cores
# pint: vector of scenarios for the covariance matrix modelling sparsity
# ncoef: number of parameters for each linear predictor
# nb: number of blocks
# param: parametrisation (MCD, 1; logM, 2)
# pint_value: percentage of intercepts for the const scenario


get_time_results <- function(nobs, dgrid,  nrun, ncores, pint = NULL, ncoef = 10, nb = 1, param = 1, pint_value = 1){
  out <- list()
  for(j in 1 : length(pint)){
    if(pint[j] != "const"){
      out[[j]] <-  time_hessian_beta(nobs, dgrid,  nrun,ncores, pint = pint[j], ncoef, nb, param)
    } else {
      out[[j]] <-  time_hessian_beta(nobs, dgrid,  nrun,ncores, pint = pint[j], ncoef, nb, param, pint_value)
    }
  }
  return(out)
}




