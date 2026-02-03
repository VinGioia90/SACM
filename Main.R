rm(list=ls())

.rs.restartR()

#####################################################
# Code for reproducing the results of the paper:    #
# "Scalable Additive Covariance Matrix Models       #
#####################################################
# This code is intended to be run in RStudio
inst_pack <- installed.packages()
if ( !require("rstudioapi") ) {
  install.packages("rstudioapi")
}

# Needed to install custom version of mgcv
if ( !require("SparseChol") ) {
  install.packages("SparseChol")
}

library(rstudioapi)
root_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(root_dir)

# Install and load specific version of mgcv (uncomment if needed)
system("rm -rf ./my_library/mgcv")
system("R CMD build mgcv")
install.packages("mgcv_9.0.tar.gz", repos = NULL, type = "source", lib = "./my_library")
library("mgcv", lib.loc="./my_library")

if(packageVersion("mgcv") != "9.0"){
  stop("Wrong version of mgcv!!")
}

# Load the required packages
source("loadPackages.R")
instload_packages()

nrun <- 10  # Set the number of runs
ncores <- 10 # Set the number of cores

###############
# SECTION 3.3 #
###############

#######################################################################
# Evaluation of the second-order derivatives w.r.t. linear predictors #
# Generates the results for Figure 1 in the paper                     #
#######################################################################
setwd("content/Section3/Comp_logM_MCD_Hessian_eta")
source("Functions_Evaluation_Hessian_eta_parLapply.R")

tic <- proc.time()

nobs <- 1000
dgrid <- seq(5, 50, by = 5)

# mcd: 2nd derivatives w.r.t. eta
param <- "mcd"
TIME_MCD_D2eta <- time_Deta(nobs, dgrid,  nrun, ncores, param = param)

# logm: 2nd derivatives w.r.t. eta
param <- "logm"
TIME_logM_D2eta <- time_Deta(nobs, dgrid,  nrun, ncores, param = param)

setwd(root_dir)
setwd("content/Section3/Results")

save(TIME_MCD_D2eta,
     file = paste0("TIME_mcd_D2eta_dgrid_min_", min(dgrid), "_max_", max(dgrid), "nobs", nobs, ".RData"))
save(TIME_logM_D2eta,
     file = paste0("TIME_logm_D2eta_dgrid_min_", min(dgrid), "_max_", max(dgrid), "nobs", nobs, ".RData"))

rm("TIME_MCD_D2eta", "TIME_logM_D2eta")
gc()

toc1 <- proc.time() - tic
print(toc1)

#######################################
# Evaluation of the overall model fit
# Generates the results for Figure 2 in the paper
#######################################
setwd(root_dir)
setwd("content/Section3/Comp_logM_MCD_Fit")
source("Functions_Evaluation_Overall_Fit_parLApply.R")

dgrid <- c(2,5,10,15,20)
nobs <- 10000
sg <- FALSE # This avoids saving the gam object (saves memory)

setwd(root_dir)
setwd("content/Section3/Results")

tic <- proc.time()
#######################
# Generation of simulated data from MCD
#######################
# Fit with MCD
sim_mcdG_mcdF <- sim_est_efs(nobs, dgrid, nrun, ncores, param1 = "mcd", param2 = "mcd", save.gam = sg,
                             expl_mean = c("x1", "x2", "x3"), expl_Theta = c("x1", "x2"), root_dir = root_dir)


save(sim_mcdG_mcdF,
     file = paste0("sim_mcdG_mcdF_nrun_", nrun, "_n_", nobs, "_d_", paste0(dgrid, collapse = "_"), ".RData"))
rm("sim_mcdG_mcdF")
gc()

# Fit with logM
sim_mcdG_logmF <- sim_est_efs(nobs, dgrid, nrun, ncores, param1 = "mcd", param2 = "logm", save.gam = sg,
                              expl_mean = c("x1", "x2", "x3"), expl_Theta = c("x1", "x2"), root_dir = root_dir)
save(sim_mcdG_logmF,
     file = paste0("sim_mcdG_logmF_nrun_", nrun, "_n_", nobs, "_d_", paste0(dgrid, collapse = "_"), ".RData"))
rm("sim_mcdG_logmF")
gc()

########################
# Generation of simulated data from logM
########################
# Fit with MCD
sim_logmG_mcdF <- sim_est_efs(nobs, dgrid, nrun, ncores, param1 = "logm", param2 = "mcd", save.gam = sg,
                              expl_mean = c("x1", "x2", "x3"), expl_Theta = c("x1", "x2"), root_dir = root_dir)
save(sim_logmG_mcdF,
     file = paste0("sim_logmG_mcdF_nrun_", nrun, "_n_", nobs, "_d_", paste0(dgrid, collapse = "_"), ".RData"))
rm("sim_logmG_mcdF")
gc()

# Fit with logM
sim_logmG_logmF <- sim_est_efs(nobs, dgrid, nrun, ncores, param1 = "logm", param2 = "logm", save.gam = sg,
                               expl_mean = c("x1", "x2", "x3"), expl_Theta = c("x1", "x2"), root_dir = root_dir)
save(sim_logmG_logmF,
     file = paste0("sim_logmG_logmF_nrun_", nrun, "_n_", nobs, "_d_", paste0(dgrid, collapse = "_"), ".RData"))
rm("sim_logmG_logmF")
gc()

toc2 <- proc.time() - tic
print(toc2)

#############
# SECTION 4 #
#############

#########################################
# Evaluation of the Hessian w.r.t. beta #
#########################################
setwd(root_dir)
setwd("content/Section4")
source("Functions_Evaluation_Hessian_beta_parLapply.R")

nobs <- 1000
dgrid <- seq(10, 120, by = 10)
ncoef <- 10
pint_type <- c("dm05", "dm1", "dm1c2", "dm2", "const") # Different parsimonious scenarios

tic <- proc.time()
##################################
# MCD: Not included in the main paper but in SM B.3
##################################
# Mean vector NOT fixed to intercepts
int.Mean <- FALSE
TIME_MCD_beta_noMeanInt <- get_time_results(nobs, dgrid,  nrun, ncores, pint = pint_type,
                                            ncoef = ncoef, nb = 1, param = 1, pint_value = 0.99, int.Mean = int.Mean)
save(TIME_MCD_beta_noMeanInt, file =
       paste0("Results/TIME_mcd_beta_d",min(dgrid),"_",max(dgrid),"_nobs",nobs,"intMean", int.Mean,".RData"))
rm("TIME_MCD_beta_noMeanInt")
gc()

# Mean vector fixed to intercepts
int.Mean <- TRUE
TIME_MCD_beta_MeanInt <- get_time_results(nobs, dgrid,  nrun, ncores, pint = pint_type,
                                          ncoef = ncoef, nb = 1, param = 1, pint_value = 0.99, int.Mean = int.Mean)
save(TIME_MCD_beta_MeanInt, file =
       paste0("Results/TIME_mcd_beta_d",min(dgrid),"_",max(dgrid),"_nobs",nobs,"intMean", int.Mean,".RData"))
rm("TIME_MCD_beta_MeanInt")
gc()



###############################
# logM
###############################
# Mean vector NOT fixed to intercepts (figure 3 in the main paper)
int.Mean <- FALSE
TIME_logM_beta_noMeanInt <- get_time_results(nobs, dgrid,  nrun,ncores,
                                             pint = pint_type, ncoef = ncoef,
                                             nb = 1, param = 2, pint_value = 0.99, int.Mean = int.Mean)

save(TIME_logM_beta_noMeanInt,
     file = paste0("Results/TIME_logm_beta_d", min(dgrid), "_", max(dgrid), "_nobs", nobs,"intMean", int.Mean,".RData"))
rm("TIME_logM_beta_noMeanInt")
gc()

# Mean vector fixed to intercepts (not included in the main paper but in SM B.3)
int.Mean <- TRUE
TIME_logM_beta_MeanInt <- get_time_results(nobs, dgrid,  nrun,ncores,
                                           pint = pint_type, ncoef = ncoef,
                                           nb = 1, param = 2, pint_value = 0.99, int.Mean = int.Mean)

save(TIME_logM_beta_MeanInt,
     file = paste0("Results/TIME_logm_beta_d", min(dgrid), "_", max(dgrid), "_nobs", nobs,"intMean", int.Mean,".RData"))
rm("TIME_logM_beta_MeanInt")
gc()

toc3 <- proc.time() - tic
print(toc3)

#############
# SECTION 5 #
#############

############################################################################
# Performance Comparison within the MCD parametrisation (FS, EFS, BAMLSS)
# This code produces the results for Table 1 in the main paper
############################################################################
setwd(root_dir)
setwd("content/Section5")
source("Functions_Evaluation_Overall_Fit_mcd_parLapply.R")

tic <- proc.time()

# FS, EFS and BAMLSS compared in up to 10 dimensions
dgrid <- c(2, 5, 10)
nobs <- 10000
sg <- FALSE # This avoids saving the gam object (saves memory)

sim_mcd_fit <- sim_est_efs_bfgs_bamlss(nobs_train = nobs, nobs_test = nobs, dgrid,  nrun, ncores, param = "mcd",
                                       expl_mean = c("x1", "x2", "x3"), expl_Theta = c("x1", "x2"), save.gam = sg,
                                       root_dir = root_dir)

save(sim_mcd_fit,
     file = paste0("Results/sim_mcd_fit_nrun_", nrun, "_n_", nobs, "_d_", paste0(dgrid, collapse = "_"), ".RData"))
rm("sim_mcd_fit")
gc()

# Above d = 10 we use only FS and EFS
dgrid <- c(15, 20)
nobs <- 10000
sg <- FALSE # This avoids saving the gam object

sim_mcd_fit_fs_efs <- sim_est_fs_efs(nobs_train = nobs, nobs_test = nobs, dgrid,  nrun, ncores, param = "mcd",
                                     expl_mean = c("x1", "x2", "x3"), expl_Theta = c("x1", "x2"), save.gam = sg,
                                     root_dir = root_dir)

save(sim_mcd_fit_fs_efs,
     file = paste0("Results/sim_mcd_fit_fs_efs_nrun_", nrun, "_n_", nobs, "_d_", paste0(dgrid, collapse = "_"), ".RData"))
rm("sim_mcd_fit_fs_efs")
gc()

toc4 <- proc.time() - tic
print(toc4)

#############
# SECTION 6 #
#############
setwd(root_dir)
setwd("content/Section6")
source("Functions_Application.R")

################################
# Application to GEFCom14 Data #
################################
source("DataPreprocessing.R") # The dataset is GEF14_data

#head(GEF14_data)
# Set the length of train set (2005 - 2010)
n_train <- which(GEF14_data$year==2011)[1]-1

d <- 24  # The test can be done with d such that d + d(d+1)/2 is divisible by 5
save.gam <- FALSE

# Mean model formula
mean_formula <- list()
for(j in 0 : (d - 1)){
  mean_formula[[j + 1]] <-
    as.formula(paste0("load_h",j," ~ load24_h",j, "+ dow + s(doy, k = 20) + s(temp95_h", j,", k = 15) + s(temp_h", j," , k = 15) +  s(progtime, k = 4)"))
}

grid_length <- 5

###############################################
# Backward effect selection on the training set
###############################################

tic <- proc.time()

#######################
# MCD parametrisation #
#######################
param <- "mcd"
for(outcome in c("residuals", "response")){
  print(outcome)
  res_mcd <- stepw_res(param = param, d = d, grid_length = grid_length, mean_model_formula = mean_formula,
                       data_train = GEF14_data[1 : n_train, ], eff_vcov = "s(doy)",
                       metric = "p",  save.gam = save.gam, outcome = outcome)
  save(res_mcd, file = paste0("Results/res_stepwise_param", param, "_d_", d, "_lstep_", grid_length, "_", outcome, ".RData"))
  
  par(mfrow = c(1, 2))
  plot(rev(unlist(res_mcd$time_fit)))
  plot(rev(unlist(res_mcd$time_fit)/unlist(res_mcd$num_iter)))
  rm("res_mcd")
  gc()
}

toc5 <- proc.time() - tic
print(toc5)


########################
# logM parametrisation #
########################
tic <- proc.time()

param <- "logm"
for(outcome in c("residuals", "response")){
  print(outcome)
  res_logm <- stepw_res(param = param, d = d, grid_length = grid_length, mean_model_formula = mean_formula,
                        data_train = GEF14_data[1 : n_train, ], eff_vcov = "s(doy)",
                        metric = "p",  save.gam = save.gam, outcome = outcome)
  save(res_logm, file = paste0("Results/res_stepwise_param", param, "_d_", d, "_lstep_", grid_length, "_", outcome, ".RData"))
  
  par(mfrow = c(1, 2))
  plot(rev(unlist(res_logm$time_fit)))
  plot(rev(unlist(res_logm$time_fit)/unlist(res_logm$num_iter)))
  rm("res_logm")
  gc()
  
}

toc6 <- proc.time() - tic
print(toc6)



###########################################################################
# Validation procedure to select the number of effects to keep in the model
# Produces data for Figure 4 in the main paper
###########################################################################

# Here, we separating in-sample ("res_h") and out-of-sample ("res_out_h") residuals
res_mat <- matrix(0, nrow(GEF14_data), d)
GEF14_data_residuals <- cbind(res_mat, res_mat, GEF14_data)
colnames(GEF14_data_residuals)[1:d] <- paste0("res_h", 0:(d-1))
colnames(GEF14_data_residuals)[(d+1):(2*d)] <- paste0("res_out_h", 0:(d-1))

# Set the rolling origin forecasting splitting
ndat <- which(GEF14_data$year == 2011)[1] - 1
ndat2 <- dim(GEF14_data)[1]
sets <- c(0, floor(seq(ndat , ndat2, length.out = 12)))

###########################################
# Preparing dataset padded with residuals #
###########################################
ncores <- 12

flag_residuals <- TRUE
if(flag_residuals){
  cl <- makePSOCKcluster(ncores)
  setDefaultCluster(cl)
  clusterExport(NULL, c("sets", "GEF14_data_residuals", "mean_formula"), envir = environment())
  clusterEvalQ(NULL, {
    library(mgcv)
  })
  
  mod <- list()
  for( j in 1:(length(sets)-1)){
    clusterExport(NULL, c("j"), envir = environment())
    
    mod[[j]] <- parLapply(NULL, 1:d, function(ii){
      gam(mean_formula[[ii]], family = gaussian,
          data = GEF14_data_residuals[(sets[1]+1):(sets[j+1]),], optimizer = "efs",
          control = list(trace = FALSE))
    })
    
    for(ii in 1:d){
      tmp_nam <- paste0("res_h", ii - 1)
      tmp_ind <- (sets[j]+1):(sets[j+1])
      tmp_mod <- mod[[j]][[ii]]
      # In-sample raw residuals
      GEF14_data_residuals[tmp_ind, tmp_nam] <- GEF14_data_residuals[tmp_ind,paste0("load_h",ii-1)] - predict(tmp_mod)[tmp_ind]
      # Out-of-sample raw residuals
      if(j < (length(sets) - 1)){
        tmp_nam2 <- paste0("res_out_h", ii-1)
        tmp_ind2 <- (sets[j+1]+1):(sets[j+2])
        GEF14_data_residuals[tmp_ind2, tmp_nam2] <- GEF14_data_residuals[tmp_ind2,paste0("load_h",ii-1)] - predict(tmp_mod, newdata = GEF14_data_residuals[tmp_ind2, ])
      }
    }
    print(j)
  }
  stopCluster(cl)
  rm(list = c("mod"))
  gc()
  
  # Now we need to inflate the marginal variance of the in-sample residuals to match the variance
  # of the out-of-sample residuals. We are using the 2011 data (last year) to compute the inflation factor.
  # The out-of-sample residuals (2011) are left unperturbed, while residuals from previous years are inflated.
  tmp <- apply(GEF14_data_residuals[GEF14_data_residuals$year == 2011, paste0("res_out_h", 0:(d-1))], 2, sd) /
    apply(GEF14_data_residuals[ GEF14_data_residuals$year == 2011, paste0("res_h", 0:(d-1))], 2, sd)
  GEF14_data_residuals[ , paste0("res_h", 0:(d-1))] <- t(t(GEF14_data_residuals[ , paste0("res_h", 0:(d-1))])*tmp)
  
  GEF14_data_residuals[GEF14_data_residuals$year == 2011, paste0("res_h", 0:(d-1))] <-  GEF14_data_residuals[GEF14_data_residuals$year == 2011, paste0("res_out_h", 0:(d-1))]
  GEF14_data_residuals[ , paste0("res_out_h", 0:(d-1))] <- NULL
}

save(GEF14_data_residuals, file = "GEF14_data_residuals.RData")

# Set the rolling origin forecasting splitting (here we do not account for the 0 at the beginning - just for implementation purposes)
ndat <- which(GEF14_data$year == 2011)[1] - 1
ndat2 <- dim(GEF14_data)[1]
sets <- floor(seq(ndat , ndat2, length.out = 12))



#######################
# MCD parametrisation #
#######################
tic <- proc.time()

ncores <- 11 # It should be 11 due to the number of sets involved in rolling origin forecasting

param <- "mcd"

# By setting a lower and upper threshold for the number of effects involved in covariance matrix
# modelling we are able to update/append further runs
low_neff_vcov <- 0
upp_neff_vcov <- 150

for(outcome in c("residuals", "response")){
  
  load(paste0("Results/res_stepwise_param", param, "_d_", d, "_lstep_", grid_length, "_", outcome, ".RData"))
  
  idx_vcov <- which((unlist(lapply(res_mcd$foo, function(x) length(x)))-d) <= (upp_neff_vcov - low_neff_vcov))
  
  if(outcome == "response"){
    cv_mcd <- cross_val(obj = res_mcd, param = param, d = d, data = GEF14_data, idx_vcov = idx_vcov,
                        sets_eval = sets, ncores = ncores, save.gam = save.gam, root_dir = root_dir)
    save(cv_mcd, file = paste0("Results/cv_res_stepwise_param", param, "_d_", d, "_lstep_",
                               grid_length, "_low_thresh_", low_neff_vcov, "_upp_thresh_", upp_neff_vcov, "_", outcome, ".RData" ))
  }
  if(outcome == "residuals"){
    cv_mcd_residuals <- cross_val(obj = res_mcd, param = param, d = d, data = GEF14_data_residuals, idx_vcov = idx_vcov,
                                  sets_eval = sets, ncores = ncores, save.gam = save.gam, root_dir = root_dir)
    save(cv_mcd_residuals, file = paste0("Results/cv_res_stepwise_param", param, "_d_", d, "_lstep_",
                                         grid_length, "_low_thresh_", low_neff_vcov, "_upp_thresh_", upp_neff_vcov, "_", outcome, ".RData" ))
  }
}

toc7 <- proc.time() - tic

########################
# logM parametrisation #
########################
tic <- proc.time()

param <- "logm"
low_neff_vcov <- 0
upp_neff_vcov <- 150

for(outcome in c("residuals", "response")){
  
  load(paste0("Results/res_stepwise_param", param, "_d_", d, "_lstep_", grid_length, "_", outcome, ".RData"))
  
  idx_vcov <- which((unlist(lapply(res_logm$foo, function(x) length(x)))-d) <= (upp_neff_vcov - low_neff_vcov))
  
  if(outcome == "response"){
    cv_logm <- cross_val(obj = res_logm, param = param, d = d, data = GEF14_data, idx_vcov = idx_vcov,
                         sets_eval = sets, ncores = ncores, save.gam = save.gam, root_dir = root_dir)
    save(cv_logm, file = paste0("Results/cv_res_stepwise_param", param, "_d_", d, "_lstep_",
                                grid_length, "_low_thresh_", low_neff_vcov, "_upp_thresh_", upp_neff_vcov, "_", outcome, ".RData" ))
  }
  if(outcome == "residuals"){
    cv_logm_residuals <- cross_val(obj = res_logm, param = param, d = d, data = GEF14_data_residuals, idx_vcov = idx_vcov,
                                   sets_eval = sets, ncores = ncores, save.gam = save.gam, root_dir = root_dir)
    save(cv_logm_residuals, file = paste0("Results/cv_res_stepwise_param", param, "_d_", d, "_lstep_",
                                          grid_length, "_low_thresh_", low_neff_vcov, "_upp_thresh_", upp_neff_vcov, "_", outcome, ".RData" ))
  }
}

toc8 <- proc.time() - tic

print(c(toc1, toc2, toc3, toc4, toc5, toc6))