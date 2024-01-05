rm(list=ls())

.rs.restartR()

#####################################################
# Code for reproducing the results of the paper:    #
# "Scalable Additive Covariance Matrix Models       #
#####################################################
library(rstudioapi)
root_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(root_dir)

# Install and load specific version of mgcv
install.packages("mgcv_9.0.tar.gz", repos = NULL, type = "source", lib = "./my_library")
library("mgcv", lib.loc="./my_library")
#library("mgcv")

# Load the needed packages
# (it might be required to do something manually)
source("loadPackages.R")
instload_packages()

nrun <- 10  # Set the number of runs
ncores <- 10 # Set the number of cores

###############
# SECTION 3.3 #
###############

#######################################################################
# Evaluation of the second-order derivatives w.r.t. linear predictors #
#######################################################################
setwd("content/Section3/Comp_logM_MCD_Hessian_eta")
source("Functions_Evaluation_Hessian_eta_parLapply.R")

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

#######################################
# Evaluation of the overall model fit #
#######################################
setwd(root_dir)
setwd("content/Section3/Comp_logM_MCD_Fit")
source("Functions_Evaluation_Overall_Fit_parLApply.R")


dgrid <- c(2,5,10,15,20)
nobs <- 10000
sg <- FALSE # This avoids saving the gam object

setwd(root_dir)
setwd("content/Section3/Results")

#######################
# Generation from MCD #
#######################
# Fit with MCD
sim_mcdG_mcdF <- sim_est_efs(nobs, dgrid, nrun, ncores, param1 = "mcd", param2 = "mcd", save.gam = sg,
                              expl_mean = c("x1", "x2", "x3"), expl_Theta = c("x1", "x2"))


save(sim_mcdG_mcdF,
    file = paste0("sim_mcdG_mcdF_nrun_", nrun, "_n_", nobs, "_d_", paste0(dgrid, collapse = "_"), ".RData"))
rm("sim_mcdG_mcdF")
gc()

# Fit with logM
sim_mcdG_logmF <- sim_est_efs(nobs, dgrid, nrun, ncores, param1 = "mcd", param2 = "logm", save.gam = sg,
                              expl_mean = c("x1", "x2", "x3"), expl_Theta = c("x1", "x2"))
save(sim_mcdG_logmF,
     file = paste0("sim_mcdG_logmF_nrun_", nrun, "_n_", nobs, "_d_", paste0(dgrid, collapse = "_"), ".RData"))
rm("sim_mcdG_logmF")
gc()

########################
# Generation from logM #
########################
# Fit with MCD
sim_logmG_mcdF <- sim_est_efs(nobs, dgrid, nrun, ncores, param1 = "logm", param2 = "mcd", save.gam = sg,
                              expl_mean = c("x1", "x2", "x3"), expl_Theta = c("x1", "x2"))
save(sim_logmG_mcdF,
     file = paste0("sim_logmG_mcdF_nrun_", nrun, "_n_", nobs, "_d_", paste0(dgrid, collapse = "_"), ".RData"))
rm("sim_logmG_mcdF")
gc()

# Fit with logM
sim_logmG_logmF <- sim_est_efs(nobs, dgrid, nrun, ncores, param1 = "logm", param2 = "logm", save.gam = sg,
                              expl_mean = c("x1", "x2", "x3"), expl_Theta = c("x1", "x2"))
save(sim_logmG_logmF,
     file = paste0("sim_logmG_logmF_nrun_", nrun, "_n_", nobs, "_d_", paste0(dgrid, collapse = "_"), ".RData"))
rm("sim_logmG_logmF")
gc()


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
dgrid <- seq(5, 100, by = 5)
ncoef <- 10
# S1 and S2 corresponds to dm05 and dm2
pint_type <- c("dm05", "dm1", "dm2", "const")

##################################
# MCD: Not included in the paper #
##################################
TIME_MCD_beta <- get_time_results(nobs, dgrid,  nrun, ncores, pint = pint_type,
                                  ncoef = ncoef, nb = 1, param = 1, pint_value = 0.99)
save(TIME_MCD_beta, file = paste0("Results/TIME_mcd_beta_d",min(dgrid),"_",max(dgrid),"_nobs",nobs,".RData"))
rm("TIME_MCD_beta")
gc()

###############################
# logM: included in the paper #
###############################
TIME_logM_beta <- get_time_results(nobs, dgrid,  nrun,ncores,
                                   pint = pint_type, ncoef = ncoef,
                                   nb = 1, param = 2, pint_value = 0.99)

save(TIME_logM_beta,
     file = paste0("Results/TIME_logm_beta_d", min(dgrid), "_", max(dgrid), "_nobs", nobs, ".RData"))

rm("TIME_logM_beta")
gc()


#############
# SECTION 6 #
#############

###################################################################################################
# Performance Comparison inside the MCD parametrisation (FS, BFGS, BFGS initialised EFS, BAMLSS)  #
###################################################################################################
setwd(root_dir)
setwd("content/Section6")
source("Functions_Evaluation_Overall_Fit_mcd_parLapply.R")


dgrid <- c(2, 5, 10)
nobs <- 10000
sg <- FALSE # This avoids saving the gam object

sim_mcd_fit <- sim_est_efs_bfgs_bamlss(nobs_train = nobs, nobs_test = nobs, dgrid,  nrun, ncores, param = "mcd",
                                       expl_mean = c("x1", "x2", "x3"), expl_Theta = c("x1", "x2"), save.gam = sg)

save(sim_mcd_fit,
     file = paste0("Results/sim_mcd_fit_nrun_", nrun, "_n_", nobs, "_d_", paste0(dgrid, collapse = "_"), ".RData"))
rm("sim_mcd_fit")
gc()


#############
# SECTION 7 #
#############
setwd(root_dir)
setwd("content/Section7")
source("Functions_Application.R")

################################
# Application to GEFCom14 Data #
################################
source("DataPreprocessing.R") # The dataset is GEF14_data

# Set the length of train set (2005 - 2010)
n_train <- which(GEF14_data$year==2011)[1]-1

d <- 24
save.gam <- FALSE

# Mean model formula
mean_formula <- list()
for(j in 0 : (d - 1)){
  mean_formula[[j + 1]] <- as.formula(paste0("load_h",j," ~ load24_h",j, "+ dow + s(doy) + s(temp95_h", j,")" ))
}

grid_length <- 5 # To change eventually

####################################
# Model selection on the train set #
####################################

#######################
# MCD parametrisation #
#######################
param <- "mcd"
res_mcd <- stepw_res(param = param, d = d, grid_length = grid_length, mean_model_formula = mean_formula,
                     data_train = GEF14_data[1 : n_train, ], eff_vcov = "s(doy)",
                     metric = "p",  save.gam = save.gam)
save(res_mcd, file = paste0("Results/res_stepwise_param", param, "_d_", d, "_lstep_", grid_length, ".RData"))
rm("res_mcd")
gc()

########################
# logM parametrisation #
########################
param <- "logm"
res_logm <- stepw_res(param = param, d = d, grid_length = grid_length, mean_model_formula = mean_formula,
                      data_train = GEF14_data[1 : n_train, ], eff_vcov = "s(doy)",
                      metric = "p",  save.gam = save.gam)
save(res_logm, file = paste0("Results/res_stepwise_param", param, "_d_", d, "_lstep_", grid_length, ".RData"))
rm("res_logm")
gc()


##############################
# Model Validation procedure #
##############################
ncores <- 11 # It should be considered 11 due to the number of sets involved in rolling origin forecasting

# Set the rolling origin forecasting splitting
ndat <- which(GEF14_data$year == 2011)[1] - 1
ndat2 <- dim(GEF14_data)[1]
sets <- floor(seq(ndat , ndat2, length.out = 12))

#######################
# MCD parametrisation #
#######################
setwd(root_dir)
setwd("content/Section7")
param <- "mcd"
load(paste0("Results/res_stepwise_param", param, "_d_", d, "_lstep_", grid_length, ".RData"))
cv_mcd <- cross_val(obj = res_mcd, param = param, d = d, data = GEF14_data, sets_eval = sets, ncores = ncores, save.gam = save.gam)
save(cv_mcd, file = paste0("Results/cv_res_stepwise_param", param, "_d_", d, "_lstep_", grid_length, ".RData"))

########################
# logM parametrisation #
########################
param <- "logm"
load(paste0("Results/res_stepwise_param", param, "_d_", d, "_lstep_", grid_length, ".RData"))
cv_logm <- cross_val(obj = res_logm, param = param, d = d, data = GEF14_data, sets_eval = sets, ncores = ncores, save.gam = save.gam)
save(cv_logm, file = paste0("Results/cv_res_stepwise_param", param, "_d_", d, "_lstep_", grid_length, ".RData"))

