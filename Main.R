#####################################################
# Code for reproducing the results of the paper:    #
# "Scalable Additive Covariance Matrix Models       #
#####################################################
library(rstudioapi)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Load the needed packages
# (it might be required to do something manually)
source("loadPackages.R")
instload_packages()

nrun <- 10  # Set the number of runs
ncores <- 5 # Set the number of cores

###############
# SECTION 3.3 #
###############

#######################################################################
# Evaluation of the second-order derivatives w.r.t. linear predictors #
#######################################################################
setwd("content/Section3/Comp_logM_MCD_Hessian_eta")
source("Functions_Evaluation_Hessian_eta.R")

nobs <- 1000
dgrid <- seq(5, 50, by = 5)

TIME_MCD_D2eta <- time_Deta(nobs, dgrid,  nrun, ncores, param = "mcd")
TIME_logM_D2eta <- time_Deta(nobs, dgrid,  nrun, ncores, param = "logm")

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("content/Section3/Results")

save(TIME_MCD_D2eta,
     file = paste0("TIME_mcd_D2eta_dgrid_min_", min(dgrid), "_max_", max(dgrid), "nobs", nobs, ".RData"))
save(TIME_logM_D2eta,
     file = paste0("TIME_logm_D2eta_dgrid_min_", min(dgrid), "_max_", max(dgrid), "nobs", nobs, ".RData"))

#######################################
# Evaluation of the overall model fit #
#######################################
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("content/Section3/Comp_logM_MCD_Fit")
source("Functions_Evaluation_Overall_Fit.R")


dgrid <- c(2,5,10,15,20)
nobs <- 10000
sg <- FALSE # This avoids saving the gam objkect

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("content/Section3/Results")

# Generation from MCD
# Fit with MCD
sim_mcdG_mcdF <- sim_est_efs(nobs, dgrid, nrun, ncores, param1 = "mcd", param2 = "mcd", save.gam = sg,
                               expl_mean = c("x1", "x2", "x3"), expl_Theta = c("x1", "x2"))
save(sim_mcdG_mcdF,
     file = paste0("sim_mcdG_mcdF_nrun_", nrun, "_n_", nobs, "_d_", paste0(dgrid, collapse = "_"), ".RData"))
rm("sim_mcdG_mcdF")

# Fit with logM
sim_mcdG_logmF <- sim_est_efs(nobs, dgrid, nrun, ncores, param1 = "mcd", param2 = "logm", save.gam = sg,
                              expl_mean = c("x1", "x2", "x3"), expl_Theta = c("x1", "x2"))
save(sim_mcdG_logmF,
     file = paste0("sim_mcdG_logmF_nrun_", nrun, "_n_", nobs, "_d_", paste0(dgrid, collapse = "_"), ".RData"))
rm("sim_mcdG_logmF")

# Generation from logM
# Fit with MCD
sim_logmG_mcdF <- sim_est_efs(nobs, dgrid, nrun, ncores, param1 = "logm", param2 = "mcd", save.gam = sg,
                             expl_mean = c("x1", "x2", "x3"), expl_Theta = c("x1", "x2"))
save(sim_logmG_mcdF,
     file = paste0("sim_logmG_mcdF_nrun_", nrun, "_n_", nobs, "_d_", paste0(dgrid, collapse = "_"), ".RData"))
rm("sim_logmG_mcdF")

# Fit with logM
sim_logmG_logmF <- sim_est_efs(nobs, dgrid, nrun, ncores, param1 = "logm", param2 = "logm", save.gam = sg,
                              expl_mean = c("x1", "x2", "x3"), expl_Theta = c("x1", "x2"))
save(sim_logmG_logmF,
     file = paste0("sim_logmG_logmF_nrun_", nrun, "_n_", nobs, "_d_", paste0(dgrid, collapse = "_"), ".RData"))
rm("sim_logmG_logmF")


#############
# SECTION 4 #
#############
# Evaluation of the Hessian w.r.t. beta
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("content/Section4")
source("Functions_Evaluation_Hessian_beta.R")

nobs <- 1000
dgrid <- seq(5, 100, by = 5)
ncoef <- 10
# S1 and S2 corresponds to dm05 and dm2
pint_type <- c("dm05", "dm1", "dm2", "const")

# MCD: Not reported in the paper
#TIME_MCD_beta <- get_time_results(nobs, dgrid,  nrun, ncores, pint = pint_type,
#                                  ncoef = ncoef, nb = 1, param = 1, pint_value = 0.99)
#save(TIME_MCD_beta, file = paste0("TIME_mcd_beta_d",min(dgrid),"_",max(dgrid),"_nobs",nobs,".RData"))


# logM: reported in the paper
TIME_logM_beta <- get_time_results(nobs, dgrid,  nrun,ncores,
                                   pint = pint_type, ncoef = ncoef,
                                   nb = 1, param = 2, pint_value = 0.99)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("content/Section4/Results")
save(TIME_logM_beta,
     file = paste0("TIME_logm_beta_d", min(dgrid), "_", max(dgrid), "_nobs", nobs, ".RData"))
