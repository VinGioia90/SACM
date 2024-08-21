###############################################################
# Code for reproducing the figures of the manuscript:         #
# "Scalable Generalized Additive Covariance Matrix Modelling" #
###############################################################
library(rstudioapi)
root_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(root_dir)

# RESTART R session before running this!!
library("mgcv", lib.loc="./my_library")

if(packageVersion("mgcv") != "9.0"){
  stop("Wrong version of mgcv!!")
}

# Load the needed packages
# (it might be required to install some packages manually)
source("loadPackages.R")
instload_packages()


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

nrun <- 10   # Set the number of runs
ncores <- 10 # Set the number of cores

scaleFUN <- function(x) sprintf("%.2f", x) # Number of decimal points set to 2

# logarithmic scale
#myscale_trans <- function(){
#  trans_new("myscale", function(x) log(x),
#            function(x) exp(x), domain = c(0, Inf))
#}

# square root scale
myscale_trans2 <- function(){
  trans_new("myscale", function(x) sqrt(x),
            function(x) x^2, domain = c(0, Inf))
}

#####################
#####################
#### SECTION 3.3 ####
#####################
#####################

#########################################################################
#########################################################################
## Evaluation of the second-order derivatives w.r.t. linear predictors ##
#########################################################################
#########################################################################
nobs <- 1000
dgrid <- seq(5, 50, by = 5)

setwd(root_dir)
setwd("content/Section3/Results")

###############################################
# Code for reproducing Figure 1 - SECTION 3.3 #
###############################################
load(file = paste0("TIME_mcd_D2eta_dgrid_min_", min(dgrid), "_max_", max(dgrid), "nobs", nobs, ".RData"))
load(file = paste0("TIME_logm_D2eta_dgrid_min_", min(dgrid), "_max_", max(dgrid), "nobs", nobs, ".RData"))

##################
setwd(root_dir)
setwd("content/Section3/Comp_logM_MCD_Hessian_eta")
source("Functions_Plots_Hessian_eta.R")               # Here there is the time_hessian function

#################################################
# MCD Computational times (overall and summary) #
#################################################
time_hessian_mcd <- time_hessian(obj = TIME_MCD_D2eta, param = "mcd",
                                 dgrid = dgrid, nrun = nrun,
                                 type = c("eff", "TMB"),
                                 type1 =  c("EFF", "AD"))


summary_time_hessian_mcd <- time_hessian_mcd$sum_res
summary_time_hessian_mcd <- summary_time_hessian_mcd[order(summary_time_hessian_mcd[,5],
                                                           summary_time_hessian_mcd[,1],
                                                           decreasing = FALSE),]
row.names(summary_time_hessian_mcd) <- NULL

all_time_hessian_mcd <- time_hessian_mcd$res
all_time_hessian_mcd <- all_time_hessian_mcd[order(all_time_hessian_mcd[,3],
                                                   all_time_hessian_mcd[,2],
                                                   decreasing = FALSE),]
row.names(all_time_hessian_mcd) <- NULL

##################################################
# logM Computational times (overall and summary) #
##################################################
time_hessian_logm <- time_hessian(obj = TIME_logM_D2eta, param = "logm",
                                        dgrid = dgrid, nrun = nrun,
                                        type = c("eff", "TMB"),
                                        type1 =  c("EFF", "AD"))

summary_time_hessian_logm <- time_hessian_logm$sum_res
summary_time_hessian_logm <- summary_time_hessian_logm[order(summary_time_hessian_logm[, 5],
                                                             summary_time_hessian_logm[, 1],
                                                             decreasing = FALSE),]
row.names(summary_time_hessian_logm) <- NULL

all_time_hessian_logm <- time_hessian_logm$res
all_time_hessian_logm <- all_time_hessian_logm[order(all_time_hessian_logm[, 3],
                                                     all_time_hessian_logm[, 2],
                                                     decreasing = FALSE),]
row.names(all_time_hessian_logm) <- NULL

################################################################
# Merge logM and MCD computational times (overall and summary) #
################################################################
summary_time_hessian <- rbind(summary_time_hessian_mcd, summary_time_hessian_logm)
summary_time_hessian$param<-factor(summary_time_hessian$param,
                                labels = c("logM", "MCD"), level = c("logm", "mcd"))

all_time_hessian <- rbind(all_time_hessian_mcd, all_time_hessian_logm)
all_time_hessian$param<-factor(all_time_hessian$param,
                               labels = c("logM", "MCD"), level = c("logm", "mcd"))


dg_sel <- dgrid #One can select a subset of the dgrid elements
all_time_hessian <- all_time_hessian[which(all_time_hessian$d %in% dg_sel),]
summary_time_hessian <- summary_time_hessian[which(summary_time_hessian$d %in% dg_sel),]

# Labels
#lab_time_logM <- with(summary_time_hessian,
#                      c(min(mean_time[param == "logM" & Type == "EFF"]),
#                        min(mean_time[param == "logM" & Type == "AD"]),
#                        max(mean_time[param == "logM" & Type == "EFF"]),
#                        max(mean_time[param == "logM" & Type == "AD"])))
#
#lab_time_MCD <- with(summary_time_hessian,
#                     c(min(mean_time[param == "MCD" & Type == "EFF"]),
#                       min(mean_time[param == "MCD" & Type == "AD"]),
#                       max(mean_time[param == "MCD" & Type == "EFF"]),
#                       max(mean_time[param == "MCD" & Type == "AD"])))



############################
# Note!!!: To set manually #
############################
lab_time_logM <- c(0, 10, 25, 50, 100, 225)
lab_time_MCD <- c(0, 0.04, 0.25, 1, 2, 3.5)


##########################################
# Plot using square root scale: Figure 1 #
##########################################
pl_Heta <- ggplot(all_time_hessian, aes(x = as.factor(d), y = time)) +
  geom_point(aes(colour = Type), size = 1, show.legend = TRUE, position = position_dodge(width = 0.3)) +
  geom_point(data = summary_time_hessian, aes(x = as.factor(d), y = mean_time, colour = Type),
             size = 2, position = position_dodge(width = 0.3), show.legend = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(data = summary_time_hessian, aes(y = mean_time, group = Type, col = Type),
            position = position_dodge(width = 0.3))+
  scale_color_manual(name = "Approach", values = c("AD" = "#00A9FF", "EFF" = "#F8766D")) +
  facet_grid2( . ~ param, scales = "free", switch = "y", independent = "all") +
  theme_bw() +
  scale_y_continuous(breaks = NULL,
                     sec.axis = sec_axis(~ . * 1, labels = scaleFUN, breaks = NULL)) +
  scale_x_discrete(breaks = dg_sel) +
  xlab("Dimension") + ylab("") +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        legend.position = "bottom", panel.spacing = unit(0.2, "lines"),
        axis.text = element_text(size = 12),
        text = element_text(size = 15),  legend.text=element_text(size=15),
        strip.text.x = element_text(size = 15)) +
  ggh4x::facetted_pos_scales(y = list(
    param == "logM" ~ scale_y_continuous(breaks = NULL, trans = myscale_trans2(),
                                         sec.axis = sec_axis(~ . * 1, labels = scaleFUN, breaks = lab_time_logM)),
    param == "MCD" ~ scale_y_continuous(breaks = NULL, trans = myscale_trans2(),
                                        sec.axis = sec_axis(~ . * 1, labels = scaleFUN, breaks = lab_time_MCD))
  ))

setwd(root_dir)
setwd("content/Section3/Plots")
ggsave("plot_TIME_hessian_eta.eps", pl_Heta, width = 30, height = 15, units = "cm")
ggsave("plot_TIME_hessian_eta.pdf", pl_Heta, width = 30, height = 15, units = "cm")


############################################################################
############################################################################
## Computational times for fitting under the MCD and logM parametrisation ##
############################################################################
############################################################################

###############################################
# Code for reproducing Figure 2 - SECTION 3.3 #
###############################################
dgrid <- c(2,5,10,15,20)
nrun <-  10
nobs <- 10000
sg <- FALSE

setwd(root_dir)
setwd("content/Section3/Comp_logM_MCD_Fit")
source("Functions_Plots_Overall_Fit.R") # Here the is the fit_time() function

######################################
# only MCD - generation in the paper #
######################################
setwd(root_dir)
setwd("content/Section3/Results")

load(paste0("sim_mcdG_mcdF_nrun_", nrun, "_n_", nobs, "_d_", paste0(dgrid, collapse="_"), ".RData"))
load(paste0("sim_mcdG_logmF_nrun_", nrun, "_n_", nobs, "_d_", paste0(dgrid, collapse = "_"), ".RData"))

# Get the computational times (overall and summary)
res_mcdG_mcdF <- fit_time(sim_mcdG_mcdF, param = "mcd", dgrid, nrun = nrun)      # Fit using the MCD parametrisation
res_mcdG_logmF <- fit_time(sim_mcdG_logmF, param = "logm", dgrid, nrun = nrun)   # Fit using the logM parametrisation
rm("sim_mcdG_mcdF", "sim_mcdG_logmF")

##########################################################################
# Merge the computational times and the iterations (overall and summary) #
##########################################################################
mean_TIME <- rbind(res_mcdG_mcdF$sum_res, res_mcdG_logmF$sum_res)
all_TIME <- rbind(res_mcdG_mcdF$res, res_mcdG_logmF$res)
mean_ITER <- rbind(res_mcdG_mcdF$sum_iter, res_mcdG_logmF$sum_iter)
all_ITER <- rbind(res_mcdG_mcdF$iter, res_mcdG_logmF$iter)


# Padding with useful information for plotting the times and iterations
colnames(mean_TIME)[2] <- colnames(all_TIME)[1] <- colnames(mean_ITER)[2] <- colnames(all_ITER)[1] <- "Value"
mean_TIME  <- cbind(mean_TIME , c(rep("Time", 2 * length(dgrid))))
colnames(mean_TIME)[6] <- "Type2"

mean_ITER  <- cbind(mean_ITER , c(rep("Iterations", 2 * length(dgrid))))
colnames(mean_ITER)[6] <- "Type2"

all_TIME  <- cbind(all_TIME , c(rep("Time", 2 * length(dgrid) * nrun)))
colnames(all_TIME)[4] <- "Type2"

all_ITER  <- cbind(all_ITER , c(rep("Iterations", 2 * length(dgrid) * nrun)))
colnames(all_ITER)[4] <- "Type2"

##########################################################
# Join the time and the iterations (overall and summary) #
##########################################################
data_time_iter <- rbind(all_TIME, all_ITER)
data_time_iter_sum <-  rbind(mean_TIME[, c(1, 2, 5, 6)], mean_ITER[, c(1, 2, 5, 6)])


data_hline <- data.frame(Type2 = c("Time","Iterations"),  # Create data for horizontal lines
                         hline = c(0, NA))
data_hline$Type2 <- factor(data_hline$Type2, levels  = c("Time", "Iterations"))


dgrid_sel <- dgrid[-1] # To select a subset of dgrid (I removed the d=2 case)
data_time_iter <- data_time_iter[which(data_time_iter$d %in% dgrid_sel),]
data_time_iter_sum <- data_time_iter_sum[which(data_time_iter_sum$d %in% dgrid_sel),]

data_time_iter$Type2 <- factor(data_time_iter$Type2, levels = c("Time", "Iterations"), labels = c("Time", "Iterations"))
data_time_iter_sum$Type2 <- factor(data_time_iter_sum$Type2, levels = c("Time", "Iterations"), labels = c("Time", "Iterations"))

# Labelling
#label_time <- with(data_time_iter_sum,
#                   c(min(Value[Type2 == "Time" & Type == "MCD"]), min(Value[Type2 == "Time" & Type == "logM"]),
#                     max(Value[Type2 == "Time" & Type == "MCD"]), max(Value[Type2 == "Time" & Type == "logM"])))

# To set manually
label_time <- c(1, 5, 25, 50, 100, 450)

##########################################
# Plot using square root scale: Figure 2 #
##########################################
pl_Fit_Time_Iter <- ggplot(data_time_iter,
                           aes(x = factor(d, labels = as.character(dgrid_sel), levels = as.character(dgrid_sel)), y = Value)) +
  geom_point(aes(colour = Type), size = 1, show.legend = TRUE, position = position_dodge(width = 0.3)) +
  geom_point(data = data_time_iter_sum,
             aes(x = factor(d, labels = as.character(dgrid_sel),levels = as.character(dgrid_sel)),
                 y = Value, colour = Type), size = 3,
             position = position_dodge(width = 0.3), show.legend = FALSE) +
  geom_line(data = data_time_iter_sum, aes(y = Value, group = Type, col = Type),
            position = position_dodge(width = 0.3)) +
  scale_color_manual(name = "Parametrisation", values = c("MCD" = "#F8766D", "logM" = "#619CFF")) +
  facet_grid2( . ~ Type2 , scales = "free", switch = "y", independent = "all") +
  theme_bw() +
  scale_y_continuous(breaks = NULL, sec.axis = sec_axis(~ . * 1,labels = scaleFUN, breaks = NULL)) +
  scale_x_discrete(breaks = dgrid_sel) +
  xlab("Dimension") + ylab("") +
  theme(panel.grid.minor = element_blank(), axis.text = element_text(size = 12), text = element_text(size = 15),
        legend.text=element_text(size=15), strip.text.x = element_text(size = 15),
        panel.grid.major = element_blank(), legend.position = "bottom", panel.spacing = unit(0.2, "lines"))+
  ggh4x::facetted_pos_scales(y = list(
    Type2 == "Time" ~ scale_y_continuous(breaks = NULL, trans = myscale_trans2(),
                                         sec.axis = sec_axis(~ . * 1, labels = scaleFUN,
                                                             breaks = label_time)),
    Type2 == "Iterations" ~ scale_y_continuous(breaks = NULL, trans = myscale_trans2(),
                                               sec.axis = sec_axis(~ . * 1,
                                                                   breaks=seq(50, 100, by = 10)))
  ))


setwd(root_dir)
setwd("content/Section3/Plots")
ggsave("plot_TIME_ITER_sqrtscale_genMCD.eps", pl_Fit_Time_Iter, width = 30, height = 15, units = "cm")
ggsave("plot_TIME_ITER_sqrtscale_genMCD.pdf", pl_Fit_Time_Iter, width = 30, height = 15, units = "cm")

##############################################################################################
# SUPPLEMENTARY MATERIAL: decide whether  including the results under logM generation        #
# !!! I think it's redundant reporting computational times and iterations                    #
# Now it corresponds to Figure B.1 of the Supplementary Material                             #
##############################################################################################
setwd(root_dir)
setwd("content/Section3/Results")

load(paste0("sim_logmG_mcdF_nrun_", nrun, "_n_", nobs, "_d_", paste0(dgrid, collapse="_"), ".RData"))
load(paste0("sim_logmG_logmF_nrun_", nrun, "_n_", nobs, "_d_", paste0(dgrid, collapse = "_"), ".RData"))


res_logmG_mcdF <- fit_time(sim_logmG_mcdF, param = "mcd", dgrid, nrun = nrun)     # MCD  fit
res_logmG_logmF <- fit_time(sim_logmG_logmF, param = "logm", dgrid, nrun = nrun)  # logM fit
rm("sim_logmG_mcdF", "sim_logmG_logmF")

##########################################################################
# Merge the computational times and the iterations (overall and summary) #
##########################################################################
mean_TIME <- rbind(res_logmG_mcdF$sum_res, res_logmG_logmF$sum_res)
all_TIME <- rbind(res_logmG_mcdF$res, res_logmG_logmF$res)
mean_ITER <- rbind(res_logmG_mcdF$sum_iter, res_logmG_logmF$sum_iter)
all_ITER <- rbind(res_logmG_mcdF$iter, res_logmG_logmF$iter)

# Padding with useful information for plotting the times and iterations
colnames(mean_TIME)[2] <- colnames(all_TIME)[1] <- colnames(mean_ITER)[2] <- colnames(all_ITER)[1] <- "Value"

mean_TIME  <- cbind(mean_TIME , c(rep("Time", 2 * length(dgrid))))
colnames(mean_TIME)[6] <- "Type2"

mean_ITER  <- cbind(mean_ITER , c(rep("Iterations", 2 * length(dgrid))))
colnames(mean_ITER)[6] <- "Type2"

all_TIME  <- cbind(all_TIME , c(rep("Time", 2 * length(dgrid) * nrun)))
colnames(all_TIME)[4] <- "Type2"

all_ITER  <- cbind(all_ITER , c(rep("Iterations", 2 * length(dgrid) * nrun)))
colnames(all_ITER)[4] <- "Type2"


##########################################################
# Join the time and the iterations (overall and summary) #
##########################################################
data_time_iter <- rbind(all_TIME, all_ITER)
data_time_iter_sum <-  rbind(mean_TIME[,c(1,2,5,6)], mean_ITER[,c(1,2,5,6)])


data_hline <- data.frame(Type2 = c("Time","Iterations"),  # Create data for horizontal lines
                         hline = c(0, NA))
data_hline$Type2 <- factor(data_hline$Type2, levels=c("Time", "Iterations"))


dgrid_sel <- dgrid[-1] # To select a subset of dgrid
data_time_iter <- data_time_iter[which(data_time_iter$d %in% dgrid_sel),]
data_time_iter_sum <- data_time_iter_sum[which(data_time_iter_sum$d %in% dgrid_sel),]

data_time_iter$Type2 <- factor(data_time_iter$Type2, levels=c("Time","Iterations"), labels=c("Time","Iterations"))
data_time_iter_sum$Type2 <- factor(data_time_iter_sum$Type2, levels=c("Time","Iterations"), labels=c("Time","Iterations"))

#label_time <- with(data_time_iter_sum,
#                   c(min(Value[Type2 == "Time" & Type == "MCD"]), min(Value[Type2 == "Time" & Type == "logM"]),
#                     max(Value[Type2 == "Time" & Type == "MCD"]), max(Value[Type2 == "Time" & Type == "logM"])))

label_time <- c(1, 5, 25, 50, 100, 450)

pl_Fit_Time_Iter <- ggplot(data_time_iter,
                           aes(x = factor(d, labels = as.character(dgrid_sel), levels = as.character(dgrid_sel)), y = Value)) +
  geom_point(aes(colour = Type), size = 1, show.legend = TRUE, position = position_dodge(width = 0.3)) +
  geom_point(data = data_time_iter_sum,
             aes(x = factor(d, labels = as.character(dgrid_sel), levels = as.character(dgrid_sel)),
                 y = Value, colour = Type), size = 3,
             position = position_dodge(width = 0.3), show.legend = FALSE) +
  geom_line(data = data_time_iter_sum, aes(y = Value, group = Type, col = Type),
            position = position_dodge(width = 0.3)) +
  scale_color_manual(name = "Parametrisation", values = c("MCD" = "#F8766D", "logM" = "#619CFF")) +
  facet_grid2( . ~ Type2 , scales = "free", switch = "y", independent = "all") +
  theme_bw() +
  scale_y_continuous(breaks = NULL, sec.axis = sec_axis(~ . * 1,labels = scaleFUN, breaks = NULL)) +
  scale_x_discrete(breaks = dgrid_sel) +
  xlab("Dimension") + ylab("") +
  theme(panel.grid.minor = element_blank(), axis.text = element_text(size = 12),  text = element_text(size = 15),
        legend.text=element_text(size=15), strip.text.x = element_text(size = 15),
        panel.grid.major = element_blank(), legend.position = "bottom", panel.spacing = unit(0.2, "lines"))+
  ggh4x::facetted_pos_scales(y = list(
    Type2 == "Time" ~ scale_y_continuous(breaks = NULL, trans = myscale_trans2(),
                                         sec.axis = sec_axis(~ . * 1, labels = scaleFUN,
                                                             breaks = label_time)),
    Type2 == "Iterations" ~ scale_y_continuous(breaks = NULL, trans = myscale_trans2(),
                                               sec.axis = sec_axis(~ . * 1,
                                                                   breaks=seq(50, 100, by = 10)))
  ))


setwd(root_dir)
setwd("content/SupplementaryMaterial/Plots")
ggsave("plot_TIME_ITER_sqrtscale_genLOGM.eps", pl_Fit_Time_Iter, width = 30, height = 15, units = "cm")
ggsave("plot_TIME_ITER_sqrtscale_genLOGM.pdf", pl_Fit_Time_Iter, width = 30, height = 15, units = "cm")

################################################################################################
################################################################################################
## SUPPLEMENTARY MATERIALS: Figure B.2 - comparison of performance MCD vs logM using logScore ##
################################################################################################
################################################################################################
setwd(root_dir)
setwd("content/Section3/Comp_logM_MCD_Fit")
source("Functions_Plots_Overall_Fit.R") #here there is the log_score() function

setwd(root_dir)
setwd("content/Section3/Results")

####################
# MCD - Generation #
####################
load(paste0("sim_mcdG_mcdF_nrun_", nrun, "_n_", nobs, "_d_", paste0(dgrid, collapse="_"), ".RData"))      # MCD Fit
load(paste0("sim_mcdG_logmF_nrun_", nrun, "_n_", nobs, "_d_", paste0(dgrid, collapse = "_"), ".RData"))   # logM fit

# Extract logScore
logS_mcdG_mcdF <- log_Score(sim_mcdG_mcdF, nrun, dgrid, nobs, param = "mcd")
colnames(logS_mcdG_mcdF) <- as.character(dgrid)
logS_mcdG_logmF <- log_Score(sim_mcdG_logmF, nrun, dgrid, nobs, param = "logm")
colnames(logS_mcdG_logmF) <- as.character(dgrid)
rm("sim_mcdG_mcdF", "sim_mcdG_logmF")

# Select a subset of dgrid
dgrid_sel <- as.character(dgrid) # Here, it is possible to remove d=2

# Merge the logScores
logS_mcd <- data.frame(logS_gen = as.vector(logS_mcdG_mcdF[,dgrid_sel]),
                       logS_nogen = as.vector(logS_mcdG_logmF[,dgrid_sel]))
# Padding with other info
logS_mcd <- cbind(logS_mcd, d = c(rep(dgrid_sel, each = nrun)))
logS_mcd <- cbind(logS_mcd, Type2 = rep("MCD - Generation"))

#####################
# logm - Generation #
#####################
load(paste0("sim_logmG_mcdF_nrun_", nrun, "_n_", nobs, "_d_", paste0(dgrid, collapse="_"), ".RData"))
load(paste0("sim_logmG_logmF_nrun_", nrun, "_n_", nobs, "_d_", paste0(dgrid, collapse = "_"), ".RData"))

logS_logmG_mcdF <- log_Score(sim_logmG_mcdF, nrun, dgrid, nobs, param = "mcd")
colnames(logS_logmG_mcdF) <- as.character(dgrid)
logS_logmG_logmF <- log_Score(sim_logmG_logmF, nrun, dgrid, nobs, param = "logm")
colnames(logS_logmG_logmF) <- as.character(dgrid)

rm("sim_logmG_mcdF", "sim_logmG_logmF")

# Select a subset of dgrid
dgrid_sel <- as.character(dgrid) # Here, it is possible to remove d=2

# Merge the logScores
logS_logm <- data.frame(logS_gen = as.vector(logS_logmG_logmF[,dgrid_sel]),
                         logS_nogen = as.vector(logS_logmG_mcdF[,dgrid_sel]))

# Padding with other info
logS_logm <- cbind(logS_logm, d = c(rep(dgrid_sel, each = nrun)))
logS_logm <- cbind(logS_logm, Type2 = rep("logM - Generation"))

# Merging the logScores under the MCD and logM generation
logS <- rbind(logS_mcd,logS_logm)

# Labelling
logS$d <- factor(logS$d, levels = dgrid_sel,
                 labels = dgrid_sel)

logS$Type2 <- factor(logS$Type2, levels = c("MCD - Generation", "logM - Generation"),
                     labels = c("MCD - Generation", "logM - Generation"))

# Setting the colours (varying with d)
col_values<- hue_pal()(length(dgrid))

# Breaks label of x and y axes
breaks_seq_MCDgen_x <- as.numeric(tapply(logS[logS$Type2 == "MCD - Generation",]$logS_gen, logS$d[logS$Type2 == "MCD - Generation"], mean))
breaks_seq_MCDgen_y <- as.numeric(tapply(logS[logS$Type2 == "MCD - Generation",]$logS_nogen, logS$d[logS$Type2 == "MCD - Generation"], mean))

# Plot under MCD generation (Note that logScores are represented in the original scale)
pl_MCD_gen <- ggplot(logS[logS$Type2 == "MCD - Generation",], aes(x = logS_gen, y = logS_nogen, color = d)) +
  geom_point(size = 2, show.legend = TRUE)+
  geom_abline(slope = 1, intercept = 0, col = "black") +
  scale_color_manual(name = "Dimension", values = col_values)+
  theme_bw() +
  facet_grid(. ~ "MCD generation") +
  xlab("LS - MCD fit") + ylab("LS - logM fit") +
  scale_y_continuous(breaks = breaks_seq_MCDgen_y,
                     sec.axis = sec_axis(~ . * 1,labels = scaleFUN, breaks = NULL),) +
  scale_x_continuous(breaks = breaks_seq_MCDgen_x,
                     sec.axis = sec_axis(~ . * 1,labels = scaleFUN, breaks = NULL),)+
  theme(panel.grid.minor = element_blank(), axis.text = element_text(size = 12),  text = element_text(size = 15),
        legend.text=element_text(size=15), strip.text.x = element_text(size = 15),
        panel.grid.major = element_blank(), legend.position = "bottom", panel.spacing = unit(0.1, "lines"))

# Breaks label of x and y axes
breaks_seq_logMgen_x <- as.numeric(tapply(logS[logS$Type2 == "logM - Generation",]$logS_gen, logS$d[logS$Type2 == "logM - Generation"], mean))
breaks_seq_logMgen_y <- as.numeric(tapply(logS[logS$Type2 == "logM - Generation",]$logS_nogen, logS$d[logS$Type2 == "logM - Generation"], mean))

# Plot under logM generation (Note that logScores are represented in the original scale)
pl_logM_gen <- ggplot(logS[logS$Type2 == "logM - Generation",], aes(x = logS_gen, y = logS_nogen, color = d)) +
  geom_point(size = 2, show.legend = TRUE)+
  geom_abline(slope = 1, intercept = 0, col = "black") +
  facet_grid(. ~ "logM generation") +
  xlab("LS - logM fit") + ylab("LS - MCD fit") +
  scale_color_manual(name = "Dimension", values = col_values)+
  theme_bw() +
  scale_y_continuous(breaks = breaks_seq_logMgen_y,
                     sec.axis = sec_axis(~ . * 1,labels = scaleFUN, breaks = NULL),) +
  scale_x_continuous(breaks=breaks_seq_logMgen_x,
                     sec.axis = sec_axis(~ . * 1,labels = scaleFUN, breaks = NULL),)+
  theme(panel.grid.minor = element_blank(), axis.text = element_text(size = 12),  text = element_text(size = 15),
        legend.text=element_text(size=15), strip.text.x = element_text(size = 15),
        panel.grid.major = element_blank(), legend.position = "bottom", panel.spacing = unit(0.1, "lines"))

# Merging the plots above
pl_logS1 <- ggarrange(pl_MCD_gen,
                      pl_logM_gen,
                      nrow = 1,
                      common.legend = TRUE,
                      legend = "bottom")

setwd(root_dir)
setwd("content/SupplementaryMaterial/Plots")
ggsave("plot_logScore1.eps", pl_logS1, width = 30, height = 15, units = "cm")
ggsave("plot_logScore1.pdf", pl_logS1, width = 30, height = 15, units = "cm")

#############################################
# Code for reproducing Figure 3 - SECTION 4 #
#############################################
setwd(root_dir)
# The function for extracting and plotiing computational times is in Section 3
setwd("content/Section3/Comp_logM_MCD_Hessian_eta")
source("Functions_Plots_Hessian_eta.R")

nobs <- 1000
dgrid <- seq(5, 100, by = 5)
step_length <- 5
pint_type <- c("dm05","dm1", "dm2","const")

##################
setwd(root_dir)
setwd("content/Section4/Results")

load(file=paste0("TIME_mcd_beta_d", min(dgrid),"_",max(dgrid),"_nobs", nobs, ".RData"))
load(file=paste0("TIME_logm_beta_d", min(dgrid),"_",max(dgrid),"_nobs", nobs, ".RData"))


# MCD: not included in the paper
# time_Hbeta_mcd <-  lapply(1 : length(pint_type),
#                             function(x) time_hessian(obj = TIME_MCD_beta[[x]], param = "mcd",
#                                                      dgrid = dgrid, nrun = nrun,
#                                                      type = c("block", "noblock"),
#                                                      type1 = c("PARS", "STD"),
#                                                      beta = TRUE))
#
#
# summary_time_hessianB_mcd <- rbind(time_Hbeta_mcd[[1]]$sum_res, time_Hbeta_mcd[[3]]$sum_res)
# row.names(summary_time_hessianB_mcd) <- NULL
#
# all_time_hessianB_mcd <- rbind(time_Hbeta_mcd[[1]]$res, time_Hbeta_mcd[[3]]$res)
# row.names(all_time_hessianB_mcd) <- NULL

# # Get the ratio of the mean and all the times
# rel_mean_time_hessianB_mcd <- data.frame(d = rep(dgrid, 2),
#                                           rel_time = summary_time_hessianB_mcd[summary_time_hessianB_mcd$Type == "STD",2]/
#                                                      summary_time_hessianB_mcd[summary_time_hessianB_mcd$Type == "PARS",2],
#                                           param = rep("mcd", 2 * length(dgrid)),
#                                           Scenario = factor(rep(c("S1", "S2"), each = length(dgrid)),
#                                                             labels = c("S1", "S2"), level = c("S1", "S2")),
#                                           Ratio = rep("STD/PARS", 2 * length(dgrid)))
#  rel_all_time_hessianB_mcd <- data.frame(d = rep(dgrid, 2),
#                                          rel_time = all_time_hessianB_mcd[all_time_hessianB_mcd$Type == "STD",1]/
#                                                     all_time_hessianB_mcd[all_time_hessianB_mcd$Type == "PARS",1],
#                                          param = rep("mcd", 2 * length(dgrid)),
#                                          Scenario = factor(rep(c("S1", "S2"), each = 2 * nrun * length(dgrid)),
#                                                          labels = c("S1", "S2"), level = c("S1", "S2")),
#                                          Ratio=rep("STD/PARS",2 * length(dgrid)))
#
#  lab_time_MCD <- c(1,max(rel_mean_time_hessianB_mcd[rel_mean_time_hessianB_mcd$Scenario == "S1",2]),
#                      max(rel_mean_time_hessianB_mcd[rel_mean_time_hessianB_mcd$Scenario == "S2",2]))
#
#
#  dg_sel <- dgrid #To select a subset of the dgrid elements
#  rel_mean_time_hessianB_MCD <- rel_mean_time_hessianB_mcd[which(rel_mean_time_hessianB_mcd$d %in% dg_sel),]
#  rel_all_time_hessianB_MCD <- rel_all_time_hessianB_mcd[which(rel_all_time_hessianB_mcd$d %in% dg_sel),]
#
#  pl_Hbeta <- ggplot(rel_all_time_hessianB_MCD, aes(x = as.factor(d), y = rel_time)) +
#    geom_point(aes(colour = Scenario), size = 1, show.legend = TRUE, position = position_dodge(width = 0.3)) +
#    geom_point(data = rel_mean_time_hessianB_mcd, aes(x = as.factor(d), y = rel_time, colour = Scenario),
#               size = 2, position = position_dodge(width = 0.3), show.legend = FALSE) +
#    geom_line(data = rel_mean_time_hessianB_mcd, aes(y = rel_time, group = Scenario, col = Scenario),
#              position = position_dodge(width = 0.3)) +
#    geom_hline(yintercept = 1, linetype = "dashed") +
#    scale_color_manual(name = "", values = c("S1" = "#00A9FF", "S2" = "#F8766D")) +
#    theme_bw() +
#    scale_y_continuous(breaks = NULL,trans = myscale_trans(),
#                       sec.axis = sec_axis(~ . * 1, labels = scaleFUN, breaks = lab_time_MCD)) +
#    scale_x_discrete(breaks = dg_sel) +
#    xlab("Dimension") + ylab("") +
#    theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
#          legend.position = "bottom", panel.spacing = unit(0.2, "lines"))
#
 # setwd(root_dir)
 # setwd("content/Section4/Plots")
 # ggsave("plot_rel_TIME_hessian_beta_MCD.eps", pl_Hbeta, width = 15, height = 15, units = "cm")
 # ggsave("plot_rel_TIME_hessian_beta_MCD.pdf", pl_Hbeta, width = 15, height = 15, units = "cm")


# logM
time_Hbeta_logm <-  lapply(1 : length(pint_type),
                           function(x) time_hessian(obj = TIME_logM_beta[[x]], param = "logm",
                                                    dgrid = dgrid, nrun = nrun,
                                                    type = c("block", "noblock"),
                                                    type1 = c("PARS", "STD"),
                                                    beta = TRUE))

# Select Scenario 1 (dm05) and Scenario 2 (dm2)
#summary_time_hessianB_logm <- rbind(time_Hbeta_logm[[1]]$sum_res, time_Hbeta_logm[[3]]$sum_res)
summary_time_hessianB_logm <- rbind(time_Hbeta_logm[[1]]$sum_res, time_Hbeta_logm[[2]]$sum_res)
row.names(summary_time_hessianB_logm) <- NULL

all_time_hessianB_logm <- rbind(time_Hbeta_logm[[1]]$res, time_Hbeta_logm[[2]]$res)
row.names(all_time_hessianB_logm) <- NULL

rel_mean_time_hessianB_logm <- data.frame(d = rep(dgrid, 2),
                                          rel_time = summary_time_hessianB_logm[summary_time_hessianB_logm$Type == "STD",2]/
                                                     summary_time_hessianB_logm[summary_time_hessianB_logm$Type == "PARS",2],
                                          param = rep("logm", 2 * length(dgrid)),
                                          Scenario = factor(rep(c("S1", "S2"), each = length(dgrid)),
                                                            labels = c("S1", "S2"), level = c("S1", "S2")),
                                          Ratio = rep("STD/PARS", 2 * length(dgrid)))


rel_all_time_hessianB_logm <- data.frame(d = rep(dgrid, 2),
                                         rel_time = all_time_hessianB_logm[all_time_hessianB_logm$Type == "STD",1]/
                                                    all_time_hessianB_logm[all_time_hessianB_logm$Type == "PARS",1],
                                         param = rep("logm",2 * length(dgrid)),
                                         Scenario = factor(rep(c("S1", "S2"), each = nrun * length(dgrid)),
                                                           labels=c("S1", "S2"), level = c("S1", "S2")),
                                         Ratio = rep("STD/PARS", 2 * length(dgrid)))

#lab_time_logM <- c(1, max(rel_mean_time_hessianB_logm[rel_mean_time_hessianB_logm$Scenario == "S1",2]),
#                   max(rel_mean_time_hessianB_logm[rel_mean_time_hessianB_logm$Scenario == "S2",2]))
lab_time_logM <- c(1, 2, 4, 8, 16, 24)


dg_sel <- dgrid #To select a subset of the dgrid elements
rel_mean_time_hessianB_logm <- rel_mean_time_hessianB_logm[which(rel_mean_time_hessianB_logm$d %in% dg_sel),]
rel_all_time_hessianB_logm <- rel_all_time_hessianB_logm[which(rel_all_time_hessianB_logm$d %in% dg_sel),]




pl_Hbeta <- ggplot(rel_all_time_hessianB_logm, aes(x = as.factor(d), y = rel_time)) +
  geom_point(aes(colour = Scenario), size = 1, show.legend = TRUE, position = position_dodge(width = 0.3)) +
  geom_point(data = rel_mean_time_hessianB_logm, aes(x = as.factor(d), y = rel_time, colour = Scenario),
             size = 2, position = position_dodge(width = 0.3), show.legend = FALSE) +
  geom_line(data = rel_mean_time_hessianB_logm, aes(y = rel_time, group = Scenario, col = Scenario),
            position = position_dodge(width = 0.3)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_color_manual(name = "Scenario", values = c("S1" = "#00A9FF", "S2" = "#F8766D")) +
  theme_bw() +
  facet_grid(. ~ "Relative times") +
  scale_y_continuous(breaks = NULL,trans = myscale_trans2(),
                     sec.axis = sec_axis(~ . * 1, labels = scaleFUN, breaks = lab_time_logM)) +
  scale_x_discrete(breaks = dg_sel) +
  xlab("Dimension") + ylab("") +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        axis.text = element_text(size = 12),  text = element_text(size = 15),
        legend.text=element_text(size=15), strip.text.x = element_text(size = 15),
        legend.position = "bottom", panel.spacing = unit(0.2, "lines"))


nspars <- data.frame(d =  rep(seq(min(dgrid), max(dgrid), by = 0.01), 2),
                     nelno0_rel = c(1/sqrt(seq(min(dgrid), max(dgrid), by = 0.01)), 1/(seq(min(dgrid), max(dgrid), by = 0.01))),
                     nelno0_abs = c(seq(min(dgrid), max(dgrid), by = 0.01) * (seq(min(dgrid), max(dgrid), by = 0.01)+1)/2 * 1/sqrt(seq(min(dgrid), max(dgrid), by = 0.01)),
                                    seq(min(dgrid), max(dgrid), by = 0.01) * (seq(min(dgrid), max(dgrid), by = 0.01)+1)/2 * 1/seq(min(dgrid), max(dgrid), by = 0.01)),
                     Scenario = factor(rep(c("S1", "S2"), each = (length(seq(min(dgrid), max(dgrid), by = 0.01)))),
                                       labels=c("S1", "S2"), level = c("S1", "S2")))
add_point <-  data.frame(d =  rep(seq(min(dgrid), max(dgrid), by = step_length), 2),
                         nelno0_rel = c(1/sqrt(seq(min(dgrid), max(dgrid), by = step_length)), 1/(seq(min(dgrid), max(dgrid), by = step_length))),
                         nelno0_abs = c(seq(min(dgrid), max(dgrid), by = step_length) * (seq(min(dgrid), max(dgrid), by = step_length)+1)/2 * 1/sqrt(seq(min(dgrid), max(dgrid), by = step_length)),
                                        seq(min(dgrid), max(dgrid), by = step_length) * (seq(min(dgrid), max(dgrid), by = step_length)+1)/2 * 1/seq(min(dgrid), max(dgrid), by = step_length)),
                         Scenario = factor(rep(c("S1", "S2"), each = length(dgrid)),
                                           labels=c("S1", "S2"), level = c("S1", "S2")))


#plot_perc_covmod_lpi <- ggplot(data.frame(nspars), aes(x = d, y = nelno0_rel)) +
#  xlab("Dimension") + ylab("") +
#  theme_bw() +
#  ylim(0, 0.5) +
#  geom_hline(yintercept = 0, linetype = "dashed") +
#  facet_grid(. ~ "logM elements modelled (%)") +  geom_point(data = add_point , aes(x = d, y = nelno0_rel, colour = Scenario), show.legend = TRUE) +
#  geom_line(aes(y = nelno0 , group = Scenario, col = Scenario), size = 0.1) +
#  scale_color_manual(name = "Scenario", values = c("S1" = "#00A9FF", "S2" = "#F8766D")) +
#  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
#        axis.text = element_text(size = 15),  text = element_text(size = 15),
#        legend.text=element_text(size=15), strip.text.x = element_text(size = 15),
#        legend.position = "bottom", panel.spacing = unit(0.2, "lines"))+
#  scale_x_continuous(breaks=add_point$d)


# 505 is the number of elements modelled for d=100 under S1
# 0.44721360 is the % number of elements modelled for d=5 under S1
# So I scaled by the ratio between the maximum of value of absolute and relative frequencies.
scale_function <- function(x){
  return (x/(505/ 0.44721360))
}
inv_scale_function <- function(x){
  return (x*(505/ 0.44721360))
}


plot_perc_covmod_lpi <- ggplot(data.frame(nspars), aes(x = d, y = nelno0_rel)) +
  xlab("Dimension") + ylab("") +
  theme_bw() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_grid(. ~ "logM elements modelled") +
  geom_point(data = add_point , aes(x = d, y = nelno0_rel, colour = Scenario), show.legend = TRUE) +
  geom_point(data = add_point , aes(x = d, y = scale_function(nelno0_abs), colour = Scenario), show.legend = TRUE) +
  geom_line(aes(y = nelno0_rel, group = Scenario, col = Scenario), size = 0.1) +
  geom_line(aes(y = scale_function(nelno0_abs), group = Scenario, col = Scenario), size = 0.1, lty = "dashed") +
  scale_color_manual(name = "Scenario", values = c("S1" = "#00A9FF", "S2" = "#F8766D")) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        axis.text = element_text(size = 12),  text = element_text(size = 15),
        legend.text=element_text(size=15), strip.text.x = element_text(size = 15),
        legend.position = "bottom", panel.spacing = unit(0.2, "lines"))+
  scale_x_continuous(breaks=add_point$d) +
  scale_y_continuous( breaks = seq(0,0.5,by=0.05), lim = c(0,0.5) ,
                     sec.axis = sec_axis(~ inv_scale_function(.), breaks = seq(0,550,50),  name = ""))


# Merging the plots above #!!! Find a way to merge the x-axis "Dimension" title !!!#
pl_comp_Hbeta <- ggarrange(plot_perc_covmod_lpi,
                           pl_Hbeta ,
                           nrow = 1,
                           common.legend = TRUE,
                           legend = "bottom")
#pl_comp_Hbeta

setwd(root_dir)
setwd("content/Section4/Plots")
ggsave("plot_rel_TIME_hessian_beta_logM_and_Sparsity.eps", pl_comp_Hbeta, width = 30, height = 15, units = "cm")
ggsave("plot_rel_TIME_hessian_beta_logM_and_Sparsity.pdf", pl_comp_Hbeta, width = 30, height = 15, units = "cm")

##############################################
# Code for reproducing Figure 4 - SECTION 6  #
##############################################
setwd(root_dir)
setwd("content/Section6")
source("Functions_Plots_Overall_Fit.R")

dgrid <- c(2, 5, 10)
nrun <- 10
nobs <- 10000


setwd(root_dir)
setwd("content/Section6/Results")
load(paste0("sim_mcd_fit_nrun_", nrun, "_n_", nobs, "_d_", paste0(dgrid, collapse = "_"), ".RData"))

LAML_extr <- LAML_extraction(sim_mcd_fit, nrun, dgrid, nobs,  param = "mcd")
#LAML_efs <- as.numeric(LAML_extr[[1]])
#LAML_bfgs <- as.numeric(LAML_extr[[2]])
#LAML_bfgsinit <- as.numeric(LAML_extr[[3]])
LAML_efs <- as.numeric(LAML_extr[[1]])
LAML_exact_efs <- as.numeric(LAML_extr[[2]])
LAML_exact_efs_initialised <- as.numeric(LAML_extr[[3]])
LAML_bfgs <- as.numeric(LAML_extr[[4]])
LAML_bfgsinit <- as.numeric(LAML_extr[[5]])


#LAML <- data.frame(LAML_efs = rep(LAML_efs,2),
#                   LAML_bfgs = c(LAML_bfgs, LAML_bfgsinit),
#                   diff_LAML_bfgs_efs = c(LAML_bfgs - LAML_efs,  LAML_bfgsinit - LAML_efs),
#                   rel_diff_LAML_bfgs_efs = c((LAML_bfgs - LAML_efs)/abs(LAML_bfgs), (LAML_bfgsinit - LAML_efs)/abs(LAML_bfgsinit)),
#                   d =  rep(factor(rep(dgrid, each = nrun ), levels = dgrid),2),
#                   init = c(rep("BFGS", length(LAML_efs)),rep("BFGS - init",length(LAML_efs))))
LAML <- data.frame(LAML_efs = rep(LAML_efs, 4),
                   LAML_other = c(LAML_exact_efs, LAML_exact_efs_initialised, LAML_bfgs, LAML_bfgsinit),
                   diff_LAML_other_efs = c(LAML_exact_efs - LAML_efs, LAML_exact_efs_initialised - LAML_efs,
                                           LAML_bfgs - LAML_efs,  LAML_bfgsinit - LAML_efs),
                   rel_diff_LAML_other_efs = c((LAML_exact_efs - LAML_efs)/abs(LAML_exact_efs),
                                               (LAML_exact_efs_initialised - LAML_efs)/abs(LAML_exact_efs),
                                               (LAML_bfgs - LAML_efs)/abs(LAML_bfgs),
                                               (LAML_bfgsinit - LAML_efs)/abs(LAML_bfgsinit)),
                   d =  rep(factor(rep(dgrid, each = nrun ), levels = dgrid),4),
                   init = c(rep("ExactFS", length(LAML_efs)), rep("ExactFS - init", length(LAML_efs)),
                            rep("BFGS", length(LAML_efs)),rep("BFGS - init",length(LAML_efs))))

TIMES <- fit_time(sim_mcd_fit, param = "mcd", dgrid, nrun)
dgrid_sel <- dgrid # To select a subset of dgrid
data_time <- TIMES$res[which(TIMES$res$d %in% dgrid_sel),]
data_time_sum <- TIMES$sum_res[which(TIMES$sum_res$d %in% dgrid_sel),]

#TIMES_sum <- data.frame(Value = c(data_time_sum$efs, data_time_sum$bamlss, data_time_sum$bfgs, data_time_sum$bfgsinit + data_time_sum$efs),
#                        d = rep(dgrid, times = 4),
#                        Type2 = rep(c("FS", "BAMLSS",  "BFGS","BFGS - init"), each=length(dgrid)))
#
#TIME_efs_bamlss_bfgs <- data.frame(Value = c(data_time$time_efs,data_time$time_bamlss, data_time$time_bfgs, data_time$time_bfgsinit + data_time$time_efs),
#                              d = rep(rep(dgrid, each = nrun),4),
#                              Type2 = rep(c("FS", "BAMLSS", "BFGS", "BFGS - init"), each=nrun*length(dgrid)))

TIMES_sum <- data.frame(Value = c(data_time_sum$efs, data_time_sum$efsExact, data_time_sum$efsExact + data_time_sum$efs,
                                  data_time_sum$bamlss, data_time_sum$bfgs, data_time_sum$bfgsinit + data_time_sum$efs),
                        d = rep(dgrid, times = 6),
                        Type2 = rep(c("FS", "ExactFS", "ExactFS - init", "BAMLSS",  "BFGS","BFGS - init"), each=length(dgrid)))

TIME_efs_bamlss_bfgs <- data.frame(Value = c(data_time$time_efs, data_time$time_exact_efs, data_time$time_exact_efs + data_time$time_efs,
                                             data_time$time_bamlss, data_time$time_bfgs, data_time$time_bfgsinit + data_time$time_efs),
                                   d = rep(rep(dgrid, each = nrun), 6),
                                   Type2 = rep(c("FS", "ExactFS", "ExactFS - init", "BAMLSS",  "BFGS", "BFGS - init"), each=nrun * length(dgrid)))

# absolute differences
pl_LAML_diff_bfgs_efs <- ggplot(LAML, aes(x = d, y = diff_LAML_other_efs, color = init)) +
  geom_point(aes(colour = as.factor(init)), size = 1, show.legend = TRUE, position = position_dodge(width = 0.3)) +
  geom_hline(yintercept = 0, col = "black", lty = "dashed") +
  theme_bw() +
  facet_grid(. ~ "(LAML(other) - LAML(EFS))") +
  xlab("Dimension") + ylab("") +
  scale_color_manual(name = "Type", values = c("ExactFS" =  "#0072B2", "ExactFS - init" =  "#00A9FF", "BFGS" = "#AE4371", "BFGS - init" = "#F8766D")) +
  theme(panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12),  text = element_text(size = 15),
        legend.text=element_text(size=15), strip.text.x = element_text(size = 15),
        panel.grid.major = element_blank(),legend.position="bottom",panel.spacing = unit(0.1, "lines"))


# Decide whether plotting a line for the mean values

#pl_LAML_rel_diff_bfgs_efs <- ggplot(LAML, aes(x = d, y = rel_diff_LAML_bfgs_efs, color = init)) +
#  geom_point(aes(colour = as.factor(init)), size = 1, show.legend = TRUE, position = position_dodge(width = 0.3)) +
#  geom_hline(yintercept = 0, col = "black", lty = "dashed") +
#  theme_bw() +
#  facet_grid(. ~ "(LAML(BFGS) - LAML(EFS)) / |LAML(BFGS)|") +
#  xlab("Dimension") + ylab("") +
#  scale_color_manual(name = "Type", values = c("BFGS" = "#00A9FF", "BFGS - init" = "#F8766D")) +
#  theme(panel.grid.minor = element_blank(),
#        axis.text = element_text(size = 12),  text = element_text(size = 15),
#        legend.text=element_text(size=15), strip.text.x = element_text(size = 15),
#        panel.grid.major = element_blank(),legend.position="bottom",panel.spacing = unit(0.1, "lines"))
pl_LAML_rel_diff_bfgs_efs <- ggplot(LAML, aes(x = d, y = rel_diff_LAML_other_efs, color = init)) +
  geom_point(aes(colour = as.factor(init)), size = 1, show.legend = TRUE, position = position_dodge(width = 0.3)) +
  geom_hline(yintercept = 0, col = "black", lty = "dashed") +
  theme_bw() +
  facet_grid(. ~ "(LAML(other) - LAML(EFS)) / |LAML(other)|") +
  xlab("Dimension") + ylab("") +
  scale_color_manual(name = "Type", values = c("ExactFS" =  "#0072B2", "ExactFS - init" =  "#00A9FF", "BFGS" = "#AE4371", "BFGS - init" = "#F8766D")) +
  theme(panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12),  text = element_text(size = 15),
        legend.text=element_text(size=15), strip.text.x = element_text(size = 15),
        panel.grid.major = element_blank(),legend.position="bottom",panel.spacing = unit(0.1, "lines"))


# Comparison EFS v s BAMLSS in terms of computational times

#TIME_efs_bamlss_bfgs$d <- factor(TIME_efs_bamlss_bfgs$d, levels=as.character(dgrid_sel), labels=as.character(dgrid_sel))
#TIME_efs_bamlss_bfgs$Type2 <- factor(TIME_efs_bamlss_bfgs$Type2, levels=c("FS","BFGS", "BFGS - init", "BAMLSS"), labels=c("FS","BFGS",  "BFGS - init", "BAMLSS"))
#TIMES_sum$d <- factor(TIMES_sum$d, levels=as.character(dgrid_sel), labels=as.character(dgrid_sel))
#TIMES_sum$Type2 <- factor(TIMES_sum$Type2, levels=c("FS", "BFGS", "BFGS - init",  "BAMLSS"), labels=c("FS", "BFGS", "BFGS - init",  "BAMLSS"))

TIME_efs_bamlss_bfgs$d <- factor(TIME_efs_bamlss_bfgs$d, levels=as.character(dgrid_sel), labels=as.character(dgrid_sel))
TIME_efs_bamlss_bfgs$Type2 <- factor(TIME_efs_bamlss_bfgs$Type2, levels=c("FS", "ExactFS", "ExactFS - init", "BFGS", "BFGS - init", "BAMLSS"),
                                     labels=c("FS", "ExactFS", "ExactFS - init", "BFGS",  "BFGS - init", "BAMLSS"))
TIMES_sum$d <- factor(TIMES_sum$d, levels=as.character(dgrid_sel), labels=as.character(dgrid_sel))
TIMES_sum$Type2 <- factor(TIMES_sum$Type2, levels=c("FS", "ExactFS", "ExactFS - init", "BFGS", "BFGS - init",  "BAMLSS"),
                          labels=c("FS", "ExactFS", "ExactFS - init", "BFGS", "BFGS - init",  "BAMLSS"))

breaks_seq_time <- c(2, 10, 25, 50, 100, 200, 400)


# pl_Fit_Time_efs_bamlss_bfgs <- ggplot(TIME_efs_bamlss_bfgs,
#                                  aes(x = as.factor(d), y = Value)) +
#   geom_point(aes(colour = Type2), size = 1, show.legend = TRUE, position = position_dodge(width = 0.1)) +
#   geom_point(data = TIMES_sum,
#              aes(x = as.factor(d),
#                  y = Value, colour = Type2), size = 3,
#              position = position_dodge(width = 0.1), show.legend = FALSE) +
#   facet_grid(. ~ "Fitting times") +
#   geom_line(data = TIMES_sum, aes(y = Value, group = Type2, col = Type2),
#             position = position_dodge(width = 0.1)) +
#   scale_color_manual(name = "Type", values = c("FS" = "#E76BF3", "BFGS - init" = "#F8766D", "BFGS" = "#00A9FF", "BAMLSS" = "#7CAE00")) +
#   theme_bw() +
#   scale_y_continuous(breaks = NULL, trans = myscale_trans2(),
#                      sec.axis = sec_axis(~ . * 1,labels = scaleFUN,
#                                          breaks = breaks_seq_time)) +
#   scale_x_discrete(breaks = dgrid_sel) +
#   xlab("Dimension") + ylab("") +
#   theme(panel.grid.minor = element_blank(),
#         axis.text = element_text(size = 12),  text = element_text(size = 15),
#         legend.text=element_text(size=15), strip.text.x = element_text(size = 15),
#         panel.grid.major = element_blank(), legend.position = "bottom", panel.spacing = unit(0.2, "lines"))

pl_Fit_Time_efs_bamlss_bfgs <- ggplot(TIME_efs_bamlss_bfgs,
                                      aes(x = as.factor(d), y = Value)) +
  geom_point(aes(colour = Type2), size = 1, show.legend = TRUE, position = position_dodge(width = 0.1)) +
  geom_point(data = TIMES_sum,
             aes(x = as.factor(d),
                 y = Value, colour = Type2), size = 3,
             position = position_dodge(width = 0.1), show.legend = FALSE) +
  facet_grid(. ~ "Fitting times") +
  geom_line(data = TIMES_sum, aes(y = Value, group = Type2, col = Type2),
            position = position_dodge(width = 0.1)) +
  scale_color_manual(name = "Type", values = c("FS" = "#E76BF3", "ExactFS" = "#0072B2", "ExactFS - init" = "#AE4371",
                                               "BFGS - init" = "#F8766D", "BFGS" = "#00A9FF", "BAMLSS" = "#7CAE00")) +
  theme_bw() +
  scale_y_continuous(breaks = NULL, trans = myscale_trans2(),
                     sec.axis = sec_axis(~ . * 1,labels = scaleFUN,
                                         breaks = breaks_seq_time)) +
  scale_x_discrete(breaks = dgrid_sel) +
  xlab("Dimension") + ylab("") +
  theme(panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12),  text = element_text(size = 15),
        legend.text=element_text(size=15), strip.text.x = element_text(size = 15),
        panel.grid.major = element_blank(), legend.position = "bottom", panel.spacing = unit(0.2, "lines"))


pl_LAML_Times <- ggarrange(pl_LAML_rel_diff_bfgs_efs,
                     pl_Fit_Time_efs_bamlss_bfgs,
                     nrow = 1,
                     common.legend = FALSE,
                     legend = "bottom")

setwd(root_dir)
setwd("content/Section6/Plots")
ggsave("plot_LAML_Times.eps", pl_LAML_Times , width = 30, height = 15, units = "cm")
ggsave("plot_LAML_Times.pdf", pl_LAML_Times , width = 30, height = 15, units = "cm")


##########################################################################################
# Code for reproducing Figure XXX - SUPPLEMENTARY MATERIAL - comparison in terms of logS #
##########################################################################################

##########################################################################################
# Decide whether including the comparison in terms of logScore on the test set in the SM #
# If so, I will improve the graphical aspects                                            #
##########################################################################################
setwd(root_dir)
setwd("content/Section6")
source("Functions_Plots_Overall_Fit.R")

dgrid <- c(2, 5, 10)
nrun <- 10
nobs <- 10000

setwd(root_dir)
setwd("content/Section6/Results")
load(paste0("sim_mcd_fit_nrun_", nrun, "_n_", nobs, "_d_", paste0(dgrid, collapse = "_"), ".RData"))

# logS_mcd_test <- log_Score_test(sim_mcd_fit, nrun, dgrid, nobs,  param = "mcd")
# logS_gen_test <- as.numeric(logS_mcd_test[[1]])
# logS_efs_test <- as.numeric(logS_mcd_test[[2]])
# logS_bfgs_test <- as.numeric(logS_mcd_test[[3]])
# logS_bfgsinit_test <- as.numeric(logS_mcd_test[[4]])
# logS_bamlss_test <- as.numeric(logS_mcd_test[[5]])

# logS_test <- data.frame(logS_efs = logS_efs_test,
#                         logS_bfgs = logS_bfgs_test,
#                         logS_bfgsinit = logS_bfgsinit_test,
#                         logS_bamlss = logS_bamlss_test,
#                         diff_logS_bfgs_efs = logS_bfgs_test - logS_efs_test,
#                         diff_logS_bfgsinit_efs = logS_bfgsinit_test - logS_efs_test,
#                         diff_logS_bamlss_efs = logS_bamlss_test - logS_efs_test,
#                         rel_diff_bfgs_efs_gen = (logS_bfgs_test - logS_efs_test)/abs(logS_gen_test),
#                         rel_diff_bfgsinit_efs_gen = (logS_bfgsinit_test - logS_efs_test)/abs(logS_gen_test),
#                         rel_diff_bamlss_efs_gen = (logS_bamlss_test - logS_efs_test)/abs(logS_gen_test),
#                         d =  factor(rep(dgrid, each = nrun ), levels = dgrid))

logS_mcd_test <- log_Score_test(sim_mcd_fit, nrun, dgrid, nobs,  param = "mcd")
logS_gen_test <- as.numeric(logS_mcd_test[[1]])
logS_efs_test <- as.numeric(logS_mcd_test[[2]])
logS_exact_efs_test <- as.numeric(logS_mcd_test[[3]])
logS_exact_efs_initialised_test <- as.numeric(logS_mcd_test[[4]])
logS_bfgs_test <- as.numeric(logS_mcd_test[[5]])
logS_bfgsinit_test <- as.numeric(logS_mcd_test[[6]])
logS_bamlss_test <- as.numeric(logS_mcd_test[[7]])

logS_test <- data.frame(logS_efs = logS_efs_test,
                        logS_exact_efs = logS_exact_efs_test,
                        logS_exact_efs_initialised = logS_exact_efs_initialised_test,
                        logS_bfgs = logS_bfgs_test,
                        logS_bfgsinit = logS_bfgsinit_test,
                        logS_bamlss = logS_bamlss_test,
                        diff_logS_Eefs_efs = logS_exact_efs_test - logS_efs_test,
                        diff_logS_EIefs_efs = logS_exact_efs_initialised_test - logS_efs_test,
                        diff_logS_bfgs_efs = logS_bfgs_test - logS_efs_test,
                        diff_logS_bfgsinit_efs = logS_bfgsinit_test - logS_efs_test,
                        diff_logS_bamlss_efs = logS_bamlss_test - logS_efs_test,
                        rel_diff_Eefs_efs_gen = (logS_exact_efs_test - logS_efs_test)/abs(logS_gen_test),
                        rel_diff_EIefs_efs_gen = (logS_exact_efs_initialised_test - logS_efs_test)/abs(logS_gen_test),
                        rel_diff_bfgs_efs_gen = (logS_bfgs_test - logS_efs_test)/abs(logS_gen_test),
                        rel_diff_bfgsinit_efs_gen = (logS_bfgsinit_test - logS_efs_test)/abs(logS_gen_test),
                        rel_diff_bamlss_efs_gen = (logS_bamlss_test - logS_efs_test)/abs(logS_gen_test),
                        d =  factor(rep(dgrid, each = nrun ), levels = dgrid))

# y axis:  (Log-score(BAMLSS) - Log-score (EFS))/(|log-Score(Generation)|)
# pl_MCD_rel_diff_bamlss_efs_gen <- ggplot(logS_test, aes(x = d, y = rel_diff_bamlss_efs_gen)) +
#   geom_point(size = 1, show.legend = TRUE)+
#   geom_hline(yintercept = 0, col = "black", lty = "dashed") +
#   theme_bw() +
#   facet_grid(. ~ "(LS(BAMLSS) - LS(FS)) / |LS(gen)|") +
#   xlab("Dimension") + ylab("") +
#   theme(panel.grid.minor = element_blank(),
#         axis.text = element_text(size = 15),  text = element_text(size = 15),
#         legend.text=element_text(size=15), strip.text.x = element_text(size = 15),
#         panel.grid.major = element_blank(),legend.position="bottom",panel.spacing = unit(0.1, "lines"))
#
# pl_MCD_rel_diff_bfgsinit_efs_gen <- ggplot(logS_test, aes(x = d, y = rel_diff_bfgsinit_efs_gen)) +
#   geom_point(size = 1, show.legend = TRUE)+
#   geom_hline(yintercept = 0, col = "black", lty = "dashed") +
#   theme_bw() +
#   facet_grid(. ~ "(LS(BFGS - init) - LS(FS)) / |LS(gen)|") +
#   xlab("Dimension") + ylab("") +
#   theme(panel.grid.minor = element_blank(),
#         axis.text = element_text(size = 12),  text = element_text(size = 15),
#         legend.text=element_text(size=15), strip.text.x = element_text(size = 15),
#         panel.grid.major = element_blank(),legend.position="bottom",panel.spacing = unit(0.1, "lines"))
#
# pl_MCD_logS_test <- ggarrange(pl_MCD_rel_diff_bamlss_efs_gen,
#                               pl_MCD_rel_diff_bfgsinit_efs_gen,
#                               nrow = 1,
#                               common.legend = FALSE,
#                               legend = "bottom")
#
# setwd(root_dir)
# setwd("content/SupplementaryMaterial/Plots")
# ggsave("pl_MCD_logS_test.eps", pl_MCD_logS_test , width = 30, height = 15, units = "cm")
# ggsave("pl_MCD_logS_test.pdf", pl_MCD_logS_test , width = 30, height = 15, units = "cm")

pl_MCD_rel_diff_bamlss_efs_gen <- ggplot(logS_test, aes(x = d, y = rel_diff_bamlss_efs_gen)) +
  geom_point(size = 1, show.legend = TRUE)+
  geom_hline(yintercept = 0, col = "black", lty = "dashed") +
  theme_bw() +
  facet_grid(. ~ "(LS(BAMLSS) - LS(FS)) / |LS(gen)|") +
  xlab("Dimension") + ylab("") +
  theme(panel.grid.minor = element_blank(),
        axis.text = element_text(size = 15),  text = element_text(size = 15),
        legend.text=element_text(size=15), strip.text.x = element_text(size = 15),
        panel.grid.major = element_blank(),legend.position="bottom",panel.spacing = unit(0.1, "lines"))

pl_MCD_rel_diff_Eefs_efs_gen <- ggplot(logS_test, aes(x = d, y = rel_diff_Eefs_efs_gen)) +
  geom_point(size = 1, show.legend = TRUE)+
  geom_hline(yintercept = 0, col = "black", lty = "dashed") +
  theme_bw() +
  facet_grid(. ~ "(LS(Exact FS) - LS(FS)) / |LS(gen)|") +
  xlab("Dimension") + ylab("") +
  theme(panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12),  text = element_text(size = 15),
        legend.text=element_text(size=15), strip.text.x = element_text(size = 15),
        panel.grid.major = element_blank(),legend.position="bottom",panel.spacing = unit(0.1, "lines"))

pl_MCD_rel_diff_EIefs_efs_gen <- ggplot(logS_test, aes(x = d, y = rel_diff_EIefs_efs_gen)) +
  geom_point(size = 1, show.legend = TRUE)+
  geom_hline(yintercept = 0, col = "black", lty = "dashed") +
  theme_bw() +
  facet_grid(. ~ "(LS(Exact FS - init) - LS(FS)) / |LS(gen)|") +
  xlab("Dimension") + ylab("") +
  theme(panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12),  text = element_text(size = 15),
        legend.text=element_text(size=15), strip.text.x = element_text(size = 15),
        panel.grid.major = element_blank(),legend.position="bottom",panel.spacing = unit(0.1, "lines"))


pl_MCD_rel_diff_bfgs_efs_gen <- ggplot(logS_test, aes(x = d, y = rel_diff_bfgs_efs_gen)) +
  geom_point(size = 1, show.legend = TRUE)+
  geom_hline(yintercept = 0, col = "black", lty = "dashed") +
  theme_bw() +
  facet_grid(. ~ "(LS(BFGS) - LS(FS)) / |LS(gen)|") +
  xlab("Dimension") + ylab("") +
  theme(panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12),  text = element_text(size = 15),
        legend.text=element_text(size=15), strip.text.x = element_text(size = 15),
        panel.grid.major = element_blank(),legend.position="bottom",panel.spacing = unit(0.1, "lines"))


pl_MCD_rel_diff_bfgsinit_efs_gen <- ggplot(logS_test, aes(x = d, y = rel_diff_bfgsinit_efs_gen)) +
  geom_point(size = 1, show.legend = TRUE)+
  geom_hline(yintercept = 0, col = "black", lty = "dashed") +
  theme_bw() +
  facet_grid(. ~ "(LS(BFGS - init) - LS(FS)) / |LS(gen)|") +
  xlab("Dimension") + ylab("") +
  theme(panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12),  text = element_text(size = 15),
        legend.text=element_text(size=15), strip.text.x = element_text(size = 15),
        panel.grid.major = element_blank(),legend.position="bottom",panel.spacing = unit(0.1, "lines"))

pl_MCD_logS_test <- ggarrange(#pl_MCD_rel_diff_bamlss_efs_gen,
  #pl_MCD_rel_diff_Eefs_efs_gen,
  #pl_MCD_rel_diff_EIefs_efs_gen,
  pl_MCD_rel_diff_bfgsinit_efs_gen,
  pl_MCD_rel_diff_bfgs_efs_gen,
  nrow = 2,
  common.legend = FALSE,
  legend = "bottom")

setwd(root_dir)
setwd("content/SupplementaryMaterial/Plots")
ggsave("pl_MCD_logS_test.eps", pl_MCD_logS_test , width = 30, height = 15, units = "cm")
ggsave("pl_MCD_logS_test.pdf", pl_MCD_logS_test , width = 30, height = 15, units = "cm")

########################################################################
# Code for reproducing Figure XXX - SECTION 7 + SUPPLEMENTARY MATERIAL #
########################################################################
setwd(root_dir)
setwd("content/Section7")
source("Functions_Plots_Application.R")
source("DataPreprocessing.R")

# Set the train set also here (reply what is done in Main_new(par_Lapply))
n_train <- which(GEF14_data$year==2011)[1]-1 # (2005 - 2010 year)

d <- 24
grid_length <- 5
name_eff <- "doy"

# Set the rolling origin forecasting splitting
ndat <- which(GEF14_data$year == 2011)[1] - 1
ndat2 <- dim(GEF14_data)[1]
sets <- floor(seq(ndat , ndat2, length.out = 12))


##############################
# Plots for validation steps #
##############################
setwd(root_dir)
setwd("content/Section7")
load("GEF14_data_residuals.RData")

low_neff_vcov <- 0
upp_neff_vcov <- 150  # Here must be setted according to the maximum value of the explored grid


flag_residuals <- TRUE
if(flag_residuals){
  param <- "mcd"
  # MCD

  load(paste0("Results/cv_res_stepwise_param", param, "_d_", d, "_lstep_",
              grid_length, "_low_thresh_", low_neff_vcov, "_upp_thresh_", upp_neff_vcov,"_residuals.RData"))
  load( paste0("Results/res_stepwise_param", param, "_d_", d, "_lstep_", grid_length, "_residuals.RData"))
  ncov_el_mcd <- sapply(1:length(res_mcd[[1]]), function(x) length(res_mcd$foo[[x]]) - d)
  ncov_el_mcd <- sort(ncov_el_mcd[ncov_el_mcd >= low_neff_vcov & ncov_el_mcd <= upp_neff_vcov], decreasing = TRUE)

  logScore_mcd_residuals <- res_perf(cv_mcd_residuals, d, GEF14_data_residuals, param, sets)

  # logM
  param <- "logm"
  load(paste0("Results/cv_res_stepwise_param", param, "_d_", d, "_lstep_",
               grid_length, "_low_thresh_", low_neff_vcov, "_upp_thresh_", upp_neff_vcov,"_residuals.RData"))
  load(paste0("Results/res_stepwise_param", param, "_d_", d, "_lstep_", grid_length, "_residuals.RData"))
  ncov_el_logm <- sapply(1:length(res_logm[[1]]), function(x) length(res_logm$foo[[x]]) - d)
  ncov_el_logm <- sort(ncov_el_logm[ncov_el_logm >= low_neff_vcov & ncov_el_logm <= upp_neff_vcov], decreasing = TRUE)

  logScore_logm_residuals <- res_perf(cv_logm_residuals, d, GEF14_data_residuals, param, sets)

  data_logScore_residuals <- data.frame("logS" = c(logScore_mcd_residuals[length(logScore_mcd_residuals):1], logScore_logm_residuals[length(logScore_logm_residuals):1]),
                              "Param" = c(rep("MCD", length(logScore_mcd_residuals)), rep("logM", length(ncov_el_logm))),
                              "Eff" = c(sort(ncov_el_mcd)[1:length(logScore_mcd_residuals)], sort(ncov_el_logm)))

  with(data_logScore_residuals[data_logScore_residuals$Param == "MCD", ],{
    plot(Eff, logS, xaxt= "n", xlab = "Number of effects (MCD)", main = "Residuals" )
    axis(1, at = Eff, labels = factor(Eff, levels = Eff))
    points(Eff[which.min(logS)], logS[which.min(logS)], col = "red", pch = 19)
  }
  )

  with(data_logScore_residuals[data_logScore_residuals$Param == "logM",],{
    plot(Eff, logS, xaxt= "n", xlab = "Number of effects (MCD)", main = "Residuals" )
    axis(1, at = Eff, labels = factor(Eff, levels = Eff))
    points(Eff[which.min(logS)], logS[which.min(logS)], col = "red", pch = 19)
  }
  )
} else {
  param <- "mcd"
  load(paste0("Results/cv_res_stepwise_param", param, "_d_", d, "_lstep_",
              grid_length, "_low_thresh_", low_neff_vcov, "_upp_thresh_", upp_neff_vcov,"_response.RData"))
  load( paste0("Results/res_stepwise_param", param, "_d_", d, "_lstep_", grid_length, "_response.RData"))
  ncov_el_mcd <- sapply(1:length(res_mcd[[1]]), function(x) length(res_mcd$foo[[x]]) - d)
  ncov_el_mcd <- sort(ncov_el_mcd[ncov_el_mcd >= low_neff_vcov & ncov_el_mcd <= upp_neff_vcov], decreasing = TRUE)
  cv_mcd <- cv_mcd_response
  logScore_mcd <- res_perf(cv_mcd, d, GEF14_data, param, sets)

  param <- "logm"
  load(paste0("Results/cv_res_stepwise_param", param, "_d_", d, "_lstep_",
              grid_length, "_low_thresh_", low_neff_vcov, "_upp_thresh_", upp_neff_vcov,"_response.RData"))
  load( paste0("Results/res_stepwise_param", param, "_d_", d, "_lstep_", grid_length, "_response.RData"))
  ncov_el_logm <- sapply(1:length(res_logm[[1]]), function(x) length(res_logm$foo[[x]]) - d)
  ncov_el_logm <- sort(ncov_el_logm[ncov_el_logm >= low_neff_vcov & ncov_el_logm <= upp_neff_vcov], decreasing = TRUE)
  cv_logm <- cv_logm_response
  logScore_logm <- res_perf(cv_logm, d, GEF14_data, param, sets)


  data_logScore <- data.frame("logS" = c(logScore_mcd[length(logScore_mcd):1], logScore_logm[length(logScore_logm):1]),
                              "Param" = c(rep("MCD", length(logScore_mcd)), rep("logM", length(ncov_el_logm))),
                              "Eff" = c(sort(ncov_el_mcd)[1:length(logScore_mcd)], sort(ncov_el_logm)))


  with(data_logScore[data_logScore$Param == "MCD", ],{
       plot(Eff, logS, xaxt= "n", xlab = "Number of effects (MCD)", main = "Response" )
       axis(1, at = Eff, labels = factor(Eff, levels = Eff))
       points(Eff[which.min(logS)], logS[which.min(logS)], col = "red", pch = 19)
  }
  )
  with(data_logScore[data_logScore$Param == "logM", ],{
    plot(Eff, logS, xaxt= "n", xlab = "Number of effects (MCD)", main = "Response" )
    axis(1, at = Eff, labels = factor(Eff, levels = Eff))
    points(Eff[which.min(logS)], logS[which.min(logS)], col = "red", pch = 19)
  }
  )

}

# Here we select the number of effects (minimum or elbow point)

setwd(root_dir)
setwd("content/Section7/Plots")


flag_residuals <- TRUE
if(flag_residuals){
  ncov_el_mcd_sel <- 65
  ncov_el_logm_sel <- 30
  data_logS_selected <- data.frame("logS" = c(logScore_mcd_residuals[which(ncov_el_mcd == ncov_el_mcd_sel)],
                                              logScore_logm_residuals[which(ncov_el_logm == ncov_el_logm_sel)]),
                                   "Param" = c("MCD","logM"),
                                   "Eff" = c(ncov_el_mcd_sel, ncov_el_logm_sel))
  pl_logS_mcd_logM <- ggplot(data_logScore_residuals,
                             aes(x = Eff, y = logS)) +
    geom_point(aes(colour = Param), size = 1, show.legend = TRUE, position = position_dodge(width = 0.3)) +
    geom_point(data = data_logS_selected, aes(x = c(ncov_el_mcd_sel, ncov_el_logm_sel), y = logS, colour = Param),
      size = 3, show.legend = FALSE, position = position_dodge(width = 0.3)) +
    geom_line(aes(y = logS, group = Param, col = Param), position = position_dodge(width = 0.3)) +
    scale_x_continuous(breaks = seq(low_neff_vcov, upp_neff_vcov, by = grid_length)) +
    scale_color_manual(name = "Parametrisation", values = c("MCD" = "#F8766D", "logM" = "#619CFF")) +
    theme_bw() +
    xlab("Number of effects") + ylab("") +
    theme(panel.grid.minor = element_blank(), axis.text.y = element_text(size = 12),
          axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5),  text = element_text(size = 15),
          legend.text=element_text(size=15), strip.text.x = element_text(size = 15),
          panel.grid.major = element_blank(), legend.position = "bottom", panel.spacing = unit(0.2, "lines"))
  pl_logS_mcd_logM
  ggsave(paste0("logS_MCDandlogM_residuals.pdf"),  plot = pl_logS_mcd_logM, width = 20, height = 20, units = "cm")
  ggsave(paste0("logS_MCDandlogM_residuals.eps"),  plot = pl_logS_mcd_logM, width = 20, height = 20, units = "cm")
} else {
  ncov_el_mcd_sel <- 80
  ncov_el_logm_sel <- 40

  data_logS_selected <- data.frame("logS" = c(logScore_mcd[which(ncov_el_mcd == ncov_el_mcd_sel)],
                                              logScore_logm[which(ncov_el_logm == ncov_el_logm_sel)]),
                                   "Param" = c("MCD","logM"),
                                   "Eff" = c(ncov_el_mcd_sel, ncov_el_logm_sel))
  pl_logS_mcd_logM <- ggplot(data_logScore,
                             aes(x = Eff, y = logS)) +
    geom_point(aes(colour = Param), size = 1, show.legend = TRUE, position = position_dodge(width = 0.3)) +
    geom_point(data = data_logS_selected, aes(x = c(ncov_el_mcd_sel, ncov_el_logm_sel), y = logS, colour = Param),
      size = 3, show.legend = FALSE, position = position_dodge(width = 0.3)) +
    geom_line(aes(y = logS, group = Param, col = Param), position = position_dodge(width = 0.3)) +
    scale_x_continuous(breaks = seq(low_neff_vcov, upp_neff_vcov, by = grid_length)) +
    scale_color_manual(name = "Parametrisation", values = c("MCD" = "#F8766D", "logM" = "#619CFF")) +
    theme_bw() +
    xlab("Number of effects") + ylab("") +
    theme(panel.grid.minor = element_blank(), axis.text.y = element_text(size = 12),
          axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5),  text = element_text(size = 15),
          legend.text=element_text(size=15), strip.text.x = element_text(size = 15),
          panel.grid.major = element_blank(), legend.position = "bottom", panel.spacing = unit(0.2, "lines"))
  pl_logS_mcd_logM
  ggsave(paste0("logS_MCDandlogM_response.pdf"),  plot = pl_logS_mcd_logM, width = 20, height = 20, units = "cm")
  ggsave(paste0("logS_MCDandlogM_response.eps"),  plot = pl_logS_mcd_logM, width = 20, height = 20, units = "cm")
}



###############################################
# Plots for showing the model selection steps #
###############################################


###############################################
# Both MCD and logM in the same plot          #
###############################################
flag_residuals <- TRUE
if(flag_residuals){
  setwd(root_dir)
  setwd("content/Section7/Results")

  # According to the previous results you must select the number of effects for
  # the MCD and logM covariance matrix model
  neff_mcd <- 65 #ncov_el_mcd_sel
  neff_logm <- 30 #ncov_el_logm_sel


  param <- "mcd"
  load( paste0("res_stepwise_param", param, "_d_", d, "_lstep_", grid_length, "_residuals.RData"))

  ncov_el_mcd <- sapply(1: length(res_mcd[[1]]), function(x) length(res_mcd$foo[[x]]) - d)

  param <- "logm"
  load( paste0("res_stepwise_param", param, "_d_", d, "_lstep_", grid_length, "_residuals.RData"))
  ncov_el_logm <- sapply(1: length(res_logm[[1]]), function(x) length(res_logm$foo[[x]]) - d)

  setwd(root_dir)
  setwd("content/Section7/Plots")

  pl_MCD_logM <- get_plots2(obj1_mcd = res_mcd,
                            obj2_logm = res_logm,
                            name_eff = name_eff,
                            d = d,
                            grid_length = grid_length,
                            neff1_mcd = neff_mcd,
                            neff2_logm = neff_logm
  )
  print(pl_MCD_logM)

  ggsave(paste0("Covmod_MCDandlogM_residuals.pdf"),  plot = pl_MCD_logM, width = 20, height = 20, units = "cm")
  ggsave(paste0("Covmod_MCDandlogM_residuals.eps"),  plot = pl_MCD_logM, width = 20, height = 20, units = "cm")
  setwd(root_dir)
  setwd("content/Section7")

} else {
  setwd(root_dir)
  setwd("content/Section7/Results")


  #According to the previous results you must select the number of effects for the MCD and logM covariance matrix model

  neff_mcd <- 80 #ncov_el_mcd_sel
  neff_logm <- 40 #ncov_el_logm_sel


  param <- "mcd"
  load( paste0("res_stepwise_param", param, "_d_", d, "_lstep_", grid_length, "_response.RData"))

  ncov_el_mcd <- sapply(1: length(res_mcd[[1]]), function(x) length(res_mcd$foo[[x]]) - d)

  param <- "logm"
  load( paste0("res_stepwise_param", param, "_d_", d, "_lstep_", grid_length, "_response.RData"))
  ncov_el_logm <- sapply(1: length(res_logm[[1]]), function(x) length(res_logm$foo[[x]]) - d)

  setwd(root_dir)
  setwd("content/Section7/Plots")

  pl_MCD_logM <- get_plots2(obj1_mcd = res_mcd,
                            obj2_logm = res_logm,
                            name_eff = name_eff,
                            d = d,
                            grid_length = grid_length,
                            neff1_mcd = neff_mcd,
                            neff2_logm =  neff_logm
  )
  print(pl_MCD_logM)

  ggsave(paste0("Covmod_MCDandlogM_response_VM.pdf"),  plot = pl_MCD_logM, width = 20, height = 20, units = "cm")
  ggsave(paste0("Covmod_MCDandlogM_response_VM.eps"),  plot = pl_MCD_logM, width = 20, height = 20, units = "cm")
  setwd(root_dir)
  setwd("content/Section7")

}


#################################################################################
# FULL MODEL SELECTION PROCESS:                                                 #
# uncomment below if you want see the evolution of the model selection process  #
# !!! You must create the folders to save the results                           #
# "content/Section7/Plots/SelectionProcess/MCD"                                 #
# content/Section7/Plots/SelectionProcess/logM#                                 #
#################################################################################

#######################
# MCD parametrisation #
#######################
# param <- "mcd"
#
# flag_residuals <- TRUE
# if(flag_residuals){
#   setwd(root_dir)
#   setwd("content/Section7/Results")
#   load( paste0("res_stepwise_param", param, "_d_", d, "_lstep_", grid_length, "_residuals.RData"))
#
#   ncov_el_mcd <- sapply(1: length(res_mcd[[1]]), function(x) length(res_mcd$foo[[x]]) - d)
#
#   setwd(root_dir)
#   setwd("content/Section7/Plots/SelectionProcess/MCD")
#
#   pl_list <- get_plots(obj = res_mcd,
#                        name_eff = name_eff,
#                        d = d,
#                       grid_length = grid_length,
#                       param = param)
#   for(j in 1:length(pl_list)){
#     ggsave(paste0("Residuals_Covmod_", param, "param_with", length(res_mcd$foo[[j+1]])-d, "Effects.pdf"),  plot=pl_list[[j]], width = 20, height = 20, units = "cm")
#   }
#   setwd(root_dir)
#   setwd("content/Section7")
# } else {
#   setwd(root_dir)
#   setwd("content/Section7/Results")
#   load( paste0("res_stepwise_param", param, "_d_", d, "_lstep_", grid_length, "_response.RData"))
#
#   ncov_el_mcd <- sapply(1: length(res_mcd[[1]]), function(x) length(res_mcd$foo[[x]]) - d)
#
#   setwd(root_dir)
#   setwd("content/Section7/Plots/SelectionProcess/MCD")
#
#   pl_list <- get_plots(obj = res_mcd,
#                        name_eff = name_eff,
#                        d = d,
#                        grid_length = grid_length,
#                        param = param)
#   for(j in 1:length(pl_list)){
#     ggsave(paste0("Response_Covmod_", param, "param_with", length(res_mcd$foo[[j+1]])-d, "Effects.pdf"),  plot=pl_list[[j]], width = 20, height = 20, units = "cm")
#   }
#   setwd(root_dir)
#   setwd("content/Section7")
# }
#
# ########################
# # logM parametrisation #
# ########################
# param <- "logm"
#
# flag_residuals <- TRUE
# if(flag_residuals){
#   setwd(root_dir)
#   setwd("content/Section7/Results")
#   load( paste0("res_stepwise_param", param, "_d_", d, "_lstep_", grid_length, "_residuals.RData"))
#
#   ncov_el_logm <- sapply(1: length(res_logm[[1]]), function(x) length(res_logm$foo[[x]]) - d)
#
#   setwd(root_dir)
#   setwd("content/Section7/Plots/SelectionProcess/logM")
#
#   pl_list <- get_plots(obj = res_logm,
#                        name_eff = name_eff,
#                        d = d,
#                        grid_length = grid_length,
#                        param = param)
#   for(j in 1:length(pl_list)){
#     ggsave(paste0("Residuals_Covmod_", param, "param_with", length(res_logm$foo[[j+1]])-d, "Effects.pdf"),  plot=pl_list[[j]], width = 20, height = 20, units = "cm")
#   }
#   setwd(root_dir)
#   setwd("content/Section7")
# } else {
#   setwd(root_dir)
#   setwd("content/Section7/Results")
#   load( paste0("res_stepwise_param", param, "_d_", d, "_lstep_", grid_length, "_response.RData"))
#
#   ncov_el_logm <- sapply(1: length(res_logm[[1]]), function(x) length(res_logm$foo[[x]]) - d)
#
#   setwd(root_dir)
#   setwd("content/Section7/Plots/SelectionProcess/logM")
#
#   pl_list <- get_plots(obj = res_logm,
#                        name_eff = name_eff,
#                        d = d,
#                        grid_length = grid_length,
#                        param = param)
#   for(j in 1:length(pl_list)){
#     ggsave(paste0("Response_Covmod_", param, "param_with", length(res_logm$foo[[j+1]])-d, "Effects.pdf"),  plot=pl_list[[j]], width = 20, height = 20, units = "cm")
#   }
#   setwd(root_dir)
#   setwd("content/Section7")
# }



############################################################################
#Computational times
############################################################################



flag_residuals <- FALSE
if(flag_residuals){
  param <- "mcd"
  setwd(root_dir)
  setwd("content/Section7/Results")
  load( paste0("res_stepwise_param", param, "_d_", d, "_lstep_", grid_length, "_residuals.RData"))


  time_mcd <- unlist(res_mcd$time_fit)/(1e9 * 60)

  param <- "logm"
  load( paste0("res_stepwise_param", param, "_d_", d, "_lstep_", grid_length, "_residuals.RData"))
  time_logm <- unlist(res_logm$time_fit)/(1e9 * 60)

  grid_d <- seq( 0, d*(d+1)/2, by = grid_length)
  data_time <- data.frame("Time" = c(time_mcd[length(time_mcd):1], time_logm[length(time_logm):1]),
                          "Param" = c(rep("MCD", length(grid_d)), rep("logM", length(grid_d))),
                          "Eff" = rep(grid_d,2))

  setwd(root_dir)
  setwd("content/Section7/Plots")


  pl_TIME_mcd_logM <- ggplot(data_time,
                             aes(x = Eff, y = Time)) + #factor(Eff, labels = as.character(Eff), levels = as.character(Eff)), y = logS)) +
    geom_point(aes(colour = Param), size = 1, show.legend = TRUE, position = position_dodge(width = 0.3)) +
    geom_line(aes(y = Time, group = Param, col = Param), position = position_dodge(width = 0.3)) +
    scale_x_continuous(breaks = grid_d) +
    scale_color_manual(name = "Parametrisation", values = c("MCD" = "#F8766D", "logM" = "#619CFF")) +
    theme_bw() +
    xlab("Number of effects") + ylab("") +
    theme(panel.grid.minor = element_blank(), axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5),
          axis.text.y = element_text(size = 12), text = element_text(size = 15),
          legend.text=element_text(size=15), strip.text.x = element_text(size = 15),
          panel.grid.major = element_blank(), legend.position = "bottom", panel.spacing = unit(0.2, "lines"))
  print(pl_TIME_mcd_logM)
  ggsave(paste0("Time_MCDandlogM_residuals.pdf"),  plot = pl_TIME_mcd_logM, width = 20, height = 20, units = "cm")
  ggsave(paste0("Time_MCDandlogM_residuals.eps"),  plot = pl_TIME_mcd_logM, width = 20, height = 20, units = "cm")

  setwd(root_dir)
  setwd("content/Section7")

} else {
  param <- "mcd"
  setwd(root_dir)
  setwd("content/Section7/Results")
  load( paste0("res_stepwise_param", param, "_d_", d, "_lstep_", grid_length, "_response.RData"))

  time_mcd <- unlist(res_mcd$time_fit)/(1e9 * 60)

  param <- "logm"
  load( paste0("res_stepwise_param", param, "_d_", d, "_lstep_", grid_length, "_response.RData"))
  time_logm <- unlist(res_logm$time_fit)/(1e9 * 60)

  grid_d <- seq( 0, d*(d+1)/2, by = grid_length)
  data_time <- data.frame("Time" = c(time_mcd[length(time_mcd):1], time_logm[length(time_logm):1]),
                          "Param" = c(rep("MCD", length(grid_d)), rep("logM", length(grid_d))),
                          "Eff" = rep(grid_d, 2))

  setwd(root_dir)
  setwd("content/Section7/Plots")


  pl_TIME_mcd_logM <- ggplot(data_time,
                             aes(x = Eff, y = Time)) + #factor(Eff, labels = as.character(Eff), levels = as.character(Eff)), y = logS)) +
    geom_point(aes(colour = Param), size = 1, show.legend = TRUE, position = position_dodge(width = 0.3)) +
    geom_line(aes(y = Time, group = Param, col = Param), position = position_dodge(width = 0.3)) +
    scale_x_continuous(breaks = grid_d) +
    scale_color_manual(name = "Parametrisation", values = c("MCD" = "#F8766D", "logM" = "#619CFF")) +
    theme_bw() +
    xlab("Number of effects") + ylab("") +
    theme(panel.grid.minor = element_blank(), axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5),
          axis.text.y = element_text(size = 12), text = element_text(size = 15),
          legend.text=element_text(size=15), strip.text.x = element_text(size = 15),
          panel.grid.major = element_blank(), legend.position = "bottom", panel.spacing = unit(0.2, "lines"))
  print(pl_TIME_mcd_logM)
  ggsave(paste0("Time_MCDandlogM_response_VM.pdf"),  plot = pl_TIME_mcd_logM, width = 20, height = 20, units = "cm")
  ggsave(paste0("Time_MCDandlogM_response_VM.eps"),  plot = pl_TIME_mcd_logM, width = 20, height = 20, units = "cm")
  setwd(root_dir)
  setwd("content/Section7")
}


#################################################
# Plot stdev and correlations of the selected model
#################################################
library(ggnewscale)
library(lubridate)

grid_d <- seq(d*(d+1)/2,  0, by = -grid_length)

# here you select to visualisize "full", "reduced" or "fixed" cases
model <- "reduced"
# MCD parametrisation
param <- "mcd"

flag_residuals <- TRUE

if(flag_residuals){
  neff_mcd <- 65
  setwd(root_dir)
  setwd("content/Section7/Results")
  param <- "mcd"
  load( paste0("res_stepwise_param", param, "_d_", d, "_lstep_", grid_length, "_residuals.RData"))

  jElem_mcd <-  which(neff_mcd == grid_d)

  if(model == "full"){
    res_mcd_final <- gam_scm(res_mcd$foo[[1]], family = mvn_scm(d = d), data = GEF14_data_residuals)
  }
  if(model == "reduced"){
    res_mcd_final <- gam_scm(res_mcd$foo[[jElem_mcd]], family = mvn_scm(d = d), data = GEF14_data_residuals)
  }
  if(model == "fixed"){
    res_mcd_final <- gam_scm(res_mcd$foo[[length(grid_d)]], family = mvn_scm(d = d), data = GEF14_data_residuals)
  }

  Sigma_pred_fitD_MCD <- predict(res_mcd_final, type = "response")
  Sigma_predMat_fitD_MCD <-  Sigma_mat(Sigma_pred_fitD_MCD[,-c(1 : d)])

  idx_min_CorrD_MCD <- which.min(unlist(lapply(1 : length(Sigma_predMat_fitD_MCD), function(x) unlist(mean(Sigma_predMat_fitD_MCD[[x]][lower.tri(Sigma_predMat_fitD_MCD[[x]])])))))
  idx_max_CorrD_MCD <- which.max(unlist(lapply(1 : length(Sigma_predMat_fitD_MCD), function(x) unlist(mean(Sigma_predMat_fitD_MCD[[x]][lower.tri(Sigma_predMat_fitD_MCD[[x]])])))))

  col_cor <- rev(colorspace::sequential_hcl(palette = "Blues 3", n = 100))
  col_var <- rev(colorspace::sequential_hcl(palette = "Red", n = 10)[1 : 5])

  label_xaxis <- c(paste0("h0", 0:9), paste0("h", 10:(d-1)))
  label_yaxis <- c(paste0("h", (d-1):10), paste0("h0", 9:0))

  days <- c("Sunday", "Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday")

  if(model == "fixed"){
    pl1_MCD <- heatmap_FitCov(Sigma_predMat_fitD_MCD[[1]], d = d, range_var = NULL, range_corr = c(0, 1), label = "h",
                              label_xaxis = label_xaxis, label_yaxis = label_yaxis,  col_var = col_var, col_cor = col_cor)
    print(pl1_MCD)
    pl2_MCD <- heatmap_FitCov(Sigma_predMat_fitD_MCD[[1]], d = d, range_var = NULL, range_corr = c(0, 1), label = "h",
                              label_xaxis = label_xaxis, label_yaxis = label_yaxis,  col_var = col_var, col_cor = col_cor)
    print(pl2_MCD)
  } else {
    pl1_MCD <- heatmap_FitCov(Sigma_predMat_fitD_MCD[[idx_min_CorrD_MCD]], d = d, range_var = NULL, range_corr = c(0, 1), label = "h",
                              label_xaxis = label_xaxis, label_yaxis = label_yaxis,  col_var = col_var, col_cor = col_cor)

    date <- as.Date(paste(GEF14_data_residuals[idx_min_CorrD_MCD, "year"], GEF14_data_residuals[idx_min_CorrD_MCD, "doy"]), format = "%Y %j")
    dow <- days[wday(date)]
    pl1_MCD <- pl1_MCD + annotate("text",  x = (d-1) + 0.25, y = (d-1) + 0.25, label = paste(dow, date), vjust=1, hjust=1, cex = 7)
    print(pl1_MCD)
    pl2_MCD <- heatmap_FitCov(Sigma_predMat_fitD_MCD[[idx_max_CorrD_MCD]], d = d, range_var = NULL, range_corr = c(0, 1), label = "h",
                              label_xaxis = label_xaxis, label_yaxis = label_yaxis,  col_var = col_var, col_cor = col_cor)
    date <- as.Date(paste(GEF14_data[idx_max_CorrD_MCD,"year"], GEF14_data[idx_max_CorrD_MCD,"doy"]), format = "%Y %j")
    dow <- days[wday(date)]
    pl2_MCD <- pl2_MCD +  annotate("text",  x = (d-1) + 0.25, y = (d-1) + 0.25, label = paste(dow, date), vjust = 1, hjust = 1, cex = 7)
    print(pl2_MCD)

  }


  setwd(root_dir)
  setwd("content/Section7/Plots")

  ggsave(paste0("Residuals_Vcormat_lowCorr_param", param, "_model_", model, ".pdf"),  plot =  pl1_MCD, width = 20, height = 20, units = "cm")
  ggsave(paste0("Residuals_Vcormat_HighCorr_param", param, "_model_", model, ".pdf"),  plot =  pl2_MCD, width = 20, height = 20, units = "cm")
  setwd(root_dir)
  setwd("content/Section7")
} else {

  neff_mcd <- 100
  setwd(root_dir)
  setwd("content/Section7/Results")
  load( paste0("res_stepwise_param", param, "_d_", d, "_lstep_", grid_length, "_response.RData"))

  jElem_mcd <-  which(neff_mcd == grid_d)


  if(model == "full"){
    res_mcd_final <- gam_scm(res_mcd$foo[[1]], family = mvn_scm(d = d), data = GEF14_data)
  }
  if(model == "reduced"){
    res_mcd_final <- gam_scm(res_mcd$foo[[jElem_mcd]], family = mvn_scm(d = d), data = GEF14_data)
  }
  if(model == "fixed"){
    res_mcd_final <- gam_scm(res_mcd$foo[[length(grid_d)]], family = mvn_scm(d = d), data = GEF14_data)
  }

  Sigma_pred_fitD_MCD <- predict(res_mcd_final, type = "response")
  Sigma_predMat_fitD_MCD <-  Sigma_mat(Sigma_pred_fitD_MCD[,-c(1 : d)])

  idx_min_CorrD_MCD <- which.min(unlist(lapply(1 : length(Sigma_predMat_fitD_MCD),
                                               function(x) unlist(mean(Sigma_predMat_fitD_MCD[[x]][lower.tri(Sigma_predMat_fitD_MCD[[x]])])))))
  idx_max_CorrD_MCD <- which.max(unlist(lapply(1 : length(Sigma_predMat_fitD_MCD),
                                               function(x) unlist(mean(Sigma_predMat_fitD_MCD[[x]][lower.tri(Sigma_predMat_fitD_MCD[[x]])])))))

  col_cor <- rev(colorspace::sequential_hcl(palette = "Blues 3", n = 100))
  col_var <- rev(colorspace::sequential_hcl(palette = "Red", n = 10)[1 : 5])

  label_xaxis <- c(paste0("h0", 0:9), paste0("h", 10:(d-1)))
  label_yaxis <- c(paste0("h", (d-1):10), paste0("h0", 9:0))

  days <- c("Sunday", "Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday")

  setwd(root_dir)
  setwd("content/Section7/Plots")


  if(model == "fixed"){
    pl1_MCD <- heatmap_FitCov(Sigma_predMat_fitD_MCD[[1]], d = d, range_var = NULL, range_corr = c(0, 1), label = "h",
                              label_xaxis = label_xaxis, label_yaxis = label_yaxis,  col_var = col_var, col_cor = col_cor)
    print(pl1_MCD)
    pl2_MCD <- heatmap_FitCov(Sigma_predMat_fitD_MCD[[1]], d = d, range_var = NULL, range_corr = c(0, 1), label = "h",
                              label_xaxis = label_xaxis, label_yaxis = label_yaxis,  col_var = col_var, col_cor = col_cor)
    print(pl2_MCD)
  } else {
    pl1_MCD <- heatmap_FitCov(Sigma_predMat_fitD_MCD[[idx_min_CorrD_MCD]], d = d, range_var = NULL, range_corr = c(0, 1), label = "h",
                              label_xaxis = label_xaxis, label_yaxis = label_yaxis,  col_var = col_var, col_cor = col_cor)

    date <- as.Date(paste(GEF14_data[idx_min_CorrD_MCD, "year"], GEF14_data[idx_min_CorrD_MCD, "doy"]), format = "%Y %j")
    dow <- days[wday(date)]
    pl1_MCD <- pl1_MCD + annotate("text",  x = (d-1) + 0.25, y = (d-1) + 0.25, label = paste(dow, date), vjust=1, hjust=1, cex = 7)
    print(pl1_MCD)
    pl2_MCD <- heatmap_FitCov(Sigma_predMat_fitD_MCD[[idx_max_CorrD_MCD]], d = d, range_var = NULL, range_corr = c(0, 1), label = "h",
                              label_xaxis = label_xaxis, label_yaxis = label_yaxis,  col_var = col_var, col_cor = col_cor)
    date <- as.Date(paste(GEF14_data[idx_max_CorrD_MCD,"year"], GEF14_data[idx_max_CorrD_MCD,"doy"]), format = "%Y %j")
    dow <- days[wday(date)]
    pl2_MCD <- pl2_MCD +  annotate("text",  x = (d-1) + 0.25, y = (d-1) + 0.25, label = paste(dow, date), vjust = 1, hjust = 1, cex = 7)
    print(pl2_MCD)

  }

  ggsave(paste0("Response_Vcormat_lowCorr_param", param, "_model_", model, ".pdf"),  plot =  pl1_MCD, width = 20, height = 20, units = "cm")
  ggsave(paste0("Response_Vcormat_highCorr_param", param, "_model_", model, ".pdf"),  plot =  pl2_MCD, width = 20, height = 20, units = "cm")
  setwd(root_dir)
  setwd("content/Section7")

}


# param <- "logm"
#
# flag_residuals <- FALSE
#
# if(flag_residuals){
#   neff_logm <- 60
#   setwd(root_dir)
#   setwd("content/Section7/Results")
#   load( paste0("res_stepwise_param", param, "_d_", d, "_lstep_", grid_length, "_residuals.RData"))
#
#   jElem_logm <-  which(neff_logm == grid_d)
#
#   if(model == "full"){
#     res_logm_final <- gam_scm(res_logm$foo[[1]], family = mvn_scm(d = d, param = "logm"), data = GEF14_data_residuals)
#   }
#   if(model == "reduced"){
#     res_logm_final <- gam_scm(res_logm$foo[[jElem_logm]], family = mvn_scm(d = d, param = "logm"), data = GEF14_data_residuals)
#   }
#   if(model == "fixed"){
#     res_logm_final <- gam_scm(res_logm$foo[[length(grid_d)]], family = mvn_scm(d = d, param = "logm"), data = GEF14_data_residuals)
#   }
#
#   Sigma_pred_fitD_logN <- predict(res_logm_final, type = "response")
#   Sigma_predMat_fitD_logM <-  Sigma_mat(Sigma_pred_fitD_logM[,-c(1 : d)])
#
#   idx_min_CorrD_logM <- which.min(unlist(lapply(1 : length(Sigma_predMat_fitD_logM), function(x) unlist(mean(Sigma_predMat_fitD_logM[[x]][lower.tri(Sigma_predMat_fitD_logM[[x]])])))))
#   idx_max_CorrD_logM <- which.max(unlist(lapply(1 : length(Sigma_predMat_fitD_logM), function(x) unlist(mean(Sigma_predMat_fitD_logM[[x]][lower.tri(Sigma_predMat_fitD_logM[[x]])])))))
#
#   col_cor <- rev(colorspace::sequential_hcl(palette = "Blues 3", n = 100))
#   col_var <- rev(colorspace::sequential_hcl(palette = "Red", n = 10)[1 : 5])
#
#   label_xaxis <- c(paste0("h0", 0:9), paste0("h", 10:(d-1)))
#   label_yaxis <- c(paste0("h", (d-1):10), paste0("h0", 9:0))
#
#   days <- c("Sunday", "Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday")
#
#   pl1_logM <- heatmap_FitCov(round(Sigma_predMat_fitD_logM[[idx_min_CorrD_MCD]],2), d = d, range_var = NULL, range_corr = c(0, 1), label = "h",
#                             label_xaxis = label_xaxis, label_yaxis = label_yaxis,  col_var = col_var, col_cor = col_cor)
#
#   date <- as.Date(paste(GEF14_data_residuals[idx_min_CorrD_logM, "year"], GEF14_data_residuals[idx_min_CorrD_logM, "doy"]), format = "%Y %j")
#   dow <- days[wday(date)]
#   pl1_logM <- pl1_logM + annotate("text",  x = (d-1) + 0.25, y = (d-1) + 0.25, label = paste(dow, date), vjust=1, hjust=1, cex = 7)
#   print(pl1_logM)
#
#   pl2_logM <- heatmap_FitCov(Sigma_predMat_fitD_logM[[idx_max_CorrD_logM]], d = d, range_var = NULL, range_corr = c(0, 1), label = "h",
#                             label_xaxis = label_xaxis, label_yaxis = label_yaxis,  col_var = col_var, col_cor = col_cor)
#   date <- as.Date(paste(GEF14_data[idx_max_CorrD_logM,"year"], GEF14_data[idx_max_CorrD_logM,"doy"]), format = "%Y %j")
#   dow <- days[wday(date)]
#   pl2_logM <- pl2_logM +  annotate("text",  x = (d-1) + 0.25, y = (d-1) + 0.25, label = paste(dow, date), vjust = 1, hjust = 1, cex = 7)
#   print(pl2_logM)
#
#   setwd(root_dir)
#   setwd("content/Section7/Plots")
#
#   ggsave(paste0("Residuals_Vcormat_lowCorr_param", param, "_model_", model, ".pdf"),  plot =  pl1_MCD, width = 20, height = 20, units = "cm")
#   ggsave(paste0("Residuals_Vcormat_lowCorr_param", param, "_model_", model, ".pdf"),  plot =  pl2_MCD, width = 20, height = 20, units = "cm")
#   setwd(root_dir)
#   setwd("content/Section7")
# } else {
#
#   neff_logm <- 60
#
#   setwd(root_dir)
#   setwd("content/Section7/Results")
#   load( paste0("res_stepwise_param", param, "_d_", d, "_lstep_", grid_length, "_response.RData"))
#
#   jElem_logm <-  which(neff_logm == grid_d)
#
#
#   if(model == "full"){
#     res_logm_final <- gam_scm(res_logm$foo[[1]], family = mvn_scm(d = d, param = "logm"), data = GEF14_data_residuals)
#   }
#   if(model == "reduced"){
#     res_logm_final <- gam_scm(res_logm$foo[[jElem_logm]], family = mvn_scm(d = d, param = "logm"), data = GEF14_data_residuals)
#   }
#   if(model == "fixed"){
#     res_logm_final <- gam_scm(res_logm$foo[[length(grid_d)]], family = mvn_scm(d = d, param = "logm"), data = GEF14_data_residuals)
#   }
#
#   Sigma_pred_fitD_logM <- predict(res_logm_final, type = "response")
#   Sigma_predMat_fitD_logM <-  Sigma_mat(Sigma_pred_fitD_logM[,-c(1 : d)])
#
#   idx_min_CorrD_logM <- which.min(unlist(lapply(1 : length(Sigma_predMat_fitD_logM),
#                                                 function(x) unlist(mean(Sigma_predMat_fitD_logM[[x]][lower.tri(Sigma_predMat_fitD_logM[[x]])])))))
#   idx_max_CorrD_logM <- which.max(unlist(lapply(1 : length(Sigma_predMat_fitD_logM),
#                                                 function(x) unlist(mean(Sigma_predMat_fitD_logM[[x]][lower.tri(Sigma_predMat_fitD_logM[[x]])])))))
#
#   col_cor <- rev(colorspace::sequential_hcl(palette = "Blues 3", n = 100))
#   col_var <- rev(colorspace::sequential_hcl(palette = "Red", n = 10)[1 : 5])
#
#   label_xaxis <- c(paste0("h0", 0:9), paste0("h", 10:(d-1)))
#   label_yaxis <- c(paste0("h", (d-1):10), paste0("h0", 9:0))
#
#   days <- c("Sunday", "Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday")
#   setwd(root_dir)
#   setwd("content/Section7/Plots")
#
#   pl1_logM <- heatmap_FitCov(round(Sigma_predMat_fitD_logM[[idx_min_CorrD_logM]],2), d = d, range_var = NULL, range_corr = c(0, 1), label = "h",
#                             label_xaxis = label_xaxis, label_yaxis = label_yaxis,  col_var = col_var, col_cor = col_cor)
#
#   date <- as.Date(paste(GEF14_data[idx_min_CorrD_logM, "year"], GEF14_data[idx_min_CorrD_logM, "doy"]), format = "%Y %j")
#   dow <- days[wday(date)]
#   pl1_logM <- pl1_logM + annotate("text",  x = (d-1) + 0.25, y = (d-1) + 0.25, label = paste(dow, date), vjust=1, hjust=1, cex = 7)
#   print(pl1_logM)
#
#   pl2_logM <- heatmap_FitCov(Sigma_predMat_fitD_logM[[idx_max_CorrD_logM]], d = d, range_var = NULL, range_corr = c(0, 1), label = "h",
#                             label_xaxis = label_xaxis, label_yaxis = label_yaxis,  col_var = col_var, col_cor = col_cor)
#   date <- as.Date(paste(GEF14_data[idx_max_CorrD_logM,"year"], GEF14_data[idx_max_CorrD_logM,"doy"]), format = "%Y %j")
#   dow <- days[wday(date)]
#   pl2_logM <- pl2_logM +  annotate("text",  x = (d-1) + 0.25, y = (d-1) + 0.25, label = paste(dow, date), vjust = 1, hjust = 1, cex = 7)
#   print(pl2_logM)
#
#   ggsave(paste0("Response_Vcormat_lowCorr_param", param, "_model_", reduced, ".pdf"),  plot =  pl1_MCD, width = 20, height = 20, units = "cm")
#   ggsave(paste0("Response_Vcormat_lowCorr_param", param, "_model_", reduced, ".pdf"),  plot =  pl2_MCD, width = 20, height = 20, units = "cm")
#   setwd(root_dir)
#   setwd("content/Section7")
# }
