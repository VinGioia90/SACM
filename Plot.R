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

# Cpp function to speed up computations
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
scaleFUN2 <- function(x) sprintf("%.2f", x/1000) # Number of decimal points set to 2

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
## Code for reproducing Figure 1                                       ##
#########################################################################
#########################################################################
nobs <- 1000
dgrid <- seq(5, 50, by = 5)

setwd(root_dir)
setwd("content/Section3/Results")

###############################################
load(file = paste0("TIME_mcd_D2eta_dgrid_min_", min(dgrid), "_max_", max(dgrid), "nobs", nobs, ".RData"))
load(file = paste0("TIME_logm_D2eta_dgrid_min_", min(dgrid), "_max_", max(dgrid), "nobs", nobs, ".RData"))

setwd(root_dir)
setwd("content/Section3/Comp_logM_MCD_Hessian_eta")
source("Functions_Plots_Hessian_eta.R") # Here there is the time_hessian function

#################################################
# MCD Computational times (overall and summary) #
#################################################
time_hessian_mcd <- time_hessian(obj = TIME_MCD_D2eta, param = "mcd",
                                 dgrid = dgrid, nrun = nrun,
                                 type = c("eff", "TMB"),
                                 type1 =  c("EFF", "AD"))


summary_time_hessian_mcd <- time_hessian_mcd$sum_res
summary_time_hessian_mcd <- summary_time_hessian_mcd[order(summary_time_hessian_mcd[, 5],
                                                           summary_time_hessian_mcd[, 1],
                                                           decreasing = FALSE),]
row.names(summary_time_hessian_mcd) <- NULL

all_time_hessian_mcd <- time_hessian_mcd$res
all_time_hessian_mcd <- all_time_hessian_mcd[order(all_time_hessian_mcd[, 3],
                                                   all_time_hessian_mcd[, 2],
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

#####################################
# Note!!!: Labels to be set by hand #
#####################################
lab_time_logM <- c(0, 1, 10, 25, 50, 100, 200)
lab_time_MCD <- c(0, 0.05, 0.25, 0.5, 1, 2, 3.5)

################################
# Figure 1 (square root scale) #
################################
pl_Heta_logM <- ggplot(all_time_hessian[all_time_hessian$param == "logM",],
                       aes(x = as.factor(d), y = time)) +
  geom_point(aes(colour = Type), size = 1, show.legend = TRUE,
             position = position_dodge(width = 0.3)) +
  geom_point(data = summary_time_hessian[summary_time_hessian$param == "logM",],
             aes(x = as.factor(d), y = mean_time, colour = Type),
             size = 2, position = position_dodge(width = 0.3), show.legend = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(data = summary_time_hessian[summary_time_hessian$param == "logM",],
            aes(y = mean_time, group = Type, col = Type),
            position = position_dodge(width = 0.3))+
  scale_color_manual(name = "Approach", values = c("AD" = "#00A9FF", "EFF" = "#F8766D")) +
  theme_bw() +
  scale_y_continuous(breaks = NULL, trans = myscale_trans2(), name = "",
                     sec.axis = sec_axis(~ . * 1, breaks = lab_time_logM, name = "")) +
  scale_x_discrete(breaks = dg_sel) +
  xlab("Dimension") +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        legend.position = "bottom",
        panel.spacing = unit(0.2, "lines"),
        axis.text = element_text(size = 12),
        text = element_text(size = 15),  legend.text=element_text(size = 15),
        strip.text = element_blank(), strip.background = element_blank(),
        axis.title.x = element_blank(),
        plot.margin=unit(c(0.5, -0.5, 1, 0), "cm"),
        axis.title.y.right = element_text(margin = margin(t = 0, r = 0, b = 0, l = 5)))#,

pl_Heta_MCD <- ggplot(all_time_hessian[all_time_hessian$param == "MCD",],
                      aes(x = as.factor(d), y = time)) +
  geom_point(aes(colour = Type), size = 1, show.legend = TRUE,
             position = position_dodge(width = 0.3)) +
  geom_point(data = summary_time_hessian[summary_time_hessian$param == "MCD",],
             aes(x = as.factor(d), y = mean_time, colour = Type),
             size = 2, position = position_dodge(width = 0.3), show.legend = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_line(data = summary_time_hessian[summary_time_hessian$param == "MCD",],
            aes(y = mean_time, group = Type, col = Type),
            position = position_dodge(width = 0.3))+
  scale_color_manual(name = "Approach", values = c("AD" = "#00A9FF", "EFF" = "#F8766D")) +
  theme_bw() +
  scale_y_continuous(breaks = NULL, trans = myscale_trans2(), name = "",
                     sec.axis = sec_axis(~ . * 1, breaks = lab_time_MCD,
                                         name = "Time (minutes)")) +
  scale_x_discrete(breaks = dg_sel) +
  xlab("Dimension") +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        legend.position = "bottom", panel.spacing = unit(0.2, "lines"),
        axis.text = element_text(size = 12),
        axis.title.x = element_blank(),
        text = element_text(size = 15),  legend.text=element_text(size = 15),
        strip.text = element_blank(), strip.background = element_blank(),
        plot.margin=unit(c(0.5, 0.10, 1, -0.5), "cm"),
        axis.title.y.right = element_text(margin = margin(t = 0, r = 0, b = 0, l = 5)))

pl_Heta <- ggarrange(pl_Heta_logM,
                     pl_Heta_MCD,
                     nrow = 1,
                     common.legend = TRUE,
                     legend = "bottom",
                     widths = c(1.0, 1.025)) +
  theme(legend.box.spacing = unit(40, "pt"))

pl_Heta <- annotate_figure(pl_Heta,
                           bottom = textGrob("Dimension",hjust = 0.7,
                                             vjust = -4, gp = gpar(cex = 1.3)))
setwd(root_dir)
setwd("content/Section3/Plots")
ggsave("plot_TIME_hessian_eta_new.eps", pl_Heta, width = 30, height = 15, units = "cm")
ggsave("plot_TIME_hessian_eta_new.pdf", pl_Heta, width = 30, height = 15, units = "cm")

############################################################################
############################################################################
## Computational times for fitting under the MCD and logM parametrisation ##
############################################################################
############################################################################

###############################################
# Code for reproducing Figure 2 - SECTION 3.3 #
###############################################
dgrid <- c(2, 5, 10, 15, 20)
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

####################################################
# Merge times and iterations (overall and summary) #
####################################################
data_time_iter <- rbind(all_TIME, all_ITER)
data_time_iter_sum <-  rbind(mean_TIME[, c(1, 2, 5, 6)], mean_ITER[, c(1, 2, 5, 6)])

data_hline <- data.frame(Type2 = c("Time","Iterations"),  # Create data for horizontal lines
                         hline = c(0, NA))
data_hline$Type2 <- factor(data_hline$Type2, levels  = c("Time", "Iterations"))

dgrid_sel <- dgrid[-1] # To select a subset of dgrid (I removed the d=2 case)
data_time_iter <- data_time_iter[which(data_time_iter$d %in% dgrid_sel),]
data_time_iter_sum <- data_time_iter_sum[which(data_time_iter_sum$d %in% dgrid_sel),]

data_time_iter$Type2 <- factor(data_time_iter$Type2,
                               levels = c("Time", "Iterations"),
                               labels = c("Time", "Iterations"))
data_time_iter_sum$Type2 <- factor(data_time_iter_sum$Type2,
                                   levels = c("Time", "Iterations"),
                                   labels = c("Time", "Iterations"))

# Label to be set manually
label_time <- c(1, 5, 25, 50, 100, 500)

################################
# Figure 2 (square root scale) #
################################
pl_Fit_Time <- ggplot(data_time_iter[data_time_iter$Type2 == "Time",],
                      aes(x = factor(d, labels = as.character(dgrid_sel),
                                     levels = as.character(dgrid_sel)), y = Value)) +
  geom_point(aes(colour = Type), size = 1, show.legend = TRUE,
             position = position_dodge(width = 0.3)) +
  geom_point(data = data_time_iter_sum[data_time_iter_sum$Type2 == "Time",],
             aes(x = factor(d, labels = as.character(dgrid_sel),
                            levels = as.character(dgrid_sel)),
                 y = Value, colour = Type), size = 3,
             position = position_dodge(width = 0.3), show.legend = FALSE) +
  geom_line(data = data_time_iter_sum[data_time_iter_sum$Type2 == "Time",],
            aes(y = Value, group = Type, col = Type),
            position = position_dodge(width = 0.3)) +
  scale_color_manual(name = "Parametrisation",
                     values = c("MCD" = "#F8766D", "logM" = "#619CFF")) +
  theme_bw() +
  scale_y_continuous(breaks = NULL, trans = myscale_trans2(),
                     sec.axis = sec_axis(~ . * 1, breaks = label_time,
                                         name = "Time (minutes)")) +
  scale_x_discrete(breaks = dgrid_sel) +
  xlab("") + ylab("")+
  theme(panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12),
        text = element_text(size = 15),
        legend.text=element_text(size = 15),
        panel.grid.major = element_blank(), legend.position = "bottom",
        panel.spacing = unit(0.2, "lines"),
        plot.margin=unit(c(0.5, 0, 1, -0.5), "cm"),
        axis.title.y.right = element_text(margin = margin(t = 0, r = 0, b = 0, l = 5)))

pl_Fit_Iter <- ggplot(data_time_iter[data_time_iter$Type2 == "Iterations",],
                   aes(x = factor(d, labels = as.character(dgrid_sel),
                                  levels = as.character(dgrid_sel)), y = Value)) +
  geom_point(aes(colour = Type), size = 1, show.legend = TRUE,
             position = position_dodge(width = 0.3)) +
  geom_point(data = data_time_iter_sum[data_time_iter_sum$Type2 == "Iterations",],
             aes(x = factor(d, labels = as.character(dgrid_sel),
                            levels = as.character(dgrid_sel)),
                 y = Value, colour = Type), size = 3,
             position = position_dodge(width = 0.3), show.legend = FALSE) +
  geom_line(data = data_time_iter_sum[data_time_iter_sum$Type2 == "Iterations",],
            aes(y = Value, group = Type, col = Type),
            position = position_dodge(width = 0.3)) +
  scale_color_manual(name = "Parametrisation",
                     values = c("MCD" = "#F8766D", "logM" = "#619CFF")) +
  theme_bw() +
  scale_y_continuous(breaks = NULL, trans = myscale_trans2(),
                     sec.axis = sec_axis(~ . * 1, breaks=seq(50, 100, by = 10),
                                         name = "Iterations")) +
  scale_x_discrete(breaks = dgrid_sel) +
  xlab("") + ylab("")+
  theme(panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12),
        text = element_text(size = 15),
        legend.text = element_text(size = 15),
        panel.grid.major = element_blank(), legend.position = "bottom",
        panel.spacing = unit(0.2, "lines"),
        plot.margin=unit(c(0.5, 0, 1, -0.5), "cm"),
        axis.title.y.right = element_text(margin = margin(t = 0, r = 0, b = 0, l = 5)))


pl_Fit_Time_Iter <- ggarrange(pl_Fit_Time,
                              pl_Fit_Iter,
                              nrow = 1,
                              common.legend = TRUE,
                              legend = "bottom",
                              widths = c(1.0, 1.0)) +
                              theme(legend.box.spacing = unit(40, "pt"))

pl_Fit_Time_Iter <- annotate_figure(pl_Fit_Time_Iter,
                                    bottom = textGrob("Dimension",hjust = 0.7,
                                                      vjust = -4, gp = gpar(cex = 1.3)))

setwd(root_dir)
setwd("content/Section3/Plots")
ggsave("plot_TIME_ITER_sqrtscale_genMCD.eps", pl_Fit_Time_Iter,
       width = 30, height = 15, units = "cm")
ggsave("plot_TIME_ITER_sqrtscale_genMCD.pdf", pl_Fit_Time_Iter,
       width = 30, height = 15, units = "cm")

################################################
################################################
#### SUPPLEMENTARY MATERIAL FOR SECTION 3.3 ####
################################################
################################################

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
data_time_iter_sum <- rbind(mean_TIME[,c(1, 2, 5, 6)], mean_ITER[,c(1, 2, 5, 6)])

data_hline <- data.frame(Type2 = c("Time","Iterations"),  # Create data for horizontal lines
                         hline = c(0, NA))
data_hline$Type2 <- factor(data_hline$Type2, levels=c("Time", "Iterations"))

dgrid_sel <- dgrid[-1] # To select a subset of dgrid
data_time_iter <- data_time_iter[which(data_time_iter$d %in% dgrid_sel),]
data_time_iter_sum <- data_time_iter_sum[which(data_time_iter_sum$d %in% dgrid_sel),]

data_time_iter$Type2 <- factor(data_time_iter$Type2,
                               levels = c("Time", "Iterations"),
                               labels = c("Time", "Iterations"))
data_time_iter_sum$Type2 <- factor(data_time_iter_sum$Type2,
                                   levels = c("Time", "Iterations"),
                                   labels = c("Time","Iterations"))

# Labels: to be set manually
label_time <- c(1, 5, 25, 50, 100, 500)

##################################
# Figure B.1 (square root scale) #
##################################
pl_Fit_Time <- ggplot(data_time_iter[data_time_iter$Type2 == "Time",],
                      aes(x = factor(d, labels = as.character(dgrid_sel),
                                     levels = as.character(dgrid_sel)), y = Value)) +
  geom_point(aes(colour = Type), size = 1, show.legend = TRUE,
             position = position_dodge(width = 0.3)) +
  geom_point(data = data_time_iter_sum[data_time_iter_sum$Type2 == "Time",],
             aes(x = factor(d, labels = as.character(dgrid_sel),
                            levels = as.character(dgrid_sel)),
                 y = Value, colour = Type), size = 3,
             position = position_dodge(width = 0.3), show.legend = FALSE) +
  geom_line(data = data_time_iter_sum[data_time_iter_sum$Type2 == "Time",],
            aes(y = Value, group = Type, col = Type),
            position = position_dodge(width = 0.3)) +
  scale_color_manual(name = "Parametrisation",
                     values = c("MCD" = "#F8766D", "logM" = "#619CFF")) +
  theme_bw() +
  scale_y_continuous(breaks = NULL, trans = myscale_trans2(),
                     sec.axis = sec_axis(~ . * 1,
                                         breaks = label_time, name = "Time (minutes)")) +
  scale_x_discrete(breaks = dgrid_sel) +
  xlab("") + ylab("")+
  theme(panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12),
        text = element_text(size = 15),
        legend.text=element_text(size=15),
        panel.grid.major = element_blank(),
        legend.position = "bottom",
        panel.spacing = unit(0.2, "lines"),
        plot.margin=unit(c(0.5,0,1,-0.5), "cm"),
        axis.title.y.right = element_text(margin = margin(t = 0, r = 0, b = 0, l = 5)))

pl_Fit_Iter <- ggplot(data_time_iter[data_time_iter$Type2 == "Iterations",],
                      aes(x = factor(d, labels = as.character(dgrid_sel),
                                     levels = as.character(dgrid_sel)), y = Value)) +
  geom_point(aes(colour = Type), size = 1, show.legend = TRUE,
             position = position_dodge(width = 0.3)) +
  geom_point(data = data_time_iter_sum[data_time_iter_sum$Type2 == "Iterations",],
             aes(x = factor(d, labels = as.character(dgrid_sel),
                            levels = as.character(dgrid_sel)),
                 y = Value, colour = Type), size = 3,
             position = position_dodge(width = 0.3), show.legend = FALSE) +
  geom_line(data = data_time_iter_sum[data_time_iter_sum$Type2 == "Iterations",],
            aes(y = Value, group = Type, col = Type),
            position = position_dodge(width = 0.3)) +
  scale_color_manual(name = "Parametrisation",
                     values = c("MCD" = "#F8766D", "logM" = "#619CFF")) +
  theme_bw() +
  scale_y_continuous(breaks = NULL, trans = myscale_trans2(),
                     sec.axis = sec_axis(~ . * 1, breaks=seq(50, 100, by = 10),
                                         name = "Iterations")) +
  scale_x_discrete(breaks = dgrid_sel) +
  xlab("") + ylab("")+
  theme(panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12),
        text = element_text(size = 15),
        legend.text=element_text(size=15),
        panel.grid.major = element_blank(),
        legend.position = "bottom",
        panel.spacing = unit(0.2, "lines"),
        plot.margin=unit(c(0.5, 0, 1, -0.5), "cm"),
        axis.title.y.right = element_text(margin = margin(t = 0, r = 0, b = 0, l = 5)))

pl_Fit_Time_Iter <- ggarrange(pl_Fit_Time,
                              pl_Fit_Iter,
                              nrow = 1,
                              common.legend = TRUE,
                              legend = "bottom",
                              widths = c(1.0, 1.0)) +
  theme(legend.box.spacing = unit(40, "pt"))

pl_Fit_Time_Iter <- annotate_figure(pl_Fit_Time_Iter,
                                    bottom = textGrob("Dimension",
                                                      hjust = 0.7, vjust = -4,
                                                      gp = gpar(cex = 1.3)))

setwd(root_dir)
setwd("content/SupplementaryMaterial/Plots")
ggsave("plot_TIME_ITER_sqrtscale_genLOGM.eps", pl_Fit_Time_Iter,
       width = 30, height = 15, units = "cm")
ggsave("plot_TIME_ITER_sqrtscale_genLOGM.pdf", pl_Fit_Time_Iter,
       width = 30, height = 15, units = "cm")

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

# Breaks labels x and y axes
breaks_seq_MCDgen_x <- as.numeric(tapply(logS[logS$Type2 == "MCD - Generation",]$logS_gen,
                                         logS$d[logS$Type2 == "MCD - Generation"], mean))
breaks_seq_MCDgen_y <- as.numeric(tapply(logS[logS$Type2 == "MCD - Generation",]$logS_nogen,
                                         logS$d[logS$Type2 == "MCD - Generation"], mean))

breaks_seq_logMgen_x <- as.numeric(tapply(logS[logS$Type2 == "logM - Generation",]$logS_nogen,
                                          logS$d[logS$Type2 == "logM - Generation"], mean))
breaks_seq_logMgen_y <- as.numeric(tapply(logS[logS$Type2 == "logM - Generation",]$logS_gen,
                                          logS$d[logS$Type2 == "logM - Generation"], mean))

############################################
# Figure B.2 (original scale in thousands) #
############################################

# Plot under MCD generation (Note that logScores are represented in the original scale)
pl_MCD_gen <- ggplot(logS[logS$Type2 == "MCD - Generation",],
                     aes(x = logS_gen, y = logS_nogen, color = d)) +
  geom_point(size = 2, show.legend = TRUE)+
  geom_abline(slope = 1, intercept = 0, col = "black") +
  scale_color_manual(name = "Dimension", values = col_values)+
  theme_bw() +
  facet_grid(. ~ "MCD generation") +
  xlab("LS - MCD fit") + ylab("LS - logM fit") +
  scale_y_continuous(breaks = breaks_seq_MCDgen_y, labels = scaleFUN2,
                     sec.axis = sec_axis(~ . * 1,
                                         labels = scaleFUN2, breaks = NULL)) +  # breaks_seq_MCDgen_y
  scale_x_continuous(breaks = breaks_seq_MCDgen_x, labels = scaleFUN2,
                     sec.axis = sec_axis(~ . * 1,
                                         labels = scaleFUN2, breaks = NULL))+  # breaks_seq_MCDgen_x
  theme(panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12),
        text = element_text(size = 15),
        legend.text = element_text(size = 15),
        strip.text.x = element_text(size = 15),
        panel.grid.major = element_blank(),
        legend.position = "bottom",
        panel.spacing = unit(0.1, "lines"),
        strip.background=element_rect(colour="black", fill="white"))


# Plot under logM generation (Note that logScores are represented in the original scale)
pl_logM_gen <- ggplot(logS[logS$Type2 == "logM - Generation",],
                      aes(x = logS_nogen, y = logS_gen, color = d)) +
  geom_point(size = 2, show.legend = TRUE)+
  geom_abline(slope = 1, intercept = 0, col = "black") +
  facet_grid(. ~ "logM generation") +
  xlab("LS - MCD fit") + ylab("LS - logM fit") +
  scale_color_manual(name = "Dimension", values = col_values)+
  theme_bw() +
  scale_y_continuous(breaks = breaks_seq_logMgen_y, labels = scaleFUN2,
                     sec.axis = sec_axis(~ . * 1,
                                        labels = scaleFUN, breaks = NULL),) + #breaks_seq_logMgen_y
  scale_x_continuous(breaks=breaks_seq_logMgen_x, labels = scaleFUN2,
                     sec.axis = sec_axis(~ . * 1,
                                         labels = scaleFUN, breaks = NULL),)+ #breaks_seq_logMgen_x
  theme(panel.grid.minor = element_blank(),
        axis.text = element_text(size = 12),
        text = element_text(size = 15),
        legend.text = element_text(size = 15),
        strip.text.x = element_text(size = 15),
        panel.grid.major = element_blank(),
        legend.position = "bottom",
        panel.spacing = unit(0.1, "lines"),
        strip.background=element_rect(colour="black", fill="white"))

# Merging the plots above
pl_logS1 <- ggarrange(pl_MCD_gen,
                      pl_logM_gen,
                      nrow = 1,
                      common.legend = TRUE,
                      legend = "bottom")

setwd(root_dir)
setwd("content/SupplementaryMaterial/Plots")
ggsave("plot_logScore1_new.eps", pl_logS1,
       width = 30, height = 15, units = "cm")
ggsave("plot_logScore1_new.pdf", pl_logS1,
       width = 30, height = 15, units = "cm")

#####################
#####################
##### SECTION 4 #####
#####################
#####################

##################################################################
##################################################################
# Code for reproducing Figure 3: Comparison Hessian w.r.t. beta ##
##################################################################
##################################################################
setwd(root_dir)
# The function for extracting and plotting computational times is in Section 3
setwd("content/Section3/Comp_logM_MCD_Hessian_eta")
source("Functions_Plots_Hessian_eta.R")

nobs <- 1000
dgrid <- seq(10, 120, by = 10)
step_length <- 10
pint_type <- c("dm05", "dm1", "dm1_c2", "dm2", "const")

##################
setwd(root_dir)
setwd("content/Section4/Results")
load(file = paste0("TIME_logm_beta_d", min(dgrid),"_",max(dgrid),"_nobs", nobs, "intMeanFALSE.RData"))

# logM - mean modelling
time_Hbeta_logm_noMeanInt <-  lapply(1 : length(pint_type),
                           function(x) time_hessian(obj = TIME_logM_beta_noMeanInt[[x]], param = "logm",
                                                    dgrid = dgrid, nrun = nrun,
                                                    type = c("block", "noblock"),
                                                    type1 = c("PARS", "STD"),
                                                    beta = TRUE))

# Select Scenario 1 (dm05) and Scenario 2 (dm1_c2)
summary_time_hessianB_logm_noMeanInt <- rbind(time_Hbeta_logm_noMeanInt[[1]]$sum_res,
                                              time_Hbeta_logm_noMeanInt[[3]]$sum_res)
row.names(summary_time_hessianB_logm_noMeanInt) <- NULL

all_time_hessianB_logm_noMeanInt <- rbind(time_Hbeta_logm_noMeanInt[[1]]$res, time_Hbeta_logm_noMeanInt[[3]]$res)
row.names(all_time_hessianB_logm_noMeanInt) <- NULL

rel_mean_time_hessianB_logm_noMeanInt <- data.frame(d = rep(dgrid, 2),
                                          rel_time = summary_time_hessianB_logm_noMeanInt[summary_time_hessianB_logm_noMeanInt$Type == "STD",2]/
                                                     summary_time_hessianB_logm_noMeanInt[summary_time_hessianB_logm_noMeanInt$Type == "PARS",2],
                                          param = rep("logm", 2 * length(dgrid)),
                                          Scenario = factor(rep(c("S1", "S2"), each = length(dgrid)),
                                                            labels = c("S1", "S2"), level = c("S1", "S2")),
                                          Ratio = rep("STD/PARS", 2 * length(dgrid)))

rel_all_time_hessianB_logm_noMeanInt <- data.frame(d = rep(dgrid, 2),
                                         rel_time = all_time_hessianB_logm_noMeanInt[all_time_hessianB_logm_noMeanInt$Type == "STD",1]/
                                                    all_time_hessianB_logm_noMeanInt[all_time_hessianB_logm_noMeanInt$Type == "PARS",1],
                                         param = rep("logm",2 * length(dgrid)),
                                         Scenario = factor(rep(c("S1", "S2"), each = nrun * length(dgrid)),
                                                           labels=c("S1", "S2"), level = c("S1", "S2")),
                                         Ratio = rep("STD/PARS", 2 * length(dgrid)))

# To set manually
lab_time_logM <- c(1, 2, 4, 8, 16, 24, 36)

dg_sel <- dgrid #To select a subset of the dgrid elements
rel_mean_time_hessianB_logm_noMeanInt <- rel_mean_time_hessianB_logm_noMeanInt[which(rel_mean_time_hessianB_logm_noMeanInt$d %in% dg_sel),]
rel_all_time_hessianB_logm_noMeanInt <- rel_all_time_hessianB_logm_noMeanInt[which(rel_all_time_hessianB_logm_noMeanInt$d %in% dg_sel),]

##################################################
# Figure 3 - right panel (logM - mean modelling) #
##################################################
pl_Hbeta <- ggplot(rel_all_time_hessianB_logm_noMeanInt,
                   aes(x = as.factor(d), y = rel_time)) +
  geom_point(aes(colour = Scenario), size = 1, show.legend = FALSE,
             position = position_dodge(width = 0.3)) +
  geom_point(data = rel_mean_time_hessianB_logm_noMeanInt,
             aes(x = as.factor(d), y = rel_time, colour = Scenario),
             size = 2, position = position_dodge(width = 0.3), show.legend = FALSE) +
  geom_line(data = rel_mean_time_hessianB_logm_noMeanInt,
            aes(y = rel_time, group = Scenario, col = Scenario),
            position = position_dodge(width = 0.3)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_color_manual(name = "Scenario",
                     values = c("S1" = "#00A9FF", "S2" = "#F8766D")) +
  theme_bw() +
  scale_y_continuous(breaks = NULL,trans = myscale_trans2(),
                     sec.axis = sec_axis(~ . * 1,
                                         breaks = lab_time_logM, name = "Relative time")) +
  scale_x_discrete(breaks = dg_sel) +
  xlab("") + ylab("") +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        axis.text = element_text(size = 12),  text = element_text(size = 15),
        axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 05, b = 0, l = 0)),
        legend.text=element_text(size=15), strip.text.x = element_text(size = 15),
        legend.position = "bottom", panel.spacing = unit(0.2, "lines"),
        legend.box.spacing = unit(10, "pt"),
        plot.margin=unit(c(0.5, 0.15, 1, -0.5), "cm"),
        axis.title.y.right = element_text(margin = margin(t = 0, r = 0, b = 0, l = 5)))


########################################
# Figure 3 - left panel (lpi modelled) #
########################################
nspars <- data.frame(d =  rep(seq(min(dgrid), max(dgrid), by = 0.01), 2),
                     nelno0_rel = c(1/sqrt(seq(min(dgrid), max(dgrid), by = 0.01)),
                                    2/(seq(min(dgrid), max(dgrid), by = 0.01))),
                     nelno0_abs = c(seq(min(dgrid), max(dgrid), by = 0.01) * (seq(min(dgrid), max(dgrid), by = 0.01)+1)/2 * 1/sqrt(seq(min(dgrid), max(dgrid), by = 0.01)),
                                    seq(min(dgrid), max(dgrid), by = 0.01) * (seq(min(dgrid), max(dgrid), by = 0.01)+1)/2 * 2/seq(min(dgrid), max(dgrid), by = 0.01)),
                     Scenario = factor(rep(c("S1", "S2"),
                                           each = (length(seq(min(dgrid), max(dgrid), by = 0.01)))),
                                       labels=c("S1", "S2"), level = c("S1", "S2")))
add_point <-  data.frame(d =  rep(seq(min(dgrid), max(dgrid), by = step_length), 2),
                         nelno0_rel = c(1/sqrt(seq(min(dgrid), max(dgrid), by = step_length)),
                                        2/(seq(min(dgrid), max(dgrid), by = step_length))),
                         nelno0_abs = c(seq(min(dgrid), max(dgrid), by = step_length) * (seq(min(dgrid), max(dgrid), by = step_length)+1)/2 * 1/sqrt(seq(min(dgrid), max(dgrid), by = step_length)),
                                        seq(min(dgrid), max(dgrid), by = step_length) * (seq(min(dgrid), max(dgrid), by = step_length)+1)/2 * 2/seq(min(dgrid), max(dgrid), by = step_length)),
                         Scenario = factor(rep(c("S1", "S2"), each = length(dgrid)),
                                           labels=c("S1", "S2"), level = c("S1", "S2")))

# Setted manually
scale_function <- function(x){
  return (x/(662.74429/ 0.3162278))
}
inv_scale_function <- function(x){
  return (x*(662.74429/ 0.3162278))
}

plot_perc_covmod_lpi <- ggplot(data.frame(nspars), aes(x = d, y = nelno0_rel)) +
  xlab("") + ylab(expression("% of linear predictors modelled")) +
  theme_bw() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(data = add_point , aes(x = d, y = nelno0_rel, colour = Scenario),
             shape = 16, size = 2.5, show.legend = F) +
  geom_point(data = add_point , aes(x = d, y = scale_function(nelno0_abs), colour = Scenario),
             shape = 17, size = 2.5, show.legend = F) +
  geom_line(aes(y = nelno0_rel, group = Scenario, col = Scenario), linewidth = 0.05) +
  geom_line(aes(y = scale_function(nelno0_abs), group = Scenario, col = Scenario),
            size = 0.05, lty = "dashed") +
  scale_color_manual(name = "Scenario",
                     values = c("S1" = "#00A9FF", "S2" = "#F8766D")) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        axis.text = element_text(size = 12),  text = element_text(size = 15),
       axis.title.y = element_text(margin = margin(t = 0, r = 05, b = 0, l = 0)),
       axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)),
        axis.title.y.right = element_text(margin = margin(t = 0, r = 0, b = 0, l = 05)),
        legend.text=element_text(size = 15), strip.text.x = element_text(size = 15),
        legend.position = "bottom", panel.spacing = unit(0.2, "lines"),
        legend.box.spacing = unit(10, "pt"), plot.margin = unit(c(0.5, 0, 1, 0), "cm"))+
  scale_x_continuous(breaks = add_point$d) +
  scale_y_continuous(breaks = seq(0, 0.35, by = 0.05), lim = c(0, 0.35),
                     sec.axis = sec_axis(~ inv_scale_function(.),
                                         breaks = seq(0, 700, 100),
                                         name="# of linear predictors modelled"))

pl_comp_Hbeta <- ggarrange(plot_perc_covmod_lpi,
                           pl_Hbeta,
                           nrow = 1,
                           common.legend = TRUE,
                           legend = "bottom",
                           widths = c(1.245, 1.0)) +
   theme(legend.box.spacing = unit(20, "pt"))

pl_comp_Hbeta <- annotate_figure(pl_comp_Hbeta,
                                 bottom = textGrob("Dimension", hjust = 0.5,
                                                   vjust = -4, gp = gpar(cex = 1.3)))

setwd(root_dir)
setwd("content/Section4/Plots")
ggsave("plot_rel_TIME_hessian_beta_logM_and_Sparsity.eps", pl_comp_Hbeta,
       width = 30, height = 15, units = "cm")
ggsave("plot_rel_TIME_hessian_beta_logM_and_Sparsity.pdf", pl_comp_Hbeta,
       width = 30, height = 15, units = "cm")

################################################
################################################
##### SUPPLEMENTARY MATERIAL FOR SECTION 4 #####
################################################
################################################

#####################################
#####################################
## Code for reproducing Figure B.3 ##
#####################################
#####################################
setwd(root_dir)
setwd("content/Section4/Results")
pint_type <- c("dm05", "dm1", "dm1_c2", "dm2", "const")
step_length <- 10
dgrid <- seq(10, 120,10)
nobs <- 1000

load(file = paste0("TIME_mcd_beta_d", min(dgrid),"_",max(dgrid),"_nobs", nobs, "intMeanTRUE.RData"))

# MCD - mean modelling
time_Hbeta_mcd_MeanInt <-  lapply(1 : length(pint_type),
                           function(x) time_hessian(obj = TIME_MCD_beta_MeanInt[[x]], param = "mcd",
                                                    dgrid = dgrid, nrun = nrun,
                                                    type = c("block", "noblock"),
                                                    type1 = c("PARS", "STD"),
                                                    beta = TRUE))

# Select Scenario 1 (dm05) and Scenario 2 (dm2)
summary_time_hessianB_mcd_MeanInt <- rbind(time_Hbeta_mcd_MeanInt[[1]]$sum_res, time_Hbeta_mcd_MeanInt[[3]]$sum_res)
row.names(summary_time_hessianB_mcd_MeanInt) <- NULL

all_time_hessianB_mcd_MeanInt <- rbind(time_Hbeta_mcd_MeanInt[[1]]$res, time_Hbeta_mcd_MeanInt[[3]]$res)
row.names(all_time_hessianB_mcd_MeanInt) <- NULL

rel_mean_time_hessianB_mcd_MeanInt <- data.frame(d = rep(dgrid, 2),
                                                 rel_time = summary_time_hessianB_mcd_MeanInt[summary_time_hessianB_mcd_MeanInt$Type == "STD",2]/
                                                 summary_time_hessianB_mcd_MeanInt[summary_time_hessianB_mcd_MeanInt$Type == "PARS",2],
                                                 param = rep("mcd", 2 * length(dgrid)),
                                                 Scenario = factor(rep(c("S1", "S2"), each = length(dgrid)),
                                                                   labels = c("S1", "S2"),
                                                                   level = c("S1", "S2")),
                                                 Ratio = rep("STD/PARS", 2 * length(dgrid)))


rel_all_time_hessianB_mcd_MeanInt <- data.frame(d = rep(dgrid, 2),
                                         rel_time = all_time_hessianB_mcd_MeanInt[all_time_hessianB_mcd_MeanInt$Type == "STD",1]/
                                         all_time_hessianB_mcd_MeanInt[all_time_hessianB_mcd_MeanInt$Type == "PARS",1],
                                         param = rep("mcd",2 * length(dgrid)),
                                         Scenario = factor(rep(c("S1", "S2"), each = nrun * length(dgrid)),
                                                           labels=c("S1", "S2"), level = c("S1", "S2")),
                                         Ratio = rep("STD/PARS", 2 * length(dgrid)))

# Setted manually
lab_time_MCD <- c(1,5,10,15, 20, 30, 50)

dg_sel <- dgrid #To select a subset of the dgrid elements
rel_mean_time_hessianB_mcd_MeanInt <- rel_mean_time_hessianB_mcd_MeanInt[which(rel_mean_time_hessianB_mcd_MeanInt$d %in% dg_sel),]
rel_all_time_hessianB_mcd_MeanInt <- rel_all_time_hessianB_mcd_MeanInt[which(rel_all_time_hessianB_mcd_MeanInt$d %in% dg_sel),]

#################################################
# Figure B.3: right side - MCD no mean modelling #
#################################################
pl_Hbeta_MeanInt <- ggplot(rel_all_time_hessianB_mcd_MeanInt,
                           aes(x = as.factor(d), y = rel_time)) +
  geom_point(aes(colour = Scenario), size = 1, show.legend = FALSE,
             position = position_dodge(width = 0.3)) +
  geom_point(data = rel_mean_time_hessianB_mcd_MeanInt,
             aes(x = as.factor(d), y = rel_time, colour = Scenario),
             size = 2, position = position_dodge(width = 0.3), show.legend = FALSE) +
  geom_line(data = rel_mean_time_hessianB_mcd_MeanInt,
            aes(y = rel_time, group = Scenario, col = Scenario),
            position = position_dodge(width = 0.3)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_color_manual(name = "Scenario",
                     values = c("S1" = "#00A9FF", "S2" = "#F8766D")) +
  theme_bw() +
  scale_y_continuous(breaks = NULL,trans = myscale_trans2(),
                     sec.axis = sec_axis(~ . * 1,  breaks = lab_time_MCD, name = "Relative time")) +
  scale_x_discrete(breaks = dg_sel) +
  xlab("") + ylab("") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text = element_text(size = 12),
        text = element_text(size = 15),
        axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 5)),
        legend.text=element_text(size = 15),
        strip.text.x = element_text(size = 15),
        plot.margin=unit(c(0.5, 0.15, 1, -0.5), "cm"),
        legend.position = "bottom", panel.spacing = unit(0.2, "lines"),
        legend.box.spacing = unit(10, "pt"))


##################
dgrid <- seq(10, 120,10)
load(file = paste0("TIME_mcd_beta_d", min(dgrid),"_",max(dgrid),"_nobs", nobs, "intMeanFALSE.RData"))

# MCD - no mean modelling
time_Hbeta_mcd_noMeanInt <-  lapply(1 : length(pint_type),
                                  function(x) time_hessian(obj = TIME_MCD_beta_noMeanInt[[x]], param = "mcd",
                                                           dgrid = dgrid, nrun = nrun,
                                                           type = c("block", "noblock"),
                                                           type1 = c("PARS", "STD"),
                                                           beta = TRUE))

# Select Scenario 1 (dm05) and Scenario 2 (dm1_c2)
summary_time_hessianB_mcd_noMeanInt <- rbind(time_Hbeta_mcd_noMeanInt[[1]]$sum_res, time_Hbeta_mcd_noMeanInt[[3]]$sum_res)
row.names(summary_time_hessianB_mcd_noMeanInt) <- NULL

all_time_hessianB_mcd_noMeanInt <- rbind(time_Hbeta_mcd_noMeanInt[[1]]$res, time_Hbeta_mcd_noMeanInt[[3]]$res)
row.names(all_time_hessianB_mcd_noMeanInt) <- NULL

rel_mean_time_hessianB_mcd_noMeanInt <- data.frame(d = rep(dgrid, 2),
                                                   rel_time = summary_time_hessianB_mcd_noMeanInt[summary_time_hessianB_mcd_noMeanInt$Type == "STD",2]/
                                                   summary_time_hessianB_mcd_noMeanInt[summary_time_hessianB_mcd_noMeanInt$Type == "PARS",2],
                                                   param = rep("mcd", 2 * length(dgrid)),
                                                   Scenario = factor(rep(c("S1", "S2"), each = length(dgrid)),
                                                                   labels = c("S1", "S2"),
                                                                   level = c("S1", "S2")),
                                                   Ratio = rep("STD/PARS", 2 * length(dgrid)))

rel_all_time_hessianB_mcd_noMeanInt <- data.frame(d = rep(dgrid, 2),
                                                  rel_time = all_time_hessianB_mcd_noMeanInt[all_time_hessianB_mcd_noMeanInt$Type == "STD",1]/
                                                  all_time_hessianB_mcd_noMeanInt[all_time_hessianB_mcd_noMeanInt$Type == "PARS",1],
                                                  param = rep("mcd",2 * length(dgrid)),
                                                  Scenario = factor(rep(c("S1", "S2"),
                                                                    each = nrun * length(dgrid)),
                                                                    labels = c("S1", "S2"),
                                                                    level = c("S1", "S2")),
                                                  Ratio = rep("STD/PARS", 2 * length(dgrid)))

lab_time_MCD <- c(0.25, 0.5, 0.75, 1,1.25,1.5)

dg_sel <- dgrid #To select a subset of the dgrid elements
rel_mean_time_hessianB_mcd_noMeanInt <- rel_mean_time_hessianB_mcd_noMeanInt[which(rel_mean_time_hessianB_mcd_noMeanInt$d %in% dg_sel),]
rel_all_time_hessianB_mcd_noMeanInt <- rel_all_time_hessianB_mcd_noMeanInt[which(rel_all_time_hessianB_mcd_noMeanInt$d %in% dg_sel),]

################################################
# Figure B.3: left side - MCD  mean modelling #
################################################
pl_Hbeta_noMeanInt <- ggplot(rel_all_time_hessianB_mcd_noMeanInt,
                             aes(x = as.factor(d), y = rel_time)) +
  geom_point(aes(colour = Scenario), size = 1, show.legend = FALSE,
             position = position_dodge(width = 0.3)) +
  geom_point(data = rel_mean_time_hessianB_mcd_noMeanInt,
             aes(x = as.factor(d), y = rel_time, colour = Scenario),
             size = 2, position = position_dodge(width = 0.3), show.legend = FALSE) +
  geom_line(data = rel_mean_time_hessianB_mcd_noMeanInt,
            aes(y = rel_time, group = Scenario, col = Scenario),
            position = position_dodge(width = 0.3)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_color_manual(name = "Scenario",
                     values = c("S1" = "#00A9FF", "S2" = "#F8766D")) +
  theme_bw() +
  scale_y_continuous(breaks = NULL,trans = myscale_trans2(),
                     sec.axis = sec_axis(~ . * 1,  breaks = lab_time_MCD)) +
  scale_x_discrete(breaks = dg_sel) +
  xlab("") + ylab("") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text = element_text(size = 12),
        text = element_text(size = 15),
        axis.title.y = element_text(margin = margin(t = 0, r = 0, b = 0, l = 04)),
        legend.text = element_text(size=15),
        strip.text.x = element_text(size = 15),
        legend.position = "bottom",
        panel.spacing = unit(0.2, "lines"),
        plot.margin=unit(c(0.5, 0.0, 1, -0.5), "cm"),
        legend.box.spacing = unit(10, "pt"))


pl_comp_Hbeta <- ggarrange(pl_Hbeta_noMeanInt ,
                           pl_Hbeta_MeanInt ,
                           nrow = 1,
                           common.legend = TRUE,
                           legend = "bottom",
                           widths = c(1.0, 1.065)) +
  theme(legend.box.spacing = unit(40, "pt"))

pl_comp_Hbeta <- annotate_figure(pl_comp_Hbeta,
                                 bottom = textGrob("Dimension", hjust = 0.7,
                                                   vjust = -4, gp = gpar(cex = 1.3)))

setwd(root_dir)
setwd("content/SupplementaryMaterial/Plots")
ggsave("plot_rel_TIME_hessian_beta_MCD_and_Sparsity.eps", pl_comp_Hbeta,
       width = 30, height = 15, units = "cm")
ggsave("plot_rel_TIME_hessian_beta_MCD_and_Sparsity.pdf", pl_comp_Hbeta,
       width = 30, height = 15, units = "cm")

#####################################
#####################################
## Code for reproducing Figure B.4 ##
#####################################
#####################################
setwd(root_dir)
setwd("content/Section4/Results")
load(file=paste0("TIME_logm_beta_d", min(dgrid),"_",max(dgrid),"_nobs", nobs, "intMeanTRUE.RData"))

# logM - no mean modelling
time_Hbeta_logm_MeanInt <-  lapply(1 : length(pint_type),
                                   function(x) time_hessian(obj = TIME_logM_beta_MeanInt[[x]], param = "logm",
                                                            dgrid = dgrid, nrun = nrun,
                                                            type = c("block", "noblock"),
                                                            type1 = c("PARS", "STD"),
                                                            beta = TRUE))

# Select Scenario 1 (dm05) and Scenario 2 (dm1_c2)
summary_time_hessianB_logm_MeanInt <- rbind(time_Hbeta_logm_MeanInt[[1]]$sum_res, time_Hbeta_logm_MeanInt[[3]]$sum_res)
row.names(summary_time_hessianB_logm_MeanInt) <- NULL

all_time_hessianB_logm_MeanInt <- rbind(time_Hbeta_logm_MeanInt[[1]]$res, time_Hbeta_logm_MeanInt[[3]]$res)
row.names(all_time_hessianB_logm_MeanInt) <- NULL

rel_mean_time_hessianB_logm_MeanInt <- data.frame(d = rep(dgrid, 2),
                                                  rel_time = summary_time_hessianB_logm_MeanInt[summary_time_hessianB_logm_MeanInt$Type == "STD", 2]/
                                                    summary_time_hessianB_logm_MeanInt[summary_time_hessianB_logm_MeanInt$Type == "PARS", 2],
                                                  param = rep("logm", 2 * length(dgrid)),
                                                  Scenario = factor(rep(c("S1", "S2"), each = length(dgrid)),
                                                                    labels = c("S1", "S2"),
                                                                    level = c("S1", "S2")),
                                                  Ratio = rep("STD/PARS", 2 * length(dgrid)))

rel_all_time_hessianB_logm_MeanInt <- data.frame(d = rep(dgrid, 2),
                                                 rel_time = all_time_hessianB_logm_MeanInt[all_time_hessianB_logm_MeanInt$Type == "STD",1]/
                                                   all_time_hessianB_logm_MeanInt[all_time_hessianB_logm_MeanInt$Type == "PARS",1],
                                                 param = rep("logm",2 * length(dgrid)),
                                                 Scenario = factor(rep(c("S1", "S2"), each = nrun * length(dgrid)),
                                                                   labels=c("S1", "S2"),
                                                                   level = c("S1", "S2")),
                                                 Ratio = rep("STD/PARS", 2 * length(dgrid)))


# Setted manually
lab_time_logM <- c(1, 2, 4, 8, 16, 32, 64)

dg_sel <- dgrid #To select a subset of the dgrid elements
rel_mean_time_hessianB_logm_MeanInt <- rel_mean_time_hessianB_logm_MeanInt[which(rel_mean_time_hessianB_logm_MeanInt$d %in% dg_sel),]
rel_all_time_hessianB_logm_MeanInt <- rel_all_time_hessianB_logm_MeanInt[which(rel_all_time_hessianB_logm_MeanInt$d %in% dg_sel),]

################
## Figure B.4 ##
################
pl_Hbeta <- ggplot(rel_all_time_hessianB_logm_MeanInt,
                   aes(x = as.factor(d), y = rel_time)) +
  geom_point(aes(colour = Scenario), size = 1, show.legend = FALSE,
             position = position_dodge(width = 0.3)) +
  geom_point(data = rel_mean_time_hessianB_logm_MeanInt,
             aes(x = as.factor(d), y = rel_time, colour = Scenario),
             size = 2, position = position_dodge(width = 0.3), show.legend = FALSE) +
  geom_line(data = rel_mean_time_hessianB_logm_MeanInt,
            aes(y = rel_time, group = Scenario, col = Scenario),
            position = position_dodge(width = 0.3)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_color_manual(name = "Scenario",
                     values = c("S1" = "#00A9FF", "S2" = "#F8766D")) +
  theme_bw() +
  scale_y_continuous(breaks = NULL,trans = myscale_trans2(),
                     sec.axis = sec_axis(~ . * 1,  breaks = lab_time_logM,
                                         name = "Relative time")) +
  scale_x_discrete(breaks = dg_sel) +
  xlab("") + ylab("") +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        axis.text = element_text(size = 12),  text = element_text(size = 15),
        axis.title.x = element_text(margin = margin(t = 15, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(margin = margin(t = 0, r = 05, b = 0, l = 0)),
        legend.text = element_text(size = 15), strip.text.x = element_text(size = 15),
        legend.position = "bottom", panel.spacing = unit(0.2, "lines"),
        legend.box.spacing = unit(10, "pt"),
        plot.margin = unit(c(0.5, 0.15, 1, -0.5), "cm"),
        axis.title.y.right = element_text(margin = margin(t = 0, r = 0, b = 0, l = 05)))

setwd(root_dir)
setwd("content/SupplementaryMaterial/Plots")
ggsave("plot_rel_TIME_hessian_beta_logM_NoMeanModelling.eps", pl_Hbeta,
       width = 15, height = 15, units = "cm")
ggsave("plot_rel_TIME_hessian_beta_logM_NoMeanModelling.pdf", pl_Hbeta,
       width = 15, height = 15, units = "cm")



#####################
#####################
##### SECTION 5 #####
#####################
#####################

################################
# Code for reproducing Table 1 #
################################
setwd(root_dir)
setwd("content/Section5") # Note: the results were saved in Section 6 (which is now Section 5)
source("Functions_Plots_Overall_Fit.R")

dgrid <- c(2, 5, 10)
dgrid2 <- c(15, 20)
nrun <- 10
nobs <- 10000

setwd(root_dir)
setwd("content/Section5/Results")
load(paste0("sim_mcd_fit_nrun_", nrun, "_n_", nobs, "_d_", paste0(dgrid, collapse = "_"), ".RData"))
load(paste0("sim_mcd_fit_fs_efs_nrun_", nrun, "_n_", nobs, "_d_", paste0(dgrid2, collapse = "_"), ".RData"))

# logS_extr <- log_Score_test(obj = sim_mcd_fit, nrun, dgrid, nobs,  param = "mcd")
# LAML_extr <- LAML_extraction(obj = sim_mcd_fit, nrun, dgrid, nobs,  param = "mcd")
# LAML_efs <- as.numeric(LAML_extr[[1]])
# LAML_exact_efs <- as.numeric(LAML_extr[[2]])
# LAML_exact_efs_initialised <- as.numeric(LAML_extr[[3]])
# LAML_bfgs <- as.numeric(LAML_extr[[4]])
# LAML_bfgsinit <- as.numeric(LAML_extr[[5]])
#
# LAML_fsVSefs <- data.frame(LAML_efs = rep(LAML_efs, 2),
#                    LAML_other = c(LAML_exact_efs, LAML_exact_efs_initialised),
#                    diff_LAML_other_efs = c(LAML_exact_efs - LAML_efs, LAML_exact_efs_initialised - LAML_efs),
#                    rel_diff_LAML_other_efs = c((LAML_exact_efs - LAML_efs)/abs(LAML_exact_efs),
#                                                (LAML_exact_efs_initialised - LAML_efs)/abs(LAML_exact_efs)),
#                    d =  rep(factor(rep(dgrid, each = nrun ), levels = dgrid),2),
#                    init = c(rep("EFS", length(LAML_efs)), rep("EFS - init", length(LAML_efs))))

# Comparison EFS v s BAMLSS in terms of computational times
TIMES <- fit_time(sim_mcd_fit, param = "mcd", dgrid, nrun)
dgrid_sel <- dgrid # To select a subset of dgrid
data_time <- TIMES$res[which(TIMES$res$d %in% dgrid_sel),]
data_time_sum <- TIMES$sum_res[which(TIMES$sum_res$d %in% dgrid_sel),]

TIMES2 <- fit_time_red(sim_mcd_fit_fs_efs, param = "mcd", dgrid2, nrun)
dgrid_sel2 <- dgrid2 # To select a subset of dgrid
data_time2 <- TIMES2$res[which(TIMES2$res$d %in% dgrid_sel2),]
data_time_sum2 <- TIMES2$sum_res[which(TIMES2$sum_res$d %in% dgrid_sel2),]


TIMES_sum_efs_bamlss <- data.frame(Value = c(data_time_sum$efs, data_time_sum$efsExact, data_time_sum$bamlss, data_time_sum2$efs, data_time_sum2$efsExact),  #data_time_sum$efsExact_initalised + data_time_sum$efs,
                                   relValue = c(rep(1, length(dgrid)), data_time_sum$efsExact/data_time_sum$efs, data_time_sum$bamlss/data_time_sum$efs, rep(1, length(dgrid2)),data_time_sum2$efsExact/data_time_sum2$efs ),
                                   d = c(rep(dgrid, times = 3), rep(dgrid2, times = 2)),
                                   Type2 = c(rep(c("FS", "EFS", "BAMLSS"), each=length(dgrid)),
                                             rep(c("FS", "EFS"), each=length(dgrid2))))

TIME_efs_bamlss <- data.frame(Value = c(data_time$time_efs, data_time$time_exact_efs, data_time$time_bamlss, data_time2$time_efs, data_time2$time_exact_efs), #data_time$time_exact_efs_initialised + data_time$time_efs,
                              relValue = c(data_time$time_efs/data_time$time_efs, data_time$time_exact_efs/data_time$time_efs, data_time$time_bamlss/data_time$time_efs,
                                           data_time2$time_efs/data_time2$time_efs, data_time2$time_exact_efs/data_time2$time_efs),
                              d = c(rep(rep(dgrid, each = nrun), 3), rep(rep(dgrid2, each = nrun), 2)),
                              Type2 = c(rep(c("FS", "EFS", "BAMLSS"), each=nrun * length(dgrid)),
                                        rep(c("FS", "EFS"), each=nrun * length(dgrid2)))) #"EFS - init",

TIME_efs_bamlss$d <- factor(TIME_efs_bamlss$d, levels=as.character(c(dgrid_sel,dgrid_sel2)), labels=as.character(c(dgrid_sel,dgrid_sel2)))
TIME_efs_bamlss$Type2 <- factor(TIME_efs_bamlss$Type2, levels=c("FS", "EFS", "BAMLSS"),
                                labels=c("FS", "EFS", "BAMLSS"))

TIMES_sum_efs_bamlss$d <- factor(TIMES_sum_efs_bamlss$d, levels=as.character(c(dgrid_sel,dgrid_sel2)), labels=as.character(c(dgrid_sel,dgrid_sel2)))
TIMES_sum_efs_bamlss$Type2 <- factor(TIMES_sum_efs_bamlss$Type2, levels = c("FS", "EFS", "BAMLSS"),
                                     labels=c("FS", "EFS", "BAMLSS"))

# Results in Table 1 extracted from
TIMES_sum_efs_bamlss
# Value  relValue  d  Type2
# 1    0.2146118  1.000000  2     FS
# 2    1.2357164  1.000000  5     FS
# 3    6.8243062  1.000000 10     FS
# 4    0.2334069  1.087577  2    EFS
# 5    1.3826317  1.118891  5    EFS
# 6    8.1075522  1.188041 10    EFS
# 7    1.0626839  4.951655  2 BAMLSS
# 8   21.1280083 17.097781  5 BAMLSS
# 9  398.8662320 58.447880 10 BAMLSS
# 10  26.9820986  1.000000 15     FS
# 11 109.0738461  1.000000 20     FS
# 12  40.9162653  1.516423 15    EFS
# 13 159.5054593  1.462362 20    EFS

with(TIME_efs_bamlss, tapply(Value, list(d, Type2), max))
# FS         EFS     BAMLSS
# 2    0.2338115   0.2708599   2.258816
# 5    1.3082293   1.4839669  34.088330
# 10   7.2579039   8.5662272 634.917791
# 15  32.9363537  47.2431301         NA
# 20 155.1233073 199.4457789         NA

with(TIME_efs_bamlss, tapply(relValue, list(d, Type2), max))
# FS      EFS   BAMLSS
# 2   1 1.238178 10.73032
# 5   1 1.161604 28.41921
# 10  1 1.275281 91.03973
# 15  1 1.771970       NA
# 20  1 1.817909       NA


# breaks_seq_time1 <- c(2, 10, 25, 50, 100, 200, 400)
# breaks_seq_time2 <- c(1, 2, 5, 10, 25, 50, 75)
#
# pl_Fit_Time_efs_bamlss <- ggplot(TIME_efs_bamlss[TIME_efs_bamlss$Type2 != "FS",], aes(x = as.factor(d), y = relValue)) +
#   geom_point(aes(colour = Type2), size = 1, show.legend = TRUE, position = position_dodge(width = 0.1)) +
#   geom_point(data = TIMES_sum_efs_bamlss[TIMES_sum_efs_bamlss$Type2 != "FS",], aes(x = as.factor(d), y = relValue, colour = Type2), size = 3, position = position_dodge(width = 0.1), show.legend = FALSE) +
#   facet_grid(. ~ "Fitting times") +
#   geom_line(data = TIMES_sum_efs_bamlss[TIMES_sum_efs_bamlss$Type2 != "FS",], aes(y = relValue, group = Type2, col = Type2), position = position_dodge(width = 0.1)) +
#   scale_color_manual(name = "", values = c("EFS" = "#F8766D", "BAMLSS" = "#0072B2")) + #"FS" = "#E76BF3", "EFS" = "#0072B2",  "BAMLSS" = "#7CAE00","EFS - init" = "#F8766D",
#   theme_bw() +
#   scale_y_continuous(breaks = NULL, trans = myscale_trans2(),
#                      sec.axis = sec_axis(~ . * 1,#labels = scaleFUN,
#                                          breaks = breaks_seq_time2)) +
#   scale_x_discrete(breaks = c(dgrid_sel,dgrid_sel2)) +
#   xlab("Dimension") + ylab("") +
#   theme(panel.grid.minor = element_blank(),
#         axis.text = element_text(size = 12),  text = element_text(size = 15),
#         legend.text=element_text(size=15), strip.text.x = element_text(size = 15),
#         panel.grid.major = element_blank(), legend.position = "bottom", panel.spacing = unit(0.2, "lines"))
# setwd(root_dir)
# setwd("content/Section6/Plots")
# ggsave("pl_Fit_Time_efs_bamlss_new.eps", pl_Fit_Time_efs_bamlss , width = 15, height = 15, units = "cm")
# ggsave("pl_Fit_Time_efs_bamlss_new.pdf", pl_Fit_Time_efs_bamlss , width = 15, height = 15, units = "cm")


# pl_Fit_Time_efs_bfgs <- ggplot(TIME_efs_bfgs, aes(x = as.factor(d), y = Value)) +
#   geom_point(aes(colour = Type2), size = 1, show.legend = TRUE, position = position_dodge(width = 0.1)) +
#   geom_point(data = TIMES_sum_efs_bfgs, aes(x = as.factor(d), y = Value, colour = Type2), size = 3,
#              position = position_dodge(width = 0.1), show.legend = FALSE) +
#   facet_grid(. ~ "Fitting times") +
#   geom_line(data = TIMES_sum_efs_bfgs, aes(y = Value, group = Type2, col = Type2), position = position_dodge(width = 0.1)) +
#   scale_color_manual(name = "", values = c("FS" = "#E76BF3",  "BFGS - init" = "#AE4371", "BFGS" = "#00A9FF", "BAMLSS" = "#7CAE00")) +
#   theme_bw() +
#   scale_y_continuous(breaks = NULL, trans = myscale_trans2(),
#                      sec.axis = sec_axis(~ . * 1,#labels = scaleFUN,
#                                          breaks = breaks_seq_time2)) +
#   scale_x_discrete(breaks = dgrid_sel) +
#   xlab("Dimension") + ylab("") +
#   theme(panel.grid.minor = element_blank(),
#         axis.text = element_text(size = 12),  text = element_text(size = 15),
#         legend.text=element_text(size=15), strip.text.x = element_text(size = 15),
#         panel.grid.major = element_blank(), legend.position = "bottom", panel.spacing = unit(0.2, "lines"))
#


#pl_LAML_Times_efs_bamlss <- ggarrange(pl_LAML_rel_diff_fs_efs,
#                     pl_Fit_Time_efs_bamlss,
#                     nrow = 1,
#                     common.legend = FALSE,
#                     legend = "bottom")

#pl_LAML_Times_efs_bfgs <- ggarrange(pl_LAML_rel_diff_fs_bfgs,# +  rremove("xlab"),
#                                    pl_Fit_Time_efs_bfgs,# + rremove("xlab"),
#                                      nrow = 1,
#                                      common.legend = FALSE,
#                                      legend = "bottom")
#annotate_figure(pl_LAML_Times_efs_bfgs, bottom = textGrob("Common x-axis", gp = gpar(cex = 1.3)))


#setwd(root_dir)
#setwd("content/SupplementaryMaterial/Plots")
#ggsave("plot_LAML_Times_efs_bamlss.eps", pl_LAML_Times_efs_bamlss,
#        width = 30, height = 15, units = "cm")
#ggsave("plot_LAML_Times_efs_bamlss.pdf", pl_LAML_Times_efs_bamlss,
#       width = 30, height = 15, units = "cm")
#
#ggsave("plot_LAML_Times_efs_bfgs.eps", pl_LAML_Times_efs_bfgs,
#       width = 30, height = 15, units = "cm")
#ggsave("plot_LAML_Times_efs_bfgs.pdf", pl_LAML_Times_efs_bfgs,
#        width = 30, height = 15, units = "cm")


##########################################################################################
# Code for reproducing Figure XXX - SUPPLEMENTARY MATERIAL - comparison in terms of logS #
##########################################################################################

##########################################################################################
# Decide whether including the comparison in terms of logScore on the test set in the SM #
# If so, I will improve the graphical aspects                                            #
##########################################################################################
# setwd(root_dir)
# setwd("content/Section6")
# source("Functions_Plots_Overall_Fit.R")
#
# dgrid <- c(2, 5, 10)
# nrun <- 10
# nobs <- 10000
#
# setwd(root_dir)
# setwd("content/Section6/Results")
# load(paste0("sim_mcd_fit_nrun_", nrun, "_n_", nobs, "_d_", paste0(dgrid, collapse = "_"), ".RData"))

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

# logS_mcd_test <- log_Score_test(sim_mcd_fit, nrun, dgrid, nobs,  param = "mcd")
# logS_gen_test <- as.numeric(logS_mcd_test[[1]])
# logS_efs_test <- as.numeric(logS_mcd_test[[2]])
# logS_exact_efs_test <- as.numeric(logS_mcd_test[[3]])
# logS_exact_efs_initialised_test <- as.numeric(logS_mcd_test[[4]])
# logS_bfgs_test <- as.numeric(logS_mcd_test[[5]])
# logS_bfgsinit_test <- as.numeric(logS_mcd_test[[6]])
# logS_bamlss_test <- as.numeric(logS_mcd_test[[7]])
#
# logS_test <- data.frame(logS_efs = logS_efs_test,
#                         logS_exact_efs = logS_exact_efs_test,
#                         logS_exact_efs_initialised = logS_exact_efs_initialised_test,
#                         logS_bfgs = logS_bfgs_test,
#                         logS_bfgsinit = logS_bfgsinit_test,
#                         logS_bamlss = logS_bamlss_test,
#                         diff_logS_Eefs_efs = logS_exact_efs_test - logS_efs_test,
#                         diff_logS_EIefs_efs = logS_exact_efs_initialised_test - logS_efs_test,
#                         diff_logS_bfgs_efs = logS_bfgs_test - logS_efs_test,
#                         diff_logS_bfgsinit_efs = logS_bfgsinit_test - logS_efs_test,
#                         diff_logS_bamlss_efs = logS_bamlss_test - logS_efs_test,
#                         rel_diff_Eefs_efs_gen = (logS_exact_efs_test - logS_efs_test)/abs(logS_gen_test),
#                         rel_diff_EIefs_efs_gen = (logS_exact_efs_initialised_test - logS_efs_test)/abs(logS_gen_test),
#                         rel_diff_bfgs_efs_gen = (logS_bfgs_test - logS_efs_test)/abs(logS_gen_test),
#                         rel_diff_bfgsinit_efs_gen = (logS_bfgsinit_test - logS_efs_test)/abs(logS_gen_test),
#                         rel_diff_bamlss_efs_gen = (logS_bamlss_test - logS_efs_test)/abs(logS_gen_test),
#                         d =  factor(rep(dgrid, each = nrun ), levels = dgrid))

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
# pl_MCD_rel_diff_Eefs_efs_gen <- ggplot(logS_test, aes(x = d, y = rel_diff_Eefs_efs_gen)) +
#   geom_point(size = 1, show.legend = TRUE)+
#   geom_hline(yintercept = 0, col = "black", lty = "dashed") +
#   theme_bw() +
#   facet_grid(. ~ "(LS(Exact FS) - LS(FS)) / |LS(gen)|") +
#   xlab("Dimension") + ylab("") +
#   theme(panel.grid.minor = element_blank(),
#         axis.text = element_text(size = 12),  text = element_text(size = 15),
#         legend.text=element_text(size=15), strip.text.x = element_text(size = 15),
#         panel.grid.major = element_blank(),legend.position="bottom",panel.spacing = unit(0.1, "lines"))
#
# pl_MCD_rel_diff_EIefs_efs_gen <- ggplot(logS_test, aes(x = d, y = rel_diff_EIefs_efs_gen)) +
#   geom_point(size = 1, show.legend = TRUE)+
#   geom_hline(yintercept = 0, col = "black", lty = "dashed") +
#   theme_bw() +
#   facet_grid(. ~ "(LS(Exact FS - init) - LS(FS)) / |LS(gen)|") +
#   xlab("Dimension") + ylab("") +
#   theme(panel.grid.minor = element_blank(),
#         axis.text = element_text(size = 12),  text = element_text(size = 15),
#         legend.text=element_text(size=15), strip.text.x = element_text(size = 15),
#         panel.grid.major = element_blank(),legend.position="bottom",panel.spacing = unit(0.1, "lines"))
#
#
# pl_MCD_rel_diff_bfgs_efs_gen <- ggplot(logS_test, aes(x = d, y = rel_diff_bfgs_efs_gen)) +
#   geom_point(size = 1, show.legend = TRUE)+
#   geom_hline(yintercept = 0, col = "black", lty = "dashed") +
#   theme_bw() +
#   facet_grid(. ~ "(LS(BFGS) - LS(FS)) / |LS(gen)|") +
#   xlab("Dimension") + ylab("") +
#   theme(panel.grid.minor = element_blank(),
#         axis.text = element_text(size = 12),  text = element_text(size = 15),
#         legend.text=element_text(size=15), strip.text.x = element_text(size = 15),
#         panel.grid.major = element_blank(),legend.position="bottom",panel.spacing = unit(0.1, "lines"))
#
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
# pl_MCD_logS_test <- ggarrange(#pl_MCD_rel_diff_bamlss_efs_gen,
#   #pl_MCD_rel_diff_Eefs_efs_gen,
#   #pl_MCD_rel_diff_EIefs_efs_gen,
#   pl_MCD_rel_diff_bfgsinit_efs_gen,
#   pl_MCD_rel_diff_bfgs_efs_gen,
#   nrow = 2,
#   common.legend = FALSE,
#   legend = "bottom")
#
# setwd(root_dir)
# setwd("content/SupplementaryMaterial/Plots")
# ggsave("pl_MCD_logS_test.eps", pl_MCD_logS_test , width = 30, height = 15, units = "cm")
# ggsave("pl_MCD_logS_test.pdf", pl_MCD_logS_test , width = 30, height = 15, units = "cm")


############################
############################
##### SECTION 6 and SM #####
############################
###########################

#################################
# Code for reproducing Figure 4 #
#################################
setwd(root_dir)
setwd("content/Section6")
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
setwd("content/Section6")
load("GEF14_data_residuals.RData")

low_neff_vcov <- 0
upp_neff_vcov <- 150  # Here must be setted according to the maximum value of the explored grid


flag_residuals <- FALSE
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
    plot(Eff, logS, xaxt= "n", xlab = "Number of effects (logM)", main = "Residuals" )
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
  #cv_mcd <- cv_mcd_response
  logScore_mcd <- res_perf(cv_mcd, d, GEF14_data, param, sets)

  param <- "logm"
  load(paste0("Results/cv_res_stepwise_param", param, "_d_", d, "_lstep_",
              grid_length, "_low_thresh_", low_neff_vcov, "_upp_thresh_", upp_neff_vcov,"_response.RData"))
  load( paste0("Results/res_stepwise_param", param, "_d_", d, "_lstep_", grid_length, "_response.RData"))
  ncov_el_logm <- sapply(1:length(res_logm[[1]]), function(x) length(res_logm$foo[[x]]) - d)
  ncov_el_logm <- sort(ncov_el_logm[ncov_el_logm >= low_neff_vcov & ncov_el_logm <= upp_neff_vcov], decreasing = TRUE)
  #cv_logm <- cv_logm_response
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
    plot(Eff, logS, xaxt= "n", xlab = "Number of effects (logM)", main = "Response" )
    axis(1, at = Eff, labels = factor(Eff, levels = Eff))
    points(Eff[which.min(logS)], logS[which.min(logS)], col = "red", pch = 19)
  }
  )

}

# Here we select the number of effects (minimum or elbow point)

setwd(root_dir)
setwd("content/Section6/Plots")


flag_residuals <- FALSE
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
    xlab("Number of effects") + ylab("LS") +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.text.y = element_text(size = 12),
          axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5),
          text = element_text(size = 15),
          legend.text=element_text(size=15),
          strip.text.x = element_text(size = 15),
          plot.margin=unit(c(0.25, 0.175, -0.2, 0.0), "cm"),
          legend.position = "bottom", panel.spacing = unit(0.2, "lines"))
  print(pl_logS_mcd_logM)
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
    xlab("Number of effects") + ylab("LS") +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.text.y = element_text(size = 12),
          axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5),
          text = element_text(size = 15),
          legend.text=element_text(size=15),
          strip.text.x = element_text(size = 15),
          plot.margin=unit(c(0.25, 0.175, -0.2, 0.0), "cm"),
          legend.position = "bottom", panel.spacing = unit(0.2, "lines"))

  print(pl_logS_mcd_logM)
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
  setwd("content/Section6/Results")

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
  setwd("content/Section6/Plots")

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
  setwd("content/Section6")

} else {
  setwd(root_dir)
  setwd("content/Section6/Results")


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
  setwd("content/Section6/Plots")

  pl_MCD_logM <- get_plots2(obj1_mcd = res_mcd,
                            obj2_logm = res_logm,
                            name_eff = name_eff,
                            d = d,
                            grid_length = grid_length,
                            neff1_mcd = neff_mcd,
                            neff2_logm =  neff_logm
  )
  print(pl_MCD_logM)

  ggsave(paste0("Covmod_MCDandlogM_response.pdf"),  plot = pl_MCD_logM, width = 20, height = 20, units = "cm")
  ggsave(paste0("Covmod_MCDandlogM_response.eps"),  plot = pl_MCD_logM, width = 20, height = 20, units = "cm")
  setwd(root_dir)
  setwd("content/Section6")

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
  setwd("content/Section6/Results")
  load( paste0("res_stepwise_param", param, "_d_", d, "_lstep_", grid_length, "_residuals.RData"))


  time_mcd <- unlist(res_mcd$time_fit)/(1e9 * 60)

  param <- "logm"
  load( paste0("res_stepwise_param", param, "_d_", d, "_lstep_", grid_length, "_residuals.RData"))
  time_logm <- unlist(res_logm$time_fit)/(1e9 * 60)

  grid_d <- seq( 0, d*(d+1)/2, by = grid_length)
  grid_d2 <- seq( 0, d*(d+1)/2, by = 2*grid_length)
  data_time <- data.frame("Time" = c(time_mcd[length(time_mcd):1], time_logm[length(time_logm):1]),
                          "Param" = c(rep("MCD", length(grid_d)), rep("logM", length(grid_d))),
                          "Eff" = rep(grid_d,2))

  setwd(root_dir)
  setwd("content/SupplementaryMaterial/Plots")


  pl_TIME_mcd_logM <- ggplot(data_time,
                             aes(x = Eff, y = Time)) + #factor(Eff, labels = as.character(Eff), levels = as.character(Eff)), y = logS)) +
    geom_point(aes(colour = Param), size = 1, show.legend = TRUE, position = position_dodge(width = 0.3)) +
    geom_line(aes(y = Time, group = Param, col = Param), position = position_dodge(width = 0.3)) +
    scale_x_continuous(breaks = grid_d2) +
    scale_y_continuous(breaks = c(0,seq(0,200,by=50)), trans = myscale_trans2(), name = "") +
    scale_color_manual(name = "Parametrisation", values = c("MCD" = "#F8766D", "logM" = "#619CFF")) +
    theme_bw() +
    xlab("Number of effects") + ylab("Time (minutes)") +
    theme(panel.grid.minor = element_blank(), axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5),
          axis.text.y = element_text(size = 12), text = element_text(size = 15),
          legend.text=element_text(size=15), strip.text.x = element_text(size = 15),
          panel.grid.major = element_blank(), legend.position = "bottom", panel.spacing = unit(0.2, "lines"))
  print(pl_TIME_mcd_logM)
  ggsave(paste0("Time_MCDandlogM_residuals.pdf"),  plot = pl_TIME_mcd_logM, width = 20, height = 20, units = "cm")
  ggsave(paste0("Time_MCDandlogM_residuals.eps"),  plot = pl_TIME_mcd_logM, width = 20, height = 20, units = "cm")

  setwd(root_dir)
  setwd("content/Section6")

} else {
  param <- "mcd"
  setwd(root_dir)
  setwd("content/Section6/Results")
  load( paste0("res_stepwise_param", param, "_d_", d, "_lstep_", grid_length, "_response.RData"))

  time_mcd <- unlist(res_mcd$time_fit)/(1e9 * 60)

  param <- "logm"
  load( paste0("res_stepwise_param", param, "_d_", d, "_lstep_", grid_length, "_response.RData"))
  time_logm <- unlist(res_logm$time_fit)/(1e9 * 60)

  grid_d <- seq( 0, d*(d+1)/2, by = grid_length)
  grid_d2 <- seq( 0, d*(d+1)/2, by = 2*grid_length)
  data_time <- data.frame("Time" = c(time_mcd[length(time_mcd):1], time_logm[length(time_logm):1]),
                          "Param" = c(rep("MCD", length(grid_d)), rep("logM", length(grid_d))),
                          "Eff" = rep(grid_d, 2))

  setwd(root_dir)
  setwd("content/SupplementaryMaterial/Plots")


  pl_TIME_mcd_logM <- ggplot(data_time,
                             aes(x = Eff, y = Time)) + #factor(Eff, labels = as.character(Eff), levels = as.character(Eff)), y = logS)) +
    geom_point(aes(colour = Param), size = 1, show.legend = TRUE, position = position_dodge(width = 0.3)) +
    geom_line(aes(y = Time, group = Param, col = Param), position = position_dodge(width = 0.3)) +
    scale_x_continuous(breaks = grid_d2) +
    scale_y_continuous(breaks = c(10,50, 100, 200), trans = myscale_trans2()) +
    scale_color_manual(name = "Parametrisation", values = c("MCD" = "#F8766D", "logM" = "#619CFF")) +
    theme_bw() +
    xlab("Number of effects") + ylab("Time (minutes)") +
    theme(panel.grid.minor = element_blank(), axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5),
          axis.text.y = element_text(size = 12), text = element_text(size = 15),
          legend.text=element_text(size=15), strip.text.x = element_text(size = 15),
          panel.grid.major = element_blank(), legend.position = "bottom", panel.spacing = unit(0.2, "lines"))
  print(pl_TIME_mcd_logM)
  ggsave(paste0("Time_MCDandlogM_response_VM.pdf"),  plot = pl_TIME_mcd_logM, width = 20, height = 20, units = "cm")
  ggsave(paste0("Time_MCDandlogM_response_VM.eps"),  plot = pl_TIME_mcd_logM, width = 20, height = 20, units = "cm")
  setwd(root_dir)
  setwd("content/Section6")
}



#################################################
# Plot stdev and correlations of the selected model
#################################################
library(ggnewscale)
library(lubridate)

grid_d <- seq(d*(d+1)/2,  0, by = -grid_length)
grid_length <- 5
########################
# MCD parametrisation  #
########################
param <- "mcd"

# here you select to visualisize "full" or "reduced"
model <- "reduced"
#flag_residuals <- FALSE

fit_mcd_static_reduced_response <- fit_model_AllData(data = GEF14_data, flag_res = FALSE,
                                                     param = "mcd", model_type = "reduced", neff_reduced = 80,
                                                     grid_length = 5, grid_d = grid_d, d = 24)

fit_mcd_static_reduced_residuals <- fit_model_AllData(data = GEF14_data_residuals, flag_res = TRUE,
                                                     param = "mcd", model_type = "reduced", neff_reduced = 65,
                                                     grid_length = 5, grid_d = grid_d, d = 24)

fit_mcd_static_full_response <- fit_model_AllData(data = GEF14_data, flag_res = FALSE,
                                                     param = "mcd", model_type = "full", neff_reduced = NULL,
                                                     grid_length = 5, grid_d = grid_d, d = 24)

fit_mcd_static_full_residuals <- fit_model_AllData(data = GEF14_data_residuals, flag_res = TRUE,
                                                     param = "mcd", model_type = "full", neff_reduced = NULL,
                                                     grid_length = 5, grid_d = grid_d, d = 24)

setwd(root_dir)
setwd("content/Section7/Results")
save(fit_mcd_static_reduced_response, file = "fit_mcd_static_reduced_response.RData")
save(fit_mcd_static_reduced_residuals, file = "fit_mcd_static_reduced_residuals.RData")
save(fit_mcd_static_full_response, file = "fit_mcd_static_full_response.RData" )
save(fit_mcd_static_full_residuals, file = "fit_mcd_static_full_residuals.RData" )

param <- "logm"
fit_logm_static_reduced_response <- fit_model_AllData(data = GEF14_data, flag_res = FALSE,
                                                     param = "logm", model_type = "reduced", neff_reduced = 40,
                                                     grid_length = 5, grid_d = grid_d, d = 24)

fit_logm_static_reduced_residuals <- fit_model_AllData(data = GEF14_data_residuals, flag_res = TRUE,
                                                     param = "logm", model_type = "reduced", neff_reduced = 30,
                                                     grid_length = 5, grid_d = grid_d, d = 24)

fit_logm_static_full_response <- fit_model_AllData(data = GEF14_data, flag_res = FALSE,
                                                   param = "logm", model_type = "full", neff_reduced = NULL,
                                                   grid_length = 5, grid_d = grid_d, d = 24)

fit_logm_static_full_residuals <- fit_model_AllData(data = GEF14_data_residuals, flag_res = TRUE,
                                                   param = "logm", model_type = "full", neff_reduced = NULL,
                                                   grid_length = 5, grid_d = grid_d, d = 24)

save(fit_logm_static_reduced_response, file = "fit_logm_static_reduced_response.RData")
save(fit_logm_static_reduced_residuals, file = "fit_logm_static_reduced_residuals.RData")
save(fit_logm_static_full_response, file = "fit_logm_static_full_response.RData" )
save(fit_logm_static_full_residuals, file = "fit_logm_static_full_residuals.RData" )


plotHM_mcd_static_reduced_response <- plot_Heatmap(data = GEF14_data, flag_res = FALSE,
                                               param = "mcd", model_type = "reduced", neff_reduced = 80,
                                               grid_length = 5, grid_d = grid_d, d = 24)

plotHM_mcd_static_reduced_residuals <- plot_Heatmap(data = GEF14_data_residuals, flag_res = TRUE,
                                                     param = "mcd", model_type = "reduced", neff_reduced = 65,
                                                     grid_length = 5, grid_d = grid_d, d = 24)

plotHM_mcd_static_full_response <- plot_Heatmap(data = GEF14_data, flag_res = FALSE,
                                                   param = "mcd", model_type = "full", neff_reduced = NULL,
                                                   grid_length = 5, grid_d = grid_d, d = 24)

plotHM_mcd_static_full_residuals <- plot_Heatmap(data = GEF14_data_residuals, flag_res = TRUE,
                                                   param = "mcd", model_type = "full", neff_reduced = NULL,
                                                   grid_length = 5, grid_d = grid_d, d = 24)

plotHM_logm_static_reduced_response <- plot_Heatmap(data = GEF14_data, flag_res = FALSE,
                                                     param = "logm", model_type = "reduced", neff_reduced = 40,
                                                     grid_length = 5, grid_d = grid_d, d = 24)

plotHM_logm_static_reduced_residuals <- plot_Heatmap(data = GEF14_data_residuals, flag_res = TRUE,
                                                     param = "logm", model_type = "reduced", neff_reduced = 30,
                                                     grid_length = 5, grid_d = grid_d, d = 24)

plotHM_logm_static_full_response <- plot_Heatmap(data = GEF14_data, flag_res = FALSE,
                                                   param = "logm", model_type = "full", neff_reduced = NULL,
                                                   grid_length = 5, grid_d = grid_d, d = 24)

plotHM_logm_static_full_residuals <- plot_Heatmap(data = GEF14_data_residuals, flag_res = TRUE,
                                                   param = "logm", model_type = "full", neff_reduced = NULL,
                                                   grid_length = 5, grid_d = grid_d, d = 24)





####################################


plotHM_mcd_logM_static_reduced_response <- plot_Heatmap2(data = GEF14_data, flag_res = FALSE,
                                                   model_type = "reduced", neff_reduced = 80,
                                                   grid_length = 5, grid_d = grid_d, d = 24)

plotHM_mcd_logM_static_reduced_residuals <- plot_Heatmap2(data = GEF14_data_residuals, flag_res = TRUE,
                                                    model_type = "reduced", neff_reduced = 65,
                                                    grid_length = 5, grid_d = grid_d, d = 24)

plotHM_mcd_logM_static_full_response <- plot_Heatmap2(data = GEF14_data, flag_res = FALSE,
                                                 model_type = "full", neff_reduced = NULL,
                                                 grid_length = 5, grid_d = grid_d, d = 24)

plotHM_mcd_logM_static_full_residuals <- plot_Heatmap2(data = GEF14_data_residuals, flag_res = TRUE,
                                                 model_type = "full", neff_reduced = NULL,
                                                 grid_length = 5, grid_d = grid_d, d = 24)


##################################### NEW

heatmap_FitCov <- function(PredCov, PredCov2 = NULL, d, range_var = NULL, range_corr = c(0, 1), label = "h",
                           label_xaxis, label_yaxis, col_var, col_cor){

  diagDf <- data.frame(
    var1 = c(paste0("h0", 1:9), paste0("h", 10:d)),
    var2 = c(paste0("h", (d):10), paste0("h0", 9:1)),
    if(is.null(PredCov2)){
      Stdev =  sqrt(diag(PredCov))
    } else {
      Stdev =  sqrt(diag(PredCov )) - sqrt(diag(PredCov2))
    }

  )

  sDf <- data.frame(X1 = rep(NA, d ^ 2), X2 = rep(NA, d ^ 2), Corr = rep(NA, d ^ 2))

  count <- 1
  for(i in 1 : d){
    for(j in 1 : d){
      if(i >= 10 & ((d - j) + 2>= 10)){
        sDf$X1[count] <- paste0(label, i)
        sDf$X2[count] <- paste0(label, (d - j) + 1)
      }
      if(i >= 10 & ((d - j) + 2 <= 10)){
        sDf$X1[count] <- paste0(label, i)
        sDf$X2[count] <- paste0(label, "0", (d - j) + 1)
      }

      if(i < 10 & ((d - j) + 2  >= 10)){
        sDf$X1[count] <- paste0(label, "0", i)
        sDf$X2[count] <- paste0(label, (d - j) + 1)
      }

      if(i < 10 & ((d - j) + 2  <= 10)){
        sDf$X1[count] <- paste0(label, "0", i)
        sDf$X2[count] <- paste0(label, "0", (d - j) + 1)
      }

      if(i == j){
        sDf$Corr[count] <- 0
      }
      if(j < i){
        sDf$Corr[count] <- NA
      }
      if(j > i){
        if(is.null(PredCov2)){
          sDf$Corr[count] <- PredCov[j, i]
        } else {
          sDf$Corr[count] <- PredCov[j, i] - PredCov2[j, i]
        }
      }
      count <- count + 1
    }
  }


  gg1 <- ggplot(sDf, aes(X1, X2)) +
    geom_tile(aes(fill = Corr)) +
    scale_fill_gradientn(colors =  col_cor, limits = range_corr, na.value = "white") +
    new_scale_fill() +
    geom_tile(data = diagDf, aes(var1, var2, fill = Stdev)) +
    scale_fill_gradientn(colors = col_var, limits = range_var) +
    theme(aspect.ratio = 1,
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x=element_text(size=10, angle = 90, vjust = 0.5, hjust=1),
          axis.text.y=element_text(size=10),
          legend.key.width  = unit(1, "lines"),
          legend.key.height = unit(1, "lines"),
          legend.text=element_text(size=12),
          legend.title=element_text(size=15))+
    scale_x_discrete(labels = label_xaxis)+
    scale_y_discrete(labels = label_yaxis)

  return(invisible(gg1))

}

plot_heat_and_traj <- function(mu_sigma, yobs, idx1, idx2, dat, nsim){

  # Extract dates
  date1 <- format(as.Date(paste(dat[idx1, "doy"], dat[idx1, "year"]), format = "%j %Y"), format="%d-%m-%Y")
  date2 <- format(as.Date(paste(dat[idx2, "doy"], dat[idx2, "year"]), format = "%j %Y"), format="%d-%m-%Y")

  # Convert to MW
  d <- ncol(yobs)
  conv_fact <- 0.1
  mu_sigma[ , 1:d] <- mu_sigma[ , 1:d] * conv_fact
  mu_sigma[ , (d+1):(2*d)] <- mu_sigma[ , (d+1):(2*d)] * conv_fact^2
  yobs <- yobs * conv_fact

  # Get covariance matrices
  Sigma_list <-  Sigma_mat(mu_sigma[,-c(1 : d)])

  # Simulate trajetories for each date
  set.seed(5151)
  COV <- Sigma_list[[idx1]]
  sdev <- sqrt(diag(COV))
  COV <- diag(sdev) %*% COV %*% diag(sdev)
  diag(COV) <- sdev^2
  X <- rmvn(nsim, rep(0, 24), COV)
  my_dat <- data.frame("hour" = rep(1:24, nsim), "y" = as.vector(t(X)))

  set.seed(5151)
  COV <- Sigma_list[[idx2]]
  sdev <- sqrt(diag(COV))
  COV <- diag(sdev) %*% COV %*% diag(sdev)
  diag(COV) <- sdev^2
  X2 <- rmvn(nsim, rep(0, 24), COV)

  # Created data frame for ggplot
  my_dat <- rbind(my_dat,
                  data.frame("hour" = rep(1:24, nsim), "y" = as.vector(t(X2))))
  my_dat$day <- rep(c(1, 2), each = 24 * nsim)
  my_dat$ID <- rep(1:(2*nsim), each = 24)
  my_dat$Day <- as.character(my_dat$day)
  my_dat$Day[my_dat$day == "1"] <- date1
  my_dat$Day[my_dat$day == "2"] <- date2
  my_dat$Day <- as.factor(my_dat$Day)

  # Heatmap for first date
  plt <- list()
  col_cor <- rev(colorspace::sequential_hcl(palette = "Blues 3", n = 100))
  col_var <- rev(colorspace::sequential_hcl(palette = "Red", n = 100)[1:90])
  label_xaxis <- c(paste0("h0", 0:9), paste0("h", 10:(d-1)))
  label_yaxis <- c(paste0("h", (d-1):10), paste0("h0", 9:0))
  days <- c("Sunday", "Monday", "Tuesday", "Wednesday", "Thursday", "Friday", "Saturday")
  range_var <- NULL #range(c(sqrt(diag(Sigma_list[[idx1]])), sqrt(diag(Sigma_list[[idx2]]))))
  pl1_MCD <- heatmap_FitCov(Sigma_list[[idx1]], d = d, range_var = range_var, range_corr = c(0, 1), label = "h",
                            label_xaxis = label_xaxis, label_yaxis = label_yaxis,  col_var = col_var, col_cor = col_cor)
  pl1_MCD <- pl1_MCD +
    annotate("text",  x = (d-1) + 0.25, y = (d-1) + 0.25, label = date1, vjust=1, hjust=1, cex = 7) +
    theme(legend.position=c(.9,.55))
  plt[[1]] <- pl1_MCD

  # Heatmap for 2nd date
  pl1_MCD <- heatmap_FitCov(Sigma_list[[idx2]], d = d, range_var = range_var, range_corr = c(0, 1), label = "h",
                            label_xaxis = label_xaxis, label_yaxis = label_yaxis,  col_var = col_var, col_cor = col_cor)
  pl1_MCD <- pl1_MCD + annotate("text",  x = (d-1) + 0.25, y = (d-1) + 0.25, label = date2,
                                vjust=1, hjust=1, cex = 7) +
    theme(legend.position=c(.9,.55))
  plt[[2]] <- pl1_MCD

  # Trajectories plot
  cols <- c("blue", "red")
  names(cols) <- as.character(c(date1, date2))
  plt[[3]] <-  ggplot(my_dat, mapping = aes(y = y, x = hour, group = ID, colour=Day, linetype = Day)) +
    geom_line(alpha = 0.3) +
    scale_color_manual(values=cols) +
    theme_bw() +
    geom_line(data = data.frame(hour = 1:24, y = (yobs[idx1, ] - mu_sigma[idx1, 1:24])),
              aes(y = y, x = hour), inherit.aes = FALSE, colour = "blue", linewidth = 1.5) +
    geom_line(data = data.frame(hour = 1:24, y = (yobs[idx2, ] - mu_sigma[idx2, 1:24])),
              aes(y = y, x = hour), inherit.aes = FALSE, colour = "red", linewidth = 1.5) +
    ylab("Residuals (GW)") + xlab("Hour") +
    theme(legend.position = c(.1,.1)) +
    guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1.5, linetype = 1))) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) #+

  return( plt )

}

library(mvnfast)

setwd(root_dir)
setwd("content/Section7/Results")

load(paste0("fit_AllData_parammcd_d_", 24, "_model_", "reduced", "_response.RData"))
mu_sigma <- predict(res_final, type = "response")

plts <- plot_heat_and_traj(mu_sigma = mu_sigma, yobs = res_final$y, idx1 = 77, idx2 = 238, dat = GEF14_data, nsim = 500)

grid.arrange(grobs = plts, ncol = 3)

