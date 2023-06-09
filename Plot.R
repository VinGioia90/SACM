######################################################
# Code for reproducing the figure of the manuscript: #
# "Scalable Additive Covariance Matrix Modelling"    #
######################################################
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
nobs <- 1000
dgrid <- seq(5, 50, by = 5)

##################
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("content/Section3/Results")


###############################################
# Code for reproducing Figure 1 - SECTION 3.3 #
###############################################
load(file = paste0("TIME_mcd_D2eta_dgrid_min_", min(dgrid), "_max_", max(dgrid), "nobs", nobs, ".RData"))
load(file = paste0("TIME_logm_D2eta_dgrid_min_", min(dgrid), "_max_", max(dgrid), "nobs", nobs, ".RData"))

##################
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("content/Section3/Comp_logM_MCD_Hessian_eta")
source("Functions_Plots_Hessian_eta.R")

# MCD
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

# logM
time_hessian_logm <- time_hessian(obj = TIME_logM_D2eta, param = "logm",
                                        dgrid = dgrid, nrun = nrun,
                                        type = c("eff", "TMB"),
                                        type1 =  c("EFF", "AD"))

summary_time_hessian_logm <- time_hessian_logm$sum_res
summary_time_hessian_logm <- summary_time_hessian_logm[order(summary_time_hessian_logm[,5],
                                                             summary_time_hessian_logm[,1],
                                                             decreasing = FALSE),]
row.names(summary_time_hessian_logm) <- NULL

all_time_hessian_logm <- time_hessian_logm$res
all_time_hessian_logm <- all_time_hessian_logm[order(all_time_hessian_logm[,3],
                                                     all_time_hessian_logm[,2],
                                                     decreasing=FALSE),]
row.names(all_time_hessian_logm) <- NULL

# Merge logM and MCD
summary_time_hessian <- rbind(summary_time_hessian_mcd, summary_time_hessian_logm)
summary_time_hessian$param<-factor(summary_time_hessian$param,
                                labels=c("logM", "MCD"), level = c("logm", "mcd"))

all_time_hessian <- rbind(all_time_hessian_mcd, all_time_hessian_logm)
all_time_hessian$param<-factor(all_time_hessian$param,
                               labels = c("logM", "MCD"), level = c("logm", "mcd"))

dg_sel <- dgrid #To select a subset of the dgrid elements
all_time_hessian <- all_time_hessian[which(all_time_hessian$d %in% dg_sel),]
summary_time_hessian <- summary_time_hessian[which(summary_time_hessian$d %in% dg_sel),]

lab_time_logM <- with(summary_time_hessian,
                      c(min(mean_time[param == "logM" & Type == "EFF"]),
                        min(mean_time[param == "logM" & Type == "AD"]),
                        max(mean_time[param == "logM" & Type == "EFF"]),
                        max(mean_time[param == "logM" & Type == "AD"])))

lab_time_MCD <- with(summary_time_hessian,
                     c(min(mean_time[param == "MCD" & Type == "EFF"]),
                       min(mean_time[param == "MCD" & Type == "AD"]),
                       max(mean_time[param == "MCD" & Type == "EFF"]),
                       max(mean_time[param == "MCD" & Type == "AD"])))

pl_Heta <- ggplot(all_time_hessian, aes(x = as.factor(d), y = time)) +
  geom_point(aes(colour = Type), size = 1, show.legend = TRUE, position = position_dodge(width = 0.3)) +
  geom_point(data = summary_time_hessian, aes(x = as.factor(d), y = mean_time, colour = Type),
             size = 2, position = position_dodge(width = 0.3), show.legend = FALSE) +
  geom_line(data = summary_time_hessian, aes(y = mean_time, group = Type, col = Type),
            position = position_dodge(width = 0.3))+
  scale_color_manual(name = "", values = c("AD" = "#00A9FF", "EFF" = "#F8766D")) +
  facet_grid2( . ~ param , scales = "free", switch = "y", independent = "all") +
  theme_bw() +
  scale_y_continuous(breaks = NULL,
                     sec.axis = sec_axis(~ . * 1, labels = scaleFUN, breaks = NULL)) +
  scale_x_discrete(breaks = dg_sel) +
  xlab("Dimension") + ylab("") +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        legend.position = "bottom", panel.spacing = unit(0.2, "lines")) +
  ggh4x::facetted_pos_scales(y = list(
    param == "logM" ~ scale_y_continuous(breaks = NULL, trans = myscale_trans(),
                                         sec.axis = sec_axis(~ . * 1, labels = scaleFUN, breaks = lab_time_logM)),
    param == "MCD" ~ scale_y_continuous(breaks = NULL, trans = myscale_trans(),
                                        sec.axis = sec_axis(~ . * 1, labels = scaleFUN, breaks = lab_time_MCD))
  ))


setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("content/Section3/Plots")
ggsave("plot_TIME_hessian_eta.eps", pl_Heta, width = 30, height = 15, units = "cm")
ggsave("plot_TIME_hessian_eta.pdf", pl_Heta, width = 30, height = 15, units = "cm")



###############################################
# Code for reproducing Figure 2 - SECTION 3.3 #
###############################################
dgrid <- c(2,5,10,15,20)
nrun <-  10
nobs <- 10000
sg <- FALSE

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("content/Section3/Comp_logM_MCD_Fit")
source("Functions_Plots_Overall_Fit.R")


# only MCD - generation in the paper (similar results by generating with the logM)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("content/Section3/Results")

load(paste0("sim_mcdG_mcdF_nrun_", nrun, "_n_", nobs, "_d_", paste0(dgrid, collapse="_"), ".RData"))
load(paste0("sim_mcdG_logmF_nrun_", nrun, "_n_", nobs, "_d_", paste0(dgrid, collapse = "_"), ".RData"))


res_mcdG_mcdF <- fit_time(sim_mcdG_mcdF, param = "mcd", dgrid, nrun = nrun)
res_mcdG_logmF <- fit_time(sim_mcdG_logmF, param = "logm", dgrid, nrun = nrun)
rm("sim_mcdG_mcdF", "sim_mcdG_logmF")

mean_TIME <- rbind(res_mcdG_mcdF$sum_res, res_mcdG_logmF$sum_res)
all_TIME <- rbind(res_mcdG_mcdF$res, res_mcdG_logmF$res)
mean_ITER <- rbind(res_mcdG_mcdF$sum_iter, res_mcdG_logmF$sum_iter)
all_ITER <- rbind(res_mcdG_mcdF$iter, res_mcdG_logmF$iter)

colnames(mean_TIME)[2] <- colnames(all_TIME)[1] <- colnames(mean_ITER)[2] <- colnames(all_ITER)[1] <- "Value"

mean_TIME  <- cbind(mean_TIME , c(rep("Time", 2 * length(dgrid))))
colnames(mean_TIME)[6] <- "Type2"

mean_ITER  <- cbind(mean_ITER , c(rep("Iterations", 2 * length(dgrid))))
colnames(mean_ITER)[6] <- "Type2"

all_TIME  <- cbind(all_TIME , c(rep("Time", 2 * length(dgrid) * nrun)))
colnames(all_TIME)[4] <- "Type2"

all_ITER  <- cbind(all_ITER , c(rep("Iterations", 2 * length(dgrid) * nrun)))
colnames(all_ITER)[4] <- "Type2"


# Join the time and the iterations (overall)
data_time_iter <- rbind(all_TIME, all_ITER)
data_time_iter_sum <-  rbind(mean_TIME[,c(1,2,5,6)], mean_ITER[,c(1,2,5,6)])


data_hline <- data.frame(Type2 = c("Time","Iterations"),  # Create data for lines
                         hline = c(0, NA))
data_hline$Type2 <- factor(data_hline$Type2, levels=c("Time", "Iterations"))


dgrid_sel <- dgrid # To select a subset f dgrid
data_time_iter <- data_time_iter[which(data_time_iter$d %in% dgrid_sel),]
data_time_iter_sum <- data_time_iter_sum[which(data_time_iter_sum$d %in% dgrid_sel),]



label_time <- with(data_time_iter_sum,
                   c(min(Value[Type2 == "Time" & Type == "MCD"]), min(Value[Type2 == "Time" & Type == "logM"]),
                     max(Value[Type2 == "Time" & Type == "MCD"]), max(Value[Type2 == "Time" & Type == "logM"])))

pl_Fit_Time_Iter <- ggplot(data_time_iter,
                           aes(x = factor(d, labels = as.character(dgrid_sel), levels = as.character(dgrid_sel)), y = Value)) +
  geom_point(aes(colour = Type), size = 1, show.legend = TRUE, position = position_dodge(width = 0.3)) +
  geom_point(data = data_time_iter_sum,
             aes(x = factor(d, labels = as.character(dgrid_sel),levels = as.character(dgrid_sel)),
                 y = Value, colour = Type), size = 3,
             position = position_dodge(width = 0.3), show.legend = FALSE) +
  geom_line(data = data_time_iter_sum, aes(y = Value, group = Type, col = Type),
            position = position_dodge(width = 0.3)) +
  scale_color_manual(name = "", values = c("MCD" = "#F8766D", "logM" = "#619CFF")) +
  facet_grid2( . ~ Type2 , scales = "free", switch = "y", independent = "all") +
  theme_bw() +
  scale_y_continuous(breaks = NULL, sec.axis = sec_axis(~ . * 1,labels = scaleFUN, breaks = NULL)) +
  scale_x_discrete(breaks = dgrid_sel) +
  xlab("Dimension") + ylab("") +
  theme(panel.grid.minor = element_blank(), axis.text = element_text(size = 9),
        panel.grid.major = element_blank(), legend.position = "bottom", panel.spacing = unit(0.2, "lines"))+
  ggh4x::facetted_pos_scales(y = list(
    Type2 == "Time" ~ scale_y_continuous(breaks = NULL, trans = myscale_trans(),
                                         sec.axis = sec_axis(~ . * 1, labels = scaleFUN,
                                                             breaks = label_time)),
    Type2 == "Iterations" ~ scale_y_continuous(breaks = NULL, trans = myscale_trans(),
                                               sec.axis = sec_axis(~ . * 1,
                                                                   breaks=seq(50, 225, by = 25)))
  ))


setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("content/Section3/Plots")
ggsave("plot_TIME_ITER.eps", pl_Fit_Time_Iter, width = 30, height = 15, units = "cm")


###############################################################
# Code for reproducing the Figure 5 - SUPPLEMENTARY MATERIALS #
###############################################################
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("content/Section3/Comp_logM_MCD_Fit")
source("Functions_Plots_Overall_Fit.R")

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("content/Section3/Results")

load(paste0("sim_mcdG_mcdF_nrun_", nrun, "_n_", nobs, "_d_", paste0(dgrid, collapse="_"), ".RData"))
load(paste0("sim_mcdG_logmF_nrun_", nrun, "_n_", nobs, "_d_", paste0(dgrid, collapse = "_"), ".RData"))

logS_mcdG_mcdF <- log_Score(obj=sim_mcdG_mcdF, nrun, dgrid, nobs,  param = "mcd", ncov = 3)
logS_mcdG_logmF <- log_Score(sim_mcdG_logmF, nrun, dgrid, nobs,  param = "logm", ncov = 3)

rm("sim_mcdG_mcdF", "sim_mcdG_logmF")

# To do: to convert into  ggplot

dgrid_sel <- dgrid #To select a subset of the grid

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("content/Section3/Plots")
pdf(file="logScore_MCDgen.pdf", width=20, height=15)
par(mfrow=c(1, length(dgrid)), pty="s")
for(j in 1: length(dgrid)){
  plot(logS_mcdG_mcdF[,j],logS_mcdG_logmF[,j], pch=16, lwd=2,
       xlab = "MCD Fit",
       ylab = "logM Fit",
       main = paste0("MCD Gen. d=", dgrid[j]))
  abline(0,1,col="red")
}
dev.off()



setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("content/Section3/Results")
load(paste0("sim_logmG_mcdF_nrun_", nrun, "_n_", nobs, "_d_", paste0(dgrid, collapse="_"), ".RData"))
load(paste0("sim_logmG_logmF_nrun_", nrun, "_n_", nobs, "_d_", paste0(dgrid, collapse = "_"), ".RData"))

logS_logmG_mcdF <- log_Score(sim_logmG_mcdF, nrun, dgrid, nobs,  param = "mcd")
logS_logmG_logmF <- log_Score(sim_logmG_logmF, nrun, dgrid, nobs,  param = "logm")

rm("sim_logmG_mcdF", "sim_logmG_logmF")


setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("content/Section3/Plots")
pdf(file="logScore_logMgen.pdf", width = 20, height = 15)
par(mfrow = c(1, length(dgrid)), pty = "s")
for(j in 1 : length(dgrid)){
  plot(logS_logmG_logmF[,j], logS_logmG_mcdF[,j], pch = 16, lwd = 2,
       xlab = "logM Fit",
       ylab = "MCD Fit",
       main = paste0("logM Gen. d=", dgrid[j]))
  abline(0, 1, col = "red")
}
dev.off()




#############################################
# Code for reproducing Figure 3 - SECTION 4 #
#############################################
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("content/Section3/Comp_logM_MCD_Hessian_eta")
source("Functions_Plots_Hessian_eta.R")

nobs <- 1000
dgrid <- seq(5, 100, by = 5)
pint_type <- c("dm05","dm1", "dm2","const")

##################
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("content/Section4/Results")

# load(file=paste0("TIME_mcd_beta_d", min(dgrid),"_",max(dgrid),"_nobs", nobs, ".RData"))
load(file=paste0("TIME_logM_beta_d", min(dgrid),"_",max(dgrid),"_nobs", nobs, ".RData"))


# MCD: not included in the paper
# time_Hbeta_mcd <-  lapply(1 : length(pint_type),
#                           function(x) time_hessian(obj = TIME_MCD_beta[[x]], param = "mcd",
#                                                    dgrid = dgrid, nrun = nrun,
#                                                    type = c("block", "noblock"),
#                                                    type1 = c("PARS", "STD"),
#                                                    beta = TRUE))
#
#
# summary_time_hessianB_mcd <- rbind(time_Hbeta_mcd[[1]]$sum_res, time_Hbeta_mcd[[3]]$sum_res)
# row.names(summary_time_hessianB_mcd) <- NULL
#
# all_time_hessianB_mcd <- rbind(time_Hbeta_mcd[[1]]$res, time_Hbeta_mcd[[3]]$res)
# row.names(all_time_hessianB_mcd) <- NULL
#
# # Get the ratio of the mean and all the times
# rel_mean_time_hessianB_mcd <- data.frame(d = rep(dgrid, 2),
#                                          rel_time = summary_time_hessianB_mcd[summary_time_hessianB_mcd$Type == "STD",2]/
#                                                     summary_time_hessianB_mcd[summary_time_hessianB_mcd$Type == "PARS",2],
#                                          param = rep("mcd", 2 * length(dgrid)),
#                                          Scenario = factor(rep(c("S1", "S2"), each = length(dgrid)),
#                                                            labels = c("S1", "S2"), level = c("S1", "S2")),
#                                          Ratio = rep("STD/PARS", 2 * length(dgrid)))
# rel_all_time_hessianB_mcd <- data.frame(d = rep(dgrid, 2),
#                                         rel_time = all_time_hessianB_mcd[all_time_hessianB_mcd$Type == "STD",1]/
#                                                    all_time_hessianB_mcd[all_time_hessianB_mcd$Type == "PARS",1],
#                                         param = rep("mcd", 2 * length(dgrid)),
#                                         Scenario = factor(rep(c("S1", "S2"), each = 2 * nrun * length(dgrid)),
#                                                         labels = c("S1", "S2"), level = c("S1", "S2")),
#                                         Ratio=rep("STD/PARS",2 * length(dgrid)))
#
# lab_time_MCD <- c(1,max(rel_mean_time_hessianB_mcd[rel_mean_time_hessianB_mcd$Scenario == "S1",2]),
#                     max(rel_mean_time_hessianB_mcd[rel_mean_time_hessianB_mcd$Scenario == "S2",2]))


# logM
time_Hbeta_logm <-  lapply(1 : length(pint_type),
                           function(x) time_hessian(obj = TIME_logM_beta[[x]], param = "logm",
                                                    dgrid = dgrid, nrun = nrun,
                                                    type = c("block", "noblock"),
                                                    type1 = c("PARS", "STD"),
                                                    beta = TRUE))

# Select Scenario 1 (dm05) and Scenario 2 (dm2)
summary_time_hessianB_logm <- rbind(time_Hbeta_logm[[1]]$sum_res, time_Hbeta_logm[[3]]$sum_res)
row.names(summary_time_hessianB_logm) <- NULL

all_time_hessianB_logm <- rbind(time_Hbeta_logm[[1]]$res, time_Hbeta_logm[[3]]$res)
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

lab_time_logM <- c(1, max(rel_mean_time_hessianB_logm[rel_mean_time_hessianB_logm$Scenario == "S1",2]),
                   max(rel_mean_time_hessianB_logm[rel_mean_time_hessianB_logm$Scenario == "S2",2]))

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
  scale_color_manual(name = "", values = c("S1" = "#00A9FF", "S2" = "#F8766D")) +
  theme_bw() +
  scale_y_continuous(breaks = NULL,trans = myscale_trans(),
                     sec.axis = sec_axis(~ . * 1, labels = scaleFUN, breaks = lab_time_logM)) +
  scale_x_discrete(breaks = dg_sel) +
  xlab("Dimension") + ylab("") +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        legend.position = "bottom", panel.spacing = unit(0.2, "lines"))

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
setwd("content/Section4/Plots")
ggsave("plot_rel_TIME_hessian_beta_logM.eps", pl_Hbeta, width = 15, height = 15, units = "cm")
ggsave("plot_rel_TIME_hessian_beta_logM.pdf", pl_Hbeta, width = 15, height = 15, units = "cm")

