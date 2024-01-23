######################################################
# Code for reproducing the figure of the manuscript: #
# "Scalable Additive Covariance Matrix Modelling"    #
######################################################
library(rstudioapi)
root_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(root_dir)

# Load the needed packages
# (it might be required to do something manually)
source("loadPackages.R")
instload_packages()

nrun <- 10   # Set the number of runs
ncores <- 10 # Set the number of cores

scaleFUN <- function(x) sprintf("%.2f", x) # Number of decimal points set to 2

# logarithmic scale
myscale_trans <- function(){
  trans_new("myscale", function(x) log(x),
            function(x) exp(x), domain = c(0, Inf))
}

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
source("Functions_Plots_Hessian_eta.R")  # Here there is the time_hessian()

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
summary_time_hessian_logm <- summary_time_hessian_logm[order(summary_time_hessian_logm[,5],
                                                             summary_time_hessian_logm[,1],
                                                             decreasing = FALSE),]
row.names(summary_time_hessian_logm) <- NULL

all_time_hessian_logm <- time_hessian_logm$res
all_time_hessian_logm <- all_time_hessian_logm[order(all_time_hessian_logm[,3],
                                                     all_time_hessian_logm[,2],
                                                     decreasing = FALSE),]
row.names(all_time_hessian_logm) <- NULL

################################################################
# Merge logM and MCD computational times (overall and summary) #
################################################################
summary_time_hessian <- rbind(summary_time_hessian_mcd, summary_time_hessian_logm)
summary_time_hessian$param<-factor(summary_time_hessian$param,
                                labels=c("logM", "MCD"), level = c("logm", "mcd"))

all_time_hessian <- rbind(all_time_hessian_mcd, all_time_hessian_logm)
all_time_hessian$param<-factor(all_time_hessian$param,
                               labels = c("logM", "MCD"), level = c("logm", "mcd"))


dg_sel <- dgrid #To select a subset of the dgrid elements
all_time_hessian <- all_time_hessian[which(all_time_hessian$d %in% dg_sel),]
summary_time_hessian <- summary_time_hessian[which(summary_time_hessian$d %in% dg_sel),]

# Labels
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

##########################################
# Plot using square root scale: Figure 1 #
##########################################
pl_Heta <- ggplot(all_time_hessian, aes(x = as.factor(d), y = time)) +
  geom_point(aes(colour = Type), size = 1, show.legend = TRUE, position = position_dodge(width = 0.3)) +
  geom_point(data = summary_time_hessian, aes(x = as.factor(d), y = mean_time, colour = Type),
             size = 2, position = position_dodge(width = 0.3), show.legend = FALSE) +
  geom_line(data = summary_time_hessian, aes(y = mean_time, group = Type, col = Type),
            position = position_dodge(width = 0.3))+
  scale_color_manual(name = "", values = c("AD" = "#00A9FF", "EFF" = "#F8766D")) +
  facet_grid2( . ~ param, scales = "free", switch = "y", independent = "all") +
  theme_bw() +
  scale_y_continuous(breaks = NULL,
                     sec.axis = sec_axis(~ . * 1, labels = scaleFUN, breaks = NULL)) +
  scale_x_discrete(breaks = dg_sel) +
  xlab("Dimension") + ylab("") +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        legend.position = "bottom", panel.spacing = unit(0.2, "lines")) +
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
label_time <- with(data_time_iter_sum,
                   c(min(Value[Type2 == "Time" & Type == "MCD"]), min(Value[Type2 == "Time" & Type == "logM"]),
                     max(Value[Type2 == "Time" & Type == "MCD"]), max(Value[Type2 == "Time" & Type == "logM"])))



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
  scale_color_manual(name = "", values = c("MCD" = "#F8766D", "logM" = "#619CFF")) +
  facet_grid2( . ~ Type2 , scales = "free", switch = "y", independent = "all") +
  theme_bw() +
  scale_y_continuous(breaks = NULL, sec.axis = sec_axis(~ . * 1,labels = scaleFUN, breaks = NULL)) +
  scale_x_discrete(breaks = dgrid_sel) +
  xlab("Dimension") + ylab("") +
  theme(panel.grid.minor = element_blank(), axis.text = element_text(size = 9),
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
ggsave("plot_TIME_ITER_sqrtscale.eps", pl_Fit_Time_Iter, width = 30, height = 15, units = "cm")
ggsave("plot_TIME_ITER_sqrtscale.pdf", pl_Fit_Time_Iter, width = 30, height = 15, units = "cm")

##############################################################################################
# SUPPLEMENTARY MATERIAL: decide whether  including the results under logM generation        #
# !!! I think it's redundant reporting computational times and iterations                    #
# Now it corresponds to Figure 6 of the Supplementary Material                               #
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

label_time <- with(data_time_iter_sum,
                   c(min(Value[Type2 == "Time" & Type == "MCD"]), min(Value[Type2 == "Time" & Type == "logM"]),

                     max(Value[Type2 == "Time" & Type == "MCD"]), max(Value[Type2 == "Time" & Type == "logM"])))

pl_Fit_Time_Iter <- ggplot(data_time_iter,
                           aes(x = factor(d, labels = as.character(dgrid_sel), levels = as.character(dgrid_sel)), y = Value)) +
  geom_point(aes(colour = Type), size = 1, show.legend = TRUE, position = position_dodge(width = 0.3)) +
  geom_point(data = data_time_iter_sum,
             aes(x = factor(d, labels = as.character(dgrid_sel), levels = as.character(dgrid_sel)),
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
    Type2 == "Time" ~ scale_y_continuous(breaks = NULL, trans = myscale_trans2(),
                                         sec.axis = sec_axis(~ . * 1, labels = scaleFUN,
                                                             breaks = label_time)),
    Type2 == "Iterations" ~ scale_y_continuous(breaks = NULL, trans = myscale_trans2(),
                                               sec.axis = sec_axis(~ . * 1,
                                                                   breaks=seq(50, 225, by = 25)))
  ))


setwd(root_dir)
setwd("content/SupplementaryMaterial/Plots")
ggsave("plot_TIME_ITER2_sqrtscale.eps", pl_Fit_Time_Iter, width = 30, height = 15, units = "cm")
ggsave("plot_TIME_ITER2_sqrtscale.pdf", pl_Fit_Time_Iter, width = 30, height = 15, units = "cm")

##############################################################################################
##############################################################################################
## SUPPLEMENTARY MATERIALS: Figure 7 - comparison of performance MCD vs logM using logScore ##
##############################################################################################
##############################################################################################
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
show_col(col_values)

# Breaks label of x and y axes
breaks_seq_MCDgen_x <- as.numeric(tapply(logS[logS$Type2 == "MCD - Generation",]$logS_gen, logS$d[logS$Type2 == "MCD - Generation"], mean))
breaks_seq_MCDgen_y <- as.numeric(tapply(logS[logS$Type2 == "MCD - Generation",]$logS_nogen, logS$d[logS$Type2 == "MCD - Generation"], mean))

# Plot under MCD generation (Note that logScores are represented in the original scale)
pl_MCD_gen <- ggplot(logS[logS$Type2 == "MCD - Generation",], aes(x = logS_gen, y = logS_nogen, color = d)) +
  geom_point(size = 2, show.legend = TRUE)+
  geom_abline(slope = 1, intercept = 0, col = "black") +
  scale_color_manual(name = "d", values = col_values)+
  theme_bw() +
  facet_grid(. ~ "MCD generation") +
  xlab("LS - MCD fit") + ylab("LS - logM fit") +
  scale_y_continuous(breaks = breaks_seq_MCDgen_y,
                     sec.axis = sec_axis(~ . * 1,labels = scaleFUN, breaks = NULL),) +
  scale_x_continuous(breaks = breaks_seq_MCDgen_x,
                     sec.axis = sec_axis(~ . * 1,labels = scaleFUN, breaks = NULL),)+
  theme(panel.grid.minor = element_blank(), axis.text = element_text(size = 9),
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
  scale_color_manual(name = "d", values = col_values)+
  theme_bw() +
  scale_y_continuous(breaks = breaks_seq_logMgen_y,
                     sec.axis = sec_axis(~ . * 1,labels = scaleFUN, breaks = NULL),) +
  scale_x_continuous(breaks=breaks_seq_logMgen_x,
                     sec.axis = sec_axis(~ . * 1,labels = scaleFUN, breaks = NULL),)+
  theme(panel.grid.minor = element_blank(), axis.text = element_text(size = 9),
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
setwd("content/Section3/Comp_logM_MCD_Hessian_eta")
source("Functions_Plots_Hessian_eta.R")

nobs <- 1000
dgrid <- seq(5, 100, by = 5)
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
  scale_y_continuous(breaks = NULL,trans = myscale_trans2(),
                     sec.axis = sec_axis(~ . * 1, labels = scaleFUN, breaks = lab_time_logM)) +
  scale_x_discrete(breaks = dg_sel) +
  xlab("Dimension") + ylab("") +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
        legend.position = "bottom", panel.spacing = unit(0.2, "lines"))

setwd(root_dir)
setwd("content/Section4/Plots")
ggsave("plot_rel_TIME_hessian_beta_logM.eps", pl_Hbeta, width = 15, height = 15, units = "cm")
ggsave("plot_rel_TIME_hessian_beta_logM.pdf", pl_Hbeta, width = 15, height = 15, units = "cm")



########################################################################
# Code for reproducing Figure XXX - SECTION 6 + SUPPLEMENTARY MATERIAL #
########################################################################
setwd(root_dir)
setwd("content/Section6")
source("Functions_Plots_Overall_Fit.R")

dgrid <- c(2, 5, 10)
nrun <- 10
nobs <- 10000

setwd(root_dir)
setwd("content/Section6/Results")
load(paste0("sim_mcd_fit_nrun_", nrun, "_n_", nobs, "_d_", paste0(dgrid, collapse = "_"), ".RData"))

logS_mcd_test <- log_Score_test(sim_mcd_fit, nrun, dgrid, nobs,  param = "mcd")
logS_gen_test <- as.numeric(logS_mcd_test[[1]])
logS_efs_test <- as.numeric(logS_mcd_test[[2]])
logS_bfgs_test <- as.numeric(logS_mcd_test[[3]])
logS_bfgsinit_test <- as.numeric(logS_mcd_test[[4]])
logS_bamlss_test <- as.numeric(logS_mcd_test[[5]])

logS_test <- data.frame(logS_efs = logS_efs_test,
                        logS_bfgs = logS_bfgs_test,
                        logS_bfgsinit = logS_bfgsinit_test,
                        logS_bamlss = logS_bamlss_test,
                        diff_logS_bfgs_efs = logS_bfgs_test - logS_efs_test,
                        diff_logS_bfgsinit_efs = logS_bfgsinit_test - logS_efs_test,
                        diff_logS_bamlss_efs = logS_bamlss_test - logS_efs_test,
                        rel_diff_bfgs_efs_gen = (logS_bfgs_test - logS_efs_test)/abs(logS_gen_test),
                        rel_diff_bfgsinit_efs_gen = (logS_bfgsinit_test - logS_efs_test)/abs(logS_gen_test),
                        rel_diff_bamlss_efs_gen = (logS_bamlss_test - logS_efs_test)/abs(logS_gen_test),
                        d =  factor(rep(dgrid, each = nrun ), levels = dgrid))

LAML_extr <- LAML_extraction(sim_mcd_fit, nrun, dgrid, nobs,  param = "mcd")
LAML_efs <- as.numeric(LAML_extr[[1]])
LAML_bfgs <- as.numeric(LAML_extr[[2]])
LAML_bfgsinit <- as.numeric(LAML_extr[[3]])

LAML <- data.frame(LAML_efs = LAML_efs,
                   LAML_bfgs = LAML_bfgs,
                   LAML_bfgsinit = LAML_bfgsinit,
                   diff_LAML_bfgs_efs = LAML_bfgs - LAML_efs,
                   diff_LAML_bfgsinit_efs = LAML_bfgsinit - LAML_efs,
                   d =  factor(rep(dgrid, each = nrun ), levels = dgrid))

col_values<- hue_pal()(length(dgrid))
show_col(col_values)

breaks_seq_MCD_efs <- aggregate(logS_test$logS_efs, list(logS_test$d), FUN = mean)[,2]
breaks_seq_MCD_bfgs <- aggregate(logS_test$logS_bfgs, list(logS_test$d), FUN = mean)[,2]
breaks_seq_MCD_bfgsinit <- aggregate(logS_test$logS_bfgsinit, list(logS_test$d), FUN = mean)[,2]
breaks_seq_MCD_bamlss <- aggregate(logS_test$logS_bamlss, list(logS_test$d), FUN = mean)[,2]

breaks_seq_LAML_efs <- aggregate(LAML$LAML_efs, list(LAML$d), FUN = mean)[,2]
breaks_seq_LAML_bfgs <- aggregate(LAML$LAML_bfgs, list(LAML$d), FUN = mean)[,2]
breaks_seq_LAML_bfgsinit <- aggregate(LAML$LAML_bfgsinit, list(LAML$d), FUN = mean)[,2]


TIMES <- fit_time(sim_mcd_fit, param = "mcd", dgrid, nrun)
rep(dgrid, times = 2)
dgrid_sel <- dgrid # To select a subset of dgrid
data_time <- TIMES$res[which(TIMES$res$d %in% dgrid_sel),]
data_time_sum <- TIMES$sum_res[which(TIMES$sum_res$d %in% dgrid_sel),]
TIMES_sum <- data.frame(Value = c(data_time_sum$efs, data_time_sum$bamlss),
                        d = rep(dgrid, times = 2),
                        Type2 = rep(c("EFS", "BAMLSS"), each=length(dgrid)))

TIME_efs_bamlss <- data.frame(Value = c(data_time$time_efs,data_time$time_bamlss),
                              d = rep(rep(dgrid, each = nrun),2),
                              Type2 = rep(c("EFS", "BAMLSS"), each=nrun*length(dgrid)))
breaks_seq_time_efs <- data_time_sum$efs
breaks_seq_time_bamlss <- data_time_sum$bamlss

###########################
# Comparisons EFS vs BFGS #
###########################
# y axis: Log-score (EFS); x axis: Log-score(BFGS))
pl_MCD_bfgs_efs <- ggplot(logS_test, aes(x = logS_bfgs, y = logS_efs, color = d)) +
  geom_point(size = 2, show.legend = TRUE)+
  geom_abline(slope = 1, intercept = 0, col = "black") +
  scale_color_manual(name = "d", values = col_values)+
  theme_bw() +
  xlab("LS - BFGS") + ylab("LS - EFS") +
  scale_y_continuous(breaks = breaks_seq_MCD_efs ,
                     sec.axis = sec_axis(~ . * 1,labels=scaleFUN, breaks = NULL),) +
  scale_x_continuous(breaks=breaks_seq_MCD_bfgs,
                     sec.axis = sec_axis(~ . * 1,labels=scaleFUN, breaks = NULL),)+
  theme(panel.grid.minor = element_blank(),axis.text=element_text(size=9),
        panel.grid.major = element_blank(),legend.position="bottom",panel.spacing = unit(0.1, "lines"))

# y axis:  Log-score(BFGS) - Log-score (EFS)
pl_MCD_diff_bfgs_efs <- ggplot(logS_test, aes(x = d, y = diff_logS_bfgs_efs)) +
  geom_point(size = 2, show.legend = TRUE)+
  geom_hline(yintercept = 0, col = "black", lty = "dashed") +
  theme_bw() +
  xlab("d") + ylab("LS(BFGS) - LS(EFS)") +
  theme(panel.grid.minor = element_blank(),axis.text=element_text(size=9),
        panel.grid.major = element_blank(),legend.position="bottom",panel.spacing = unit(0.1, "lines"))

# y axis:  (Log-score(BFGS) - Log-score (EFS))/(|log-Score(Generation)|)
pl_MCD_rel_diff_bfgs_efs_gen <- ggplot(logS_test, aes(x = d, y = rel_diff_bfgs_efs_gen)) +
  geom_point(size = 2, show.legend = TRUE)+
  geom_hline(yintercept = 0, col = "black", lty = "dashed") +
  theme_bw() +
  xlab("d") + ylab("(LS(BFGS) - LS(EFS))/(|LS(gen)|)") +
  theme(panel.grid.minor = element_blank(),axis.text=element_text(size=9),
        panel.grid.major = element_blank(),legend.position="bottom",panel.spacing = unit(0.1, "lines"))


# y axis: LAML (EFS); x axis: LAML(BFGS))
pl_LAML_bfgs_efs <- ggplot(LAML, aes(x = LAML_bfgs, y = LAML_efs, color = d)) +
  geom_point(size = 2, show.legend = TRUE)+
  geom_abline(slope = 1, intercept = 0, col = "black") +
  scale_color_manual(name = "d", values = col_values)+
  theme_bw() +
  xlab("LAML - BFGS") + ylab("LAML - EFS") +
  scale_y_continuous(breaks = breaks_seq_LAML_efs ,
                     sec.axis = sec_axis(~ . * 1,labels=scaleFUN, breaks = NULL),) +
  scale_x_continuous(breaks=breaks_seq_LAML_bfgs,
                     sec.axis = sec_axis(~ . * 1,labels=scaleFUN, breaks = NULL),)+
  theme(panel.grid.minor = element_blank(),axis.text=element_text(size=9),
        panel.grid.major = element_blank(),legend.position="bottom",panel.spacing = unit(0.1, "lines"))

# y axis:  LAML(BFGS) - LAML (EFS)
pl_LAML_diff_bfgs_efs <- ggplot(LAML, aes(x = d, y = diff_LAML_bfgs_efs)) +
  geom_point(size = 2, show.legend = TRUE)+
  geom_hline(yintercept = 0, col = "black", lty = "dashed") +
  theme_bw() +
  xlab("d") + ylab("LAML(BFGS) - LAML(EFS)") +
  theme(panel.grid.minor = element_blank(),axis.text=element_text(size=9),
        panel.grid.major = element_blank(),legend.position="bottom",panel.spacing = unit(0.1, "lines"))

##################################################
# Comparisons EFS vs BFGS (initialised with EFS) #
##################################################
# y axis: Log-score (EFS); x axis: Log-score(BFGS init)
pl_MCD_bfgsinit_efs <- ggplot(logS_test, aes(x = logS_bfgsinit, y = logS_efs, color = d)) +
  geom_point(size = 2, show.legend = TRUE)+
  geom_abline(slope = 1, intercept = 0, col = "black") +
  scale_color_manual(name = "d", values = col_values)+
  theme_bw() +
  xlab("LS - BFGS (init.)") + ylab("LS - EFS") +
  scale_y_continuous(breaks = breaks_seq_MCD_efs,
                     sec.axis = sec_axis(~ . * 1,labels=scaleFUN, breaks = NULL),) +
  scale_x_continuous(breaks=breaks_seq_MCD_bfgsinit,
                     sec.axis = sec_axis(~ . * 1,labels=scaleFUN, breaks = NULL),)+
  theme(panel.grid.minor = element_blank(),axis.text=element_text(size=9),
        panel.grid.major = element_blank(),legend.position="bottom",panel.spacing = unit(0.1, "lines"))

# y axis:  Log-score(BFGS init) - Log-score (EFS)
pl_MCD_diff_bfgsinit_efs <- ggplot(logS_test, aes(x = d, y = diff_logS_bfgsinit_efs)) +
  geom_point(size = 2, show.legend = TRUE)+
  geom_hline(yintercept = 0, col = "black", lty = "dashed") +
  theme_bw() +
  xlab("d") + ylab("LS(init. BFGS) - LS(EFS)") +
  theme(panel.grid.minor = element_blank(),axis.text=element_text(size=9),
        panel.grid.major = element_blank(),legend.position="bottom",panel.spacing = unit(0.1, "lines"))

# y axis:  (Log-score(BFGS init) - Log-score (EFS))/(|log-Score(Generation)|)
pl_MCD_rel_diff_bfgsinit_efs_gen <- ggplot(logS_test, aes(x = d, y = rel_diff_bfgsinit_efs_gen)) +
  geom_point(size = 2, show.legend = TRUE)+
  geom_hline(yintercept = 0, col = "black", lty = "dashed") +
  theme_bw() +
  xlab("d") + ylab("(LS(BFGS init) - LS(EFS))/(|LS(gen)|)") +
  theme(panel.grid.minor = element_blank(),axis.text=element_text(size=9),
        panel.grid.major = element_blank(),legend.position="bottom",panel.spacing = unit(0.1, "lines"))

# y axis: LAML (EFS); x axis: LAML(BFGS init))
pl_LAML_bfgsinit_efs <- ggplot(LAML, aes(x = LAML_bfgsinit, y = LAML_efs, color = d)) +
  geom_point(size = 2, show.legend = TRUE)+
  geom_abline(slope = 1, intercept = 0, col = "black") +
  scale_color_manual(name = "d", values = col_values)+
  theme_bw() +
  xlab("LAML - BFGS (init.)") + ylab("LAML - EFS") +
  scale_y_continuous(breaks = breaks_seq_LAML_efs ,
                     sec.axis = sec_axis(~ . * 1,labels=scaleFUN, breaks = NULL),) +
  scale_x_continuous(breaks=breaks_seq_LAML_bfgsinit,
                     sec.axis = sec_axis(~ . * 1,labels=scaleFUN, breaks = NULL),)+
  theme(panel.grid.minor = element_blank(),axis.text=element_text(size=9),
        panel.grid.major = element_blank(),legend.position="bottom",panel.spacing = unit(0.1, "lines"))

# y axis:  LAML(BFGS init.) - LAML (EFS)
pl_LAML_diff_bfgsinit_efs <- ggplot(LAML, aes(x = d, y = diff_LAML_bfgsinit_efs)) +
  geom_point(size = 2, show.legend = TRUE)+
  geom_hline(yintercept = 0, col = "black", lty = "dashed") +
  theme_bw() +
  xlab("d") + ylab("LAML(BFGS init.) - LAML(EFS)") +
  theme(panel.grid.minor = element_blank(),axis.text=element_text(size=9),
        panel.grid.major = element_blank(),legend.position="bottom",panel.spacing = unit(0.1, "lines"))


#############################
# Comparisons EFS vs BAMLSS #
#############################
# y axis: Log-score (EFS); x axis: Log-score(BAMLSS)
pl_MCD_bamlss_efs <- ggplot(logS_test, aes(x = logS_bamlss, y = logS_efs, color = d)) +
  geom_point(size = 2, show.legend = TRUE)+
  geom_abline(slope = 1, intercept = 0, col = "black") +
  scale_color_manual(name = "d", values = col_values)+
  theme_bw() +
  xlab("LS - BAMLSS") + ylab("LS - EFS") +
  scale_y_continuous(breaks = breaks_seq_MCD_efs,
                     sec.axis = sec_axis(~ . * 1,labels=scaleFUN, breaks = NULL),) +
  scale_x_continuous(breaks=breaks_seq_MCD_bamlss,
                     sec.axis = sec_axis(~ . * 1,labels=scaleFUN, breaks = NULL),)+
  theme(panel.grid.minor = element_blank(),axis.text=element_text(size=9),
        panel.grid.major = element_blank(),legend.position="bottom",panel.spacing = unit(0.1, "lines"))


pl_MCD_diff_bamlss_efs <- ggplot(logS_test, aes(x = d, y = diff_logS_bamlss_efs)) +
  geom_point(size = 2, show.legend = TRUE)+
  geom_hline(yintercept = 0, col = "black", lty = "dashed") +
  theme_bw() +
  xlab("d") + ylab("LS(BAMLSS) - LS(EFS)") +
  theme(panel.grid.minor = element_blank(),axis.text=element_text(size=9),
        panel.grid.major = element_blank(),legend.position="bottom",panel.spacing = unit(0.1, "lines"))


# y axis:  (Log-score(BAMLSS) - Log-score (EFS))/(|log-Score(Generation)|)
pl_MCD_rel_diff_bamlss_efs_gen <- ggplot(logS_test, aes(x = d, y = rel_diff_bamlss_efs_gen)) +
  geom_point(size = 2, show.legend = TRUE)+
  geom_hline(yintercept = 0, col = "black", lty = "dashed") +
  theme_bw() +
  xlab("d") + ylab("(LS(BAMLSS) - LS(EFS))/(|LS(gen)|)") +
  theme(panel.grid.minor = element_blank(),axis.text=element_text(size=9),
        panel.grid.major = element_blank(),legend.position="bottom",panel.spacing = unit(0.1, "lines"))


# Comparison EFS v s BAMLSS in terms of computational times
TIME_efs_bamlss$d <- factor(TIME_efs_bamlss$d, levels=as.character(dgrid_sel), labels=as.character(dgrid_sel))
TIME_efs_bamlss$Type2 <- factor(TIME_efs_bamlss$Type2, levels=c("EFS", "BAMLSS"), labels=c("EFS", "BAMLSS"))
TIMES_sum$d <- factor(TIMES_sum$d, levels=as.character(dgrid_sel), labels=as.character(dgrid_sel))
TIMES_sum$Type2 <- factor(TIMES_sum$Type2, levels=c("EFS", "BAMLSS"), labels=c("EFS", "BAMLSS"))

pl_Fit_Time_efs_bamlss <- ggplot(TIME_efs_bamlss,
                                  aes(x = d, y = Value)) +
  geom_point(aes(colour = Type2), size = 1, show.legend = TRUE, position = position_dodge(width = 0.3)) +
  geom_point(data = TIMES_sum,
             aes(x = d,
                 y = Value, colour = Type2), size = 3,
             position = position_dodge(width = 0.3), show.legend = FALSE) +
  geom_line(data = TIMES_sum, aes(y = Value, group = Type2, col = Type2),
            position = position_dodge(width = 0.3)) +
  scale_color_manual(name = "Type", values = c("EFS" = "#F8766D", "BAMLSS" = "#619CFF")) +
  theme_bw() +
  scale_y_continuous(breaks = NULL, trans = myscale_trans2(),
                     sec.axis = sec_axis(~ . * 1,labels = scaleFUN,
                                         breaks = c(breaks_seq_time_efs,breaks_seq_time_bamlss))) +
  scale_x_discrete(breaks = dgrid_sel) +
  xlab("Dimension") + ylab("Fitting times") +
  theme(panel.grid.minor = element_blank(), axis.text = element_text(size = 9),
        panel.grid.major = element_blank(), legend.position = "bottom", panel.spacing = unit(0.2, "lines"))


# Save the results
setwd(root_dir)
setwd("content/Section6/Plots")
ggsave("plot_logS_comparison_efs_bfgs_test.eps", pl_MCD_bfgs_efs, width = 15, height = 15, units = "cm")
ggsave("plot_logS_comparison_efs_bfgs_test.pdf", pl_MCD_bfgs_efs, width = 15, height = 15, units = "cm")
ggsave("plot_logS_difference_efs_bfgs_test.eps", pl_MCD_diff_bfgs_efs, width = 15, height = 15, units = "cm")
ggsave("plot_logS_difference_efs_bfgs_test.pdf", pl_MCD_diff_bfgs_efs, width = 15, height = 15, units = "cm")
ggsave("plot_logS_rel_difference_bfgs_efs_gen_test.eps", pl_MCD_rel_diff_bfgs_efs_gen, width = 15, height = 15, units = "cm")
ggsave("plot_logS_rel_difference_bfgs_efs_gen_test.pdf", pl_MCD_rel_diff_bfgs_efs_gen, width = 15, height = 15, units = "cm")
ggsave("plot_LAML_efs_bfgs.eps", pl_LAML_bfgs_efs, width = 15, height = 15, units = "cm")
ggsave("plot_LAML_efs_bfgs.pdf", pl_LAML_bfgs_efs, width = 15, height = 15, units = "cm")
ggsave("plot_LAML_difference_bfgs_efs.eps", pl_LAML_diff_bfgs_efs, width = 15, height = 15, units = "cm")
ggsave("plot_LAML_difference_bfgs_efs.pdf", pl_LAML_diff_bfgs_efs, width = 15, height = 15, units = "cm")

ggsave("plot_logS_comparison_efs_bfgsinit_test.eps", pl_MCD_bfgsinit_efs, width = 15, height = 15, units = "cm")
ggsave("plot_logS_comparison_efs_bfgsinit_test.pdf", pl_MCD_bfgsinit_efs, width = 15, height = 15, units = "cm")
ggsave("plot_logS_difference_efs_bfgsinit_test.eps", pl_MCD_diff_bfgsinit_efs, width = 15, height = 15, units = "cm")
ggsave("plot_logS_difference_efs_bfgsinit_test.pdf", pl_MCD_diff_bfgsinit_efs, width = 15, height = 15, units = "cm")
ggsave("plot_logS_rel_difference_efs_bfgsinit_gen_test.eps", pl_MCD_rel_diff_bfgsinit_efs_gen, width = 15, height = 15, units = "cm")
ggsave("plot_logS_rel_difference_efs_bfgsinit_gen_test.pdf", pl_MCD_rel_diff_bfgsinit_efs_gen, width = 15, height = 15, units = "cm")
ggsave("plot_LAML_efs_bfgsinit.eps", pl_LAML_bfgsinit_efs, width = 15, height = 15, units = "cm")
ggsave("plot_LAML_efs_bfgsinit.pdf", pl_LAML_bfgsinit_efs, width = 15, height = 15, units = "cm")
ggsave("plot_LAML_difference_bfgsinit_efs.eps", pl_LAML_diff_bfgsinit_efs, width = 15, height = 15, units = "cm")
ggsave("plot_LAML_difference_bfgsinit_efs.pdf", pl_LAML_diff_bfgsinit_efs, width = 15, height = 15, units = "cm")

ggsave("plot_logS_comparison_efs_bamlss_test.eps", pl_MCD_bamlss_efs, width = 15, height = 15, units = "cm")
ggsave("plot_logS_comparison_efs_bamlss_test.pdf", pl_MCD_bamlss_efs, width = 15, height = 15, units = "cm")
ggsave("plot_logS_difference_efs_bamlss_test.eps", pl_MCD_diff_bamlss_efs, width = 15, height = 15, units = "cm")
ggsave("plot_logS_difference_efs_bamlss_test.pdf", pl_MCD_diff_bamlss_efs, width = 15, height = 15, units = "cm")
ggsave("plot_logS_rel_difference_efs_bamlss_gen_test.eps", pl_MCD_rel_diff_bamlss_efs_gen, width = 15, height = 15, units = "cm")
ggsave("plot_logS_rel_difference_efs_bamlss_gen_test.pdf", pl_MCD_rel_diff_bamlss_efs_gen, width = 15, height = 15, units = "cm")
ggsave("ComputationalTimes_efs_bamlss.eps", pl_Fit_Time_efs_bamlss, width = 15, height = 15, units = "cm")
ggsave("ComputationalTimes_efs_bamlss.pdf", pl_Fit_Time_efs_bamlss, width = 15, height = 15, units = "cm")


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


###############################################
# Plots for showing the model selection steps #
###############################################

#######################
# MCD parametrisation #
#######################
param <- "mcd"

setwd(root_dir)
setwd("content/Section7/Results")
load( paste0("res_stepwise_param", param, "_d_", d, "_lstep_", grid_length, ".RData"))

ncov_el_mcd <- sapply(1: length(res_mcd[[1]]), function(x) length(res_mcd$foo[[x]]) - d)

setwd(root_dir)
setwd("content/Section7/Plots")

pl_list <- get_plots(obj = res_mcd,
                     name_eff = name_eff,
                     d = d,
                     grid_length = grid_length,
                     param = param)

for(j in 1:length(pl_list)){
  ggsave(paste0("Covmod_", param, "param_with", length(res_mcd$foo[[j+1]])-d, "Effects.pdf"),  plot=pl_list[[j]], width = 20, height = 20, units = "cm")
}

########################
# logM parametrisation #
########################
param <- "logm"
setwd(root_dir)
setwd("content/Section7/Results")
load( paste0("res_stepwise_param", param, "_d_", d, "_lstep_", grid_length, ".RData"))

ncov_el_logm <- sapply(1: length(res_logm[[1]]), function(x) length(res_logm$foo[[x]]) - d)

setwd(root_dir)
setwd("content/Section7/Plots")

pl_list <- get_plots(obj = res_logm,
                     name_eff = name_eff,
                     d = d,
                     grid_length = grid_length,
                     param = param)

for(j in 1:length(pl_list)){
  ggsave(paste0("Covmod_", param, "param_with", length(res_logm$foo[[j+1]]), "Effects.pdf"),  plot=pl_list[[j]], width = 20, height = 20, units = "cm")
}


##############################
# Plots for validation steps #
##############################
setwd(root_dir)
setwd("content/Section7")

param <- "mcd"
low_neff_vcov <- 0
upp_neff_vcov <- 150

load(paste0("Results/cv_res_stepwise_param", param, "_d_", d, "_lstep_",
            grid_length, "_low_thresh_", low_neff_vcov, "_upp_thresh_", upp_neff_vcov,   ".RData"))


logScore <- res_perf(cv_mcd, d, GEF14_data, param, sets)

load( paste0("Results/res_stepwise_param", param, "_d_", d, "_lstep_", grid_length, ".RData"))
ncov_el_mcd <- sapply(1:length(res_mcd[[1]]), function(x) length(res_mcd$foo[[x]]) - d)
ncov_el_mcd <- sort(ncov_el_mcd[ncov_el_mcd >= low_neff_vcov & ncov_el_mcd <= upp_neff_vcov], decreasing = TRUE)


setwd(root_dir)
setwd("content/Section7/Plots")
png(paste0("LogScore_param", param ,".png"))
plot(ncov_el_mcd, logScore, xaxt= "n", xlab = "Number of effects (MCD)" )
axis(1, at = ncov_el_mcd, labels = factor(ncov_el_mcd, levels = ncov_el_mcd))
points(ncov_el_mcd[which.min(logScore)], logScore[which.min(logScore)], col = "red",  bg=25)
dev.off()

########################
# logM parametrisation #
########################
setwd(root_dir)
setwd("content/Section7")

param <- "logm"
low_neff_vcov <- 0
upp_neff_vcov <- 150
load(paste0("Results/cv_res_stepwise_param", param, "_d_", d, "_lstep_",
            grid_length, "_low_thresh_", low_neff_vcov, "_upp_thresh_", upp_neff_vcov,   ".RData"))

logScore <- res_perf(cv_logm, d, data, param, sets_eval)

load( paste0("Results/res_stepwise_param", param, "_d_", d, "_lstep_", grid_length, ".RData"))
ncov_el_logm <- sapply(1:length(res_logm[[1]]), function(x) length(res_logm$foo[[x]]) - d)
ncov_el_logm <- sort(ncov_el_logm[ncov_el_logm >= low_neff_vcov & ncov_el_logm <= upp_neff_vcov], decreasing = TRUE)


setwd(root_dir)
setwd("content/Section7/Plots")
png(paste0("LogScore_param", param ,".png"))
plot(ncov_el_logm, logScore, xaxt= "n", xlab = "Number of effects (logM)" )
axis(1, at = ncov_el_logm, labels = factor(ncov_el_logm, levels = ncov_el_logm))
points(ncov_el_logm[which.min(logScore)], logScore[which.min(logScore)], col = "red",  bg=25)
dev.off()


