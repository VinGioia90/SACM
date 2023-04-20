#######################
# Some useful package #
#######################
library(rstudioapi)
library(Gmisc, quietly = TRUE)
library(glue)
library(htmlTable)
library(grid)
library(magrittr)

library(ggpubr)
library(gridExtra)
library(ggplot2)
library(lattice)


library(plyr)
library(dplyr)
library(magrittr)
library(tidyr)
library(scales)
library(geojsonio)
library(rgeos)
library(maptools)
library(gdata)
library(ggh4x)




################################
# Visualization of the results #
################################
# This function allows to summarise the computational times
summary_time <- function(obj, param = c("mcd", "logm"), dgrid1 = NULL, nrun){
  if(param == "mcd") param2 <- "MCD" else param2 <- "logM"

  str_type <- "obj[[x]]$gen[[z]]$time"
  data_time <-  data.frame(unlist(lapply(1 : length(dgrid),
                                         function(z) unlist(lapply(1 : nrun, function(x) obj[[x]]$gen[[z]]$time / 1e9)))),
                           rep(dgrid,each = nrun))
  colnames(data_time) <- c("time","d")
  data <- aggregate(data_time$time, list(data_time$d), FUN = mean)
  data_min <- aggregate(data_time$time, list(data_time$d), FUN = min)
  data_max <- aggregate(data_time$time, list(data_time$d), FUN = max)

  out <- cbind(data,data_min[,2], data_max[,2])

  res <- as.data.frame(out)
  res <- cbind(res, rep(param2, each = length(dgrid)))
  colnames(res) <- c("d", "mean_time", "min_time", "max_time", "Type")
  return(data.frame(res))
}

# This function extracts all the computational times
all_time <- function(obj, param = c("mcd", "logm"), dgrid, nrun){
  if(param == "mcd") param2 <- "MCD" else param2 <- "logM"

  str_type <- "obj[[x]]$gen[[z]]$time"
  data_time <-  data.frame(unlist(lapply(1 : length(dgrid),
                                         function(z) unlist(lapply(1 : nrun, function(x) obj[[x]]$gen[[z]]$time/1e9)))),
                           rep(dgrid,each = nrun))
  colnames(data_time) <- c("time","d")

  res <- data_time

  res <- cbind(res, rep(param2, each = length(dgrid)))
  colnames(res) <- c("time", "d", "Type")
  return(data.frame(res))
}

# Function to extract the number of iterations
ITER_extract <- function(obj, dgrid, nrun){
  return(unlist(lapply(1 : length(dgrid),
                       function(z) unlist(lapply(1 : nrun,
                                                 function(x) obj[[x]]$gen[[z]]$inner)))))
}

# Auxiliary functions for the plot
scaleFUN <- function(x) sprintf("%.2f", x)

myscale_trans <- function()
{
  trans_new("myscale", function(x) sqrt(x),
            function(x) x^2, domain = c(0, Inf))
}


##########################################
# Code for reproducibility of Figure 2.7 #
##########################################
dgrid <- c(2,3,5)
nrun <-  7
ncores <- 7
nobs <- 5000

sg <- FALSE


load(paste0("sim_mcdG_mcdF_nrun_", nrun, "_n_", nobs, "_d_", paste0(dgrid, collapse="_"), ".RData"))
load(paste0("sim_mcdG_logmF_nrun_", nrun, "_n_", nobs, "_d_", paste0(dgrid, collapse = "_"), ".RData"))

###############################################################
# Times
mean_TIME_MCDg_MCDf <- summary_time(sim_mcdG_mcdF, param = "mcd", dgrid, nrun = nrun)
TIME_MCDg_MCDf <- all_time(sim_mcdG_mcdF, param = "mcd", dgrid, nrun = nrun)

mean_TIME_MCDg_logMf <- summary_time(sim_mcdG_logmF, param = "logm", dgrid, nrun = nrun)
TIME_MCDg_logMf <- all_time(sim_mcdG_logmF, param = "logm", dgrid, nrun = nrun)

mean_TIME <- rbind(mean_TIME_MCDg_MCDf, mean_TIME_MCDg_logMf)
TIME_overall <- rbind(TIME_MCDg_MCDf, TIME_MCDg_logMf)

data_time_mean <- cbind(mean_TIME, Type1 = c(rep("S1",length(dgrid))))
colnames(data_time_mean)[2] <- "Value"
data_time_mean  <- cbind(data_time_mean , c(rep("Time",6)))
colnames(data_time_mean)[7] <- "Type2"


###############################################################
# Iterations
ITER_MCDg_MCDf <- ITER_extract(sim_mcdG_mcdF, dgrid, nrun)
ITER_MCDg_logMf <- ITER_extract(sim_mcdG_logmF, dgrid, nrun)
ITER_overall <- c(ITER_MCDg_MCDf, ITER_MCDg_logMf)

TYPE <- factor(c(rep("MCD", length(dgrid) * nrun),rep("logM", length(dgrid) * nrun)),
               levels = c("MCD", "logM"))

d <- factor(rep(c(dgrid,dgrid), each = nrun))
Scen <- rep("S1", 2 * length(dgrid) * nrun)
data_iter = data.frame(ITER_overall, d, TYPE, Scen)
colnames(data_iter) <- c("Value","d","Type","Type1")


###################################################################
# Join the time and the iterations (overall)
data_time_iter <- cbind(TIME_overall, Type1 = rep("S1", 2 * length(dgrid) * nrun))
colnames(data_time_iter)[1] <- "Value"

data_time_iter <- rbind(data_time_iter, data_iter)
data_time_iter <- cbind(data_time_iter, c(rep("Time", 2 * length(dgrid)*nrun), rep("Iterations",2 * length(dgrid) * nrun)))
colnames(data_time_iter)[5] <- "Type2"

###################################################################
# Join the time and the iterations (average)
data_iter_mean <- aggregate(data_iter$Value, list(data_iter$d, data_iter$Type), FUN = mean)
data_iter_mean <- cbind(data_iter_mean, rep("S1", 2 * length(dgrid)), rep("Iterations",2*length(dgrid)))
colnames(data_iter_mean) <- c("d", "Type", "Value", "Type1", "Type2")


data_time_mean <- rbind(data_time_mean[,c(1,2,5,6,7)], data_iter_mean)
data_time_mean$Type2 <- factor(data_time_mean$Type2, levels = c("Time", "Iterations"))
data_time_iter$Type2 <- factor(data_time_iter$Type2, levels=c("Time", "Iterations"))

data_hline <- data.frame(Type2 = c("Time","Iterations"),  # Create data for lines
                         hline = c(0, NA))
data_hline$Type2 <- factor(data_hline$Type2, levels=c("Time", "Iterations"))

p1 <- ggplot(data_time_iter[which(data_time_iter$Type1 == "S1"),], aes(x = as.factor(d), y = Value)) +
      geom_point(aes(colour = Type),size = 1,
                 show.legend = TRUE,
                 position = position_dodge(width = 0.3)) +
      geom_point(data = data_time_mean[which(data_time_mean$Type1=="S1"),],
                 aes(x = as.factor(d), y = Value, colour = Type), size = 3,
                 position = position_dodge(width = 0.3), show.legend = FALSE) +
      geom_hline(data = data_hline,
                 aes(yintercept = hline),
                 linetype = "dashed") +
      geom_line(data = data_time_mean[which(data_time_mean$Type1 == "S1"),],
                aes(y = Value, group = Type, col = Type),
                position = position_dodge(width = 0.3)) +
      scale_color_manual(name = "", values = c("MCD" = "#F8766D", "logM" = "#619CFF")) +
      facet_grid2(c("Type2","Type1") , scales="free_y",switch = "y")+
      theme_bw() +
      scale_y_continuous(breaks=NULL,
                         sec.axis = sec_axis(~ . * 1,labels=scaleFUN, breaks=NULL)) +
      scale_x_discrete(breaks = dgrid) +
      xlab("Dimension") +
      ylab("") +
      theme(panel.grid.minor = element_blank(),
            axis.text = element_text(size = 9),
            panel.grid.major = element_blank(),
            legend.position = "bottom",
            panel.spacing = unit(0.1, "lines"))+
      ggh4x::facetted_pos_scales(y = list(Type2 == "Time" ~ scale_y_continuous(breaks = NULL,
                                                                               trans = myscale_trans(),
                                                                               sec.axis = sec_axis(~ . * 1,
                                                                                                   labels = scaleFUN,
                                                                                                   breaks = c(min(data_time_mean[data_time_mean$Type2 == "Time", "Value"]),
                                                                                                              max(data_time_mean[data_time_mean$Type2 == "Time", "Value"])))),
                                          Type2 == "Iterations" ~ scale_y_continuous(breaks = NULL,
                                                                                     trans = myscale_trans(),
                                                                                     sec.axis = sec_axis(~ . * 1,
                                                                                                         breaks=c(70,80,90,100)))
  ))
p1


ggsave("plot_TIME_ITER.eps",ggarrange(p1,ggplot() + theme_void(), p2, nrow=1,align = "h",widths = c(1.01,0.001,1),common.legend = TRUE, legend="bottom") +
         theme(plot.margin = margin(-0.1,-0.1,-0.1,-0.1, "cm")),width=20, height=15, units="cm")
