#####################################################################
# Code for reproducibility of the results of Section 2.3            #
# Plot computational comparison adopting the Hessian block strategy #
#####################################################################
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



library(rstudioapi)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

################################
# Plots of computational times #
################################
summary_time <- function(obj, dgrid1=NULL, dgrid2=NULL,nrun){
  type <- c("block","noblock")
  type1 <- c("PARS",  "STD")

  out <- lapply(1:length(type), function(jj){
    str_type <- paste0("obj[[x]]$time$thessian_",type[[jj]])
    data_time <-  data.frame(cbind(unlist(lapply(1:nrun, function(x) eval(parse(text=str_type))/1e9)),rep(dgrid,times=nrun)))
    colnames(data_time) <- c("time","d")
    data <- aggregate(data_time$time, list(data_time$d), FUN=mean)
    data_min <- aggregate(data_time$time, list(data_time$d), FUN=min)
    data_max <- aggregate(data_time$time, list(data_time$d), FUN=max)
    data <- cbind(data,data_min[,2],data_max[,2])
    return(data)
  })


  res <- matrix(0, length(dgrid)*length(out), 4)
  for(j in 1:length(out)){
    res[c((length(dgrid)*(j-1)+1):(length(dgrid)*j)),] <- as.matrix(out[[j]])
  }
  res <- as.data.frame(res)
  res <- cbind(res, rep(type1,each=length(dgrid)))
  colnames(res) <- c("d","mean_time","min_time","max_time","Type")
  return(data.frame(res))
}


###################################
# Extract the computational times #
###################################
all_time <- function(obj, dgrid1=NULL, dgrid2=NULL,nrun){
  type <- c("block","noblock")
  type1 <- c("PARS",  "STD")

  out <- lapply(1:length(type), function(jj){
    str_type <- paste0("obj[[x]]$time$thessian_",type[[jj]])
    data_time <-  data.frame(cbind(unlist(lapply(1:nrun, function(x) eval(parse(text=str_type))/1e9)),rep(dgrid,times=nrun)))
    colnames(data_time) <- c("time","d")
    return(data_time)
  })

  res <- out[[1]]
  for(j in 2:length(out)){
    res <- rbind(res, out[[j]])
  }

 res <- cbind(res, rep(type1,each=length(dgrid)*nrun))

  colnames(res) <- c("time","d", "Type")
  return(data.frame(res))
}

myscale_trans <- function(){
  trans_new("myscale", function(x) sqrt(x),
            function(x) x^2, domain = c(0, Inf))
}

scaleFUN <- function(x) sprintf("%.2f", x)

###################################
# Code for reproducing Figure 2.3 #
###################################
setwd("C:/Users/Gioia/Desktop/PhD Thesis/Code_Reproducibility/Chapter2/Section2_3")
load("TIME_mcd_beta.RData")
load("TIME_logm_beta.RData")


dgrid1 <-c(2,3,5,10,15,20)
nrun <- 1
pint_type <- c("dm05","dm1", "dm2","const")

time_mcd_dmod <-  lapply(1:length(pint_type), function(x)
                              all_time(obj=TIME_MCD_beta[[x]], dgrid1=dgrid1, nrun=nrun))
time_logm_dmod <-  lapply(1:length(pint_type), function(x)
                          all_time(obj=TIME_logM_beta[[x]], dgrid1=dgrid1, nrun=nrun))


mean_time_mcd_dmod <-  lapply(1:length(pint_type), function(x)
                         summary_time(obj=TIME_MCD_beta[[x]], dgrid1=dgrid1, nrun=nrun))
mean_time_logm_dmod <-  lapply(1:length(pint_type), function(x)
                          summary_time(obj=TIME_logM_beta[[x]], dgrid1=dgrid1, nrun=nrun))



load("TIME_mcd_beta_large.RData")
load("TIME_logm_beta_large.RData")
dgrid2 <-c(25,50,75)


time_mcd_dlarge <-  lapply(1:length(pint_type), function(x)
                           all_time(obj=TIME_MCD_beta[[x]], dgrid1=dgrid2, nrun=nrun))
time_logm_dlarge <-  lapply(1:length(pint_type), function(x)
                            all_time(obj=TIME_logM_beta[[x]], dgrid1=dgrid2, nrun=nrun))


mean_time_mcd_dlarge <-  lapply(1:length(pint_type), function(x)
                                summary_time(obj=TIME_MCD_beta[[x]], dgrid1=dgrid2, nrun=nrun))

mean_time_logm_dlarge <-  lapply(1:length(pint_type), function(x)
                                 summary_time(obj=TIME_logM_beta[[x]], dgrid1=dgrid2, nrun=nrun))

time_mcd <- list()
time_logm <- list()

for(j in 1:4){
  time_mcd[[j]]<-rbind(time_mcd_dmod[[j]][1:(nrun*length(dgrid1)),],
                       time_mcd_dlarge[[j]][1:(nrun*length(dgrid2)),],
                       time_mcd_dmod[[j]][(nrun*length(dgrid1) + 1):(2*nrun*length(dgrid1)),],
                       time_mcd_dlarge[[j]][(nrun*length(dgrid2) + 1):(2*nrun*length(dgrid2)),])
  time_logm[[j]]<-rbind(time_logm_dmod[[j]][1:(nrun*length(dgrid1)),],
                        time_logm_dlarge[[j]][1:(nrun*length(dgrid2)),],
                        time_logm_dmod[[j]][(nrun*length(dgrid1) + 1):(2*nrun*length(dgrid1)),],
                        time_logm_dlarge[[j]][(nrun*length(dgrid2) + 1):(2*nrun*length(dgrid2)),])

}

Sc1_time_mcd <- time_mcd[[1]]
Sc1_time_mcd <- cbind(Sc1_time_mcd, Scenario = rep("Scenario 1",nrow(time_mcd[[1]])),
                      param = rep("MCD", nrow(time_mcd[[1]])))

Sc2_time_mcd <- time_mcd[[3]]
Sc2_time_mcd <- cbind(Sc2_time_mcd, Scenario = rep("Scenario 2",nrow(time_mcd[[1]])),
                      param = rep("MCD", nrow(time_mcd[[1]])))

Sc1_time_logm <- time_logm[[1]]
Sc1_time_logm <- cbind(Sc1_time_logm, Scenario = rep("Scenario 1",nrow(time_logm[[1]])),
                       param = rep("logM", nrow(time_mcd[[1]])))

Sc2_time_logm <- time_logm[[3]]
Sc2_time_logm <- cbind(Sc2_time_logm, Scenario = rep("Scenario 2",nrow(time_logm[[3]])),
                       param = rep("logM", nrow(time_mcd[[1]])))

Sc1_time <- rbind(Sc1_time_mcd,Sc1_time_logm)
Sc2_time <- rbind(Sc2_time_mcd,Sc2_time_logm)

mean_time_mcd <- list()
mean_time_logm <- list()

for(j in 1:4){
  mean_time_mcd[[j]]<-rbind(mean_time_mcd_dmod[[j]][1:length(dgrid1),],
                            mean_time_mcd_dlarge[[j]][1:length(dgrid2),],
                            mean_time_mcd_dmod[[j]][(length(dgrid1)+1):(2*length(dgrid1)),],
                            mean_time_mcd_dlarge[[j]][(length(dgrid2)+1):(2*length(dgrid2)),])
  mean_time_logm[[j]]<-rbind(mean_time_logm_dmod[[j]][1:length(dgrid1),],
                             mean_time_logm_dlarge[[j]][1:length(dgrid2),],
                             mean_time_logm_dmod[[j]][(length(dgrid1)+1):(2*length(dgrid1)),],
                             mean_time_logm_dlarge[[j]][(length(dgrid2)+1):(2*length(dgrid2)),])
}

Sc1_mean_time_mcd <- mean_time_mcd[[1]]
Sc1_mean_time_mcd <- cbind(Sc1_mean_time_mcd,
                           Scenario = rep("Scenario 1",nrow(mean_time_mcd[[1]])),
                           param = rep("MCD", nrow(mean_time_mcd[[1]])))

Sc2_mean_time_mcd <- mean_time_mcd[[3]]
Sc2_mean_time_mcd <- cbind(Sc2_mean_time_mcd,
                           Scenario = rep("Scenario 2",nrow(mean_time_mcd[[3]])),
                           param = rep("MCD", nrow(mean_time_mcd[[3]])))


Sc1_mean_time_logm <- mean_time_logm[[1]]
Sc1_mean_time_logm <- cbind(Sc1_mean_time_logm,
                            Scenario = rep("Scenario 1",nrow(mean_time_logm[[1]])),
                            param = rep("logM", nrow(mean_time_logm[[1]])))

Sc2_mean_time_logm <- mean_time_logm[[3]]
Sc2_mean_time_logm <- cbind(Sc2_mean_time_logm,
                            Scenario = rep("Scenario 2",nrow(mean_time_logm[[3]])),
                            param = rep("logM", nrow(mean_time_logm[[3]])))



Sc1_mean_time <- rbind(Sc1_mean_time_mcd,Sc1_mean_time_logm)
Sc2_mean_time <- rbind(Sc2_mean_time_mcd,Sc2_mean_time_logm)

ll_dgrid <- length(dgrid1) + length(dgrid2)
label_SC1 <- Sc1_mean_time[c(ll_dgrid,2*ll_dgrid,3*ll_dgrid,4*ll_dgrid),"mean_time"]
label_SC2 <- Sc2_mean_time[c(ll_dgrid,2*ll_dgrid,3*ll_dgrid,4*ll_dgrid), "mean_time"]

#Plot Scenario 1
p1 <- ggplot(Sc1_time, aes(x=as.factor(d), y=time)) +
             geom_point(aes(colour=Type),size=1,
                        show.legend=TRUE,
                        position=position_dodge(width=0.3))+
  geom_point(data=Sc1_mean_time, aes(x=as.factor(d), y=mean_time, colour=Type),size=3,
             position=position_dodge(width=0.3),show.legend=FALSE)+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_line(data=Sc1_mean_time, aes(y=mean_time, group=Type, col=Type),position=position_dodge(width=0.3))+
  scale_color_manual(name="",values = c("PARS" = "#F8766D", "STD" = "#619CFF"))+
  facet_grid2(c("param","Scenario") , scales="free_y",switch = "y")+
  theme_bw() +
  scale_y_continuous(breaks=NULL,
                     sec.axis = sec_axis(~ . * 1,labels=scaleFUN, breaks=NULL),
  ) +
  scale_x_discrete(breaks = c(dgrid1,dgrid2)) + xlab("Dimension") +
  ylab("") +
  theme(panel.grid.minor = element_blank(),axis.text=element_text(size=9), panel.grid.major = element_blank(),legend.position="bottom",panel.spacing = unit(0.1, "lines")) +
  ggh4x::facetted_pos_scales(y = list(
    param == "MCD" ~ scale_y_continuous(breaks=NULL, trans = myscale_trans(), sec.axis = sec_axis(~ . * 1,labels=scaleFUN,breaks=label_SC1[1:2])),
    param == "logM" ~ scale_y_continuous(breaks=NULL, trans = myscale_trans(),sec.axis = sec_axis(~ . * 1,labels=scaleFUN,breaks=label_SC1[3:4]))
  ))
p1

#Plot Scenario 2
p2 <- ggplot(Sc2_time, aes(x=as.factor(d), y=time)) +
  geom_point(aes(colour=Type),size=1,
             show.legend=TRUE,
             position=position_dodge(width=0.3))+
  geom_point(data=Sc2_mean_time, aes(x=as.factor(d), y=mean_time, colour=Type),size=3,
             position=position_dodge(width=0.3),show.legend=FALSE)+
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_line(data=Sc2_mean_time, aes(y=mean_time, group=Type, col=Type),position=position_dodge(width=0.3))+
  scale_color_manual(name="",values = c("PARS" = "#F8766D", "STD" = "#619CFF")) +
  facet_grid2(c("param","Scenario") , scales="free_y",switch = "y")+
  theme_bw() +
  scale_y_continuous(breaks=NULL,
                     sec.axis = sec_axis(~ . * 1,labels=scaleFUN, breaks=NULL), #  , breaks = unique(ax_lab2[[j]]), labels=scaleFUN),
  ) +
  scale_x_discrete(breaks = c(dgrid1,dgrid2)) + xlab("Dimension") +
  ylab("") +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(),strip.text.y = element_blank(),legend.position="bottom",panel.spacing = unit(0.1, "lines")) +
  ggh4x::facetted_pos_scales(y = list(
    param == "MCD" ~ scale_y_continuous(breaks=NULL, trans = myscale_trans(), sec.axis = sec_axis(~ . * 1,labels=scaleFUN,breaks=label_SC2[1:2])),
    param == "logM" ~ scale_y_continuous(breaks=NULL, trans = myscale_trans(),sec.axis = sec_axis(~ . * 1,labels=scaleFUN,breaks=label_SC2[3:4]))
  ))
p2

ggsave("plot_TIME_blocks_dm05dm2.eps",ggarrange(p1,ggplot() + theme_void(), p2, nrow=1,align = "h",widths = c(1.10,0.0000001,1),common.legend=TRUE, legend="bottom") +
         theme(plot.margin = margin(-0.1,-0.1,-0.1,-0.1, "cm")),width=20, height=15, units="cm")



