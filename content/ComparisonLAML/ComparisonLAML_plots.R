###########################################################################################################################################
# Code for reproducibility of the results of Section 2.4:                                                                                 #
# Visualising the computational time for comparing the old and the new formulation of the term involved in LAML-derivative based approach #
# and the time for fomputing the third derivatives is also included                                                                       #
###########################################################################################################################################
rm(list = ls())
setwd("C:/Users/Gioia/Desktop/PhD Thesis/Code_Reproducibility")
source("loadPackages.R")
setwd("C:/Users/Gioia/Desktop/PhD Thesis/Code_Reproducibility/Chapter2/Section2_4")


#######################################
# Visualising the computational times #
#######################################
summary_time <- function(obj,
                         dgrid1=NULL, dgrid2=NULL,nrun){
  type <- c("mcd_eta","complete","new") 
  type1 <- c("D3",  "LAML - STD", "LAML - EFF")
  
  out <- lapply(1:length(type), function(jj){
    str_type <- paste0("obj[[x]]$time$td3_",type[[jj]])
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

myscale_trans <- function()
{
  trans_new("myscale", function(x) sqrt(x),
            function(x) x^2, domain = c(0, Inf))
}

scaleFUN <- function(x) sprintf("%.2f", x)

###################################
# Extract the computational times #
###################################
all_time <- function(obj, dgrid1=NULL, dgrid2=NULL,nrun){
  type <- c("mcd_eta","complete","new") 
  type1 <- c("D3",  "LAML - STD", "LAML - EFF")
  
  out <- lapply(1:length(type), function(jj){
    str_type <- paste0("obj[[x]]$time$td3_",type[[jj]])
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


###################################
# Code for reproducing Figure 2.5 #
###################################
setwd("C:/Users/Gioia/Desktop/PhD Thesis/Code_Reproducibility/Chapter2/Section2_4")
load("TIME_d3.RData")
dgrid <-c(2,3,5,10,15,20)
nrun <- 7

time_d3_dHess_drho <-  all_time(obj=TIME_d3, dgrid1=dgrid, nrun=7)
mean_time_d3_dHess_drho <-  summary_time(obj=TIME_d3, dgrid1=dgrid, nrun=7)

label_mean_time <- mean_time_d3_dHess_drho[c(6,12,18),"mean_time"]

p1 <- ggplot(time_d3_dHess_drho, aes(x=as.factor(d), y=time)) + 
  geom_point(aes(colour=Type),size=1,
             show.legend=TRUE,
             position=position_dodge(width=0.3))+      
  geom_point(data=mean_time_d3_dHess_drho , aes(x=as.factor(d), y=mean_time, colour=Type),size=3,
             position=position_dodge(width=0.3),show.legend=FALSE)+      
  geom_hline(yintercept = 0, linetype="dashed") +
  geom_line(data=mean_time_d3_dHess_drho, aes(y=mean_time, group=Type, col=Type),position=position_dodge(width=0.3))+
  scale_color_manual(name="",values = c("D3" = "#00BA38", "LAML - STD" = "#619CFF", "LAML - EFF" = "#F8766D"))+
  theme_bw() +
  scale_y_continuous(breaks=NULL,  trans = myscale_trans(),
                     sec.axis = sec_axis(~ . * 1,labels=scaleFUN, breaks=label_mean_time), 
  ) +  
  scale_x_discrete(breaks = c(2,3,5,10,15,20)) + xlab("Dimension") +
  ylab("") +
  theme(panel.grid.minor = element_blank(),axis.text=element_text(size=9), panel.grid.major = element_blank(),legend.position="bottom",panel.spacing = unit(0.1, "lines")) 

p1

ggsave("plot_TIME_d3.eps",p1 + 
         theme(plot.margin = margin(-0.01,-0.01,-0.01,-0.01, "cm")),width=20, height=15, units="cm")




