library(bamlss)
library(mvnchol)
library(SCM)
seed <- 13
n <- 50000 #1e6
d <- 3

set.seed(seed)
data <- matrix(0, n, d)
x0 <- runif(n, 0, 1)
sx0 <- sort(x0)

fmean <- function(x, m1, m2, m3, m4) (m1+0.2) * x^11 * (m2 * (1 - x))^6 + m3 * (m4 * x)^3 * (1 - x)^10
fvcov <- function(x, c1, c2, c3, beta0c) beta0c + c1 * sin(2 * pi * (x + c2)) + c3 * cos(2 * pi * (x+c2))

c1 <- runif(d * (d + 1)/2, -0.5, 0.5)
c2 <- runif(d * (d + 1)/2, -0.5, 0.5)
c3 <- runif(d * (d + 1)/2, -0.25, 0.25)
beta0c <- runif(d * (d + 1)/2, -0.25, 0.25)

m1 <- runif(d, 0, 0.5)
m2 <- runif(d, 9, 11)
m3 <- runif(d, 9, 11)
m4 <- runif(d, 9, 11)

mu <- list()
Sigma<- list()
Theta <- list()

for(i in 1:n){
  mu[[i]] <- rep(0,d)
  for(j in 1:d) mu[[i]][j] <- fmean(sx0[i], m1[j], m2[j], m3[j], m4[j])

  Theta[[i]] <- matrix(0,d,d)
  Theta[[i]][1,1] <- fvcov(sx0[i], c1[1], c2[1], c3[1], beta0c[1])
  count <- d + 1
  for(j in 2:d){
    Theta[[i]][j,j] <- fvcov(sx0[i], c1[j], c2[j], c3[j], beta0c[j])
    for(k in 1 : (j - 1)){
      Theta[[i]][j,k] <- Theta[[i]][k,j] <- fvcov(sx0[i], c1[count], c2[count], c3[count], beta0c[count])
      count <- count + 1
    }
  }
    D <- exp(0.5*diag(Theta[[i]]))
    Dm2 <-diag(D^(-2))

    T <- matrix(0, d, d)
    for(j in 2:d){
      for(k in 1:(j-1)){
        T[j,k] <- Theta[[i]][j,k]
      }
    }
    diag(T) <- rep(1,d)
    iSigma <- t(T)%*%Dm2%*%T
    Sigma[[i]] <- solve(iSigma)
    L <- solve(T)
    u <- mvnfast::rmvn(1, rep(0, d), diag(rep(1, d)))
    data[i,] <- t(mu[[i]] + L %*% t(u*D))

}

colnames(data) <-  paste0("y_", 1 : d)

expl <- "sx0"

mformula <- function(d=2, expl){

  meanv_f <- paste0("y_", 1:(d-1), sep = "|", collapse="")
  meanv_l <- paste0("y_",d , "~ s(", expl,", k=10)")
  mean_foo <- formula(paste0(meanv_f, meanv_l))
  labelTh <- rep(0, d*(d+1)/2)

  for(j in (d+1):(d+d*(d+1)/2)){
    labelTh[j-d] <- SCM:::labTh(d,j)
  }


  covx_foo_f <- covx_foo_l <- c()
  covx_foo_f <- paste0( labelTh[1:(d*(d+1)/2 - 1)], sep = "|", collapse = "" )
  covx_foo_l <- paste0( labelTh[d*(d+1)/2], "~ s(", expl,")")
  covx_foo <- formula(paste0(covx_foo_f , covx_foo_l))

  return(c(mean_foo,covx_foo))
}


foo <- mformula(d=3, expl="sx0")
old <- Sys.time()
g1 <- gam_scm(foo, family = mvn_scm(d = 3), data = data.frame(data))#,aGam = list( control = list(trace=TRUE)))
new <- Sys.time() - old
#Time difference of 1.007471 mins

setwd("C:/Users/Gioia/Desktop/ScalableAdditiveCovarianceMatrixModels")
save(g1, file="fit_efs.RData")
load("fit_efs.RData")


old <- Sys.time()
g2 <- gam_scm(foo, family = mvn_scm(d = 3), data = data.frame(data), optimizer="bfgs")
new <- Sys.time() - old
#Time difference of 7.369316 mins

save(g2, file="fit_bfgs.RData")
load("fit_bfgs.RData")


g1$l
#[1] -205305.6
g1$aic
#[1] 410783.3

g2$aic
#[1] 410783.2
g2$l
#[1] -205305.5



f <- make_formula(y_1 | y_2 | y_3 ~ s(sx0,k=10) | s(sx0) | s(sx0))
b1 <- bamlss(f, family = mvnchol_bamlss(k = 3, type="modified"), data = data.frame(data), sampler=FALSE)
#AICc 410782.8 logPost -205631. logLik -205305. edf 86.031 eps 0.0000 iteration  10
#elapsed time:  5.36min

save(b1, file="fit_bamlss.RData")
load("fit_bamlss.RData")



set.seed(456)
n <- nrow(data)
batch_ids <- lapply(1:500, function(...) sample(n, size = 1000)) # 100/10000
#batch_ids <- c("nobs" = 1000, "nbatch" = 5)
length(unique(unlist(batch_ids)))


b1_batch <- bamlss(f, data = data.frame(data), family = mvnchol_bamlss(k = 3, type="modified"),
                   sampler = FALSE, optimizer = opt_bbfit,
                   nu = 0.1, always = TRUE, AIC = TRUE,
                   batch_ids = batch_ids) #nbatch = 0.5)#, #)
#elapsed time: 10.73min

save(b1_batch, file="fit_bamlss_batch500.RData")
load("fit_bamlss_batch500.RData")

b1_resampling <- bamlss(f, data = data.frame(data), family = mvnchol_bamlss(k = 3, type="modified"),
                   sampler = FALSE, optimizer = opt_bbfit,
                   slice = TRUE, AIC = TRUE,
                   batch_ids = batch_ids) #nbatch = 0.5)#, #)
#elapsed time: 11.18min
save(b1_resampling, file="fit_bamlss_batchR500.RData")
load("fit_bamlss_batchR500.RData")


plot(g1, page=1, scale=FALSE,se=FALSE)
plot(g2, page=1, scale=FALSE,se=FALSE)

par(mfrow=c(3,3))
plot(b1)
plot(b1_batch)
plot(b1_resampling)

par(mfrow=c(3,3))
pathplot(b1_batch, name="mu1.s.s(sx0).b")
pathplot(b1_batch, name="mu1.p", add=TRUE)
pathplot(b1_batch, name="mu2.s.s(sx0).b")
pathplot(b1_batch, name="mu2.p", add=TRUE)
pathplot(b1_batch, name="mu3.s.s(sx0).b")
pathplot(b1_batch, name="mu3.p", add=TRUE)
pathplot(b1_batch, name="innov1.s.s(sx0).b")
pathplot(b1_batch, name="innov1.p", add=TRUE)
pathplot(b1_batch, name="innov2.s.s(sx0).b")
pathplot(b1_batch, name="innov2.p", add=TRUE)
pathplot(b1_batch, name="innov3.s.s(sx0).b")
pathplot(b1_batch, name="innov3.p", add=TRUE)
pathplot(b1_batch, name="phi12.s.s(sx0).b")
pathplot(b1_batch, name="phi12.p", add=TRUE)
pathplot(b1_batch, name="phi13.s.s(sx0).b")
pathplot(b1_batch, name="phi13.p", add=TRUE)
pathplot(b1_batch, name="phi23.s.s(sx0).b")
pathplot(b1_batch, name="phi23.p", add=TRUE)

par(mfrow=c(3,3))
pathplot(b1_resampling, name="mu1.s.s(sx0).b")
pathplot(b1_resampling, name="mu1.p", add=TRUE)
pathplot(b1_resampling, name="mu2.s.s(sx0).b")
pathplot(b1_resampling, name="mu2.p", add=TRUE)
pathplot(b1_resampling, name="mu3.s.s(sx0).b")
pathplot(b1_resampling, name="mu3.p", add=TRUE)
pathplot(b1_resampling, name="innov1.s.s(sx0).b")
pathplot(b1_resampling, name="innov1.p", add=TRUE)
pathplot(b1_resampling, name="innov2.s.s(sx0).b")
pathplot(b1_resampling, name="innov2.p", add=TRUE)
pathplot(b1_resampling, name="innov3.s.s(sx0).b")
pathplot(b1_resampling, name="innov3.p", add=TRUE)
pathplot(b1_resampling, name="phi12.s.s(sx0).b")
pathplot(b1_resampling, name="phi12.p", add=TRUE)
pathplot(b1_resampling, name="phi13.s.s(sx0).b")
pathplot(b1_resampling, name="phi13.p", add=TRUE)
pathplot(b1_resampling, name="phi23.s.s(sx0).b")
pathplot(b1_resampling, name="phi23.p", add=TRUE)


mu1<-unlist(lapply(1:n, function(x) mu[[x]][1]))
mu2<-unlist(lapply(1:n, function(x) mu[[x]][2]))
mu3<-unlist(lapply(1:n, function(x) mu[[x]][3]))
Th11<-unlist(lapply(1:n, function(x) Theta[[x]][1,1]))
Th22<-unlist(lapply(1:n, function(x) Theta[[x]][2,2]))
Th33<-unlist(lapply(1:n, function(x) Theta[[x]][3,3]))
Th12<-unlist(lapply(1:n, function(x) Theta[[x]][1,2]))
Th13<-unlist(lapply(1:n, function(x) Theta[[x]][1,3]))
Th23<-unlist(lapply(1:n, function(x) Theta[[x]][2,3]))

pg1 <- predict(g1)
pg2 <- predict(g2)
pb1 <- predict(b1)
pb1_batch <- predict(b1_batch)
pb1_resampling <- predict(b1_resampling)

max(abs(pg1-pg2))
#[1] 0.002499082

max(abs(pg1-cbind(pb1$mu1,pb1$mu2,pb1$mu3,pb1$innov1, pb1$innov2, pb1$innov3, -pb1$phi12, -pb1$phi13, -pb1$phi23)))
#[1] 0.02103329


par(mfrow=c(1,1))
plot(mu1, ylim=c(-1,18))
lines(pg1[,1], col="red")
lines(pg2[,1], col="blue")
lines(pb1$mu1, col="green")
lines(pb1_batch$mu1, col="violet")
lines(pb1_resampling$mu1, col="purple")


par(mfrow=c(1,1))
plot(mu2, ylim=c(-1,11))
lines(pg1[,2], col="red")
lines(pg2[,2], col="blue")
lines(pb1$mu2, col="green")
lines(pb1_batch$mu2, col="violet")
lines(pb1_resampling$mu2, col="purple")


par(mfrow=c(1,1))
plot(mu3, ylim=c(-1,9))
lines(pg1[,3], col="red")
lines(pg2[,3], col="blue")
lines(pb1$mu3, col="green")
lines(pb1_batch$mu3, col="violet")
lines(pb1_resampling$mu3, col="purple")





par(mfrow=c(1,1))
plot(Th11, ylim=c(-1,1))
lines(pg1[,4], col="red")
lines(pg2[,4], col="blue")
lines(pb1$innov1, col="green")
lines(pb1_batch$innov1, col="violet")
lines(pb1_resampling$innov1, col="purple")

par(mfrow=c(1,1))
plot(Th22, ylim=c(-0.75,0.5))
lines(pg1[,5], col="red")
lines(pg2[,5], col="blue")
lines(pb1$innov2, col="green")
lines(pb1_batch$innov2, col="violet")
lines(pb1_resampling$innov2, col="purple")

par(mfrow=c(1,1))
plot(Th33, ylim=c(-0.5,0.5))
lines(pg1[,6], col="red")
lines(pg2[,6], col="blue")
lines(pb1$innov3, col="green")
lines(pb1_batch$innov3, col="violet")
lines(pb1_resampling$innov3, col="purple")


par(mfrow=c(1,1))
plot(Th12, ylim=c(-0.75,0.5))
lines(pg1[,7], col="red")
lines(pg2[,7], col="blue")
lines(-pb1$phi12, col="green")
lines(-pb1_batch$phi12, col="violet")
lines(-pb1_resampling$phi12, col="purple")

par(mfrow=c(1,1))
plot(Th13, ylim=c(-0.75,0.5))
lines(pg1[,8], col="red")
lines(pg2[,8], col="blue")
lines(-pb1$phi13, col="green")
lines(-pb1_batch$phi13, col="violet")
lines(-pb1_resampling$phi13, col="purple")

par(mfrow=c(1,1))
plot(Th23, ylim=c(-0.5,0.75))
lines(pg1[,9], col="red")
lines(pg2[,9], col="blue")
lines(-pb1$phi23, col="green")
lines(-pb1_batch$phi23, col="violet")
lines(-pb1_resampling$phi23, col="purple")





set.seed(456)
n <- nrow(data)
batch_ids <- lapply(1:100, function(...) sample(n, size = 1000))
#batch_ids <- c("nobs" = 1000, "nbatch" = 5)
length(unique(unlist(batch_ids)))


b2_batch <- bamlss(f, data = data.frame(data), family = mvnchol_bamlss(k = 3, type="modified"),
                   sampler = FALSE, optimizer = opt_bbfit,
                   nu = 0.1, always = TRUE, AIC = TRUE,
                   batch_ids = batch_ids) #nbatch = 0.5)#, #)
#elapsed time:  2.22min
save(b2_batch, file="fit_bamlss_batch100.RData")
load("fit_bamlss_batch100.RData")

b2_resampling <- bamlss(f, data = data.frame(data), family = mvnchol_bamlss(k = 3, type="modified"),
                        sampler = FALSE, optimizer = opt_bbfit,
                        slice = TRUE, AIC = TRUE,
                        batch_ids = batch_ids) #nbatch = 0.5)#, #)
#elapsed time:  2.52min
save(b2_resampling, file="fit_bamlss_batchR100.RData")
load("fit_bamlss_batchR100.RData")


plot(g1, page=1, scale=FALSE,se=FALSE)
plot(g2, page=1, scale=FALSE,se=FALSE)

par(mfrow=c(3,3))
plot(b2_batch)
plot(b2_resampling)

par(mfrow=c(3,3))
pathplot(b2_batch, name="mu1.s.s(sx0).b")
pathplot(b2_batch, name="mu1.p", add=TRUE)
pathplot(b2_batch, name="mu2.s.s(sx0).b")
pathplot(b2_batch, name="mu2.p", add=TRUE)
pathplot(b2_batch, name="mu3.s.s(sx0).b")
pathplot(b2_batch, name="mu3.p", add=TRUE)
pathplot(b2_batch, name="innov1.s.s(sx0).b")
pathplot(b2_batch, name="innov1.p", add=TRUE)
pathplot(b2_batch, name="innov2.s.s(sx0).b")
pathplot(b2_batch, name="innov2.p", add=TRUE)
pathplot(b2_batch, name="innov3.s.s(sx0).b")
pathplot(b2_batch, name="innov3.p", add=TRUE)
pathplot(b2_batch, name="phi12.s.s(sx0).b")
pathplot(b2_batch, name="phi12.p", add=TRUE)
pathplot(b2_batch, name="phi13.s.s(sx0).b")
pathplot(b2_batch, name="phi13.p", add=TRUE)
pathplot(b2_batch, name="phi23.s.s(sx0).b")
pathplot(b2_batch, name="phi23.p", add=TRUE)

par(mfrow=c(3,3))
pathplot(b2_resampling, name="mu1.s.s(sx0).b")
pathplot(b2_resampling, name="mu1.p", add=TRUE)
pathplot(b2_resampling, name="mu2.s.s(sx0).b")
pathplot(b2_resampling, name="mu2.p", add=TRUE)
pathplot(b2_resampling, name="mu3.s.s(sx0).b")
pathplot(b2_resampling, name="mu3.p", add=TRUE)
pathplot(b2_resampling, name="innov1.s.s(sx0).b")
pathplot(b2_resampling, name="innov1.p", add=TRUE)
pathplot(b2_resampling, name="innov2.s.s(sx0).b")
pathplot(b2_resampling, name="innov2.p", add=TRUE)
pathplot(b2_resampling, name="innov3.s.s(sx0).b")
pathplot(b2_resampling, name="innov3.p", add=TRUE)
pathplot(b2_resampling, name="phi12.s.s(sx0).b")
pathplot(b2_resampling, name="phi12.p", add=TRUE)
pathplot(b2_resampling, name="phi13.s.s(sx0).b")
pathplot(b2_resampling, name="phi13.p", add=TRUE)
pathplot(b2_resampling, name="phi23.s.s(sx0).b")
pathplot(b2_resampling, name="phi23.p", add=TRUE)


pb2_batch <- predict(b2_batch)
pb2_resampling <- predict(b2_resampling)


par(mfrow=c(1,1))
plot(mu1, ylim=c(-1,18))
lines(pg1[,1], col="red")
lines(pg2[,1], col="blue")
lines(pb1$mu1, col="green")
lines(pb1_batch$mu1, col="violet")
lines(pb1_resampling$mu1, col="purple")


par(mfrow=c(1,1))
plot(mu2, ylim=c(-1,11))
lines(pg1[,2], col="red")
lines(pg2[,2], col="blue")
lines(pb1$mu2, col="green")
lines(pb2_batch$mu2, col="violet")
lines(pb2_resampling$mu2, col="purple")


par(mfrow=c(1,1))
plot(mu3, ylim=c(-1,9))
lines(pg1[,3], col="red")
lines(pg2[,3], col="blue")
lines(pb1$mu3, col="green")
lines(pb2_batch$mu3, col="violet")
lines(pb2_resampling$mu3, col="purple")


par(mfrow=c(1,1))
plot(Th11, ylim=c(-1,1))
lines(pg1[,4], col="red")
lines(pg2[,4], col="blue")
lines(pb1$innov1, col="green")
lines(pb2_batch$innov1, col="violet")
lines(pb2_resampling$innov1, col="purple")


par(mfrow=c(1,1))
plot(Th22, ylim=c(-0.75,0.5))
lines(pg1[,5], col="red")
lines(pg2[,5], col="blue")
lines(pb1$innov2, col="green")
lines(pb2_batch$innov2, col="violet")
lines(pb2_resampling$innov2, col="purple")

par(mfrow=c(1,1))
plot(Th33, ylim=c(-0.5,0.5))
lines(pg1[,6], col="red")
lines(pg2[,6], col="blue")
lines(pb1$innov3, col="green")
lines(pb2_batch$innov3, col="violet")
lines(pb2_resampling$innov3, col="purple")


par(mfrow=c(1,1))
plot(Th12, ylim=c(-0.75,0.5))
lines(pg1[,7], col="red")
lines(pg2[,7], col="blue")
lines(-pb1$phi12, col="green")
lines(-pb2_batch$phi12, col="violet")
lines(-pb2_resampling$phi12, col="purple")

par(mfrow=c(1,1))
plot(Th13, ylim=c(-0.75,0.5))
lines(pg1[,8], col="red")
lines(pg2[,8], col="blue")
lines(-pb1$phi13, col="green")
lines(-pb2_batch$phi13, col="violet")
lines(-pb2_resampling$phi13, col="purple")

par(mfrow=c(1,1))
plot(Th23, ylim=c(-0.5,0.75))
lines(pg1[,9], col="red")
lines(pg2[,9], col="blue")
lines(-pb1$phi23, col="green")
lines(-pb2_batch$phi23, col="violet")
lines(-pb2_resampling$phi23, col="purple")




