# nohup R CMD BATCH oubmbmabc.r >/dev/null &
rm(list=ls())
library(TreeSim)
library(EasyABC)
library(coda)
library(Sim.DiffProc)
library(MCMCpack)
library(ape)
library(abc)
source("/Users/TerryLai/Dropbox/TerryLai/R_code/abc_V2/oubmbmabc.r")
source("/Users/TerryLai/Dropbox/TerryLai/R_code/abc_V2/ououbmabc.r")
source("/Users/TerryLai/Dropbox/TerryLai/R_code/abc_V2/oubmcirabc.r")
source("/Users/TerryLai/Dropbox/TerryLai/R_code/abc_V2/ououcirabc.r")


setwd("/Users/TerryLai/Dropbox/TerryLai/R_code/abc_V2")
system("ls")
library(abc)
library(TreeSim)
library(ape)
n<-20
numbsim<-1;lambda<-2.0;mu<-0.5;frac<-0.6;age<-2
tree<-sim.bd.taxa.age(n=n,numbsim=1,lambda=lambda,mu=mu,frac=frac,age=age,mrca=TRUE)[[1]]
tree<-reorder(tree,"postorder")
tree$edge
plot(tree)
nodelabels()
tiplabels()



models<-rep(c("oubmbm","ououbm","oubmcir","ououcir"), each=50000)
full.stat.sim.sum <- array(0, dim=c(200000,6))
dim(full.stat.sim.sum)
raw.sum.stat <- array(0, dim=c(4,6))
load("test_oubmbm.Rdata")
dim(sim.sum.stat)
full.stat.sim.sum[1:50000,]=sim.sum.stat
raw.sum.stat[1,]<-c(raw.sum.stat.y,raw.sum.stat.x1,raw.sum.stat.x2)
load("test_ououbm.Rdata")
full.stat.sim.sum[50001:100000,]=sim.sum.stat
raw.sum.stat[2,]<-c(raw.sum.stat.y,raw.sum.stat.x1,raw.sum.stat.x2)
load("test_oubmcir.Rdata")
full.stat.sim.sum[100001:150000,]=sim.sum.stat
raw.sum.stat[3,]<-c(raw.sum.stat.y,raw.sum.stat.x1,raw.sum.stat.x2)
load("test_ououcir.Rdata")
full.stat.sim.sum[150001:200000,]=sim.sum.stat
raw.sum.stat[4,]<-c(raw.sum.stat.y,raw.sum.stat.x1,raw.sum.stat.x2)
colnames(full.stat.sim.sum)<-c("y.mean","y.sd","x1.mean","x1.sd","x2.mean","x2.sd")
params <- c("y.mean","y.sd","x1.mean","x1.sd","x2.mean","x2.sd")
head(full.stat.sim.sum)
tail(full.stat.sim.sum)
cv.modsel.rej<-cv4postpr(models, full.stat.sim.sum, nval=5, tols=0.01, method="rejection")
s.rej<-summary(cv.modsel.rej)
plot(cv.modsel.rej, names.arg=c("oubmbm","ououbm","oubmcir","ououcir"))
par(mfrow=c(2,3))
for(i in 1:6){
  boxplot(full.stat.sim.sum[,i]~models, main=params[i])
}
# cv.modsel.log<-cv4postpr(models, full.stat.sim.sum, nval=5, tols=0.10, method="mnlogistic")
# s.log<-summary(cv.modsel.log)
# cv.modsel.neu<-cv4postpr(models, full.stat.sim.sum, nval=5, tols=0.01, method="neuralnet")
# s.neu<-summary(cv.modsel.neu)
raw.sum.stat <- c(raw.sum.stat.y, raw.sum.stat.x1, raw.sum.stat.x2)

### model seiection by BF
modsel<-postpr(target=raw.sum.stat, index=models, sumstat=full.stat.sim.sum, tol=0.05, method="mnlogistic")
summary(modsel)

### Goodness of fit
res.gfit.oubmbm <- gfit(target=raw.sum.stat, sumstat=full.stat.sim.sum[models=="oubmbm",], statistic=mean, nb.replicate = 100);summary(res.gfit.oubmbm)

res.gfit.ououbm <- gfit(target=raw.sum.stat, sumstat=full.stat.sim.sum[models=="ououbm",], statistic=mean, nb.replicate = 100);summary(res.gfit.ououbm)

res.gfit.oubmcir <- gfit(target=raw.sum.stat, sumstat=full.stat.sim.sum[models=="oubmcir",], statistic=mean, nb.replicate = 100);summary(res.gfit.oubmcir)

res.gfit.ououcir <- gfit(target=raw.sum.stat, sumstat=full.stat.sim.sum[models=="ououcir",], statistic=mean, nb.replicate = 100);summary(res.gfit.ououcir)
#gfitpca(target=raw.sum.stat, sumstat=full.stat.sim.sum, index=models, cprob=0.1)
### Posterior predictive checks ---  currently Not running
labels <- c("y.mean","y.sd","x1.mean","x1.sd","x2.mean","x2.sd")
par(mfrow=c(2,3))
for(i in c(1:6)){
  hist(full.stat.sim.sum[,i], breaks=10000, xlab=labels[i], main="")
  abline(v=raw.sum.stat[i], col=2)
}
### Cross-validation
stat.oubmbm.sim <- subset(full.stat.sim.sum, subset=models=="oubmbm")
cv.res.rej.oubmbm <- cv4abc(full.stat.sim.sum, stat.oubmbm.sim, nval=10, tols=c(0.005,0.01,0.05),method="rejection")
summary(cv.res.rej.oubmbm)

stat.ououbm.sim <- subset(full.stat.sim.sum, subset=models=="ououbm")
cv.res.rej.ououbm <- cv4abc(full.stat.sim.sum, stat.ououbm.sim, nval=10, tols=c(0.005,0.01,0.05),method="rejection")
summary(cv.res.rej.ououbm)

stat.oubmcir.sim <- subset(full.stat.sim.sum, subset=models=="oubmcir")
cv.res.rej.oubmcir <- cv4abc(full.stat.sim.sum, stat.oubmcir.sim, nval=10, tols=c(0.005,0.01,0.05),method="rejection")
summary(cv.res.rej.oubmcir)

stat.ououcir.sim <- subset(full.stat.sim.sum, subset=models=="ououcir")
cv.res.rej.ououcir <- cv4abc(full.stat.sim.sum, stat.ououcir.sim, nval=10, tols=c(0.005,0.01,0.05),method="rejection")
summary(cv.res.rej.ououcir)

### Cross-validation for parameter inference
par(mfrow=c(2,3))
plot(cv.res.rej.oubmbm,caption = c(expression(mu[y]),expression(sigma[y]),expression(mu[x[1]]),expression(sigma[x[1]]),expression(mu[x[2]]),expression(sigma[x[2]])))

plot(cv.res.rej.ououbm,caption = c(expression(mu[y]),expression(sigma[y]),expression(mu[x[1]]),expression(sigma[x[1]]),expression(mu[x[2]]),expression(sigma[x[2]])))

plot(cv.res.rej.oubmcir,caption = c(expression(mu[y]),expression(sigma[y]),expression(mu[x[1]]),expression(sigma[x[1]]),expression(mu[x[2]]),expression(sigma[x[2]])))

plot(cv.res.rej.ououcir,caption = c(expression(mu[y]),expression(sigma[y]),expression(mu[x[1]]),expression(sigma[x[1]]),expression(mu[x[2]]),expression(sigma[x[2]])))

res <- abc(target = raw.sum.stat,param= ,sumstat=full.stat.sim.sum, tol=0.05, method="rejection")
