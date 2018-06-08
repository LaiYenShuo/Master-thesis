# nohup R CMD BATCH oubmbmabc.r >/dev/null &
#rm(list=ls())
library(TreeSim)
library(EasyABC)
library(coda)
library(Sim.DiffProc)
library(MCMCpack)
library(ape)
library(abc)
source("~/Dropbox/TerryLai/R_code/abcmcmc/ouest_Tony.R")
source("~/Dropbox/TerryLai/R_code/abc_infor_prior/oubmbmabc.r")
source("~/Dropbox/TerryLai/R_code/abc_infor_prior/ououbmabc.r")
source("~/Dropbox/TerryLai/R_code/abc_infor_prior/oubmcirabc.r")
source("~/Dropbox/TerryLai/R_code/abc_infor_prior/ououcirabc.r")
regboundfcn<-function(olssum=olssum){
  bdd<-  c(olssum$coefficients[,1]-4*olssum$coefficients[,2],olssum$coefficients[,1]+4*olssum$coefficients[,2])
  bdd<-bdd[c(1,4,2,5,3,6)]  
  return(bdd)
}
olssum <- summary(lm(resptrait~predtrait1+predtrait2))
regbound<-regboundfcn(olssum=olssum)


sims <- 50000
# hyper paramters
alpha.y.rate <-5# assume exponential
alpha.x.rate <- 5# assume exponential
theta.x.mean <- 0# assume normal
theta.x.sd  <- 1
tau.rate <- 3 # assume invgamma
sigma.x.rate <-2# assume invgamma
alpha.tau.rate <- 5# assume exponential
theta.tau.df <- 2 
sigma.tau.rate <- 2

b0.min=regbound[1]
b0.max=regbound[2]
b1.min=regbound[3]
b1.max=regbound[4]
b2.min=regbound[5]
b2.max=regbound[6]

prior.model.params=c(alpha.y.rate, alpha.x.rate, theta.x.mean, theta.x.sd, tau.rate, sigma.x.rate, alpha.tau.rate, theta.tau.df, sigma.tau.rate)
names(prior.model.params)<-c("alpha.y.rate", "alpha.x.rate", "theta.x.mean", "theta.x.sd", "tau.rate" , "sigma.x.rate", "alpha.tau.rate", "theta.tau.df", "sigma.tau.rate")
prior.reg.params=c(b0.min, b0.max, b1.min, b1.max, b2.min, b2.max)
###############
### oubmbm ###
##############
prior.params.oubmbm <- oubmbmprior(prior.model.params = prior.model.params, prior.reg.params = prior.reg.params)
sim.oubmbm.trait<-array(0,c(dim(dataset)[1],3,sims))
model.params.oubmbm<-array(0,c(3,sims))
rownames(model.params.oubmbm)<-c("alpha.y","sigma.x","tau")
reg.params.oubmbm<-array(0,c(3,sims))
row.names(reg.params.oubmbm)<-c("b0", "b1", "b2")
y.sum.stat.oubmbm<-array(0,c(2,sims))
rownames(y.sum.stat.oubmbm)<-c("y.mean","y.sd")
x1.sum.stat.oubmbm<-array(0,c(2,sims))
rownames(x1.sum.stat.oubmbm)<-c("x1.mean","x1.sd")
x2.sum.stat.oubmbm<-array(0,c(2,sims))
rownames(x2.sum.stat.oubmbm)<-c("x2.mean","x2.sd")

for(simIndex in 1:sims){
  if(simIndex %%5000==0){print(paste("oubmbm",simIndex,sep=""))}
  prior.params<-oubmbmprior(prior.model.params = prior.model.params, prior.reg.params = prior.reg.params)
  model.params.oubmbm[,simIndex]<-prior.params$model.params#for record only
  reg.params.oubmbm[,simIndex]<-prior.params$reg.params#for record only
  
  sim.trait <-oubmbmmodel(model.params=prior.params$model.params,reg.params=prior.params$reg.params,root=root,tree=tree)
  sim.oubmbm.trait[,1,simIndex]<-sim.trait$y
  sim.oubmbm.trait[,2,simIndex]<-sim.trait$x1
  sim.oubmbm.trait[,3,simIndex]<-sim.trait$x2
  y.sum.stat.oubmbm[,simIndex]<- sum.stat(trait=sim.trait$y,tree=tree)
  x1.sum.stat.oubmbm[,simIndex]<- sum.stat(trait=sim.trait$x1,tree=tree)
  x2.sum.stat.oubmbm[,simIndex]<- sum.stat(trait=sim.trait$x2,tree=tree)
}# end of loop

sum.stat.oubmbm <- cbind(t(y.sum.stat.oubmbm),t(x1.sum.stat.oubmbm),t(x2.sum.stat.oubmbm))
oubmbm.par.sim <- cbind(t(model.params.oubmbm),t(reg.params.oubmbm))
###############
### ououbm ###
##############
prior.params.ououbm <- ououbmprior(prior.model.params = prior.model.params, prior.reg.params = prior.reg.params)
sim.ououbm.trait<-array(0,c(dim(dataset)[1],5,sims))
model.params.ououbm<-array(0,c(5,sims))
rownames(model.params.ououbm)<-c("alpha.y","alpha.x","theta.x","sigma.x","tau")
reg.params.ououbm<-array(0,c(3,sims))
row.names(reg.params.ououbm)<-c("b0", "b1", "b2")
y.sum.stat.ououbm<-array(0,c(2,sims))
rownames(y.sum.stat.ououbm)<-c("y.mean","y.sd")
x1.sum.stat.ououbm<-array(0,c(2,sims))
rownames(x1.sum.stat.ououbm)<-c("x1.mean","x1.sd")
x2.sum.stat.ououbm<-array(0,c(2,sims))
rownames(x2.sum.stat.ououbm)<-c("x2.mean","x2.sd")
# rownames(post.model.params.ououbm)<-c("alpha.y","alpha.x","theta.x","sigma.x","tau")
# row.names(post.reg.params.ououbm)<-c("b0", "b1", "b2")

prior.params.ououbm <- ououbmprior(prior.model.params=prior.model.params,prior.reg.params=prior.reg.params)

for(simIndex in 1:sims){
  if(simIndex %%5000==0){print(paste("ououbm",simIndex,sep=""))}
  prior.params <- ououbmprior(prior.model.params=prior.model.params,prior.reg.params=prior.reg.params)
  model.params.ououbm[,simIndex]<-prior.params$model.params#for record only
  reg.params.ououbm[,simIndex]<-prior.params$reg.params#for record only
  
  sim.trait <-ououbmmodel(model.params=prior.params$model.params,reg.params=prior.params$reg.params,root=root,tree=tree)
  sim.ououbm.trait[,1,simIndex]<-sim.trait$y
  sim.ououbm.trait[,2,simIndex]<-sim.trait$x1
  sim.ououbm.trait[,3,simIndex]<-sim.trait$x2
  y.sum.stat.ououbm[,simIndex]<- sum.stat(trait=sim.trait$y,tree=tree)
  x1.sum.stat.ououbm[,simIndex]<- sum.stat(trait=sim.trait$x1,tree=tree)
  x2.sum.stat.ououbm[,simIndex]<- sum.stat(trait=sim.trait$x2,tree=tree)
}#end of loop
sum.stat.ououbm <- cbind(t(y.sum.stat.ououbm),t(x1.sum.stat.ououbm),t(x2.sum.stat.ououbm))
ououbm.par.sim <- cbind(t(model.params.ououbm),t(reg.params.ououbm))

###############
### oubmcir ###
##############
prior.params.oubmcir <- oubmcirprior(prior.model.params = prior.model.params, prior.reg.params = prior.reg.params)
sim.oubmcir.trait<-array(0,c(dim(dataset)[1],3,sims))
model.params.oubmcir<-array(0,c(5,sims))
rownames(model.params.oubmcir)<-c("alpha.y","sigma.x","alpha.tau","theta.tau","sigma.tau")
reg.params.oubmcir<-array(0,c(3,sims))
row.names(reg.params.oubmcir)<-c("b0", "b1", "b2")
y.sum.stat.oubmcir<-array(0,c(2,sims))
rownames(y.sum.stat.oubmcir)<-c("y.mean","y.sd")
x1.sum.stat.oubmcir<-array(0,c(2,sims))
rownames(x1.sum.stat.oubmcir)<-c("x1.mean","x1.sd")
x2.sum.stat.oubmcir<-array(0,c(2,sims))
rownames(x2.sum.stat.oubmcir)<-c("x2.mean","x2.sd")

prior.params.oubmcir <- oubmcirprior(prior.model.params = prior.model.params, prior.reg.params = prior.reg.params)
for(simIndex in 1:sims){
  if(simIndex %%5000==0){print(paste("oubmcir",simIndex,sep=""))}
  prior.params<-oubmcirprior(prior.model.params = prior.model.params, prior.reg.params = prior.reg.params)
  model.params.oubmcir[,simIndex]<-prior.params$model.params#for record only
  reg.params.oubmcir[,simIndex]<-prior.params$reg.params#for record only
  
  sim.trait <-oubmcirmodel(model.params=prior.params$model.params,reg.params=prior.params$reg.params,root=root,tree=tree)
  sim.oubmcir.trait[,1,simIndex]<-sim.trait$y
  sim.oubmcir.trait[,2,simIndex]<-sim.trait$x1
  sim.oubmcir.trait[,3,simIndex]<-sim.trait$x2
  y.sum.stat.oubmcir[,simIndex]<- sum.stat(trait=sim.trait$y,tree=tree)
  x1.sum.stat.oubmcir[,simIndex]<- sum.stat(trait=sim.trait$x1,tree=tree)
  x2.sum.stat.oubmcir[,simIndex]<- sum.stat(trait=sim.trait$x2,tree=tree)
}# end of loop

sum.stat.oubmcir <- cbind(t(y.sum.stat.oubmcir),t(x1.sum.stat.oubmcir),t(x2.sum.stat.oubmcir))
oubmcir.par.sim <- cbind(t(model.params.oubmcir),t(reg.params.oubmcir))

###############
### ououcir ###
##############
prior.params.ououcir <- ououcirprior(prior.model.params = prior.model.params, prior.reg.params = prior.reg.params)
sim.ououcir.trait<-array(0,c(dim(dataset)[1],3,sims))
model.params.ououcir<-array(0,c(7,sims))
rownames(model.params.ououcir)<-c("alpha.y", "alpha.x", "theta.x", "sigma.x","alpha.tau","theta.tau","sigma.tau")
reg.params.ououcir<-array(0,c(3,sims))
row.names(reg.params.ououcir)<-c("b0", "b1", "b2")
y.sum.stat.ououcir<-array(0,c(2,sims))
rownames(y.sum.stat.ououcir)<-c("y.mean","y.sd")
x1.sum.stat.ououcir<-array(0,c(2,sims))
rownames(x1.sum.stat.ououcir)<-c("x1.mean","x1.sd")
x2.sum.stat.ououcir<-array(0,c(2,sims))
rownames(x2.sum.stat.ououcir)<-c("x2.mean","x2.sd")

prior.params.ououcir <- ououcirprior(prior.model.params = prior.model.params, prior.reg.params = prior.reg.params)
for(simIndex in 1:sims){
  if(simIndex %%5000==0){print(paste("ououcir",simIndex,sep=""))}
  prior.params<-ououcirprior(prior.model.params = prior.model.params, prior.reg.params = prior.reg.params)
  model.params.ououcir[,simIndex]<-prior.params$model.params#for record only
  reg.params.ououcir[,simIndex]<-prior.params$reg.params#for record only
  
  sim.trait <-ououcirmodel(model.params=prior.params$model.params,reg.params=prior.params$reg.params,root=root,tree=tree)
  sim.ououcir.trait[,1,simIndex]<-sim.trait$y
  sim.ououcir.trait[,2,simIndex]<-sim.trait$x1
  sim.ououcir.trait[,3,simIndex]<-sim.trait$x2
  y.sum.stat.ououcir[,simIndex]<- sum.stat(trait=sim.trait$y,tree=tree)
  x1.sum.stat.ououcir[,simIndex]<- sum.stat(trait=sim.trait$x1,tree=tree)
  x2.sum.stat.ououcir[,simIndex]<- sum.stat(trait=sim.trait$x2,tree=tree)
}# end of loop

sum.stat.ououcir <- cbind(t(y.sum.stat.ououcir),t(x1.sum.stat.ououcir),t(x2.sum.stat.ououcir))
ououcir.par.sim <- cbind(t(model.params.ououcir),t(reg.params.ououcir))


###############
### Analysis ##
###############
models<-rep(c("oubmbm","ououbm","oubmcir","ououcir"), each=sims)
full.sum.stat <- rbind(sum.stat.oubmbm, sum.stat.ououbm, sum.stat.oubmcir, sum.stat.ououcir)

#################
### rejection ###
#################

### posterior
rej.result.oubmbm <- abc(target = raw.sum.stat,param = data.frame(oubmbm.par.sim), sumstat = sum.stat.oubmbm,tol = 0.05,method="rejection")
rej.result.ououbm <- abc(target = raw.sum.stat,param = data.frame(ououbm.par.sim), sumstat = sum.stat.ououbm,tol = 0.05,method="rejection")
rej.result.oubmcir <- abc(target = raw.sum.stat,param = data.frame(oubmcir.par.sim), sumstat = sum.stat.oubmcir,tol = 0.05,method="rejection")
rej.result.ououcir <- abc(target = raw.sum.stat,param = data.frame(ououcir.par.sim), sumstat = sum.stat.ououcir,tol = 0.05,method="rejection")

post.oubmbm <- as.data.frame(rej.result.oubmbm$unadj.values)
post.ououbm <- as.data.frame(rej.result.ououbm$unadj.values)
post.oubmcir <- as.data.frame(rej.result.oubmcir$unadj.values)
post.ououcir <- as.data.frame(rej.result.ououcir$unadj.values)

#######################
### model selection ###
#######################
modsel <- postpr(target=raw.sum.stat,models,full.sum.stat, tol=0.05,method="rejection")



#######################
### goodness of fit ###
#######################

res.gfit.oubmbm=gfit(target=raw.sum.stat,sumstat = full.sum.stat[models=="oubmbm",],statistic = mean, nb.replicate = 100)
#summary(res.gfit.oubmbm)
res.gfit.ououbm=gfit(target=raw.sum.stat,sumstat = full.sum.stat[models=="ououbm",],statistic = mean, nb.replicate = 100)
#summary(res.gfit.ououbm)
res.gfit.oubmcir=gfit(target=raw.sum.stat,sumstat = full.sum.stat[models=="oubmcir",],statistic = mean, nb.replicate = 100)
#summary(res.gfit.oubmcir)
res.gfit.ououcir=gfit(target=raw.sum.stat,sumstat = full.sum.stat[models=="ououcir",],statistic = mean, nb.replicate = 100)
#summary(res.gfit.ououcir)

########################
### Cross vaildation ###
########################
### ABC ###
cv.result.oubmbm <- cv4abc(param=oubmbm.par.sim, sumstat=sum.stat.oubmbm, abc.out=rej.result.oubmbm, nval=5, tols=c(0.1,0.2,0.3))
#summary(cv.result.oubmbm)

cv.result.ououbm <- cv4abc(param=ououbm.par.sim, sumstat=sum.stat.ououbm, abc.out=rej.result.ououbm, nval=5, tols=c(0.1,0.2,0.3))
#summary(cv.result.ououbm)

cv.result.oubmcir <- cv4abc(param=oubmcir.par.sim, sumstat=sum.stat.oubmcir, abc.out=rej.result.oubmcir, nval=5, tols=c(0.1,0.2,0.3))
#summary(cv.result.oubmcir)

cv.result.ououcir <- cv4abc(param=ououcir.par.sim, sumstat=sum.stat.ououcir, abc.out=rej.result.ououcir, nval=5, tols=c(0.1,0.2,0.3))
#summary(cv.result.ououcir)

### Cross vaildation model selection
cv.modsel <- cv4postpr(models, full.sum.stat, nval=5, tol=0.01, method="rejection")
#summary(cv.modsel)

############
### PLoT ###
############

#path="/Users/TerryLai/Dropbox/TerryLai/treetraitdata/webster.purvis_foram/webster_jpeg/"
#######################
#Plot OUBMBM Posterior#
#######################
# oubmbm.params.names <- c("alpha.y","sigma.x","tau","b0","b1","b2") 
# oubmbm.xlab.names <- c(expression(alpha[y]),expression(sigma[x]^2),expression(tau),expression(b[0]),expression(b[1]),expression(b[2]))
# 
# for (i in 1:length(post.oubmbm)) {
#   jpeg(filename=paste("OUBMBM_",oubmbm.params.names[i],".jpeg",sep=""), width=400, height=400, units= "mm", res=500)
#   truehist(post.oubmbm[,i],col="gray",main="Histogram of posterior parameters",xlab=oubmbm.xlab.names[i])
#   lines(density(post.oubmbm[,i]),col="blue")
#   dev.off()
# }
#######################
#Plot OUOUBM Posterior#
#######################
# ououbm.params.names <- c("alpha.y","alpha.x","theta.x","sigma.x","tau","bo","b1","b2")
# ououbm.xlab.names <- c(expression(alpha[y]),expression(alpha[x]),expression(theta[x]),expression(sigma[x]^2),expression(tau),expression(b[0]),expression(b[1]),expression(b[2]))
# 
# for (i in 1:length(post.ououbm)){
#   jpeg(filename = paste("OUOUBM_",ououbm.params.names[i],".jpeg",sep=""), width=400, height=400, units= "mm", res=500)
#   truehist(post.ououbm[,i],col="gray",main="Histogram of posterior parameters",xlab=ououbm.xlab.names[i])
#   lines(density(post.ououbm[,i]),col="blue")
#   dev.off()
# }
########################
#Plot OUBMCIR Posterior#
########################
# oubmcir.params.names <- c("alpha.y","sigma.x","alpha.tau","theta.tau","sigma.tau","bo","b1","b2")
# oubmcir.xlab.names <- c(expression(alpha[y]),expression(sigma[x]^2),expression(alpha[tau]),expression(theta[tau]),expression(sigma[tau]),expression(b[0]),expression(b[1]),expression(b[2]))
# 
# for(i in 1:length(post.oubmcir)){
#   jpeg(filename = paste("OUBMCIR_",oubmcir.params.names[i],".jpeg",sep=""), width=400, height=400, units= "mm", res=500)
#   truehist(post.oubmcir[,i],col="gray",main="Histogram of posterior parameters",xlab=oubmcir.xlab.names[i])
#   lines(density(post.oubmcir[,i]),col="blue")
#   dev.off()
# }
########################
#Plot OUOUCIR Posterior#
######################## 
# ououcir.params.names <- c("alpha.y","alpha.x","theta.x","sigma.x","alpha.tau","theta.tau","sigma.tau","bo","b1","b2")
# ououcir.xlab.names <- c(expression(alpha[y]),expression(alpha[x]),expression(theta[x]),expression(sigma[x]^2),expression(alpha[tau]),expression(theta[tau]),expression(sigma[tau]),expression(b[0]),expression(b[1]),expression(b[2]))
# for(i in 1:length(post.ououcir)){
#   jpeg(filename = paste("OUOUCIR_",ououcir.params.names[i],".jpeg",sep=""), width=400, height=400, units= "mm", res=500)
#   truehist(post.ououcir[,i],col="gray",main="Histogram of posterior parameters",xlab=ououcir.xlab.names[i])
#   lines(density(post.ououcir[,i]),col="blue")
#   dev.off()
# }
########################
### sum.stat.boxplot ###
########################

# boxplot.names <- c("y.mean","y.sd","x1.mean","x1.sd","x2.mean","x2.sd")
# boxplot.mains <- c("Mean of Y","Standard Deviation of Y","Mean of X1","Standard Deviation of X1","Mean of X2","Standard Deviation of X2")
# for(i in 1:length(raw.sum.stat)){
#   jpeg(filename = paste("boxplot_",i,".jpeg",sep=""),width=400,height=400,units="mm",res=500)
#   boxplot(full.sum.stat[,boxplot.names[i]]~models,main=boxplot.mains[i],col="gray")
#   abline(h=raw.sum.stat[i],col="blue")
#   dev.off()
# }

############################
### plot goodness of fit ###
############################
# {jpeg(filename = paste("gfit.oubmbm.jpeg",sep=""),width=400,height=400,units="mm",res=500)
#   plot(res.gfit.oubmbm,col="gray")
#   dev.off()
#   jpeg(filename = paste("gfit.oubmcir.jpeg",sep=""),width=400,height=400,units="mm",res=500)
#   plot(res.gfit.oubmcir,col="gray")
#   dev.off()
#   jpeg(filename = paste("gfit.ououbm.jpeg",sep=""),width=400,height=400,units="mm",res=500)
#   plot(res.gfit.ououbm,col="gray")
#   dev.off()
#   jpeg(filename = paste("gfit.ououcir.jpeg",sep=""),width=400,height=400,units="mm",res=500)
#   plot(res.gfit.ououcir,col="gray")
#   dev.off()}

########################
### Confusion matrix ###
########################

# {jpeg(filename = paste("Confusion_matrix.jpeg",sep=""),width=400,height=400,units="mm",res=500)
#   plot(cv.modsel, names.arg=c("OUBMBM", "OUBMCIR", "OUOUBM", "OUOUCIR"))
#   dev.off()}

