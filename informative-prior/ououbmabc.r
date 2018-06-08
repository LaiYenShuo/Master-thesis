# rm(list=ls())
# #install.packages("TreeSim")
library(TreeSim)
library(EasyABC)
library(coda)
library(Sim.DiffProc)
library(MCMCpack)
library(ape)
library(abc)
#oubmbmintegrand<-function(s, alpha.y=alpha.y){
#    alpha.y*exp(alpha.y*s)*rnorm(n=1,mean=0,sd=s)
#    }

ououbmmodel<-function(model.params,reg.params,root=root,tree=tree){
  alpha.y<<-model.params[1]
  alpha.x<-model.params[2]
  theta.x<-model.params[3]
  sigma.x<-model.params[4]
  tau<-model.params[5]

  b0<-reg.params[1]
  b1<-reg.params[2]
  b2<-reg.params[3]

  n<-Ntip(tree)
  x1nodestates<-array(0,c(2*n-1))
  x1nodestates[n+1]<-root
  x2nodestates<-array(0,c(2*n-1))
  x2nodestates[n+1]<-root
  optimnodestates<-array(0,c(2*n-1))
  optimnodestates[n+1]<-root
  sigmasqnodestates<-array(0,c(2*n-1))
  sigmasqnodestates[n+1]<-root
  ynodestates<-array(0,c(2*n-1))
  ynodestates[n+1]<-root

  N<-dim(tree$edge)[1]
  anc<-tree$edge[,1]
  des<-tree$edge[,2]
  treelength<-tree$edge.length

  for(index in N:1){

    x1ou.mean<-x1nodestates[anc[index]]*exp(-alpha.x*treelength[index]) + theta.x*(1-exp(-alpha.x*treelength[index]))
    x1ou.sd<-sqrt((sigma.x^2/(2*alpha.x))*(1-exp(-2*alpha.x*treelength[index])))
    x1nodestates[des[index]]<-rnorm(n=1,mean=x1ou.mean,sd=x1ou.sd)

    x2ou.mean<-x2nodestates[anc[index]]*exp(-alpha.x*treelength[index]) + theta.x*(1-exp(-alpha.x*treelength[index]))
    x2ou.sd<-sqrt((sigma.x^2/(2*alpha.x))*(1-exp(-2*alpha.x*treelength[index])))
    x2nodestates[des[index]]<-rnorm(n=1,mean=x2ou.mean,sd=x2ou.sd)

    optimnodestates[des[index]]<- b0 + b1*x1nodestates[des[index]] + b2*x2nodestates[des[index]]

    sigmasqnodestates[des[index]]<-rnorm(n=1,mean=sigmasqnodestates[anc[index]],sd=sqrt(tau*treelength[index]))
    sigma.sq.theta<- b1^2*sigma.x^2 + b2^2*sigma.x^2
    #inttheta<-integrate(oubmbmintegrand,lower=0 ,upper=treelength[index], alpha.y=true.alpha.y)

    A1<-(alpha.y*optimnodestates[des[index]]/ (alpha.y-alpha.x)) *(exp((alpha.y-alpha.x)*treelength[index]) -1)
    theta1<-0
    A2<- theta1*(exp(alpha.y*treelength[index])-1) - (alpha.y*theta1/(alpha.y-alpha.x))*(exp((alpha.y-alpha.x)*treelength[index]) -1)
    A3<-1/(sigma.sq.theta*alpha.y^2)*rnorm(1,mean=0,sd=sqrt((sigma.sq.theta*alpha.y^2)^3 *treelength[index]^3 /3))
    INTtime<- A1+A2+A3
    INT1<-exp(-alpha.y*treelength[index])*INTtime

    fexpr<-expression(alpha.y*exp(alpha.y*t)*w)
    res<-st.int(fexpr,type="ito",M=1,lower=0,upper=treelength[index])
    INT2<-median(res$X)

    ynodestates[des[index]]<-ynodestates[anc[index]] + INT1 + INT2
    }

  simtrait<-ynodestates[1:n]
  return(list(y=simtrait,x1=x1nodestates[1:n],x2=x2nodestates[1:n]))
  #return(c(mean(simtrait),sd(simtrait)))
  }

#print(model(model.params=c(5,1,1,3,2),reg.params=c(1,1,1),root=root,tree=tree))

ououbmprior<-function(prior.model.params=prior.model.params,prior.reg.params=prior.reg.params){
  alpha.y.rate <-prior.model.params["alpha.y.rate"]
  alpha.x.rate <- prior.model.params["alpha.x.rate"]
  theta.x.mean <-prior.model.params["theta.x.mean"]
  theta.x.sd<-prior.model.params["theta.x.sd"]
  sigma.x.rate <-prior.model.params["sigma.x.rate"]
  tau.rate <- prior.model.params["tau.rate"]


  alpha.y<-rexp(n=1, rate=alpha.y.rate)
  alpha.x<-rexp(n=1,rate=alpha.x.rate)
  theta.x<-rnorm(n=1,mean=theta.x.mean,sd=theta.x.sd)
  sigma.x<-rexp(n=1, rate=sigma.x.rate)
  tau<-rexp(n=1, rate=tau.rate)

  b0.min<-prior.reg.params[1]
  b0.max<-prior.reg.params[2]
  b1.min<-prior.reg.params[3]
  b1.max<-prior.reg.params[4]
  b2.min<-prior.reg.params[5]
  b2.max<-prior.reg.params[6]
  b0<-runif(n=1, min=b0.min, max=b0.max)
  b1<-runif(n=1, min=b1.min, max=b1.max)
  b2<-runif(n=1, min=b2.min, max=b2.max)

  model.params<-c(alpha.y, alpha.x, theta.x, sigma.x, tau)
  reg.params<- c(b0, b1, b2)

  return(list(model.params=model.params, reg.params=reg.params))
}

sum.stat<-function(trait=trait,tree=tree){
  pic.trait<-pic(x=trait,phy=tree)
  return(c(mean(pic.trait),sd(pic.trait)))
}


sum.stat.distance<-function(raw.sum.stat=raw.sum.stat,sim.sum.stat=sim.sum.stat){
  raw.sum.stat <- raw.sum.stat/median(raw.sum.stat-median(raw.sum.stat))
  sim.sum.stat <- sim.sum.stat/median(sim.sum.stat-median(sim.sum.stat))
  return(sum((raw.sum.stat-sim.sum.stat)^2))
}

# ## main
# 
# n<-20
# numbsim<-1;lambda<-2.0;mu<-0.5;frac<-0.6;age<-2
# tree<-sim.bd.taxa.age(n=n,numbsim=1,lambda=lambda,mu=mu,frac=frac,age=age,mrca=TRUE)[[1]]
# tree<-reorder(tree,"postorder")
# tree$edge
# plot(tree)
# nodelabels()
# tiplabels()
# 
# root<-0
# true.alpha.y<-1
# true.alpha.x<-0.2
# true.theta.x<-0.5
# true.sigma.x<-2
# true.tau<- 1
# 
# true.b0 <- 0
# true.b1 <- 1
# true.b2 <- 0.2
# 
# #hyper parameters
# alpha.y.rate <- 5
# alpha.x.rate <- 5
# theta.x.mean <- 0
# theta.x.sd  <- 1
# sigma.x.shape <-2
# sigma.x.scale <-1
# tau.shape <- 3
# tau.scale <- 1
# b0.min = -5
# b0.max = 5
# b1.min = -5
# b1.max = 5
# b2.min = -5
# b2.max = 5
# 
# prior.model.params=c(alpha.y.rate, alpha.x.rate, theta.x.mean, theta.x.sd, sigma.x.shape, sigma.x.scale, tau.shape, tau.scale)
# names(prior.model.params)<-c("alpha.y.rate", "alpha.x.rate", "theta.x.mean","theta.x.sd","sigma.x.shape", "sigma.x.scale", "tau.shape", "tau.scale")
# prior.reg.params=c(b0.min, b0.max, b1.min, b1.max, b2.min, b2.max)
# 
# prior.params <- ououbmprior(prior.model.params=prior.model.params,prior.reg.params=prior.reg.params)
# 
# true.trait<-ououbmmodel(model.params=c(true.alpha.y, true.alpha.x, true.theta.x, true.sigma.x, true.tau),reg.params=c(true.b0,true.b1,true.b2),root=root,tree=tree)
# 
# sim.trait <-ououbmmodel(model.params=prior.params$model.params,reg.params=prior.params$reg.params,root=root,tree=tree)
# true.trait
# 
# 
# raw.sum.stat.y<-sum.stat(trait=true.trait$y,tree=tree)
# raw.sum.stat.x1<-sum.stat(trait=true.trait$x1,tree=tree)
# raw.sum.stat.x2<-sum.stat(trait=true.trait$x2,tree=tree)
# 
# 
# sims<-50000
# sim.ououbm.trait<-array(0,c(n,5,sims))
# model.params.array<-array(0,c(5,sims))
# rownames(model.params.array)<-c("alpha.y","alpha.x","theta.x","sigma.x","tau")
# reg.params.array<-array(0,c(3,sims))
# row.names(reg.params.array)<-c("b0", "b1", "b2")
# y.sum.stat.array<-array(0,c(2,sims))
# rownames(y.sum.stat.array)<-c("y.mean","y.sd")
# x1.sum.stat.array<-array(0,c(2,sims))
# rownames(x1.sum.stat.array)<-c("x1.mean","x1.sd")
# x2.sum.stat.array<-array(0,c(2,sims))
# rownames(x2.sum.stat.array)<-c("x2.mean","x2.sd")
# post.model.params.array<-array(0,c(5,sims))
# rownames(post.model.params.array)<-c("alpha.y","alpha.x","theta.x","sigma.x","tau")
# post.reg.params.array<-array(0,c(3,sims))
# row.names(post.reg.params.array)<-c("b0", "b1", "b2")
# 
# prior.params <- ououbmprior(prior.model.params=prior.model.params,prior.reg.params=prior.reg.params)
# post.model.params.array[,1]<-prior.params$model.params
# post.reg.params.array[,1]<-prior.params$reg.params
# sum.stat.distance.array<-array(0,c(sims))
# for(simIndex in 1:sims){
#   if(simIndex %%1000==0){print(simIndex)}
#   prior.params <- ououbmprior(prior.model.params=prior.model.params,prior.reg.params=prior.reg.params)
#   model.params.array[,simIndex]<-prior.params$model.params#for record only
#   reg.params.array[,simIndex]<-prior.params$reg.params#for record only
# 
#   sim.trait <-ououbmmodel(model.params=prior.params$model.params,reg.params=prior.params$reg.params,root=root,tree=tree)
#   sim.ououbm.trait[,1,simIndex]<-sim.trait$y
#   sim.ououbm.trait[,2,simIndex]<-sim.trait$x1
#   sim.ououbm.trait[,3,simIndex]<-sim.trait$x2
#   y.sum.stat.array[,simIndex]<- sum.stat(trait=sim.trait$y,tree=tree)
#   x1.sum.stat.array[,simIndex]<- sum.stat(trait=sim.trait$x1,tree=tree)
#   x2.sum.stat.array[,simIndex]<- sum.stat(trait=sim.trait$x2,tree=tree)
#   }#end of loop
# 
# 
# sim.sum.stat <- cbind(t(y.sum.stat.array),t(x1.sum.stat.array),t(x2.sum.stat.array))
# head(sim.sum.stat)
# ououbm.par.sim <- cbind(t(model.params.array),t(reg.params.array))
# setwd("/Users/TerryLai/Dropbox/TerryLai/R_code/abc_V2")
# #setwd("~/Dropbox/FCU/Teaching/Mentoring/2017Spring/TerryLai/R_code/abc_V2/")
# save.image("test_ououbm.RData")
# 
# #ououbm.post.sim <- cbind(t(post.model.params.array),t(post.reg.params.array))
# rej <- abc(target=c(raw.sum.stat.y,raw.sum.stat.x1,raw.sum.stat.x2), param=ououbm.par.sim, sumstat=sim.sum.stat, tol=0.05, method="rejection"  )
# rej$unadj.values
# rej$ss
# rej$transf
# rej$logit.bounds
# rej$numparam
# rej$numstat
# rej$names
# 
# lin <- abc(target=c(raw.sum.stat.y,raw.sum.stat.x1,raw.sum.stat.x2), param=ououbm.par.sim, sumstat=sim.sum.stat, tol=0.2,hcorr=FALSE, method="loclinear", transf=c(rep("none",rej$numparam)))
# ls(lin)
# ls(lin)
# lin[c(1:18)]
# linhc <- abc(target=c(raw.sum.stat.y,raw.sum.stat.x1,raw.sum.stat.x2), param=ououbm.par.sim, sumstat=sim.sum.stat, tol=0.2, method="loclinear", transf=c(rep("none",rej$numparam)))
# 
# ls(lin)
# ls(linhc)
# linhc[c(1:18)]
# head(lin$adj.values);head(linhc$adj.values)
# ## posterior summaries
# ##
# linhc$adj.values
# linsum <- summary(linhc, intvl=0.9)
# linsum
# ## compare with the rejection sampling
# summary(linhc, unadj=TRUE, intvl=0.9)
# 
# 
# hist(linhc, breaks=30, caption=c( expression(alpha[y]),  expression(alpha[x]),expression(theta[x]),
#   expression(sigma[x]),expression(tau),expression(b[0]),expression(b[1]),
#   expression(b[2])))
# #  ,file = "~/Desktop/linhc")
# 
