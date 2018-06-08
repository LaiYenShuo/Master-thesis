library(TreeSim)
library(EasyABC)
library(coda)
library(Sim.DiffProc)
library(MCMCpack)
library(ape)
library(abc)
oubmcirmodel<-function(model.params,reg.params,root=root,tree=tree){
  alpha.y<-model.params[1]
  sigma.x<-model.params[2]
  alpha.tau<-model.params[3]
  theta.tau<-model.params[4]
  sigma.tau<-model.params[5]
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
  sigmasqnodestates[n+1]<-1
  ynodestates<-array(0,c(2*n-1))
  ynodestates[n+1]<-root

  N<-dim(tree$edge)[1]
  anc<-tree$edge[,1]
  des<-tree$edge[,2]
  treelength<-tree$edge.length
  for(index in N:1){

    x1nodestates[des[index]]<-rnorm(n=1,mean=x1nodestates[anc[index]],sd= sqrt(sigma.x^2*treelength[index]))
    x2nodestates[des[index]]<-rnorm(n=1,mean=x2nodestates[anc[index]],sd= sqrt(sigma.x^2*treelength[index]))
    optimnodestates[des[index]]<- b0 + b1*x1nodestates[des[index]] + b2*x2nodestates[des[index]]

    c<- sigma.tau^2*(1-exp(-alpha.tau*treelength[index]))/(4*alpha.tau)
    k<- (4*theta.tau*alpha.tau)/sigma.tau^2
    lambda<-4*alpha.tau*exp(-alpha.tau*treelength[index])/(sigma.tau^2*(1-exp(-alpha.tau*treelength[index])))*sigmasqnodestates[anc[index]]
    tmp = rchisq(n=1, df=k, ncp = lambda)
    sig_u <- c*tmp
    sigmasqnodestates[des[index]]<-sig_u

    sigma.sq.theta<- b1^2*sigma.x^2 + b2^2*sigma.x^2
    #inttheta<-integrate(oubmbmintegrand,lower=0 ,upper=treelength[index], alpha.y=true.alpha.y)

    INT1var<-treelength[index]*exp(2*alpha.y*treelength[index]) - 2*(exp(2*alpha.y*treelength[index])-exp(alpha.y*treelength[index])) + alpha.y/2*(exp(2*alpha.y*treelength[index])-1)
    INT1<-exp(-alpha.y*treelength[index])* rnorm(n=1,mean=0,sd=sqrt(abs(INT1var)))

    a <- rnorm(n=1, mean=0, sd=sqrt(theta.tau^2*(exp(2*alpha.y*treelength[index])-1)/(2*alpha.y)))
    b <- rnorm(n=1, mean=0, sd=sqrt(((sigmasqnodestates[des[index]]-theta.tau)^2/(2*(alpha.y-alpha.tau)))*(exp(2*(alpha.y-alpha.tau)*treelength[index])-1)))
n_t <- 20
n_s <- 20
    outer.int.sum=0
    for(outer.index in 1:n_t){
      inner.int.sum = 0
      for(inner.index in 1:n_s){
        c<- sigma.tau^2*(1-exp(-alpha.tau*(inner.index/n_s)))/(4*alpha.tau)
        k<- (4*theta.tau*alpha.tau)/sigma.tau^2
        lambda<- 4*alpha.tau*exp(-alpha.tau*(inner.index/n_s))/(sigma.tau^2*(1-exp(-alpha.tau*(inner.index/n_s))))*sigmasqnodestates[des[index]]
        tmp = rchisq(n=1, df=k, ncp = lambda)
        sig_u <- c*tmp
        inner.int.sum  <-  inner.int.sum + exp(alpha.tau*(inner.index/n_s))*rnorm(n=1,mean=0, sd=sqrt(1/n_s))*sqrt(sig_u)
      }
      outer.int.sum <- outer.int.sum + exp((alpha.y-alpha.tau)*(outer.index/n_t))*inner.int.sum*rnorm(n=1,mean=0, sd=sqrt(1/n_t))
    }
    outer.int.sum
    c <- sigma.tau*outer.int.sum

    INTsigdWy <- exp(-alpha.y*treelength[index])*(a + b + c)

    ynodestates[des[index]]<-ynodestates[anc[index]] + INT1 + INTsigdWy
  }
  simtrait<-ynodestates[1:n]
  return(list(y=simtrait,x1=x1nodestates[1:n],x2=x2nodestates[1:n]))
  }

oubmcirprior<-function(prior.model.params=prior.model.params,prior.reg.params=prior.reg.params){
  alpha.y.min <-prior.model.params["alpha.y.min"]# assume uniform
  alpha.y.max <-prior.model.params["alpha.y.max"] 
  sigma.x.min <-prior.model.params["sigma.x.min"]# 
  sigma.x.max <-prior.model.params["sigma.x.max"]
  alpha.tau.min <- prior.model.params["alpha.tau.min"]# 
  alpha.tau.max <- prior.model.params["alpha.tau.max"]
  theta.tau.a <- prior.model.params["theta.tau.a"] #assume uniform
  theta.tau.b <- prior.model.params["theta.tau.b"]
  sigma.tau.min <- prior.model.params["sigma.tau.min"]# 
  sigma.tau.max <- prior.model.params["sigma.tau.max"]

  alpha.y<-runif(n=1, min = alpha.y.min, max = alpha.y.max)
  sigma.x<-runif(n=1, min = sigma.x.min, max = sigma.x.max)
  alpha.tau <- runif(n=1, min = alpha.tau.min, max = alpha.tau.max)
  theta.tau <- runif(n=1, min = theta.tau.a, max=theta.tau.b)
  sigma.tau<-runif(n=1, min = sigma.tau.min, max = sigma.tau.max)

  b0.min<-prior.reg.params[1]
  b0.max<-prior.reg.params[2]
  b1.min<-prior.reg.params[3]
  b1.max<-prior.reg.params[4]
  b2.min<-prior.reg.params[5]
  b2.max<-prior.reg.params[6]
  b0<-runif(n=1, min=b0.min, max=b0.max)
  b1<-runif(n=1, min=b1.min, max=b1.max)
  b2<-runif(n=1, min=b2.min, max=b2.max)

  model.params<-c(alpha.y, sigma.x, alpha.tau, theta.tau, sigma.tau)
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
# ### main
# n<-20
# numbsim<-1;lambda<-2.0;mu<-0.5;frac<-0.6;age<-2
# tree<-sim.bd.taxa.age(n=n,numbsim=1,lambda=lambda,mu=mu,frac=frac,age=age,mrca=TRUE)[[1]]
# tree<-reorder(tree,"postorder")
# tree$edge
# plot(tree)
# nodelabels()
# tiplabels()
# 
# root<-1
# true.alpha.y<-1
# true.sigma.x<-2
# true.alpha.tau<-0.5
# true.theta.tau<-3
# true.sigma.tau<- 1
# n_s = 10
# n_t = 10
# true.b0 <- 0
# true.b1 <- 1
# true.b2 <- 0.2
# # hyper paramters
# alpha.y.rate <-5# assume exponential
# sigma.x.shape <-2# assume invgamma
# sigma.x.scale <-1
# alpha.tau.rate <- 5# assume exponential
# theta.tau.a <- 0 #assume uniform
# theta.tau.b <- 100
# sigma.tau.shape <- 2# assume inv gamma
# sigma.tau.scale <- 1
# b0.min=-5
# b0.max=5
# b1.min=-5
# b1.max=5
# b2.min=-5
# b2.max=5
# 
# prior.model.params=c(alpha.y.rate, sigma.x.shape, sigma.x.scale, alpha.tau.rate, theta.tau.a, theta.tau.b, sigma.tau.shape, sigma.tau.scale)
# names(prior.model.params)<-c("alpha.y.rate", "sigma.x.shape", "sigma.x.scale", "alpha.tau.rate", "theta.tau.a", "theta.tau.b", "sigma.tau.shape", "sigma.tau.scale")
# prior.reg.params=c(b0.min, b0.max, b1.min, b1.max, b2.min, b2.max)
# prior.params <- oubmcirprior(prior.model.params = prior.model.params, prior.reg.params=prior.reg.params)
# 
# true.trait <- oubmcirmodel(model.params = c(true.alpha.y, true.sigma.x, true.alpha.tau, true.theta.tau, true.sigma.tau),reg.params = c(true.b0,true.b1,true.b2), root=root, tree=tree)
# sim.trait <- oubmcirmodel(model.params = prior.params$model.params, reg.params = prior.params$reg.params, root=root, tree=tree)
# 
# raw.sum.stat.y <- sum.stat(trait = true.trait$y, tree=tree)
# raw.sum.stat.x1 <- sum.stat(trait = true.trait$x1, tree=tree)
# raw.sum.stat.x2 <- sum.stat(trait = true.trait$x2, tree=tree)
# 
# sims=50000
# sim.oubmcir.trait<-array(0,c(n,3,sims))
# model.params.array<-array(0,c(5,sims))
# rownames(model.params.array)<-c("alpha.y","sigma.x","alpha.tau","theta.tau","sigma.tau")
# reg.params.array<-array(0,c(3,sims))
# row.names(reg.params.array)<-c("b0", "b1", "b2")
# y.sum.stat.array<-array(0,c(2,sims))
# rownames(y.sum.stat.array)<-c("y.mean","y.sd")
# x1.sum.stat.array<-array(0,c(2,sims))
# rownames(x1.sum.stat.array)<-c("x1.mean","x1.sd")
# x2.sum.stat.array<-array(0,c(2,sims))
# rownames(x2.sum.stat.array)<-c("x2.mean","x2.sd")
# 
# 
# prior.params <- oubmcirprior(prior.model.params = prior.model.params, prior.reg.params = prior.reg.params)
# sum.stat.distance.array<-array(0,c(sims))
# for(simIndex in 1:sims){
#   if(simIndex %%1000==0){print(simIndex)}
#   prior.params<-oubmcirprior(prior.model.params = prior.model.params, prior.reg.params = prior.reg.params)
#   model.params.array[,simIndex]<-prior.params$model.params#for record only
#   reg.params.array[,simIndex]<-prior.params$reg.params#for record only
# 
#   sim.trait <-oubmcirmodel(model.params=prior.params$model.params,reg.params=prior.params$reg.params,root=root,tree=tree)
#   sim.oubmcir.trait[,1,simIndex]<-sim.trait$y
#   sim.oubmcir.trait[,2,simIndex]<-sim.trait$x1
#   sim.oubmcir.trait[,3,simIndex]<-sim.trait$x2
#   y.sum.stat.array[,simIndex]<- sum.stat(trait=sim.trait$y,tree=tree)
#   x1.sum.stat.array[,simIndex]<- sum.stat(trait=sim.trait$x1,tree=tree)
#   x2.sum.stat.array[,simIndex]<- sum.stat(trait=sim.trait$x2,tree=tree)
# }# end of loop
# 
# ### Use abc package
# sim.sum.stat <- cbind(t(y.sum.stat.array),t(x1.sum.stat.array),t(x2.sum.stat.array))
# oubmcir.par.sim <- cbind(t(model.params.array),t(reg.params.array))
# setwd("/Users/TerryLai/Dropbox/TerryLai/R_code/abc_V2")
# #setwd("~/Dropbox/FCU/Teaching/Mentoring/2017Spring/TerryLai/R_code/abc_V2/")
# save.image("test_oubmcir.RData")
# ### The rejection alogortim
# rej <- abc(target=c(raw.sum.stat.y,raw.sum.stat.x1,raw.sum.stat.x2), param=oubmcir.par.sim, sumstat=sim.sum.stat, tol=0.05, method="rejection"  )
# ls(rej)
# rej[1:12]
# ## ABC with local linear regression correction without/with correction
# ## for heteroscedasticity
# ##
# 
# lin <- abc(target=c(raw.sum.stat.y,raw.sum.stat.x1,raw.sum.stat.x2), param=oubmcir.par.sim, sumstat=sim.sum.stat, tol=0.2,hcorr=FALSE, method="loclinear",transf = c("log","log","log","log","log","log","log","log"))
# ls(lin)
# lin[c(1:18)]
# linhc <- abc(target=c(raw.sum.stat.y,raw.sum.stat.x1,raw.sum.stat.x2), param=oubmcir.par.sim, sumstat=sim.sum.stat, tol=0.2, method="loclinear", transf=c("log","log","log","log","log","log","log","log"))
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
# ## posterior histograms
# ##
# ## or send histgrams to pdf file
# hist(linhc, breaks=30, caption=c(
#   expression(alpha[y]),
#   expression(sigma[x]),
#   expression(alpha[tau]),
#   expression(theta[tau]),
#   expression(sigma[tau]),
#   expression(b[0]),
#   expression(b[1]),
#   expression(b[2])))
# 
