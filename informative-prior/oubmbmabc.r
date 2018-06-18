# rm(list=ls())
# #install.packages("TreeSim")
library(TreeSim)
library(EasyABC)
library(coda)
library(Sim.DiffProc)
library(MCMCpack)
library(ape)
library(abc)

oubmbmmodel<-function(model.params,reg.params,root=root,tree=tree){
  alpha.y<<-model.params[1]
  sigma.x<-model.params[2]
  tau<-model.params[3]
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
    x1nodestates[des[index]]<-rnorm(n=1,mean=x1nodestates[anc[index]],sd= sqrt(sigma.x^2*treelength[index]))
    x2nodestates[des[index]]<-rnorm(n=1,mean=x2nodestates[anc[index]],sd= sqrt(sigma.x^2*treelength[index]))
    optimnodestates[des[index]]<- b0 + b1*x1nodestates[des[index]] + b2*x2nodestates[des[index]]
    sigmasqnodestates[des[index]]<-rnorm(n=1,mean=sigmasqnodestates[anc[index]],sd=sqrt(tau*treelength[index]))
    sigma.sq.theta<- b1^2*sigma.x^2 + b2^2*sigma.x^2
    
    INT1mean<- optimnodestates[des[index]]*exp(alpha.y*treelength[index])  - optimnodestates[anc[index]]
    INT1var<-  (exp(2*alpha.y*treelength[index])-1)/ (2*alpha.y)
    
    INT1<-exp(-alpha.y*treelength[index])* rnorm(n=1,mean= INT1mean, sd=sqrt(abs(INT1var)))
      
    fexpr<-expression(alpha.y*exp(alpha.y*t)*w)
    res<-st.int(fexpr,type="ito",M=1,lower=0,upper=treelength[index])
    INT2<-median(res$X)
    ynodestates[des[index]]<-ynodestates[anc[index]] + INT1 + INT2
  }
  
  simtrait<-ynodestates[1:n]
  
  return(list(y=simtrait,x1=x1nodestates[1:n],x2=x2nodestates[1:n]))

}

oubmbmprior<-function(prior.model.params=prior.model.params,prior.reg.params=prior.reg.params){
  alpha.y.rate <-prior.model.params["alpha.y.rate"]
  sigma.x.rate <-prior.model.params["sigma.x.rate"]
  tau.rate <- prior.model.params["tau.rate"]
  
  alpha.y<-rexp(n=1, rate=alpha.y.rate)
  sigma.x<-rexp(n=1, rate = sigma.x.rate)
  tau<-rexp(n=1, rate = tau.rate)
  
  b0.min<-prior.reg.params[1]
  b0.max<-prior.reg.params[2]
  b1.min<-prior.reg.params[3]
  b1.max<-prior.reg.params[4]
  b2.min<-prior.reg.params[5]
  b2.max<-prior.reg.params[6]
  b0<-runif(n=1, min=b0.min, max=b0.max)
  b1<-runif(n=1, min=b1.min, max=b1.max)
  b2<-runif(n=1, min=b2.min, max=b2.max)
  
  model.params<-c(alpha.y, sigma.x, tau)
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

