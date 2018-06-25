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
