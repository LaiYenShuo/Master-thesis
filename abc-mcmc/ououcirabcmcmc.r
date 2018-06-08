# rm(list=ls())
# #install.packages("TreeSim")
library(TreeSim)
library(EasyABC)
library(coda)
library(Sim.DiffProc)
library(MCMCpack)
library(ape)
library(abc)

ououcirmodel<-function(model.params,reg.params,root=root,tree=tree){
  alpha.y<-model.params[1]
  alpha.x<-model.params[2]
  theta.x<-model.params[3]
  sigma.x<-model.params[4]
  alpha.tau<-model.params[5]
  theta.tau<-model.params[6]
  sigma.tau<-model.params[7]
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
    
    c<- sigma.tau^2*(1-exp(-alpha.tau*treelength[index]))/(4*alpha.tau)
    k<- (4*theta.tau*alpha.tau)/sigma.tau^2
    lambda<-4*alpha.tau*exp(-alpha.tau*treelength[index])/(sigma.tau^2*(1-exp(-alpha.tau*treelength[index])))*sigmasqnodestates[anc[index]]
    tmp = rchisq(n=1, df=k, ncp = lambda)
    sig_u <- c*tmp
    sigmasqnodestates[des[index]]<-sig_u
    
    sigma.sq.theta<- b1^2*sigma.x^2 + b2^2*sigma.x^2
    
    A1<-  (alpha.y*optimnodestates[des[index]]/ (alpha.y-alpha.x)) *(exp((alpha.y-alpha.x)*treelength[index]) -1)
    theta1<-0
    A2<- theta1*(exp(alpha.y*treelength[index])-1) - (alpha.y*theta1/(alpha.y-alpha.x))*(exp((alpha.y-alpha.x)*treelength[index]) -1)
    A3<-1/(sigma.sq.theta*alpha.y^2)*rnorm(1,mean=0,sd=sqrt((sigma.sq.theta*alpha.y^2)^3 *treelength[index]^3 /3))
    INTtime<- A1+A2+A3
    INT1<-exp(-alpha.y*treelength[index])*INTtime
    
    a <- rnorm(n=1, mean=0, sd=sqrt(theta.tau^2*(exp(2*alpha.y*treelength[index])-1)/(2*alpha.y)))
    b <- rnorm(n=1, mean=0, sd=sqrt(((sigmasqnodestates[des[index]]-theta.tau)^2/(2*(alpha.y-alpha.tau)))*(exp(2*(alpha.y-alpha.tau)*treelength[index])-1)))
    # n_t=30
    # n_s=30
    n_t=25
    n_s=25
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

ououcirprior <- function(prior.model.params=prior.model.params,prior.reg.params=prior.reg.params){
  alpha.y.min <- prior.model.params["alpha.y.min"]
  alpha.y.max <- prior.model.params["alpha.y.max"]
  alpha.x.min <- prior.model.params["alpha.x.min"]
  alpha.x.max <- prior.model.params["alpha.x.max"]
  theta.x.min <-prior.model.params["theta.x.min"]
  theta.x.max <-prior.model.params["theta.x.max"]
  sigma.x.min <-prior.model.params["sigma.x.min"]
  sigma.x.max <-prior.model.params["sigma.x.max"]
  alpha.tau.min <- prior.model.params["alpha.tau.min"]# 
  alpha.tau.max <- prior.model.params["alpha.tau.max"]
  theta.tau.min <- prior.model.params["theta.tau.min"] #assume uniform
  theta.tau.max <- prior.model.params["theta.tau.max"]
  sigma.tau.min <- prior.model.params["sigma.tau.min"]# 
  sigma.tau.max <- prior.model.params["sigma.tau.max"]
  
  alpha.y<-runif(n=1, min = alpha.y.min, max = alpha.y.max)
  alpha.x<-runif(n=1, min = alpha.x.min, max = alpha.x.max)
  theta.x<-runif(n=1, min = theta.x.min, max = theta.x.max)
  sigma.x<-runif(n=1, min = sigma.x.min, max = sigma.x.max)
  alpha.tau <- runif(n=1, min = alpha.tau.min, max = alpha.tau.max)
  theta.tau <- runif(n=1, min=theta.tau.min, max=theta.tau.max)
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
  
  model.params<-c(alpha.y, alpha.x, theta.x, sigma.x, alpha.tau, theta.tau, sigma.tau)
  reg.params<- c(b0, b1, b2)
  
  return(list(model.params=model.params, reg.params=reg.params))
  
}

d.ououcirprior <- function(param=param,regbound=regbound,sigup=sigup){
  alpha.y<-param[1]
  alpha.x<-param[2]
  theta.x<-param[3]
  sigma.x<-param[4]
  alpha.tau<-param[5]
  theta.tau<-param[6]
  sigma.tau<-param[7]
  b0<-param[8]
  b1<-param[9]
  b2<-param[10]
  prior.alpha.y <- dexp(alpha.y, rate=1/10, log=T)
  prior.alpha.x <- dexp(alpha.x, rate=1/10, log=T)
  prior.theta.x <- dnorm(theta.x, mean=mean(regbound[1],regbound[2]), sd=(regbound[2]-regbound[1])/3, log=T)
  prior.sigma.x <- dexp(sigma.x, rate=1/sigup, log=T)
  prior.alpha.tau <- dexp(alpha.tau, rate=10, log=T)
  prior.theta.tau <- dchisq(theta.tau, df=2, ncp=0, log=T)
  prior.sigma.tau <- dexp(sigma.tau, rate=sigup, log=T)
  
  prior.b0 <- dnorm(b0, mean=mean(regbound[1],regbound[2]), sd=(regbound[2]-regbound[1])/3, log=T)
  prior.b1 <- dunif(b1, min=regbound[3], max=regbound[4], log=T)
  prior.b2 <- dunif(b2, min=regbound[5], max=regbound[6], log=T)
  return(sum(prior.alpha.y,prior.alpha.x,prior.theta.x,prior.sigma.x,prior.alpha.tau,prior.theta.tau,prior.sigma.tau,prior.b0,prior.b1,prior.b2))
}

ououcirproposal <- function(param=param,regbound=regbound,sigup=sigup){
  alpha.y<-param[1]
  alpha.x<-param[2]
  theta.x<-param[3]
  sigma.x<-param[4]
  alpha.tau<-param[5]
  theta.tau<-param[6]
  sigma.tau<-param[7]
  b0<-param[8]
  b1<-param[9]
  b2<-param[10]
  
  prop.alpha.y<-rnorm(1,mean=alpha.y,sd=sqrt((10-0)^2/12)/4)
  while (prop.alpha.y < 0){
    prop.alpha.y<-rnorm(1,mean=alpha.y,sd=sqrt((10-0)^2/12)/4)
  }
  prop.alpha.x<-rnorm(1,mean=alpha.x,sd=sqrt((10-0)^2/12)/4)
  while (prop.alpha.x < 0){
    prop.alpha.x<-rnorm(1,mean=alpha.x,sd=sqrt((10-0)^2/12)/4)
  }
  
  prop.theta.x <- rnorm(1,mean=theta.x,sd=sqrt((regbound[2]-regbound[1])^2/12)/4)
  prop.sigma.x <- rnorm(1,mean=sigma.x,sd=sqrt(sigup^2/12)/4)
  while (prop.sigma.x < 0) {
    prop.sigma.x <- rnorm(1,mean=sigma.x,sd=sqrt(sigup^2/12)/4)
  }
  prop.alpha.tau<-rnorm(1,mean=alpha.tau,sd=sqrt((10-0)^2/12)/4)
  while (prop.alpha.tau < 0){
    prop.alpha.tau<-rnorm(1,mean=alpha.tau,sd=sqrt((10-0)^2/12)/4)
  }
  prop.theta.tau <- rnorm(1,mean=theta.tau,sd=sqrt(sigup^2/12)/4)
  while (prop.theta.tau < 0) {
    prop.theta.tau <- rnorm(1,mean=theta.tau,sd=sqrt(sigup^2/12)/4)
  }
  prop.sigma.tau <- rnorm(1,mean=sigma.tau,sd=sqrt(sigup^2/12)/4)
  while (prop.sigma.tau < 0) {
    prop.sigma.tau <- rnorm(1,mean=sigma.tau,sd=sqrt(sigup^2/12)/4)
  }
  prop.b0 <- rnorm(1,mean=b0,sd=sqrt((regbound[2]-regbound[1])^2/12)/4)
  prop.b1 <- rnorm(1,mean=b1,sd=sqrt((regbound[4]-regbound[3])^2/12)/4)
  prop.b2 <- rnorm(1,mean=b2,sd=sqrt((regbound[6]-regbound[5])^2/12)/4)
  
  return(c(prop.alpha.y,prop.alpha.x,prop.theta.x,prop.sigma.x,prop.alpha.tau,prop.theta.tau,prop.sigma.tau,prop.b0,prop.b1,prop.b2))
}

ououcir_abc_MCMC <- function(startvalue=startvalue, iterations=iterations, y=y,x1=x1,x2=x2, regbound=regbound,sigup=sigup, tree=tree, rootvalue=rootvalue, errorbound=errorbound){
  
  #y<- resptrait
  #x1<-predtrait1
  #x2<-predtrait2
  
  S0y <- sum.stat(trait = y, tree = tree)
  S0x1 <- sum.stat(trait = x1, tree = tree)
  S0x2 <- sum.stat(trait = x2, tree = tree)
  S0<-c(S0y,S0x1,S0x2)
  #print(S0)
  distS0S1 <- c()
  chain = array(dim=c(iterations+1,10))
  chain[1,] = startvalue
  for (i in 1:iterations){
    if(i %%5000==0){print(paste("ououcir",i,sep=""))}
    #print(i)
    #    Sys.sleep(1)
    proposal <- ououcirproposal(chain[i,],regbound=regbound, sigup=sigup)
    # print(proposal)
    simD <- ououcirmodel(model.params=proposal[1:7],reg.params=proposal[8:10],root=rootvalue,tree=tree)
    
    # print(y)
    # print(simD$y)
    # print("-------------------------")
    # print(x1)
    # print(simD$x1)
    # print("-------------------------")
    # print(x2)
    # print(simD$x2)
    # print("-------------------------")
    
    S1y <- sum.stat(trait = simD$y, tree=tree)
    S1x1 <- sum.stat(trait = simD$x1, tree=tree)
    S1x2 <- sum.stat(trait = simD$x2, tree=tree)
    S1<-c(S1y,S1x1,S1x2)  
    # print(S0)
    # print(S1)
    # print(dist(rbind(S0,S1)))
    # cat("\n\n")
    # 
    distS0S1 <- c(distS0S1,dist(rbind(S0,S1)))
    if(dist(rbind(S0,S1))<errorbound){
      priorratio <- exp(d.ououcirprior(proposal, regbound=regbound,sigup=sigup)-d.ououcirprior(chain[i,],regbound=regbound,sigup=sigup))
      h=min(1,priorratio)
      if(runif(1)<h){
        chain[i+1,]=proposal
        #S0 <- S1
        #print("updated proposal successful")
        #print(c(i,proposal))
      }else{
        chain[i+1,]=chain[i,]
      }
    }else{
      chain[i+1,] = chain[i,]
    }
  }
  return(list(chain=chain,distS0S1=distS0S1))
}

sum.stat<-function(trait=trait,tree=tree){
  pic.trait<-pic(x=trait,phy=tree)
  #print("pic mean and pic sd")
  #print(c(mean(pic.trait),sd(pic.trait)))
  return(c(mean(pic.trait),sd(pic.trait)))
}

sum.stat.distance<-function(raw.sum.stat=raw.sum.stat,sim.sum.stat=sim.sum.stat){
  raw.sum.stat <- raw.sum.stat/median(raw.sum.stat-median(raw.sum.stat))
  sim.sum.stat <- sim.sum.stat/median(sim.sum.stat-median(sim.sum.stat))
  return(sqrt(sum((raw.sum.stat-sim.sum.stat)^2)))
}

regboundfcn<-function(olssum=olssum){
  bdd<-  c(olssum$coefficients[,1]-4*olssum$coefficients[,2],olssum$coefficients[,1]+4*olssum$coefficients[,2])
  bdd<-bdd[c(1,4,2,5,3,6)]  
  return(bdd)
}

# sigsqx.upper<-function(resptrait=resptrait,predtrait1=predtrait1,predtrait2=predtrait2){
#   return(max(diff(range(resptrait)),diff(range(predtrait1)),diff(range(predtrait2)))/4)
# }

sigsqx.upper<-function(predtrait1=predtrait1,predtrait2=predtrait2){
  return(max(diff(range(predtrait1)),diff(range(predtrait2)))/4)
}