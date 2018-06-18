#rm(list=ls())
library(TreeSim)
library(EasyABC)
library(coda)
library(Sim.DiffProc)
library(MCMCpack)
library(ape)
library(abc)

#source("/Users/Administrator/Dropbox/TerryLai/R_code/abcmcmc/ouest_Tony.R")
source("~/Dropbox/TerryLai/R_code/abcmcmc/ouest_Tony.R")
abcmcmc <- function(tree=tree, resptrait=resptrait, predtrait1=predtrait1,predtrait2=predtrait2, iterations=iterations, errorbound=errorbound, model=model){
  if(model == "oubmbm"){
    
    source("~/Dropbox/TerryLai/R_code/abcmcmc/oubmbmabcmcmc.r")
    print(model)
    estOU <-getOUest(Y=resptrait,tree=tree)
    olssum <- summary(lm(resptrait~predtrait1+predtrait2))
    sigup<-sigsqx.upper(predtrait1=predtrait1,predtrait2=predtrait2)
    regbound<-regboundfcn(olssum=olssum)
    alpha.y.min = 0.67*estOU$alpha
    alpha.y.max = 1.33*estOU$alpha
    sigma.x.min <- 0.5*min(var(predtrait1),var(predtrait2))
    sigma.x.max <- 1.5*max(var(predtrait1),var(predtrait2))
    tau.min <- 0
    tau.max <- sigup
    b0.min=regbound[1]
    b0.max=regbound[2]
    b1.min=regbound[3]
    b1.max=regbound[4]
    b2.min=regbound[5]
    b2.max=regbound[6]
    
    prior.model.params=c(alpha.y.min,alpha.y.max, sigma.x.min, sigma.x.max, tau.min, tau.max)
    
    names(prior.model.params)<-c("alpha.y.min","alpha.y.max", "sigma.x.min", "sigma.x.max", "tau.min", "tau.max") 
    prior.reg.params=c(b0.min, b0.max, b1.min, b1.max, b2.min, b2.max)
    prior.params <- oubmbmprior(prior.model.params=prior.model.params,prior.reg.params=prior.reg.params)
    
    startvalue <- c(unlist(prior.params))
    
    rootresp <- ace(resptrait, tree, type="continuous", method="REML")
    ls(rootresp)
    rootvalue<-rootresp$ace[1]
    
    oubmbmresult <- oubmbm_abc_MCMC(startvalue=startvalue, iterations=iterations, y=resptrait,x1=predtrait1,x2=predtrait2, regbound=regbound,sigup=sigup, tree=tree, rootvalue=rootvalue, errorbound=errorbound)
    return(oubmbmresult=oubmbmresult)
  }
  
  if(model == "ououbm"){
    source("~/Dropbox/TerryLai/R_code/abcmcmc/ououbmabcmcmc.r")
    print(model)
    estOU <-getOUest(Y=resptrait,tree=tree)
    estOUpred1 <- getOUest(Y=predtrait1, tree=tree)
    estOUpred2 <- getOUest(Y=predtrait2, tree=tree) 
    
    olssum <- summary(lm(resptrait~predtrait1+predtrait2))

    sigup<-sigsqx.upper(predtrait1=predtrait1,predtrait2=predtrait2)
    regbound<-regboundfcn(olssum=olssum)
    
    names(regbound)
    
    alpha.y.min <- 0.67*estOU$alpha #0.67*oufitresp$opt$alpha
    alpha.y.max <- 1.33*estOU$alpha #1.33*oufitresp$opt$alpha
    alpha.x.min <- 0.67*min(estOUpred1$alpha,estOUpred2$alpha)
    alpha.x.max <- 1.33*max(estOUpred1$alpha,estOUpred2$alpha)
    theta.x.min <- 0.67*min(estOUpred1$mu,estOUpred2$mu)
    theta.x.max <- 1.33*max(estOUpred1$mu,estOUpred2$mu) 
    sigma.x.min <- 0.5*min(var(predtrait1),var(predtrait2))
    sigma.x.max <- 1.5*max(var(predtrait1),var(predtrait2))
    tau.min <- 0
    tau.max <- sigup
    b0.min=regbound[1]
    b0.max=regbound[2]
    b1.min=regbound[3]
    b1.max=regbound[4]
    b2.min=regbound[5]
    b2.max=regbound[6]
    # 
    prior.model.params=c(alpha.y.min,alpha.y.max, alpha.x.min, alpha.x.max, theta.x.min, theta.x.max, sigma.x.min, sigma.x.max, tau.min, tau.max)
    
    names(prior.model.params)<-c("alpha.y.min","alpha.y.max", "alpha.x.min", "alpha.x.max", "theta.x.min", "theta.x.max", "sigma.x.min", "sigma.x.max", "tau.min", "tau.max") 
    prior.reg.params=c(b0.min, b0.max, b1.min, b1.max, b2.min, b2.max)
    prior.params <- ououbmprior(prior.model.params=prior.model.params,prior.reg.params=prior.reg.params)
    
    startvalue <- c(unlist(prior.params))
    rootresp <- ace(resptrait, tree, type="continuous", method="REML")
    ls(rootresp)
    rootvalue<-rootresp$ace[1]
    ououbmresult <- ououbm_abc_MCMC(startvalue=startvalue, iterations=iterations, y=resptrait,x1=predtrait1,x2=predtrait2, regbound=regbound,sigup=sigup, tree=tree, rootvalue=rootvalue, errorbound=errorbound)
    return(ououbmresult=ououbmresult)
  }
  if(model == "oubmcir"){

    source("~/Dropbox/TerryLai/R_code/abcmcmc/oubmcirabcmcmc.r")
    print(model)
    estOU <-getOUest(Y=resptrait,tree=tree)
    olssum <- summary(lm(resptrait~predtrait1+predtrait2))
    #sigup<-sigsqx.upper(resptrait=resptrait,predtrait1=predtrait1,predtrait2=predtrait2)
    sigup<-sigsqx.upper(predtrait1=predtrait1,predtrait2=predtrait2)
    regbound<-regboundfcn(olssum=olssum)
    
    alpha.y.min = 0.67*estOU$alpha #0.67*oufitresp$opt$alpha
    alpha.y.max = 1.33*estOU$alpha #1.33*oufitresp$opt$alpha
    sigma.x.min <- 0.5*min(var(predtrait1),var(predtrait2))
    sigma.x.max <- 1.5*max(var(predtrait1),var(predtrait2))
    alpha.tau.min <- 0.67*estOU$alpha
    alpha.tau.max <- 1.33*estOU$alpha
    theta.tau.min <- 0.67*estOU$mu
    theta.tau.max <- 1.33*estOU$mu #sigup#diff(range(resptrait))/4#sd(resptrait)
    sigma.tau.min <- 0.5*min(var(predtrait1),var(predtrait2))#0.67*17#0.67*estOU$sigsq
    sigma.tau.max <- 1.5*max(var(predtrait1),var(predtrait2))#1.33*17#1.33*estOU$sigsq
    
    b0.min=regbound[1]
    b0.max=regbound[2]
    b1.min=regbound[3]
    b1.max=regbound[4]
    b2.min=regbound[5]
    b2.max=regbound[6]
    # 
    prior.model.params=c(alpha.y.min,alpha.y.max, sigma.x.min, sigma.x.max, alpha.tau.min, alpha.tau.max,theta.tau.min,theta.tau.max,sigma.tau.min,sigma.tau.max)
    names(prior.model.params)<-c("alpha.y.min","alpha.y.max", "sigma.x.min", "sigma.x.max", "alpha.tau.min", "alpha.tau.max","theta.tau.min","theta.tau.max","sigma.tau.min","sigma.tau.max") 
    prior.reg.params=c(b0.min, b0.max, b1.min, b1.max, b2.min, b2.max)
    prior.params <- oubmcirprior(prior.model.params=prior.model.params,prior.reg.params=prior.reg.params)
    
    startvalue <- c(unlist(prior.params))
    rootresp <- ace(resptrait, tree, type="continuous", method="REML")
    rootvalue<-rootresp$ace[1]
    oubmcirresult <- oubmcir_abc_MCMC(startvalue=startvalue, iterations=iterations, y=resptrait,x1=predtrait1,x2=predtrait2, regbound=regbound,sigup=sigup, tree=tree, rootvalue=rootvalue, errorbound=errorbound)
    return(oubmcirresult=oubmcirresult)
  }
  if(model == "ououcir"){

    source("~/Dropbox/TerryLai/R_code/abcmcmc/ououcirabcmcmc.r")
    print(model)
    estOU <-getOUest(Y=resptrait,tree=tree)
    estOUpred1 <- getOUest(Y=predtrait1, tree=tree)
    estOUpred2 <- getOUest(Y=predtrait2, tree=tree)
    olssum <- summary(lm(resptrait~predtrait1+predtrait2))
    #sigup<-sigsqx.upper(resptrait=resptrait,predtrait1=predtrait1,predtrait2=predtrait2)
    sigup<-sigsqx.upper(predtrait1=predtrait1,predtrait2=predtrait2)
    regbound<-regboundfcn(olssum=olssum)
    
    alpha.y.min <- 0.67*estOU$alpha #0.67*oufitresp$opt$alpha
    alpha.y.max <- 1.33*estOU$alpha #1.33*oufitresp$opt$alpha
    alpha.x.min <- 0.67*min(estOUpred1$alpha,estOUpred2$alpha)
    alpha.x.max <- 1.33*max(estOUpred1$alpha,estOUpred2$alpha)
    theta.x.min <- 0.67*min(estOUpred1$mu,estOUpred2$mu)
    theta.x.max <- 1.33*max(estOUpred1$mu,estOUpred2$mu) 
    sigma.x.min <- 0.5*min(var(predtrait1),var(predtrait2))
    sigma.x.max <- 1.5*max(var(predtrait1),var(predtrait2))
    alpha.tau.min <- 0.67*estOU$alpha
    alpha.tau.max <- 1.33*estOU$alpha
    theta.tau.min <- 0.67*estOU$mu
    theta.tau.max <- 1.33*estOU$mu #sigup#diff(range(resptrait))/4#sd(resptrait)
    sigma.tau.min <- 0.5*min(var(predtrait1),var(predtrait2))
    sigma.tau.max <- 1.5*max(var(predtrait1),var(predtrait2))
    
    b0.min=regbound[1]
    b0.max=regbound[2]
    b1.min=regbound[3]
    b1.max=regbound[4]
    b2.min=regbound[5]
    b2.max=regbound[6]
    # 
    prior.model.params=c(alpha.y.min,alpha.y.max,alpha.x.min,alpha.x.max,theta.x.min,theta.x.max, sigma.x.min, sigma.x.max, alpha.tau.min, alpha.tau.max,theta.tau.min,theta.tau.max,sigma.tau.min,sigma.tau.max)
    names(prior.model.params)<-c("alpha.y.min","alpha.y.max","alpha.x.min","alpha.x.max","theta.x.min","theta.x.max", "sigma.x.min", "sigma.x.max", "alpha.tau.min", "alpha.tau.max","theta.tau.min","theta.tau.max","sigma.tau.min","sigma.tau.max") 
    prior.reg.params=c(b0.min, b0.max, b1.min, b1.max, b2.min, b2.max)
    prior.params <- ououcirprior(prior.model.params=prior.model.params,prior.reg.params=prior.reg.params)
    startvalue <- c(unlist(prior.params))
    rootresp <- ace(resptrait, tree, type="continuous", method="REML")
    rootvalue<-rootresp$ace[1]
    ououcirresult <- ououcir_abc_MCMC(startvalue=startvalue, iterations=iterations, y=resptrait,x1=predtrait1,x2=predtrait2, regbound=regbound,sigup=sigup, tree=tree, rootvalue=rootvalue, errorbound=errorbound)
    return(ououcirresult=ououcirresult)
  }
}

###### main 
# tree<-read.tree("~/Dropbox/TerryLai/treetraitdata/bonine_lizard/tree.phy")
# trait<-read.csv("~/Dropbox/TerryLai/treetraitdata/bonine_lizard/trait.csv")
# resptrait<-trait$bodymass
# predtrait1<-trait$snoutlength
# predtrait2<-trait$thighmuscle

oubmbmresult <- abcmcmc(tree = tree, resptrait = resptrait, predtrait1 = predtrait1, predtrait2 = predtrait2, iterations = iterations, errorbound = errorbound, model = "oubmbm")
ououbmresult <- abcmcmc(tree = tree, resptrait = resptrait, predtrait1 = predtrait1, predtrait2 = predtrait2, iterations = iterations, errorbound = errorbound, model = "ououbm")
oubmcirresult <- abcmcmc(tree = tree, resptrait = resptrait, predtrait1 = predtrait1, predtrait2 = predtrait2, iterations = iterations, errorbound = errorbound, model = "oubmcir")
ououcirresult <- abcmcmc(tree = tree, resptrait = resptrait, predtrait1 = predtrait1, predtrait2 = predtrait2, iterations = iterations, errorbound = errorbound, model = "ououcir")
