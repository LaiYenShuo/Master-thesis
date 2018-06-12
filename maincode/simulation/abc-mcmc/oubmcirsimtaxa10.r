#### Simulation ABC-MCMC
#### tree, root, parameters, taxa.size, startvalue, S0 <- truetrait, distS0S1, errorbound, 

rm(list=ls())
source("~/Dropbox/TerryLai/R_code/abcmcmc/oubmcirabcmcmc.r")
#source("/Users/ms022/Dropbox/TerryLai/R_code/abcmcmc/oubmcirabcmcmc.r")
#source("C:/Users/ms022/Dropbox/TerryLai/R_code/abcmcmc/oubmcirabcmcmc.r")
#source("~/Dropbox/TerryLai/R_code/abcmcmc/ouest_Tony.R")
#source("/Users/ms022/Dropbox/TerryLai/R_code/abcmcmc/ouest_Tony.R")

root<-0
#true params
true.alpha.y<-0.15
true.sigma.x<-1
true.alpha.tau<-0.20
true.theta.tau<-30
true.sigma.tau<- 0.5


true.b0 <- 0
true.b1 <- 0.5
true.b2 <- 0.5

true.param.array<-c(true.alpha.y,true.sigma.x,true.alpha.tau,true.theta.tau,true.sigma.tau,true.b0,true.b1,true.b2)
#hyper parameters
alpha.y.min <- 0
alpha.y.max <- 0.3
sigma.x.min <- 0
sigma.x.max <- 2
alpha.tau.min <- 0
alpha.tau.max <- 0.4
theta.tau.min <- 0 #assume uniform
theta.tau.max <- 60
sigma.tau.min <- 0# assume inv gamma
sigma.tau.max <- 1
b0.min=-1
b0.max=1
b1.min=0
b1.max=1
b2.min=0
b2.max=1

regbound <- c(b0.min, b0.max, b1.min, b1.max, b2.min, b2.max)

prior.model.params=c(alpha.y.min,alpha.y.max, sigma.x.min, sigma.x.max, alpha.tau.min,alpha.tau.max, theta.tau.min, theta.tau.max, sigma.tau.min, sigma.tau.max)
names(prior.model.params)<-c("alpha.y.min","alpha.y.max", "sigma.x.min", "sigma.x.max", "alpha.tau.min","alpha.tau.max", "theta.tau.min", "theta.tau.max", "sigma.tau.min", "sigma.tau.max")
prior.reg.params=c(b0.min, b0.max, b1.min, b1.max, b2.min, b2.max)

numbsim<-1;lambda<-2.0;mu<-0.5;frac<-0.6;age<-2
taxa.size.array<-c(10)#,20,50,100)

for(taxa.size.Index in 1:length(taxa.size.array)){
  #taxa.size.Index = 1
  n<-taxa.size.array[taxa.size.Index]
  print(paste("taxa",n,sep=""))
  #sim.oubmcir.trait<-array(0,c(n,3,sims))
  tree<-sim.bd.taxa.age(n=n,numbsim=1,lambda=lambda,mu=mu,frac=frac,age=age,mrca=TRUE)[[1]]
  tree<-reorder(tree,"postorder")
  # tree$edge
  # plot(tree)
  # nodelabels()
  # tiplabels()
  true.trait<-oubmcirmodel(model.params=c(true.alpha.y,true.sigma.x,true.alpha.tau,true.theta.tau,true.sigma.tau),reg.params=c(true.b0,true.b1,true.b2),root=root,tree=tree)
  assign(paste("true.trait.taxa",taxa.size.array[taxa.size.Index],sep=""),true.trait)
  # for(simIndex in 1:sims){
  #   if(simIndex %%10==0){print(simIndex)}
  #   prior.params <- oubmcirprior(prior.model.params=prior.model.params,prior.reg.params=prior.reg.params)
  #   model.params.array[,simIndex]<-prior.params$model.params#for record only
  #   reg.params.array[,simIndex]<-prior.params$reg.params#for record only
  #   sim.trait <-oubmcirmodel(model.params=prior.params$model.params,reg.params=prior.params$reg.params,root=root,tree=tree)
  #   sim.oubmcir.trait[,1,simIndex]<-sim.trait$y
  #   sim.oubmcir.trait[,2,simIndex]<-sim.trait$x1
  #   sim.oubmcir.trait[,3,simIndex]<-sim.trait$x2
  #   y.sum.stat.array[,simIndex]<- sum.stat(trait=sim.trait$y,tree=tree)
  #   x1.sum.stat.array[,simIndex]<- sum.stat(trait=sim.trait$x1,tree=tree)
  #   x2.sum.stat.array[,simIndex]<- sum.stat(trait=sim.trait$x2,tree=tree)
  ####  run abc-mcmc    
  startvalue <- true.param.array #c(model.params.array[,simIndex],reg.params.array[,simIndex])
  errorbound = 100
  iterations = 50000
  sigup <- sigsqx.upper(predtrait1 = true.trait$x1, predtrait2 = true.trait$x2)
  
  oubmcirresult10 <- oubmcir_abc_MCMC(startvalue = startvalue, iterations = iterations, y = true.trait$y, x1=true.trait$x1, x2=true.trait$x2, regbound = regbound, sigup = sigup, tree=tree, rootvalue = root, errorbound = errorbound)
  #m.param[simIndex,] <-apply(result$chain, 2, mean)
  #apply(result$chain,2,mean)
  setwd("/Users/terrylai/Documents/simulation_mcmc/")
  save.image(paste("oubmcirsimstaxa",taxa.size.array[taxa.size.Index], ".RData", sep=""))
  
}# end of taxasize
