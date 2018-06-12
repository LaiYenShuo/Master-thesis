#### Simulation ABC-MCMC
#### tree, root, parameters, taxa.size, startvalue, S0 <- truetrait, distS0S1, errorbound, 

rm(list=ls())
source("~/Dropbox/TerryLai/R_code/abcmcmc/oubmbmabcmcmc.r")

root<-0
#true params
true.alpha.y<-0.15
true.sigma.x<-1
true.tau<- 0.35

true.b0 <- 0
true.b1 <- 0.5
true.b2 <- 0.5

true.param.array<-c(true.alpha.y,true.sigma.x,true.tau,true.b0,true.b1,true.b2)

#hyper parameters
alpha.y.min <- 0
alpha.y.max <- 0.3
sigma.x.min <- 0
sigma.x.max <- 2
tau.min <- 0
tau.max <- 0.7
b0.min=-1
b0.max=1
b1.min=0
b1.max=1
b2.min=0
b2.max=1
regbound <- c(b0.min, b0.max, b1.min, b1.max, b2.min, b2.max)

prior.model.params=c(alpha.y.min,alpha.y.max, sigma.x.min, sigma.x.max, tau.min, tau.max)
names(prior.model.params)<-c("alpha.y.min","alpha.y.max", "sigma.x.min", "sigma.x.max", "tau.min", "tau.max") 
prior.reg.params=c(b0.min, b0.max, b1.min, b1.max, b2.min, b2.max)

numbsim<-1;lambda<-2.0;mu<-0.5;frac<-0.6;age<-2
taxa.size.array<-c(10)#,20,50,100)

for(taxa.size.Index in 1:length(taxa.size.array)){
  n<-taxa.size.array[taxa.size.Index]
  print(paste("taxa",n,sep=""))
  tree<-sim.bd.taxa.age(n=n,numbsim=1,lambda=lambda,mu=mu,frac=frac,age=age,mrca=TRUE)[[1]]
  tree<-reorder(tree,"postorder")
  true.trait<-oubmbmmodel(model.params=c(true.alpha.y,true.sigma.x,true.tau),reg.params=c(true.b0,true.b1,true.b2),root=root,tree=tree)
  assign(paste("true.trait.taxa",taxa.size.array[taxa.size.Index],sep=""),true.trait)
####  run abc-mcmc    
startvalue <- true.param.array #c(model.params.array[,simIndex],reg.params.array[,simIndex])
errorbound = 100
iterations = 50000
sigup <- sigsqx.upper(predtrait1 = true.trait$x1, predtrait2 = true.trait$x2)

oubmbmresult10 <- oubmbm_abc_MCMC(startvalue = startvalue, iterations = iterations, y = true.trait$y, x1=true.trait$x1, x2=true.trait$x2, regbound = regbound, sigup = sigup, tree=tree, rootvalue = root, errorbound = errorbound)

setwd("/Users/terrylai/Documents/simulation_mcmc/")
save.image(paste("oubmbmsimstaxa",taxa.size.array[taxa.size.Index], ".RData", sep=""))

}# end of taxasize
#colnames(m.param) <- c("alpha.y", "sigma.x", "tau", "b0", "b1", "b2")
