# rm(list=ls())
# #install.packages("TreeSim")
library(TreeSim)
library(EasyABC)
library(coda)
library(Sim.DiffProc)
library(MCMCpack)
library(ape)
library(abc)
#install.packages("MCMCpack")

#oubmbmintegrand<-function(s, alpha.y=alpha.y){
#    alpha.y*exp(alpha.y*s)*rnorm(n=1,mean=0,sd=s)
#    }
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
    #inttheta<-integrate(oubmbmintegrand,lower=0 ,upper=treelength[index], alpha.y=true.alpha.y)
    
    #INT1var<-treelength[index]*exp(2*alpha.y*treelength[index]) - 2*(exp(2*alpha.y*treelength[index])-exp(alpha.y*treelength[index])) + alpha.y/2*(exp(2*alpha.y*treelength[index])-1)
    #print(alpha.y)
    INT1mean<- optimnodestates[des[index]]*exp(alpha.y*treelength[index])  - optimnodestates[anc[index]]
    INT1var<-  (exp(2*alpha.y*treelength[index])-1)/ (2*alpha.y)
    
    INT1<-exp(-alpha.y*treelength[index])* rnorm(n=1,mean= INT1mean, sd=sqrt(abs(INT1var)))
    #INT1<-exp(-alpha.y*treelength[index])* rnorm(n=1,mean=0,sd=sqrt(sigma.sq.theta*abs(INT1var)))
    
    fexpr<-expression(alpha.y*exp(alpha.y*t)*w)
    #print(alpha.y)
    res<-st.int(fexpr,type="ito",M=1,lower=0,upper=treelength[index])
    INT2<-median(res$X)
    ynodestates[des[index]]<-ynodestates[anc[index]] + INT1 + INT2
  }
  
  simtrait<-ynodestates[1:n]
  
  return(list(y=simtrait,x1=x1nodestates[1:n],x2=x2nodestates[1:n]))
  #return(c(mean(simtrait),sd(simtrait)))
}

oubmbmprior<-function(prior.model.params=prior.model.params,prior.reg.params=prior.reg.params){
  alpha.y.rate <-prior.model.params["alpha.y.rate"]
  sigma.x.rate <-prior.model.params["sigma.x.rate"]
  tau.rate <- prior.model.params["tau.rate"]
  # sigma.x.shape <-prior.model.params["sigma.x.shape"]
  # sigma.x.scale <-prior.model.params["sigma.x.scale"]
  # tau.shape <- prior.model.params["tau.shape"]
  # tau.scale <- prior.model.params["tau.scale"]
  
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

#sum.stat.distance(raw.sum.stat =raw.sum.stat.y, sim.sum.stat = y.sum.stat.array[,simIndex] )
#raw.sum.stat.y/median(raw.sum.stat.y-median(raw.sum.stat.y))

sum.stat.distance<-function(raw.sum.stat=raw.sum.stat,sim.sum.stat=sim.sum.stat){
  raw.sum.stat <- raw.sum.stat/median(raw.sum.stat-median(raw.sum.stat))
  sim.sum.stat <- sim.sum.stat/median(sim.sum.stat-median(sim.sum.stat))
  return(sum((raw.sum.stat-sim.sum.stat)^2))
}

# #main
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
# true.sigma.x<-2
# true.tau<- 1
# 
# true.b0 <- 0
# true.b1 <- 1
# true.b2 <- 0.2
# 
# #hyper parameters
# alpha.y.rate <- 5
# sigma.x.shape <-2
# sigma.x.scale <-1
# tau.shape <- 3
# tau.scale <- 1
# b0.min=-5
# b0.max=5
# b1.min=-5
# b1.max=5
# b2.min=-5
# b2.max=5
# 
# prior.model.params=c(alpha.y.rate, sigma.x.shape, sigma.x.scale, tau.shape, tau.scale)
# names(prior.model.params)<-c("alpha.y.rate", "sigma.x.shape", "sigma.x.scale", "tau.shape", "tau.scale")
# prior.reg.params=c(b0.min, b0.max, b1.min, b1.max, b2.min, b2.max)
# # alpha.y <- rexp(n=1, rate=5);alpha.y
# # sigma.x <- rinvgamma(n=1, shape=2, scale=1);sigma.x
# # tau <- rinvgamma(n=1, shape=3, scale=1);tau
# #
# # b0 <- runif(n=1, min=-5, max=5);b0
# # b1 <- runif(n=1, min=-5, max=5);b1
# # b2 <- runif(n=1, min=-5, max=5);b2
# prior.params <- oubmbmprior(prior.model.params=prior.model.params,prior.reg.params=prior.reg.params)
# 
# true.trait<-oubmbmmodel(model.params=c(true.alpha.y,true.sigma.x,true.tau),reg.params=c(true.b0,true.b1,true.b2),root=root,tree=tree)
# sim.trait <-oubmbmmodel(model.params=prior.params$model.params,reg.params=prior.params$reg.params,root=root,tree=tree)
# true.trait
# 
# raw.sum.stat.y<-sum.stat(trait=true.trait$y,tree=tree)
# raw.sum.stat.x1<-sum.stat(trait=true.trait$x1,tree=tree)
# raw.sum.stat.x2<-sum.stat(trait=true.trait$x2,tree=tree)
# 
# sims<-50000
# sim.oubmbm.trait<-array(0,c(n,3,sims))
# model.params.array<-array(0,c(3,sims))
# rownames(model.params.array)<-c("alpha.y","sigma.x","tau")
# reg.params.array<-array(0,c(3,sims))
# row.names(reg.params.array)<-c("b0", "b1", "b2")
# y.sum.stat.array<-array(0,c(2,sims))
# rownames(y.sum.stat.array)<-c("y.mean","y.sd")
# x1.sum.stat.array<-array(0,c(2,sims))
# rownames(x1.sum.stat.array)<-c("x1.mean","x1.sd")
# x2.sum.stat.array<-array(0,c(2,sims))
# rownames(x2.sum.stat.array)<-c("x2.mean","x2.sd")
# 
# prior.params <- oubmbmprior(prior.model.params=prior.model.params,prior.reg.params=prior.reg.params)
# sum.stat.distance.array<-array(0,c(sims))
# for(simIndex in 1:sims){
#   if(simIndex %%1000==0){print(simIndex)}
#   prior.params <- oubmbmprior(prior.model.params=prior.model.params,prior.reg.params=prior.reg.params)
#   model.params.array[,simIndex]<-prior.params$model.params#for record only
#   reg.params.array[,simIndex]<-prior.params$reg.params#for record only
# 
#   sim.trait <-oubmbmmodel(model.params=prior.params$model.params,reg.params=prior.params$reg.params,root=root,tree=tree)
#   sim.oubmbm.trait[,1,simIndex]<-sim.trait$y
#   sim.oubmbm.trait[,2,simIndex]<-sim.trait$x1
#   sim.oubmbm.trait[,3,simIndex]<-sim.trait$x2
#   y.sum.stat.array[,simIndex]<- sum.stat(trait=sim.trait$y,tree=tree)
#   x1.sum.stat.array[,simIndex]<- sum.stat(trait=sim.trait$x1,tree=tree)
#   x2.sum.stat.array[,simIndex]<- sum.stat(trait=sim.trait$x2,tree=tree)
# 
#   #y.dist<-sum.stat.distance(raw.sum.stat =raw.sum.stat.y, sim.sum.stat = y.sum.stat.array[,simIndex] )
#   #x1.dist<-sum.stat.distance(raw.sum.stat =raw.sum.stat.x1, sim.sum.stat = x1.sum.stat.array[,simIndex] )
#   #x2.dist<-sum.stat.distance(raw.sum.stat =raw.sum.stat.x2, sim.sum.stat = x2.sum.stat.array[,simIndex] )
#   #sum.stat.distance.array[simIndex] <- sqrt(y.dist + x1.dist + x2.dist)
#   # sum.stat.distance.array[simIndex] <- sum.stat.distance(raw.sum.stat =c(raw.sum.stat.y,raw.sum.stat.x1,raw.sum.stat.x2), sim.sum.stat = c(y.sum.stat.array[,simIndex],x1.sum.stat.array[,simIndex],x2.sum.stat.array[,simIndex]))
#   #
#   # if(sum.stat.distance.array[simIndex] <= threhold){
#   #   post.model.params.array[,simIndex]<-prior.params$model.params
#   #   post.reg.params.array[,simIndex]<-prior.params$reg.params
#   # }else{
#   #   post.model.params.array[,simIndex]<-post.model.params.array[,simIndex-1]
#   #   post.reg.params.array[,simIndex]<-post.reg.params.array[,simIndex-1]
#   #   }
#   }#end of loop
# 
# ### Use abc package
# sim.sum.stat <- cbind(t(y.sum.stat.array),t(x1.sum.stat.array),t(x2.sum.stat.array))
# oubmbm.par.sim <- cbind(t(model.params.array),t(reg.params.array))
# setwd("/Users/TerryLai/Dropbox/TerryLai/R_code/abc_V2")
# #setwd("~/Dropbox/FCU/Teaching/Mentoring/2017Spring/TerryLai/R_code/abc_V2/")
# save.image("test_oubmbm.RData")
# 
# ### The rejection alogoritm
# rej <- abc(target=c(raw.sum.stat.y,raw.sum.stat.x1,raw.sum.stat.x2), param=oubmbm.par.sim, sumstat=sim.sum.stat, tol=0.05, method="rejection"  )
# ls(rej)
# rej$unadj.values
# rej$ss
# rej$transf
# rej$logit.bounds
# rej$numparam
# rej$numstat
# rej$names
# ## ABC with local linear regression correction without/with correction
# ## for heteroscedasticity
# ##
# 
# lin <- abc(target=c(raw.sum.stat.y,raw.sum.stat.x1,raw.sum.stat.x2), param=oubmbm.par.sim, sumstat=sim.sum.stat, tol=0.2,hcorr=FALSE, method="loclinear", transf=c("none","log","none","log","none","log"))
# ls(lin)
# lin[c(1:18)]
# linhc <- abc(target=c(raw.sum.stat.y,raw.sum.stat.x1,raw.sum.stat.x2), param=oubmbm.par.sim, sumstat=sim.sum.stat, tol=0.2, method="loclinear", transf=c("none","log","none","log","none","log"))
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
#   expression(tau),
#   expression(b[0]),
#   expression(b[1]),
#   expression(b[2])),file = "~/Desktop/linhc")
# 
# 
# (linhc$adj.values[,4])
# 
# 
# #setwd("~/Dropbox/FCU/Teaching/Mentoring/2017Spring/TerryLai/R_code/abc_V2/")
# #save.image("test1.RData")
# ## diagnotic plots: compare the 2 'abc' objects:"localinear",
# ##"localinear" with correction for heteroscedasticity
# ##
# lin$numparam
# lin$adj.values
# plot(lin, param = oubmbm.par.sim)## warning
# plot(linhc, param = oubmbm.par.sim)## warning
# ## example illustrates how to add "true" parameter values to a plot
# ##
# oubmbm.par.sim
# postmod <- c(post.model.params.array[match(max(post.model.params.array[,2]), post.model.params.array[,2]),1], post.reg.params.array[match(max(post.reg.params.array[,2]),post.reg.params.array[,2]),1])
# plot(linhc, param = oubmbm.par.sim, true=postmod)
# ## artficial example to show how to use the logit transformations
# ##
# myp <- data.frame(par1=runif(n=1000, min=-1, max=1),par2=rnorm(n=1000), par3=runif(n=1000, min=0, max=2))
# mys <- myp + rnorm(n=1000, sd=0.1)
# myt <- c(0,0,1.5)
# lin2 <- abc(target=myt, param=myp, sumstat=mys, tol=.1, method="loclinear", transf=c("logit", "none", "logit"), logit.bounds= rbind(c(-1,1), c(NA,NA), c(0,2)))
# summary(lin2)
# 
# ##cv4abc 交叉驗證
# ## cv4abc() calls abc(). Here we show two ways for the supplying
# ## arguments of abc(). 1st way: passing arguments directly. In this
# ## example only 'param', 'sumstat', 'tol', and 'method', while default
# ## values are used for the other arguments.
# ## Number of eval. should be much more greater in realistic settings
# cv.rej <- cv4abc(param = oubmbm.par.sim, sumstat = sim.sum.stat,abc.out = rej, nval=5, tols=c(.1,.2,.3), methods = "rejection")
# ls(cv.rej)
# ## abc.out is NULL
# ## 2nd way: first creating an object of class 'abc', and then using it
# ## to pass its arguments to abc().
# ##
# lin <- abc(target=c(raw.sum.stat.y,raw.sum.stat.x1,raw.sum.stat.x2), param=oubmbm.par.sim, sumstat=sim.sum.stat, tol=0.2,hcorr=FALSE, method="loclinear", transf=c("none","log","none","log","none","log"))
# 
# cv.lin <- cv4abc(param = oubmbm.par.sim, sumstat = sim.sum.stat,abc.out=lin, nval=5, tols=c(.1,.2,.3))
# ls(cv.lin)
# cv.lin$true
# cv.lin$cvsamples
# cv.lin$estim
# cv.lin$seed
# cv.linhc <- cv4abc(param = oubmbm.par.sim, sumstat = sim.sum.stat,abc.out=linhc, nval=5, tols=c(.1,.2,.3))
# ## using the plot method. Different tolerance levels are plotted with
# ## different heat.colors. Smaller the tolerance levels correspond to
# ## more red points.
# ## !!! consider using the argument 'exclude' (plot.cv4abc) to supress
# ## the plotting of any outliers that mask readibility !!!
# plot(cv.lin, log=c("xy","xy","xy","xy","xy","xy"),caption=c(
#   expression(alpha[y]),
#   expression(sigma[x]),
#   expression(tau),
#   expression(b[0]),
#   expression(b[1]),
#   expression(b[2])))
# ## comparing with the rejection sampling
# plot(cv.rej, log=c("","xy","","xy","","xy"),caption=c(
#   expression(alpha[y]),
#   expression(sigma[x]),
#   expression(tau),
#   expression(b[0]),
#   expression(b[1]),
#   expression(b[2])))
# ## or printing results directly to a postscript file...
# plot(cv.lin, log=c("xy","xy","xy","xy","xy","xy"),caption=c(
#   expression(alpha[y]),
#   expression(sigma[x]),
#   expression(tau),
#   expression(b[0]),
#   expression(b[1]),
#   expression(b[2])),file="~/Desktop/CVrej", postscript = FALSE)
# ## using the summary method to calculate the prediction error
# summary(cv.lin)
# ## campare with rejection sampling
# summary(cv.rej)
# 
# ## cv4postpr Leave-one-our cross validation for model selection ABC
# require(abc.data)
# data(human)
# ###Reduce the sample size of the simulations to reduce the running time.
# ###Do not do that with your own data!
# ss<-c(1:1000,50001:51000,100001:101000)
# cv.modsel <- cv4postpr(models[ss], stat.3pops.sim[ss,], nval=5, tols=c(.05,.1), method="rejection")
# summary(cv.modsel)
# plot(cv.modsel, names.arg=c("Bottleneck", "Constant", "Exponential"))
# 
# 
# ## expected.deviance Expected deviance
# expected.deviance(target = c(raw.sum.stat.y,raw.sum.stat.x1,raw.sum.stat.x2), postsumstat = oubmbm.post.sim)
# ## warning
# 
# ## gift Goodness of fit
# oubmbmb.gfit <- gfit(target = c(raw.sum.stat.y,raw.sum.stat.x1,raw.sum.stat.x2), sumstat = sim.sum.stat, nb.replicate = 10, tol=0.01, statistic = mean)
# oubmbmb.gfit$dist.sim
# oubmbmb.gfit$dist.obs
# ## Plot the distribution of the null statistic and indicate where is the
# ## observed value.
# plot(oubmbmb.gfit, main="Histogram under H0")
# summary(oubmbmb.gfit)
# ## gfitpca Goodness of fit with principal component analysis
# oubmbm.gfitpca <- gfitpca(target=c(raw.sum.stat.y,raw.sum.stat.x1,raw.sum.stat.x2), sumstat=sim.sum.stat, cprob=0.1)
# ## no index we have one model can not model selection
# 
# ## hist.abc Posterior histograms
# hist(lin, unadj = FALSE, true = NULL, file = "~/Desktop/histabc_lin",
#      postscript = FALSE, onefile = TRUE, ask =!is.null(deviceIsInteractive()), col.hist = "grey", col.true = "red",caption = c(expression(alpha[y]),
#  expression(sigma[x]), expression(tau),expression(b[0]),expression(b[1]),expression(b[2])))
# 
# hist(linhc, unadj = FALSE, true = NULL, file = "~/Desktop/histabc_linhc",
#      postscript = FALSE, onefile = TRUE, ask =!is.null(deviceIsInteractive()), col.hist = "grey", col.true = "red",caption = c(expression(alpha[y]),
# expression(sigma[x]), expression(tau),expression(b[0]),expression(b[1]),expression(b[2])))
# ## plot.abc Diagnostic plots for ABC
# ## S3 method for class 'abc'
# plot(lin, param=oubmbm.par.sim, subsample = sims, true = NULL, file = "~/Desktop/plotabc_lin",postscript = FALSE, onefile = TRUE, ask =
#        !is.null(deviceIsInteractive()),caption = c(
#  expression(alpha[y]),expression(sigma[x]),expression(tau),expression(b[0]),expression(b[1]),expression(b[2])))
# plot(linhc, param=oubmbm.par.sim, subsample = sims, true = NULL, file = "~/Desktop/plotabc_linhc",postscript = FALSE, onefile = TRUE, ask =
#   !is.null(deviceIsInteractive()),caption = c(expression(alpha[y]), expression(sigma[x]),expression(tau),expression(b[0]),expression(b[1]),expression(b[2])))
# ## plot.cv4abc Cross-vaildation plots for ABC
# ## S3 method for class 'cv4abc'
# plot(cv.lin, exclude = NULL, log = NULL, file = "~/Desktop/cvlin",
#      postscript = FALSE, onefile = TRUE, ask =
#        !is.null(deviceIsInteractive()), caption = c(
#          expression(alpha[y]),
#          expression(sigma[x]),
#          expression(tau),
#          expression(b[0]),
#          expression(b[1]),
#          expression(b[2])))
# plot(cv.linhc, exclude = NULL, log = NULL, file = "~/Desktop/cvlinhc",
#      postscript = FALSE, onefile = TRUE, ask =
#        !is.null(deviceIsInteractive()), caption = c(
#          expression(alpha[y]),
#          expression(sigma[x]),
#          expression(tau),
#          expression(b[0]),
#          expression(b[1]),
#          expression(b[2])))
# ## plot.cv4postpr Barplot of model misclassification
# ## no cv4postpr
# 
# ##plot.gfit Goodness-of-fit plot for ABC
# plot(oubmbmb.gfit,  breaks="Freedman-Diaconis")
# 
# ## postpr Estimating posterior model probabilities
# ## Model selection with Approximate Bayesian Computation (ABC).
# ## one model do not use
# ## an artifical example
# ss <- cbind(runif(1000),rt(1000,df=20))
# postpr(target=c(3), index=c(rep("norm",500),rep("t",500)),
#        sumstat=ss[,1], tol=.1, method="rejection")
# ## human demographic history
# require(abc.data)
# data(human)
# ## five R objects are loaded. See ?human and vignette("abc") for details.
# ## illustrate the informativeness of two summary statistics: mean and
# ## variance of Tajima's D
# par(mfcol = c(1,3))
# boxplot(stat.3pops.sim[,"pi"]~models, main="Mean nucleotide diversity")
# boxplot(stat.3pops.sim[,"TajD.m"]~models, main="Mean Tajima's D")
# boxplot(stat.3pops.sim[,"TajD.v"]~models, main="Var in Tajima's D")
# ## model selection with ABC for the European population
# modsel.it <- postpr(stat.voight["italian",], models, stat.3pops.sim, tol=.05, method="mnlogistic")
# summary(modsel.it)
# 
# ## summary.abc Summaries of posterior samples generated by ABC algortithms
# ## Calculates simple summaries of posterior samples: the minimum and maximum, the weighted mean, median, mode, and credible intervals.
# summary(rej, unadj = FALSE, intvl = .95, print = TRUE,
#         digits = max(3, getOption("digits")-3))
# summary(lin, unadj = FALSE, intvl = .95, print = TRUE,
#         digits = max(3, getOption("digits")-3))
# summary(linhc, unadj = FALSE, intvl = .95, print = TRUE,
#         digits = max(3, getOption("digits")-3))
# 
# ## summary.cv4abc Calculates the cross-validation prediction error
# ## This function calculates the prediction error from an object of class "cv4abc" for each parameter and tolerance level.
# 
# summary(cv.rej, print = TRUE, digits = max(3, getOption("digits")-3))
# summary(cv.lin, print = TRUE, digits = max(3, getOption("digits")-3))
# summary(cv.linhc, print = TRUE, digits = max(3, getOption("digits")-3))
# ##summary.cv4postpr Confusion matrix and misclassification probabilities of model
# ## no cv4postpr
# 
# ## summary.gfit Calculates the p-value of the goodness-of-fit test
# summary(oubmbmb.gfit)
# ##summary.postpr Posterior model probabilities and Bayes factors
# ## no postpr
