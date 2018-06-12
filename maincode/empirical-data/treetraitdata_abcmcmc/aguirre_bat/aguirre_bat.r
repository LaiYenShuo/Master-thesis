rm(list=ls())
setwd("~/Dropbox/TerryLai/treetraitdata/aguirre_bat/")
library(phangorn)

Gmatrix<-read.table("treeJSTOR-66.txt")
Gmatrix<-Gmatrix+t(Gmatrix)
diag(Gmatrix)<-1
Dmatrix<-2*(1-Gmatrix)
tree<-upgma(Dmatrix)

#tree<-read.tree("~/Dropbox/TerryLai/treetraitdata/aguirre_bat/intree")
trait <- as.data.frame(read.table(file="traitJSTOR-66.txt"))
colnames(trait) <- c("mass","head_height","head_length","max_bite_force")
resptrait<-trait$mass
predtrait1<-trait$head_height
predtrait2<-trait$head_length
iterations = 50000
errorbound = 100

source("~/Dropbox/TerryLai/R_code/empirical_abc_mcmc.r")

setwd("~/Documents/empirical_abcmcmc/")
save.image("aguirre_bat.RData")
#oubmcirresult$distS0S1
#oubmbmresult$distS0S1
# oubmcirresult <- abcmcmc(tree = tree, resptrait = resptrait, predtrait1 = predtrait1, predtrait2 = predtrait2, iterations = iterations, errorbound = errorbound, model = "oubmcir")
# ououcirresult <- abcmcmc(tree = tree, resptrait = resptrait, predtrait1 = predtrait1, predtrait2 = predtrait2, iterations = iterations, errorbound = errorbound, model = "ououcir")