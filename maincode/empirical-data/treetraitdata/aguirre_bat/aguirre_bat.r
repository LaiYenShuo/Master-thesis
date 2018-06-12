rm(list=ls())
setwd("~/Dropbox/TerryLai/treetraitdata/aguirre_bat/")
library(phangorn)
sum.stat<-function(trait=trait,tree=tree){
  pic.trait<-pic(x=trait,phy=tree)
  return(c(mean(pic.trait),sd(pic.trait)))
}

Gmatrix<-read.table("treeJSTOR-66.txt")
Gmatrix<-Gmatrix+t(Gmatrix)
diag(Gmatrix)<-1
Dmatrix<-2*(1-Gmatrix)
tree<-upgma(Dmatrix)

#tree<-read.tree("~/Dropbox/TerryLai/treetraitdata/aguirre_bat/intree")
tree$root.edge<-0
root <- 0
#plot(tree)

dataset <- as.data.frame(read.table(file="traitJSTOR-66.txt"))
colnames(dataset) <- c("mass","head_height","head_length","max_bite_force")

resptrait<-dataset$mass
predtrait1<-dataset$head_height
predtrait2<-dataset$head_length

raw.sum.stat.y <- sum.stat(trait = dataset$mass, tree=tree)
raw.sum.stat.x1 <- sum.stat(trait = dataset$head_height, tree=tree)
raw.sum.stat.x2 <- sum.stat(trait = dataset$head_length, tree=tree)
raw.sum.stat<-c(raw.sum.stat.y,raw.sum.stat.x1,raw.sum.stat.x2)

source("~/Dropbox/TerryLai/R_code/empirical_analysis.r")
setwd("~/Documents/empiricalresult/")
save.image("aguirre_bat.RData")