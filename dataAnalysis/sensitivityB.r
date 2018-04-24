rm(list=ls())
load("HDIdata.rda")

y=array(NA,c(164,6,3)) 
nams=HDIdata[1:164,2]

u=unique(HDIdata[,1])
for(j in 1:6) {
for(h in 1:3) {
y[,j,h]=HDIdata[HDIdata[,1]==u[j],h+2]}}

source("codeKvariable.r")

library(snipEM)
library(mvtnorm)
set.seed(12345)


load("rl.rda")
Bs=c(10,50,100,200,500)
rlB=list()
inits=rl[[7]]
inits$k=rep(3,6)
iter=1
for(B in Bs) {
    rlB[[iter]]=rlm(y,0.3,4,inits=inits,hits=50,B=B)
iter=iter+1
}

rlB2=list()
inits$k=rep(4,6)
iter=1
for(B in Bs) {
    rlB2[[iter]]=rlm(y,0.3,4,inits=inits,hits=50,B=B)
iter=iter+1
}


save.image(file="B.rda")
