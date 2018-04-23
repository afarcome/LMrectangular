rm(list=ls())
load("HDIdata.rda")

y=array(NA,c(164,6,3)) 
nams=HDIdata[1:164,2]

u=unique(HDIdata[,1])
for(j in 1:6) {
for(h in 1:3) {
y[,j,h]=HDIdata[HDIdata[,1]==u[j],h+2]}}

source("codeKvariable.r")

x=c(0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5)

library(snipEM)
library(mvtnorm)
set.seed(12345)

rl=list()

inits0=rlm.fixed(y,c(4,4,4,3,3,3))
inits0$k=rep(4,6)

rl[[1]]=rlm(y,0,4,inits=inits0)

for(whl in 2:length(x)) {
    inits=rl[[whl-1]]
    inits$k=rep(3,6)
    rl[[whl]]=rlm(y,x[whl],4,inits=inits,hits=50)
}

liks=sapply(1:length(x),function(j) rl[[j]]$lik)
crit=diff(liks)/(mean(liks)*diff(x))

save.image(file="rl.rda")
