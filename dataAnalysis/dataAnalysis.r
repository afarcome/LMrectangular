rm(list=ls())
load("HDIdata.rda")

y=array(NA,c(164,6,3)) 
nams=HDIdata[1:164,2]

u=unique(HDIdata[,1])
for(j in 1:6) {
for(h in 1:3) {
y[,j,h]=HDIdata[HDIdata[,1]==u[j],h+2]}}

source("codeKvariable.r")

x=c(0,0.02,0.05,0.07,0.1,0.15,0.2,0.25,0.3,0.35,0.375,0.4,0.425,0.45,0.5)

set.seed(1234)
rl0=rlm(y,0,4,cmax=2)
rl=list()
l=length(x)
rl[[1]]=rl0
for(j in 2:l) {
rl[[j]]=rlm(y,x[j],4,inits=rl[[j-1]],cmax=2)}

liks=sapply(1:l,function(j) rl[[j]]$lik)
crit=diff(liks)/(mean(liks)*diff(x))


save.image(file="rl.rda")
