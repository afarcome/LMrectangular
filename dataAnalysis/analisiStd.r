rm(list=ls())
load("HDIdata.rda")

y=array(NA,c(164,9,3)) 
nams=HDIdata[1:164,2]

u=unique(HDIdata[,1])
for(j in 1:9) {
for(h in 1:3) {
y[,j,h]=HDIdata[HDIdata[,1]==u[j],h+2]}}
#y[,,1]=log(y[,,1])

y=y[,-c(6:8),]
source("codeKvariable.r")

yo=y
yo[,,1]=(y[,,1]-apply(y[,,1],2,mean))
yo[,,2]=(y[,,2]-apply(y[,,2],2,mean))
yo[,,3]=(y[,,3]-apply(y[,,3],2,mean))

x=c(0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.5)

rl0=rlm(y,0,4)
rl=list()
l=length(x)
rl[[1]]=rl0
for(j in 2:l) {
rl[[j]]=rlm(y,x[j],4,inits=rl[[j-1]])}

save.image()

liks=sapply(1:l,function(j) rl[[j]]$lik)
jpeg("screeplot.jpg",units="in",width=7,height=7,res=1000)
plot(x,liks,type="l")
dev.off()

save.image(file="rl.rda")
