library(foreach)
library(doSNOW)

source("codeKfixedLOG.r")
source("codeKvariable.r")
source("ncpu.r")

cl=makeCluster(ncpu, type = "SOCK")
registerDoSNOW(cl)

res=foreach(iter = 1:1000)%dopar%{
library(snipEM)
library(mvtnorm)
set.seed(iter+100)
pi=array(NA,c(kmax,kmax))
PI=array(NA,c(kmax,kmax,kmax,kmax))
pi[1,1]=PI[,1,,1]=1
for(j in 2:kmax) {
pi[j,1:j]=1/j}
for(h in 2:kmax) {
for(j in 1:kmax) {
if(j>1) {
PI[j,h,1:j,1:h]=1
if(h>1 & j>1) {
diag(PI[j,h,1:j,1:h])=pers+j
PI[j,h,1:j,1:h]=PI[j,h,1:j,1:h]/apply(PI[j,h,1:j,1:h],1,sum)
}}
if(j==1) {
PI[j,h,1,1:h]=1
PI[j,h,1,1]=pers+j
PI[j,h,1:j,1:h]=PI[j,h,1:j,1:h]/sum(PI[j,h,1:j,1:h])
}}}

inits=list()
inits$pi=pi
inits$PI=PI

xi=array(1,c(kmax,kmax,3)) 
sigma=xi
for(j in 1:kmax) {
for(u in 1:r) {
xi[j,1:j,u]=seq(0,sepp*j,length=j)}}
inits$xi=xi
inits$sigma=sigma
## generate data ###

y=array(NA,c(n,Ti,3))
prob=pi[k[1],1:k[1]]
u=sample(k[1],n,replace=TRUE,prob=prob)
y[,1,1]=rnorm(n,xi[k[1],u,1])
y[,1,2]=rnorm(n,xi[k[1],u,2])
y[,1,3]=rnorm(n,xi[k[1],u,3])
for(ti in 2:Ti) {
prob=PI[k[ti-1],k[ti],u,1:k[ti]]
for(i in 1:n) {
u[i]=sample(k[ti],1,replace=TRUE,prob=prob[i,])}
y[,ti,1]=rnorm(n,xi[k[ti],u,1])
y[,ti,2]=rnorm(n,xi[k[ti],u,2])
y[,ti,3]=rnorm(n,xi[k[ti],u,3])
}

## estimate models

rl2=try(rlm.fixed(y,k))
rlmin=try(rlm.fixed(y,k=rep(min(k),Ti)))
rlmax=try(rlm.fixed(y,k=rep(max(k),Ti)))
rl=try(rlm(y,lambda=0.1,kmax,cmax=1))
res=list()

# eval 
res$rl=rl
res$rlfixed=rl2
res$rlmin=rlmin
res$rlmax=rlmax
res$ktrue=k

return(res)}

stopCluster(cl)

sp=ifelse(sepp>2.5,3,2)
save(res,file=paste("s",sp,"n",n,"TT",Ti,"r",r,"k",paste(k,collapse=""),".rda",sep=""))




