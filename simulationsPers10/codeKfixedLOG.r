library(snipEM)
library(mvtnorm)

likco.fixed=function(xi,sigma,pi,PI,k,kmax,n,Ti) {

qu=array(NA,c(n,Ti,kmax))

for(j in 1:k[1])
{qu[,1,j]=dmvnorm(y[,1,],xi[k[1],j,],diag(sigma[k[1],j,]^2),log=TRUE)+log(pi[k[1],j])}

for(t in 2:Ti) {
for(j in 1:k[t]) {
if(k[t-1]>1) {
qu[,t,j]=(qu[,t-1,1]+log(PI[k[t-1],k[t],1,j]))
for(d in 2:k[t-1]) {
qu[,t,j]=apply(cbind(qu[,t,j],(qu[,t-1,d]+log(PI[k[t-1],k[t],d,j]))),1,sumlog)}
qu[,t,j]=qu[,t,j]+dmvnorm(y[,t,],xi[k[t],j,],diag(sigma[k[t],j,]^2),log=TRUE)}


if(k[t-1]==1) {
qu[,t,j]=dmvnorm(y[,t,],xi[k[t],j,],diag(sigma[k[t],j,]^2),log=TRUE)+qu[,t-1,1]+log(PI[k[t-1],k[t],1,j])
}

}}

if(k[Ti]>1) {
liks=apply(qu[,Ti,1:k[Ti]],1,sumlog)}
if(k[Ti]==1) {liks=qu[,Ti,1]}
res=sum(liks)

qub=qu
qub[,Ti,]=0
for(t in (Ti-1):1) {
for(c in 1:k[t]) {
qub[,t,c]=dmvnorm(y[,t+1,],xi[k[t+1],1,],diag(sigma[k[t+1],1,]^2),log=TRUE)+qub[,t+1,1]+log(PI[k[t],k[t+1],c,1])
if(k[t+1]>1) {
for(j in 2:k[t+1]) {
jnk=dmvnorm(y[,t+1,],xi[k[t+1],j,],diag(sigma[k[t+1],j,]^2),log=TRUE)+qub[,t+1,j]+log(PI[k[t],k[t+1],c,j])
qub[,t,c]=apply(cbind(qub[,t,c],jnk),1,sumlog)
}}}}

return(list(lik=res,liks=liks,qub=qub,qu=qu))}

rlm.fixed=function(y,k,tol=1e-4,maxit=Inf,inits=NULL,verbose=FALSE,debug=FALSE,exceed=NULL) {

n=nrow(y)
Ti=ncol(y)
kmax=max(k)
r=dim(y)[3]

if(all(k==1)) {

xi=apply(y,3,mean)
sigma=apply(y,3,sd)
lik=dmvnorm(matrix(y,ncol=3),xi,diag(sigma^2),log=TRUE)
aic=-2*lik+2*2*r
bic=-2*lik+log(n)*2*r
return(list(V=array(1,c(n,Ti,1)),pi=1,PI=1,xi=xi,sigma=sigma,lik=lik,aic=aic,bic=bic))

}


xi=array(NA,c(kmax,kmax,3)) 
sigma=xi

pi=array(NA,c(kmax,kmax))
PI=array(NA,c(kmax,kmax,kmax,kmax))
## inits ## 

if(!is.null(inits)) {
xi=inits$xi
sigma=inits$sigma
PI=inits$PI
pi=inits$pi
}

if(is.null(inits)) {
for(j in 1:kmax) {
for(u in 1:r) {
km=kmeans(as.vector(y[,,u]),j)
xi[j,1:j,u]=km$centers
sigma[j,1:j,u]=sqrt(km$withinss/km$size)}
if(j>1) {
for(l in 1:r) {
o=order(xi[j,1:j,l])
xi[j,1:j,l]=xi[j,1:j,l][o]
sigma[j,1:j,l]=sigma[j,1:j,l][o]}
}}}

if(is.null(inits$PI)) {

pi[1,1]=PI[,1,,1]=1
for(j in 2:kmax) {
pi[j,1:j]=1/j}
for(h in 2:kmax) {
for(j in 1:kmax) {
if(j>1) {
PI[j,h,1:j,1:h]=1
if(h>1 & j>1) {
diag(PI[j,h,1:j,1:h])=8+j
PI[j,h,1:j,1:h]=PI[j,h,1:j,1:h]/apply(PI[j,h,1:j,1:h],1,sum)
}}
if(j==1) {
PI[j,h,1,1:h]=1
PI[j,h,1,1]=8+j
PI[j,h,1:j,1:h]=PI[j,h,1:j,1:h]/sum(PI[j,h,1:j,1:h])
}}}
}
## E-step ##

lst=likco.fixed(xi,sigma,pi,PI,k,kmax,n,Ti)

qu=lst$qu
lik=lst$lik
liks=lst$liks
qub=lst$qub
V=qu+qub
for(j in 1:kmax) {
V[,,j]=sweep(V[,,j],1,liks)
}
V=exp(V)

Z=array(0,c(kmax,kmax,kmax,kmax))
for(j in 1:kmax) {
for(h in 1:kmax) {
trans=k[-Ti]==j & k[-1]==h 
if(any(trans)) {
wtrans=which(trans)
for(c in 1:j) {
for(d in 1:h) {
sw=matrix(y[,wtrans+1,],ncol=r)
ap=matrix(dmvnorm(sw,xi[h,d,],diag(sigma[h,d,]^2),log=TRUE),ncol=length(wtrans),byrow=FALSE)
Z[j,h,c,d]=sumlog(log(PI[j,h,c,d])+sweep(qu[,wtrans,c]+qub[,wtrans+1,d]+ap,1,liks))}}
}}}

lkold=lik*2 
iters=1
if(is.null(exceed)) {exceed=Inf}
while(lik-lkold>tol & iters<maxit & lik<exceed) {

#if(verbose) {print("## update initial probabilities")}

if(k[1]>1) {
pi[k[1],1:k[1]]=apply(V[,1,1:k[1]],2,sum)
pi[k[1],]=pi[k[1],]/sum(pi[k[1],1:k[1]])
}
if(debug) {
lkpartial=likco.fixed(xi,sigma,pi,PI,k,kmax,n,Ti)$lik
if(lik-lkpartial>0) {print("problems with pi!!!")}}

#if(verbose) {print("## update hidden transitions")}

for(j in 1:kmax) {
for(h in 1:kmax) {
trans=k[-Ti]==j & k[-1]==h 
if(any(trans)) {
if(j>1 & h>1) {
PI[j,h,1:j,1:h]=exp(sweep(Z[j,h,1:j,1:h],1,apply(Z[j,h,1:j,1:h],1,sumlog)))}
if(j==1) {
PI[j,h,1,1:h]=exp(Z[j,h,1,1:h]-sumlog(Z[j,h,1,1:h]))}
}}}
if(debug) {
lkpartold=lkpartial
lkpartial=likco.fixed(xi,sigma,pi,PI,k,kmax,n,Ti)$lik
if(lkpartold-lkpartial>0) {print(paste("problems with PI at iteration",iters,"!!!"))}}

#if(verbose) {print("## update xi & sigma")}

for(ks in unique(k)) {
for(l in 1:r) {
for(d in 1:ks) {
xi[ks,d,l]=sum(V[,which(k==ks),d]*y[,which(k==ks),l])/sum(V[,which(k==ks),d])
sigma[ks,d,l]=sqrt(sum(V[,which(k==ks),d]*(y[,which(k==ks),l]-xi[ks,d,l])^2)/(sum(V[,which(k==ks),d])))
}}
if(ks>1) {
o=order(xi[ks,1:ks,1])
xi[ks,1:ks,]=xi[ks,1:ks,][o,]
sigma[ks,1:ks,]=sigma[ks,1:ks,][o,]}
}

if(debug) {
lkpartold=lkpartial
lkpartial=likco.fixed(xi,sigma,pi,PI,k,kmax,n,Ti)$lik
if(lkpartold-lkpartial>0) {print(paste("problems with xi or sigma at iteration",iters,"!!!"))}}

if(verbose) {print("## E-step")}
lst=likco.fixed(xi,sigma,pi,PI,k,kmax,n,Ti)

lkold=lik 

qu=lst$qu
lik=lst$lik
liks=lst$liks
qub=lst$qub
V=qu+qub
for(j in 1:kmax) {
V[,,j]=sweep(V[,,j],1,liks)
}
V=exp(V)

Z=array(0,c(kmax,kmax,kmax,kmax))
for(j in 1:kmax) {
for(h in 1:kmax) {
trans=k[-Ti]==j & k[-1]==h 
if(any(trans)) {
wtrans=which(trans)
for(c in 1:j) {
for(d in 1:h) {
sw=matrix(y[,wtrans+1,],ncol=r)
ap=matrix(dmvnorm(sw,xi[h,d,],diag(sigma[h,d,]^2),log=TRUE),ncol=length(wtrans),byrow=FALSE)
Z[j,h,c,d]=sumlog(log(PI[j,h,c,d])+sweep(qu[,wtrans,c]+qub[,wtrans+1,d]+ap,1,liks))}}
}}}



## drop empty clusters ## 

#for(t in 1:Ti) {
#if(k[t]>1) {
#if(any(apply(V[,t,1:k[t]],2,sum)==0)) {
#print("DROPPED EMPTY CLUSTER!!!") 
#k[t]=k[t]-1}}}

if(verbose) {print("lk diff")}
if(verbose) {print(c(lkold,lik,lik-lkold))}
iters=iters+1
}


np=sum(2*r*unique(k))+k[1]-1

for(j in 1:kmax) {
for(h in 1:kmax) {
trans=k[-Ti]==j & k[-1]==h 
if(any(trans)) {
np=np+j*(h-1)
}}}

aic=-2*lik+2*np
bic=-2*lik+log(n)*np

return(list(V=V,pi=pi,PI=PI,xi=xi,sigma=sigma,lik=lik,aic=aic,bic=bic,liks=liks))}

library(compiler)

rlm.fixed=cmpfun(rlm.fixed)


