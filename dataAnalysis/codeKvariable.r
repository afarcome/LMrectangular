library(snipEM)
library(mvtnorm)
source("codeKfixedLOG.r")
source("chooseD.r")
obj=function(x,lambda,xi,sigma,k,kmax,n,Ti,pi,PI,done) {

if(k[1]>1) {
pi[k[1],2:k[1]]=exp(x[1:(k[1]-1)])/(1+sum(exp(x[1:(k[1]-1)])))
pi[k[1],1]=1-sum(pi[k[1],2:k[1]])
if(any(pi[k[1],1:k[1]]==0)) {
w0=which(pi[k[1],]==0)
pi[k[1],w0]=1e-323
pi[k[1],1:k[1]]=pi[k[1],1:k[1]]/sum(pi[k[1],1:k[1]])}
x=x[-c(1:(k[1]-1))]}
for(itd in 1:nrow(done)) {
j=done[itd,1]
h=done[itd,2]
for(ku in 1:j) {
jnk=x[(ku-1)*(h-1)+1:(h-1)]
PI[j,h,ku,2:h]=exp(jnk)/(1+sum(exp(jnk)))
PI[j,h,ku,1]=1-sum(PI[j,h,ku,2:h])
if(any(PI[j,h,ku,1:h]==0)) {
w0=which(PI[j,h,ku,1:h]==0)
PI[j,h,ku,w0]=1e-323
PI[j,h,ku,1:h]=PI[j,h,ku,1:h]/sum(PI[j,h,ku,1:h])}}
x=x[-c(1:(j*(h-1)))]
}
PI[PI<1e-323]=1e-323
-forw(xi,sigma,pi,PI,k,kmax,n,Ti,lambda)}

forw=function(xi,sigma,pi,PI,k,kmax,n,Ti,lambda) {

qu=array(NA,c(n,Ti,kmax))

for(j in 1:k[1])
{qu[,1,j]=dmvnorm(y[,1,],xi[k[1],j,],diag(sigma[k[1],j,]^2),log=TRUE)+log(pi[k[1],j])}

for(t in 2:Ti) {
for(j in 1:k[t]) {
if(k[t-1]>1) {
qu[,t,j]=(qu[,t-1,1]+log(PI[k[t-1],k[t],1,j]))
for(d in 2:k[t-1]) {
    cb=cbind(qu[,t,j],(qu[,t-1,d]+log(PI[k[t-1],k[t],d,j])))
    cb[cb < -743.7]=-743.7
qu[,t,j]=apply(cb,1,sumlog)}
qu[,t,j]=qu[,t,j]+dmvnorm(y[,t,],xi[k[t],j,],diag(sigma[k[t],j,]^2),log=TRUE)}


if(k[t-1]==1) {
qu[,t,j]=dmvnorm(y[,t,],xi[k[t],j,],diag(sigma[k[t],j,]^2),log=TRUE)+qu[,t-1,1]+log(PI[k[t-1],k[t],1,j])
}

}}

    if(k[Ti]>1) {
        cb=qu[,Ti,1:k[Ti]]
        cb[cb < -743.7]=-743.7
liks=apply(cb,1,sumlog)}
if(k[Ti]==1) {liks=qu[,Ti,1]}
lik=sum(liks)
jnk=pi[k[1],1:k[1]]
jnk2=jnk[jnk>1e-15]
lik=lik+lambda*n*sum(jnk2*log(jnk2))
for(j in 2:Ti) {
jnk=jnk%*%PI[k[j-1],k[j],1:k[j-1],1:k[j]]
jnk2=jnk[jnk>1e-15]
lik=lik+lambda*n*sum(jnk2*log(jnk2))
}
return(lik)}

likco=function(xi,sigma,pi,PI,k,kmax,n,Ti) {

qu=array(NA,c(n,Ti,kmax))

for(j in 1:k[1])
{qu[,1,j]=dmvnorm(y[,1,],xi[k[1],j,],diag(sigma[k[1],j,]^2),log=TRUE)+log(pi[k[1],j])}

for(t in 2:Ti) {
for(j in 1:k[t]) {
if(k[t-1]>1) {
qu[,t,j]=(qu[,t-1,1]+log(PI[k[t-1],k[t],1,j]))
for(d in 2:k[t-1]) {
    cb=cbind(qu[,t,j],(qu[,t-1,d]+log(PI[k[t-1],k[t],d,j])))
    cb[cb< -743.7]=-743.7
qu[,t,j]=apply(cb,1,sumlog)}
qu[,t,j]=qu[,t,j]+dmvnorm(y[,t,],xi[k[t],j,],diag(sigma[k[t],j,]^2),log=TRUE)}


if(k[t-1]==1) {
qu[,t,j]=dmvnorm(y[,t,],xi[k[t],j,],diag(sigma[k[t],j,]^2),log=TRUE)+qu[,t-1,1]+log(PI[k[t-1],k[t],1,j])
}

}}

    if(k[Ti]>1) {
        cb=qu[,Ti,1:k[Ti]]
        cb[cb < -743.7]=-743.7
liks=apply(cb,1,sumlog)}
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
    cb=cbind(qub[,t,c],jnk)
    cb[cb < -743.7]=-743.7
qub[,t,c]=apply(cb,1,sumlog)
}}}}

return(list(lik=res,liks=liks,qub=qub,qu=qu))}

rlm=function(y,lambda,kmax,B=200,Di=NULL,cmax=1,tol=1e-4,maxit=Inf,hits=20,inits=NULL,verbose=FALSE,debug=FALSE) {

n=nrow(y)
Ti=ncol(y)
r=dim(y)[3]

xi=array(NA,c(kmax,kmax,r)) 
sigma=xi

pi=array(NA,c(kmax,kmax))
PI=array(NA,c(kmax,kmax,kmax,kmax))
## inits ## 

if(!is.null(inits)) {
xi=inits$xi
sigma=inits$sigma
PI=inits$PI
pi=inits$pi
k=inits$k
}

if(is.null(inits)) {
library(mclust)
k=rep(round(kmax/2),Ti)
for(j in 1:Ti) {
k[j]=max(min(kmax,Mclust(y[,j,],modelNames="VVI")$G),2)
}
for(j in 1:kmax) {
for(u in 1:r) {
km=kmeans(as.vector(y[,,u]),j,iter.max=1000,nstart=10)
xi[j,1:j,u]=km$centers
sigma[j,1:j,u]=sqrt(km$withinss/km$size)}
if(j>1) {
for(l in 1:r) {
o=order(xi[j,1:j,l])
xi[j,1:j,l]=xi[j,1:j,l][o]
sigma[j,1:j,l]=sigma[j,1:j,l][o]}
}}
}

if(is.null(inits$PI)) {

pi[1,1]=PI[,1,,1]=1
for(j in 2:kmax) {
km=kmeans(y[,1,],j,nstart=10,iter.max=1000)
o=order(km$center[,1])
pi[j,1:j]=prop.table(table(km$cluster))[o]}
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

lst=likco(xi,sigma,pi,PI,k,kmax,n,Ti)

qu=lst$qu
lik=lst$lik
liks=lst$liks
qub=lst$qub
V=qu+qub
for(j in 1:kmax) {
V[,,j]=sweep(V[,,j],1,liks)
}
V=exp(V)
jnk=pi[k[1],1:k[1]]
jnk2=jnk[jnk>1e-15]
lik=lik+lambda*n*sum(jnk2*log(jnk2))
for(j in 2:Ti) {
jnk=jnk%*%PI[k[j-1],k[j],1:k[j-1],1:k[j]]
jnk2=jnk[jnk>1e-15]
lik=lik+lambda*n*sum(jnk2*log(jnk2))
}

lkold=lik*2 
iters=1
    flagK=TRUE
    Z=array(0,c(kmax,kmax,kmax,kmax))
while(lik-lkold>tol & iters<maxit) {

if(all(k==1)) {

xi=apply(y,3,mean)
sigma=apply(y,3,sd)
lkold=lik
lik=dmvnorm(matrix(y,ncol=r),xi,diag(sigma^2),log=TRUE)

}

if(any(k>1)) {
    if(verbose) {print("## update pi & PI")}

    if(lambda>0) {
done=NULL
initM=log(pi[k[1],2:k[1]]/pi[k[1],1])
for(j in 1:kmax) {
for(h in 1:kmax) {
trans=k[-Ti]==j & k[-1]==h 
if(any(trans)) {
if(debug) {print(PI[j,h,,])}
done=rbind(done,c(j,h))
for(ku in 1:j) {
initM=c(initM,log(PI[j,h,ku,2:h]/PI[j,h,ku,1]))
}}}}

op=try(optim(initM,obj,lambda=lambda,xi=xi, sigma=sigma, k=k, kmax=kmax, n=n, Ti=Ti,
pi=pi, PI=PI, done=done, control=list(maxit=hits)),silent=!debug)
x=initM
if(!inherits(op,"try-error")){x=op$par}

pi[k[1],2:k[1]]=exp(x[1:(k[1]-1)])/(1+sum(exp(x[1:(k[1]-1)])))
pi[k[1],1]=1-sum(pi[k[1],2:k[1]])
if(any(pi[k[1],1:k[1]]==0)) {
w0=which(pi[k[1],]==0)
pi[k[1],w0]=1e-16
pi[k[1],1:k[1]]=pi[k[1],1:k[1]]/sum(pi[k[1],1:k[1]])}
x=x[-c(1:(k[1]-1))]
for(itd in 1:nrow(done)) {
j=done[itd,1]
h=done[itd,2]
for(ku in 1:j) {
jnk=x[(ku-1)*(h-1)+1:(h-1)]
PI[j,h,ku,2:h]=exp(jnk)/(1+sum(exp(jnk)))
PI[j,h,ku,1]=1-sum(PI[j,h,ku,2:h])
if(any(PI[j,h,ku,1:h]==0)) {
w0=which(PI[j,h,ku,1:h]==0)
PI[j,h,ku,w0]=1e-16
PI[j,h,ku,1:h]=PI[j,h,ku,1:h]/sum(PI[j,h,ku,1:h])}}
x=x[-c(1:(j*(h-1)))]
}
    }


    if(lambda==0) {

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

if(k[1]>1) {
pi[k[1],1:k[1]]=apply(V[,1,1:k[1]],2,sum)
pi[k[1],]=pi[k[1],]/sum(pi[k[1],1:k[1]])
}
for(j in 1:kmax) {
for(h in 1:kmax) {
trans=k[-Ti]==j & k[-1]==h 
if(any(trans)) {
if(j>1 & h>1) {
PI[j,h,1:j,1:h]=exp(sweep(Z[j,h,1:j,1:h],1,apply(Z[j,h,1:j,1:h],1,sumlog)))}
if(j==1) {
PI[j,h,1,1:h]=exp(Z[j,h,1,1:h]-sumlog(Z[j,h,1,1:h]))}
}}}

    }
if(verbose) {print("## update xi & sigma")}

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

#if(debug) {
#lkpartold=lkpartial
#lkpartial=likco(xi,sigma,pi,PI,k,kmax,n,Ti)$lik
#if(lkpartold-lkpartial>0) {print(paste("problems with xi or sigma at iteration",iters,"!!!"))}}

if(verbose) {print("## E-step")}
lst=try(likco(xi,sigma,pi,PI,k,kmax,n,Ti),silent=!debug)

if(inherits(lst,"try-error")) {break}

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
jnk=pi[k[1],1:k[1]]
jnk2=jnk[jnk>1e-15]
lik=lik+lambda*n*sum(jnk2*log(jnk2))
for(j in 2:Ti) {
jnk=jnk%*%PI[k[j-1],k[j],1:k[j-1],1:k[j]]
jnk2=jnk[jnk>1e-15]
lik=lik+lambda*n*sum(jnk2*log(jnk2))
}}

## MM part ### 
if(verbose) {print("## update k")}
if(!flagK) {B=round(B/2)}
flagK=FALSE
    if(B>0) {
        Di2=Di
        if(is.null(Di)) {
        current=list()
        current$pi=pi
        current$PI=PI
        current$xi=xi
        current$sigma=sigma
        current$k=k
       Di2=chooseD(y,lambda,kmax,verbose=FALSE,tol=tol,inits=current)}
for(b in 1:B) {
ti=sample(Ti,1) 
u=sample(c(-1,1),1)
inits=list()
inits$xi=xi
inits$sigma=sigma
inits$PI=PI
inits$pi=pi
kc=k
if(k[ti]==1) {kc[ti]=kc[ti]+1}
if(k[ti]==kmax) {kc[ti]=kc[ti]-1}
if(k[ti]>1 & k[ti]<kmax) {kc[ti]=kc[ti]+u}
rl=try(rlm.fixed(y,kc,maxit=cmax,inits=inits,exceed=sum(liks)),silent=!debug)
if(!inherits(rl,"try-error")) {
jnk=rl$pi[kc[1],1:kc[1]]
jnk2=jnk[jnk>1e-15]
likc=rl$lik+lambda*n*sum(jnk2*log(jnk2))
for(j in 2:Ti) {
jnk=jnk%*%rl$PI[kc[j-1],kc[j],1:kc[j-1],1:kc[j]]
jnk2=jnk[jnk>1e-15]
likc=likc+lambda*n*sum(jnk2*log(jnk2))
}

pb=min(c(exp(-log(b+1)*(lik-likc)/Di2),1))
if(runif(1)<pb) {k=kc
flagK=TRUE
xi=rl$xi
sigma=rl$sigma
PI=rl$PI
pi=rl$pi
lik=likc
liks=rl$liks
V=rl$V
if(verbose) {print("kandidate accepted!")
             print(paste("k=",k,sep="",collapse=","))}
}}}}

if(verbose) {print(paste("iters:",iters,sep=" "))}
if(verbose) {print("lk diff")}
if(verbose) {print(c(lkold,lik,lik-lkold))}
iters=iters+1
}

if(all(k==1)) {

xi=apply(y,3,mean)
sigma=apply(y,3,sd)
lik=dmvnorm(matrix(y,ncol=r),xi,diag(sigma^2),log=TRUE)
aic=-2*lik+2*2*r
bic=-2*lik+log(n)*2*r

}

if(any(k>1)) {
np=sum(2*r*unique(k))+k[1]-1

for(j in 1:kmax) {
for(h in 1:kmax) {
trans=k[-Ti]==j & k[-1]==h 
if(any(trans)) {
np=np+j*(h-1)
}}}

aic=-2*lik+2*np
bic=-2*lik+log(n)*np}

return(list(V=V,pi=pi,PI=PI,xi=xi,sigma=sigma,lik=lik,aic=aic,bic=bic,k=k))}

library(compiler)

rlm=cmpfun(rlm)

