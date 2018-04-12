library(snipEM)
library(mvtnorm)
source("codeKfixedLOG.r")

obj=function(x,lambda,xi,sigma,k,kmax,n,Ti,pi,PI,done) {

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
qu[,t,j]=apply(cbind(qu[,t,j],(qu[,t-1,d]+log(PI[k[t-1],k[t],d,j]))),1,sumlog)}
qu[,t,j]=qu[,t,j]+dmvnorm(y[,t,],xi[k[t],j,],diag(sigma[k[t],j,]^2),log=TRUE)}


if(k[t-1]==1) {
qu[,t,j]=dmvnorm(y[,t,],xi[k[t],j,],diag(sigma[k[t],j,]^2),log=TRUE)+qu[,t-1,1]+log(PI[k[t-1],k[t],1,j])
}

}}

if(k[Ti]>1) {
liks=apply(qu[,Ti,1:k[Ti]],1,sumlog)}
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

chooseD=function(y,lambda,kmax,B=1,cmax=1,maxit=1,inits=NULL,verbose=FALSE,debug=FALSE,tol=1e-4) {

n=nrow(y)
Ti=ncol(y)
r=dim(y)[3]

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
k=inits$k
}

if(is.null(inits)) {
k=rep(round(kmax/2),Ti)
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
}}
}

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
output=tol
for(iters in 1:maxit) {

if(all(k==1)) {

xi=apply(y,3,mean)
sigma=apply(y,3,sd)
lkold=lik
lik=dmvnorm(matrix(y,ncol=3),xi,diag(sigma^2),log=TRUE)

}
    
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
}

## MM part ### 
if(verbose) {print("## update k")}
if(B>0) {
for(b in 1:B) {
for(ti in 1:Ti) {
for(u in c(-1,1)) {
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
output=max(c(output,abs(lik-likc)))
if(verbose) print(output)
}}}}}

if(verbose) {print(paste("iters:",iters,sep=" "))}}

return(1/(ncol(y)*output))}

library(compiler)

chooseD=cmpfun(chooseD)

