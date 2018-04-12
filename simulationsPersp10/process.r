bigres=matrix(NA,length(dir(pattern="*.rda")),12)

k=1
for(I in dir(pattern="*.rda")) {

load(I)
sa=sapply(1:1000,function(j) {
jnk=res[[j]]$rl
ifelse(!inherits(jnk,"try-error"),
sum(abs(res[[j]]$rl$k-res[[1]]$ktrue)),NA)})
bigres[k,1]=as.numeric(substr(I,4,6))
bigres[k,3]=ifelse(as.numeric(substr(I,2,2))==2,0,1)
    bigres[k,2]=as.numeric(substr(I,9,9))
    Ti=bigres[k,2]
    bigres[k,1:Ti+3]=res[[1]]$ktrue
bigres[k,12]=mean(sa==0,na.rm=TRUE)
k=k+1
}
colnames(bigres)=c("n","T","highlysep","k1","k2","k3","k4","k5","k6","k7","k8","pCorrect")

bigres3=matrix(NA,length(dir(pattern="k4")),14)

k=1
for(I in dir(pattern="k4")) {

load(I)

bigres3[k,1]=as.numeric(substr(I,4,6))
bigres3[k,2]=as.numeric(substr(I,9,9))
bigres3[k,3]=ifelse(as.numeric(substr(I,2,2))==2,0,1)
    Ti=bigres3[k,2]
    bigres3[k,1:Ti+3]=res[[1]]$ktrue
# MSE xi2

xi=array(1,c(4,4,3)) 
    sigma=xi
    sepp=ifelse(as.numeric(substr(I,2,2))==3,4,2.5)
for(j in 1:4) {
for(u in 1:3) {
xi[j,1:j,u]=seq(0,sepp*j,length=j)}}
       
bigres3[k,12]=mean(sapply(1:1000, function(j) {
jnk=res[[j]]$rl
ifelse(!is.null(jnk) & !inherits(jnk,"try-error"), ifelse(any(jnk$k==4),
sum((res[[j]]$rl$xi[4,1:4,]-xi[4,1:4,])^2,na.rm=TRUE),NA),NA)}),na.rm=TRUE)
  
bigres3[k,13]=mean(sapply(1:1000, function(j) {
jnk=res[[j]]$rlmax
ifelse(!is.null(jnk) & !inherits(jnk,"try-error"),
sum((res[[j]]$rlmax$xi[4,1:4,]-xi[4,1:4,])^2,na.rm=TRUE),NA)}),na.rm=TRUE)

    bigres3[k,14]=mean(sapply(1:1000, function(j) {
jnk=res[[j]]$rlfixed
ifelse(!is.null(jnk)  & !inherits(jnk,"try-error"),
sum((res[[j]]$rlfixed$xi[4,1:4,]-xi[4,1:4,])^2,na.rm=TRUE),NA)}),na.rm=TRUE)


k=k+1
}

colnames(bigres3)=c("n","T","highlysep","k1","k2","k3","k4","k5","k6","k7","k8",paste("msexi4",c("est","max","fix")))

rm(res)
save.image(file="resultsPI10.RData")
