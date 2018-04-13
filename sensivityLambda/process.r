bigres=matrix(0,length(dir(pattern="*.rda")),7)

k=1
Id=dir(pattern="*.rda")
for(I in Id) {
load(I)
for(i in 1:6) {
sa=sapply(1:1000,function(j) {
jnk=res[[j]][[i]]
ifelse(!inherits(jnk,"try-error"),
       sum(abs(jnk$k-res[[1]]$ktrue)),NA)})
bigres[k,i]=mean(sa==0,na.rm=TRUE)
}
lambdas=c(0.01,0.05,0.1,0.2,0.35,0.5)


for(j in 1:1000) {
    liks=sapply(res[[j]][1:6],function(x) x$lik)
    crit=diff(liks)/(mean(liks)*diff(lambdas))
        w=which.min(crit)+1
        bigres[k,7]=bigres[k,7]+ifelse(sum(abs(res[[j]][[w]]$k-res[[1]]$ktrue))==0,1/1000,0)}
        k=k+1}

rm(res)
save.image(file="L.RData")
