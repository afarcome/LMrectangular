bigres=matrix(NA,length(dir(pattern="*.rda"))+1,5)

k=1
Id=dir(pattern="*.rda")
for(I in Id) {
load(I)
for(i in 1:5) {
sa=sapply(1:1000,function(j) {
jnk=res[[j]][[i]]
ifelse(!inherits(jnk,"try-error"),
       sum(abs(jnk$k-res[[1]]$ktrue)),NA)})
its=sapply(1:1000,function(j) res[[j]][[i]]$iters)
bigres[k,i]=mean(sa==0,na.rm=TRUE)
}
k=k+1}

rm(res)
save.image(file="B.RData")
