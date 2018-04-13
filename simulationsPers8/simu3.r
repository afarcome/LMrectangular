library(foreach)
library(doSNOW)

r=3
ktot=c(4,4,3,3,2,2,4,4)
source("codeKfixedLOG.r")
source("codeKvariable.r")
kmax=4

n=c(100,200)
Ti=c(4,6,8)
pers=8
sepp=c(2.5,4)


ex=expand.grid(n,Ti,sepp)
nex=nrow(ex)
for(j in 1:nex) {
    n=ex[j,1]
    Ti=ex[j,2]
    sepp=ex[j,3]
    k=ktot[1:Ti]
source("doSimu.r")
}
