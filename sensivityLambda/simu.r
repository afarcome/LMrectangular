library(foreach)
library(doSNOW)

r=3
k=c(2,2,2,3)
source("codeKfixedLOG.r")
source("codeKvariable.r")
kmax=4

n=100
Ti=4
pers=12
sepp=2.5

source("doSimu.r")
pers=10
source("doSimu.r")
pers=8
source("doSimu.r")
