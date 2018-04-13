load("resultsPI8.RData")

                                        # 2 & 3

colnames(bigres2)[3]="s"
colnames(bigres2)[4:11]=paste("$k_",1:8,"$",sep="")
w=which(bigres2[,3]==0)
bigres2[w,3]=2.5
bigres2[w,3]=4
bigres2=bigres2[,-14]
w=which(is.na(bigres2[,17]))
bigres2[w,17]=bigres2[w,16]
bigres2=bigres2[,-16]
library(xtable)
print(xtable(bigres2,digits=c(0,0,0,1,rep(0,8),rep(3,6))),include.rownames=FALSE)
