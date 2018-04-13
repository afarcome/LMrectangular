Files in this folder can be used to reproduce the sensitivity analysis
with respect to the penalty parameter lambda. 

Simulations are based on parallel computing. First of all, edit ncpu.r to set the number of CPUs to be used (ncpu = 16). 
Note: ncpu must be at least 2. 

Then, in R, one shall execute 

source("simu.r")

After execution as above, one shall run: 

source("process.r")

to generate the following data frame: 

bigres: 3 by 7 matrix where bigres[a,b] gives the proportion of times
lambda[b] lead to the correct selection of latent configuration
when persistence was pers[a]; with lambda=c(0.01,0.05,0.1,0.2,0.35,0.5) if b<=6 and optimal lambda if b=7; and for
pers=c(8,10,12). The remaining parameters are set as
k=c(2,2,2,3), n=100, separation s=2.5.


