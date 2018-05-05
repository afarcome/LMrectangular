Files in this folder can be used to reproduce the data analysis.

Data are in R workspace HDIdata.rda, which contains a data frame with

> names(HDIdata)

[1] "year"           "Country"        "per.capita.GNI" "Life_Exp" "Exp_Educ"       "Pop"

where pop is Population and the other variable names are
auto-explicative. 

Results are in R workspace rl.rda, which can be reproduced (it takes
about 1.5 hours) by executing 

source("dataAnalysis.r")

This replaces file rl.rda with an identical one. 

PLEASE NOTE THAT THE FILE ABOVE IS SET UP SIMPLY FOR PURPOSES OF REPRODUCING 
THE DATA ANALYSIS IN THE PAPER SWIFTLY. WITH A NEW DATA SET A MULTISTART APPROACH WITH A FEW DIFFERENT INITIAL STARTING SOLUTIONS (POSSIBLY WELL CHOSEN) 
ARE RECOMMENDED. THE HOT STARTING / PRETRAINING inits0 IS NEVERTHELESS A GOOD STARTING POINT IN GENERAL. 

The results workspace rl.rda contains:

y: array with data arranged for analysis with rlm or rlm.fixed; where
y[,1,] is GNI; y[,2,] is life expectancy and y[,3,] is expected
education.

x: vector of lambda values used (0 to 0.5)

nams: vector of country names, where nams[i] is the name of the
country whose values are recorded in y[i,,]

rl: list of results of function rlm, where rl[[i]] corresponds to a
penalty parameter x[i]

crit: (lik[i]-lik[i-1])/(mean(lik)(x[i]-x[i-1])), where x[which.n(crit)] is the optimal lambda as defined in the paper

A few notes: 

-- in this reproducible example we have used simple deterministic starting solutions for demonstrative purposes. The same (seemingly optimal) solution is obtained after trying several different starting solutions. 

-- these data set actually makes a good example of basically all problems that can occurr with rlm and rlm.fixed: 

a) latent transitions are unlikely, making some PI[,,i,j] close to zero when i!=j or, worse, PI[,,i,i] close to 1. This might lead the code to numerical issues which can be avoided using a high tolerance for convergence 

b) due to the low sample size (LM models are more common for larger n), local optima are abundant for this data set. Refer also to 
Bartolucci and Farcomeni (2015) - Biometrics and Farcomeni (2017) - Biometrical Journal on this point. 

c) probably due to the large number of clusters in the face of a low sample size, the likelihood is not very steep as a function of PI. This results in certain values of PI not being updated much at the M step (when lambda>0), unless parameter hits is slightly increased from its default. 

Secondly, 

source("sensitivityB.r")

reproduces the sensitivity analysis with respect to tuning parameter B, with 
possible values Bs=c(10,50,100,200,500).

As result a list is obtained, named rlB, where

rlB[[i]] is the outcome when lambda=0 and B=Bs[i].

