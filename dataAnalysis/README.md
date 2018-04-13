Files in this folder can be used to reproduce the data analysis.

Data are in R workspace HDIdata.rda, which contains a data frame with

> names(HDIdata)

[1] "year"           "Country"        "per.capita.GNI" "Life_Exp" "Exp_Educ"       "Pop"

where pop is Population and the other variable names are
auto-explicative. 

Results are in R workspace rl.rda, which can be reproduced (it takes
about 1 hour) by executing 

source("dataAnalysis.r")

This replaces file rl.rda with an identical one.

The results workspace rl.rda contains:

y: array with data arranged for analysis with rlm or rlm.fixed; where
y[,1,] is GNI; y[,2,] is life expectancy and y[,3,] is expected
education.

x: vector of lambda values used (0 to 0.5)

nams: vector of country names, where nams[i] is the name of the
country whose values are recorded in y[i,,]

rl: list of results of function rlm, where rl[[i]] corresponds to a
penalty parameter x[i]

crit: (lik[i]-lik[i-1])/(mean(lik)(x[i]-x[i-1])), where x[min(crit)+1] is the optimal lambda as defined in the paper



