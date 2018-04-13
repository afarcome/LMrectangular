# LMrectangular
Code, data analysis and simulation studies for "Rectangular latent Markov models for time-specific   clustering, with an analysis of the well being of                           nations"

This project reports R code for implementing the methods and reproducing
data analysis and simulation studies in Anderson et al. "Rectangular
latent Markov models for time-specific   clustering, with an analysis
of the well being of                           nations"

The following libraries shall be installed in R before using code from
this project:

library(snipEM)

library(mvtnorm)

library(mclust)

library(compiler)

library(foreach)

library(doSNOW)

library(xtable)

Code for the general methods, with a readme giving explanations, is in folder mainCode. 

To reproduce the main simulation studies in the paper, please read the
readme files and execute code in folders: simulationsPers12, simulationsPers10, simulationsPers8. 

To reproduce the sensitivity analyses with respect to B and lambda, please refer to folders sensitivityB and sensitivityLambda

To download the data and/or reproduce the data analysis, please refer ti folder dataAnalysis

