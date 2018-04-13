Files in this folder can be used to reproduce the main simulation study. All cases with persistence 8 are generated.

Simulations are based on parallel computing. First of all, edit ncpu.r to set the number of CPUs to be used (ncpu = 16). 
Note: ncpu must be at least 2. 

Then, in R, one shall execute 

source("simu.r")

source("simu2.r")

source("simu3.r")

source("simu4.r")

each file corresponds to all settings for a given configuration of latent groups k1,k2,..,kT

Note: this is time consuming and generates about 2.5Gb of R workspaces. 

After execution as above, one shall run: 

source("process.r")

source("table2.r")

to generate the following data frames: 

bigres: subset of Table 1 in the main paper, reporting on the proportion of times the correct configuration is selected by
the proposed method

bigres2: MSE for the centroid when k=2 and when k=3 for different methods 

bigres3: MSE for the centroid when k=4 

