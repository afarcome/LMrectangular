Files in this folder provide the main functions for fitting rectangular latent Markov models. 

source("codeKvariable.r")

provides the following main functions: 

rlm.fixed(y, k, tol = 1e-04, maxit = Inf, inits = NULL, verbose = FALSE, 
    debug = FALSE, exceed = NULL) 

which can be used to fit a rectangular latent Markov model for a *fixed* configuration of latent states k. In function rlm.fixed the arguments are: 

y: n by T by r array of continuous measurements, where n: sample size,
T: number of measurement occasions, r: number of outcomes

k: vector of T integers indicating the number of groups at each
measurement occasion

tol: tolerance for convergence in the objective function (default: 1e-4). Note: in case you get a "Not Summable!" error, consider decreasing tolerance or scaling data appropriately. 

maxit: maximal number of iterations of the EM (default: Inf)

inits: list with the initial solutions, whose elements are
xi, sigma, PI, pi. xi and sigma must be kmax by kmax by r arrays (where
kmax=max(k)); pi a kmax by kmax matrix, PI a kmax by kmax by kmax by
kmax array. Inits default to NULL, which leads to a deterministic initial
solution based on k-means separately at each time occasion.

verbose: if TRUE, some summary information at each iteration is
provided. Defaults to FALSE. 

debug: if TRUE, some more information, mostly for internal checking of
errors, is provided. Defaults to FALSE

exceed: if lik>exceed, stop iterations. Mostly for internal
use. Defaults to NULL, mapping to Inf. 

The rlm.fixed function returns a list with elements:

V: n by T by kmax array, where V[i,t,d] is the posterior probability
that the i-th subject is in latent state d at time t. 

pi: kmax by kmax matrix, where pi[k[1],d] gives the prior probability
of being in latent state d at time T=1. Note that for j different that
k[1] the initial solution is returned, which can be ignored. 

PI: kmax by kmax by kmax by kmax array, where PI[a,b,c,d] gives the
prior probability that a subject would switch from state c at time t-1
to state d at time t, provided that k[t-1]=a and k[t]=b. Note that if
no k[t-1]=a and k[t]=b, PI[a,b,,] corresponds to the initial solution,
which can be ignored. 

xi: kmax by kmax by r array, where xi[a,b,c] gives the b-th estimated
centroid for the c-th outcome when there are a groups. As above, if
all(k!=a), xi[a,,] corresponds to the initial solution. 

sigma: kmax by kmax by r array, where sigma[a,b,c] gives the 
standard deviation for the c-th outcome in the b-th group
when there are a groups. As above, if
all(k!=a), sigma[a,,] corresponds to the initial solution. 

lik: objective function at convergence

aic: Akaike Information Criterion

bic: Bayesian Information Criterion

liks: individual contributions to the likelihood


The other main function is 

rlm(y, lambda, kmax, B = 200, Di = NULL, cmax = 1, tol = 1e-04, 
    maxit = Inf, inits = NULL, verbose = FALSE, debug = FALSE) 

which can be used to fit a rectangular latent Markov models by optimizing a penalized likelihood as described in the paper, in order to simultaneously estimate 
the parameters and the sequence of number of latent groups. 

In function rlm the arguments are as above, plus: 

lambda: penalty parameter

kmax: maximum number of groups to be considered at each time occasion

B: number of candidate solutions to be generated at each time occasion

Di: parameter D for acceptance probability. Defaults to NULL, which
gives the optimal acceptance probability. Set a different value at
your own risk. 

cmax: number of candidate EM iterations-1. Defaults to 1, that is, no
candidate EM iterations: only the current initial candidate solution
is considered. If cmax>1, the current initial candidate solution is
improved through cmax-1 candidate EM iterations. 

The function returns a list with elements as that of rlm.fixed, plus

k: final optimal configuration of the number of latent groups. 


