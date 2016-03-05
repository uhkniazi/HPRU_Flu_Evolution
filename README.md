# HPRU_Flu_Evolution
Viral passage experiment NGS data to check for virus evolution

# header.R
contains the functions and variables used in other scripts

## Utility functions
###logit = function(p) log(p/(1-p))
###logit.inv = function(p) {exp(p)/(exp(p)+1) }


# NOTE: There is a problem with getPosteriorPredict of this function and needs fixing. See the getPosteriorPredict function in getSequenceParameters function.
## getSequenceVariance
### ARGS:  
seq = named vector A, T, G and C, giving the counts for each nucleotide.
prior = uniform prior by default, but can be changed
size = default 1, but can be increased to 1000 to represent read depth.
### DESC:
the function models the variable seq, as a multinomial variable with a uniform dirichlet prior, and uses a conjugate model to calculate the posterior dirichlet parameters. These parameters are used to get the posterior predictive values for the sequence using a multinomial model, where a sample is generated to calculate the variance for each component of the multinomial data. The variance is added together to get the total variance for the sequence.  
### RETS:
numeric variable with the total variance for the sequence seq.  

Contains internal functions to calculate:  
#### getAlpha
get the conjugate posterior alpha for dirichlet distribution
#### getPosterior
get the random sample from the dirichlet posterior and use the sample mean to parameterize the posterior predictive multinomial model.
#### getPosteriorPredict
get the multinomial simulated sample using dirichlet posterior parameters 
  
## f_getSeq
### ARGS:
pile = object returned from Rsamtools:pileup function  
### DESC:
splits each sequence position into a component of a list, where each component includes a count of nucleotides and the counts of those nucleotides at that position. This is used to generate a sequence matrix (rows = unique number of positions, columns = A, T, G, C). 
### RETS:
matrix with the counts of A, T, G, C at each position of the sequence.  

## getSequenceParameters
### ARGS:
ivSeq = named integer vector A, T, G and C with number of nucleotides at each position. Typically this will be the output from the function f_getSeq  
cRefBase = single character that represents the reference base  
prior = default is jefferey's dirichlet prior i.e. 1/2 for each alpha.  
iSize = simulation size, default = 1000  
### DESC:
This function requires R package LearnBayes to use the function rdirichlet. It uses ivSeq and prior to get the alpha for the dirichlet posterior. Calculates dirichlet variance and adds this vector (other possibilities include perhaps l1 or l2 norm, try these later). Calculates the posterior mean by simulating from rdirichlet, and saves the posterior mean 'theta' for the reference base. Finally the posterior for the non-refrence base can be a proportion very close to zero (as mutations are rare events), so we can calculate a rate instead of a proportion. We do this by modelling the components of the posterior dirichlet alpha as independent gamma variables with parameters (shape=alpha, rate=1) (see getPosteriorGamma for more details below). We calculate two rates, one for the reference base (which is not really required, and may be removed to speed up things) and a rate for non-reference base (sum of 3 non-reference base individual components of alpha).  
### RETS:
A named integer vector with components  
theta = proportion of times reference base was seen.  
var = total dirichlet variance.  
lambda.base = the reference base rate per 1000  
lambda.other = non-reference base rate per 1000  
A, T, G and C = the actual number of times the bases were seen, basically similar to input ivSeq.  

### internal functions  
### getAlpha = function(seq, prior=c(A=1/2, T=1/2, G=1/2, C=1/2))  
Calculates posterior alpha for dirichlet posterior  

### getDirichletVariance = function(alpha)
analytically calculates dirichlet variance vector, see Gelman 2008 Bayesian Data Analysis appendix section for detailed formula.  

### getPosterior = function(alpha, n=1000)
Returns the matrix of 1000 siumulations from the rdirichlet. May replace this for an analytic calculation, to increase speed.  

### getPosteriorPredict = function(theta, n=1)
For each vector of the simulation from dirichlet posterior using getPosterior, it simulates a multinomial variable. Not currently used.  

### getPosteriorGamma = function(alpha.scale, base, n=1000, prior)
[see gelman P 583 and bayesian computations with R page 66]. The function first substracts the prior (which was a dirichlet prior) from the alpha.scale (we could just use ivSeq instead to speed things up). It adds the components of the non-reference base from the alpha to create one component and the other component is the reference base (base in the function argument). It then adds a Jefferey's prior of 0.5 for the alpha while beta=1 for the Gamma distribution. This is because the likelihood is a poisson rate, and the poisson lambda is ditributed as a gamma variable. We convert the rate to per 1000, simulate 1000 from rgamma.


# bam_variance.R [NOTE: there is a minor problem with this script as getSequenceVariance needs fixing for calculation of getPosteiorPredict]. Use next script bam_mutation_probability.R instead.
Loads a list of sorted indexed bamfiles in the given folder, calculates the variance vector for each bam file. The vector for the variance of a bam file, has components equal to the areas of the sequence that have coverage in the bam file, and the variance at that position. Various plots are produced as output.  

# bam_mutation_probability.R
[add details here after deciding on which plots and parameters to report]

