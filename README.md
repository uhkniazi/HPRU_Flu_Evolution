# HPRU_Flu_Evolution
Viral passage experiment NGS data to check for virus evolution

# header.R
contains the functions and variables used in other scripts

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


# bam_variance.R
Loads a list of sorted indexed bamfiles in the given folder, calculates the variance vector for each bam file. The vector for the variance of a bam file, has components equal to the areas of the sequence that have coverage in the bam file, and the variance at that position. Various plots are produced as output.  


