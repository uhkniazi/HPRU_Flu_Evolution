# Name: header.R
# Auth: u.niazi@imperial.ac.uk
# Date: 27/01/2016
# Desc: header file with global variables and functions




############## functions
getSequenceVariance = function(seq, prior=c(A=1, T=1, G=1, C=1), size=1){
  ## internal functions
  # get alpha values for dirichlet posterior
  getAlpha = function(seq, prior=c(A=1, T=1, G=1, C=1)){
    a = letterFrequency(seq, letters='ATGC', OR=0, as.prob=F)
    alpha = a + prior
    return(alpha)
  }
  # get posterior theta from posterior dirichlet
  getPosterior = function(alpha, n=1000){
    if(!require(LearnBayes)) stop('Package LearnBayes required')
    p = rdirichlet(n, alpha)
    colnames(p) = names(alpha)
    m = colMeans(p)
    return(m)
  }
  # get posterior predictive values using multinomial random sample
  getPosteriorPredict = function(theta, n=1){
    return(t(rmultinom(1000, n, theta)))
  }
  
  a = getAlpha(DNAString(as.character(seq)), prior)
  p = getPosterior(a)
  r = getPosteriorPredict(p, size)
  # get variance by adding variance of each binomial component of the posterior predictive data
  return(sum(apply(r, 2, var)))
}