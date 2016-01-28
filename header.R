# Name: header.R
# Auth: u.niazi@imperial.ac.uk
# Date: 27/01/2016
# Desc: header file with global variables and functions




############## functions
getSequenceVariance = function(seq, prior=c(A=1, T=1, G=1, C=1), size=1){
  ## internal functions
  # get alpha values for dirichlet posterior
  getAlpha = function(seq, prior=c(A=1, T=1, G=1, C=1)){
    #a = letterFrequency(seq, letters='ATGC', OR=0, as.prob=F)
    alpha = seq + prior #a + prior
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
  
  #a = getAlpha(DNAString(as.character(seq)), prior)
  a = getAlpha(seq, prior)
  p = getPosterior(a)
  r = getPosteriorPredict(p, size)
  # get variance by adding variance of each binomial component of the posterior predictive data
  return(sum(apply(r, 2, var)))
}


f_getSeq = function(pile){
  # split the pile up data frame on positions
  df = data.frame(P=pile$pos, N=pile$nucleotide, C=pile$count)
  l = split(df, df$P)
  seq = array(0, dim=c(length(l), 4), dimnames = list(names(l), c('A', 'T', 'G', 'C')))
  for (i in 1:nrow(seq)){
    s2 = seq[i,]
    # sequence at current nucleotide position
    x = l[[i]]
    s = x$C
    names(s) = x$N
    # match the positions of names
    p = match(names(s), names(s2))
    # if a non standard residue is present then the position will not match
    # and will have an NA there
    p2 = !is.na(p)
    p = p[p2]
    s2[p] = s[p2]
    seq[i,] = s2
  }
  return(seq)
}

