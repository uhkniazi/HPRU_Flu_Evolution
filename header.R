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



getSequenceParameters = function(ivSeq, cRefBase, prior=c(A=1/2, T=1/2, G=1/2, C=1/2), iSize=1000){
  if(!require(LearnBayes)) stop('R Package LearnBayes required')
  ## internal functions
  # get alpha values for dirichlet posterior
  getAlpha = function(seq, prior=c(A=1/2, T=1/2, G=1/2, C=1/2)){
    #a = letterFrequency(seq, letters='ATGC', OR=0, as.prob=F)
    alpha = seq + prior #a + prior
    return(alpha)
  }
  # calculate dirichlet variance 
  getDirichletVariance = function(alpha){
    al0 = sum(alpha)
    denom = al0^2 * (al0+1)
    ret = sapply(alpha, function(x){
      return(x * (al0 - x) / denom)
    })
    return(ret)
  }
  # get posterior theta from posterior dirichlet
  getPosterior = function(alpha, n=1000){
    p = rdirichlet(n, alpha)
    colnames(p) = names(alpha)
    #m = colMeans(p)
    return(p)
  }
  # posterior preditive distribution based on the theta from dirichlet posterior
  # set n to adjust size of sum alpha
  getPosteriorPredict = function(theta, n=1){
    ret = t(apply(theta, 1, function(x) rmultinom(1, n, x)))
    colnames(ret) = colnames(theta)
    return(ret)
  }
  # posterior gamma, a component of the dirichlet distribution
  ## see gelman P 583 and bayesian computations with R page 66
  getPosteriorGamma = function(alpha.scale, base, n=1000, prior){
    # adjust alpha by removing dirichlet prior
    alpha.scale = alpha.scale - prior
    i = which(names(alpha.scale) == base)
    alpha.new = c(alpha.scale[i], sum(alpha.scale[-i]))
    names(alpha.new) = c('Base', 'Other')
    ## set a non-informative jeffery's prior for gamma distributed rate
    jef.prior = c(alpha=0.5, beta=1)
    # adjust the new alpha for gamma posterior
    alpha.new = alpha.new + jef.prior['alpha']
    # convert the alpha to rate per 1000
    alpha.new = (alpha.new/sum(alpha.new)) * 1000
    rg = sapply(seq_along(alpha.new), function(x) {
      return(rgamma(n, alpha.new[x], jef.prior['beta']))
    })
    colnames(rg) = names(alpha.new)
    return(rg)
  }
  ###### processing steps
  ## get posterior values
  a = getAlpha(ivSeq, prior)
  # get posterior dirichlet variance
  var = getDirichletVariance(a)
  # get posterior via simulation
  p = getPosterior(a, iSize)
  #   r = getPosteriorPredict(colMeans(p), 1000)
  #   # get variance by adding variance of each binomial component of the posterior predictive data
  #   (sum(apply(r, 2, var)))
  var = sum(var)
  theta = colMeans(p)[cRefBase]
  if (is.na(theta)) {rate=c('Base'= NA, 'Other'=NA) } else {
    rate = round(colMeans(getPosteriorGamma(a, cRefBase, iSize, prior)), 2)}
  return(c(theta=theta, var=var, lambda.base=rate['Base'], lambda.other=rate['Other'], ivSeq))
}


logit = function(p) log(p/(1-p))
logit.inv = function(p) {exp(p)/(exp(p)+1) }


