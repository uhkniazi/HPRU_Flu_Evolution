# temp.R

library(LearnBayes)
library(Biostrings)
seq.1 = DNAString('AAAAAAAAAATGCC')
#letterFrequency(seq, 'ATCG', as.prob=T)
#source('~/Dropbox/Home/Data/R/My_Libraries/NGS_functions.R')

#alpha = letterFrequency(seq.1, letters='ATGC', OR=0, as.prob=F)
getAlpha = function(seq, base){
  a = letterFrequency(seq, letters='ATGC', OR=0, as.prob=F)
  i = which(names(a) == base)
  alpha = c(a[i], sum(a[-i]))
  names(alpha) = c(base, 'Other')
  return(alpha + 1)
}

getBase.prob = function(alpha, base){
  p = rdirichlet(1000, alpha) 
  i = which(names(alpha) == base)
  return(mean(p[,i] > p[,-i]))
}

getColumn.var = function(theta){
  #theta = min(c(theta, 1-theta))
  return(var(rbinom(1000, 1, theta)))
}

seq.1 = DNAString('AAAAAAAAAATGCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAAAAAAAAAA')
seq.2 = DNAString('AAAATGCCCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA')
seq.3 = DNAString('AAAAAAAATTTTTTTCCCCCCCCGGGGGGGGG')
seq.4 = DNAString('AAAAAAAACCCCCCCCC')
s = list(seq.1, seq.2, seq.3, seq.4)
s = DNAStringSet(list(seq.1, seq.2, seq.3, seq.4))
ref = c('A', 'A', 'T', 'C')


v = sapply(1:4, function(x) {
  a = getAlpha(seq = DNAString(as.character(s[x])), base = ref[x])
  p = getBase.prob(a, ref[x])
  getColumn.var(p)
})


v2 = sapply(1:4, function(x){
seq = s[x]
a = letterFrequency(seq, letters='ATGC', OR=0, as.prob=F)
alpha = a+1
p = rdirichlet(1000, alpha)
colnames(p) = colnames(alpha)
head(p)
m = colMeans(p)
m
rmultinom(1, 100, m)
r = t(rmultinom(1000, 1, m))
#r = p
apply(r, 2, var)
sum((apply(r, 2, var)))
#sum(rowSums(cov(r)))
})


getAlpha = function(seq, prior=c(A=1, T=1, G=1, C=1)){
  a = letterFrequency(seq, letters='ATGC', OR=0, as.prob=F)
  alpha = a + prior
  return(alpha)
}

getPosterior = function(alpha, n=1000){
  if(!require(LearnBayes)) stop('Package LearnBayes required')
  p = rdirichlet(n, alpha)
  colnames(p) = names(alpha)
  m = colMeans(p)
  return(m)
}

getPosteriorPredict = function(theta, n=1){
  return(t(rmultinom(1000, n, theta)))
}

sum(apply(r, 2, var))


theta = rdirichlet(10, alpha)
theta
attach(election.2008)
prob.Obama=function(j)
{
  p=rdirichlet(5000,
               500*c(M.pct[j],O.pct[j],100-M.pct[j]-O.pct[j])/100+1)
  mean(p[,2]>p[,1])
}

Obama.win.probs=sapply(1:51,prob.Obama)
j = 4
v = rep(NA, length.out=length(Obama.win.probs))
v = sapply(1:51, function(j) var(rbinom(1000, 1, Obama.win.probs[j])))
m = sapply(1:51, function(j) mean(rbinom(1000, 1, Obama.win.probs[j])))
 
w1 = rbinom(1000, 1, Obama.win.probs[j])
mean(w1); var(w1)
alpha = c(1, 1) # uniform prior
p = rdirichlet(1000, alpha)
mean(p[,1] > p[,2])


## generate some sequences
p = c(10, 10, 10, 10)
s1 = paste0(sample(c('A', 'T', 'G', 'C'), 30, replace = T, prob = p/sum(p) ), collapse = '')
getSequenceVariance(s1)



