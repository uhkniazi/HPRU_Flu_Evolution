# temp.R

library(LearnBayes)
library(Biostrings)
seq.1 = DNAString('AAAAAAAAAATGCC')
#letterFrequency(seq, 'ATCG', as.prob=T)
source('~/Dropbox/Home/Data/R/My_Libraries/NGS_functions.R')

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
p = c(10, 1, 5, 0)
s1 = paste0(sample(c('A', 'T', 'G', 'C'), 30, replace = T, prob = p/sum(p) ), collapse = '')
getSequenceVariance(s1)


## loading bam files
csBamfile = file.choose()
oGAbam = readGAlignments(csBamfile, param = ScanBamParam(what = 'seq'))

oGRbam = reduce(as(oGAbam, 'GRanges'))
strand(oGRbam) = '*'
oGRbam = reduce(oGRbam)

oPile = pileup(csBamfile, scanBamParam = ScanBamParam(which = as(oGAbam[1], 'GRanges')), 
               pileupParam = PileupParam(distinguish_strands = F))

# process chunks at a time
oGRbam.bin = f_split_GRanges_by_length(oGRbam, 1000)

oPile = pileup(csBamfile, scanBamParam = ScanBamParam(which = oGRbam.bin[1]), 
               pileupParam = PileupParam(distinguish_strands = F))

seq = f_getSeq(oPile)

v = apply(seq, 1, getSequenceVariance)



df = data.frame(P=oPile$pos, N=oPile$nucleotide, C=oPile$count)
l = split(df, df$P)
x = l[['9']]
s2 = c(A=0, T=0, G=0, C=0)
s = x$C
names(s) = x$N
s = c(A=1, T=10, X=11)
s2 = c(A=0, T=0, G=0, C=0)
i = match(names(s), names(s2))
i2 = is.na(i)
i = i[i2]
s2[i] = s[i2]

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



