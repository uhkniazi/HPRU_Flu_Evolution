logit = function(p) log(p/(1-p))
logit.inv = function(p) exp(p)/(exp(p)+1) 

mFile = matrix(NA, nrow=width(refseq), ncol=2, dimnames=list(1:width(refseq), c('theta', 'var')))

m = match(rownames(file17), rownames(mFile))

mFile[m,] = file17
mFile17 = mFile

mFile = matrix(NA, nrow=width(refseq), ncol=2, dimnames=list(1:width(refseq), c('theta', 'var')))
m = match(rownames(file15), rownames(mFile))
mFile[m,] = file15
mFile15 = mFile

mTheta = cbind(theta15=mFile15[,'theta'], theta17=mFile17[,'theta'])

m15 = na.omit(mFile15)

quantile(m15[,'var'], 0:10/10)
var.groups = cut(m15[,'var'], breaks = quantile(m15[,'var'], 0:10/10), labels = 0:9)

stripchart(logit(m15[,'theta']) ~ var.groups, method='jitter', pch=20)
table(var.groups)

x = as.numeric(rownames(m15[var.groups != '9',]))
plot(x, logit(m15[var.groups != '9', 'theta']), pch=20, cex=0.5)
plot(rownames(m15), logit(m15[, 'theta']), pch=20, cex=0.5)

m17 = na.omit(mFile17)

quantile(m17[,'var'], 0:10/10)
var.groups = cut(m17[,'var'], breaks = quantile(m17[,'var'], 0:10/10), labels = 0:9)

stripchart((m17[,'theta']) ~ var.groups, method='jitter', pch=20)
table(var.groups)

x = as.numeric(rownames(m17[var.groups != '9',]))
plot(x, logit(m17[var.groups != '9', 'theta']), pch=20, cex=0.5)
plot(rownames(m17), logit(m17[, 'theta']), pch=20, cex=0.5)

x = as.numeric(rownames(m15[var.groups != '9' & m15[,'theta'] < 0.2,]))

seqcheck = refseq[[1]][x[10]:x[length(x)]]
writeXStringSet(DNAStringSet(seqcheck), 'Results/seq.fasta')


file13 = na.omit(t(file13))
head(file13)

file11 = f_getMutations(file.choose(), refseq)
file11 = na.omit(t(file11))
head(file11)

file12 = f_getMutations(file.choose(), refseq)
file12 = na.omit(t(file12))
head(file12)

file15 = f_getMutations(file.choose(), refseq)
file15 = na.omit(t(file15))
head(file15)

file17 = f_getMutations(file.choose(), refseq)
file17 = na.omit(t(file17))
head(file17)

file19 = f_getMutations(file.choose(), refseq)
file19 = na.omit(t(file19))
head(file19)

file16 = f_getMutations(file.choose(), refseq)
file16 = na.omit(t(file16))
head(file16)


par(mfrow=c(2,1))
plot(log(file13[,1]), pch=20, cex=0.5, main='file13', ylab='theta' )
plot(log(file13[,2]), pch=20, cex=0.5, main='file13', ylab='log var' )

par(mfrow=c(2,1))
plot(log(file11[,1]), pch=20, cex=0.5, main='file11', ylab='theta' )
plot(log(file11[,2]), pch=20, cex=0.5, main='file11', ylab='log var' )

par(mfrow=c(1,1))
plot(1:nrow(file13), log(file13[,1]), pch=20, cex=0.5)
points(1:nrow(file13), log(file11[,1]), pch=20, cex=0.5, col='red')
points(1:nrow(file12), log(file12[,1]), pch=20, cex=0.5, col='blue')
points(1:nrow(file15), log(file15[,1]), pch=20, cex=0.5, col='green')

plot(1:nrow(file19), log(file19[,1]), pch=20, cex=0.5)
points(1:nrow(file17), log(file17[,1]), pch=20, cex=0.5, col='blue')

plot(1:nrow(file15), log(file15[,1]), pch=20, cex=0.5, col='grey')
points(1:nrow(file17), log(file17[,1]), pch=20, cex=0.5, col='blue')

plot(1:nrow(file15), (file15[,1]), pch=20, cex=0.5, col='grey')
points(1:nrow(file17), (file17[,1]), pch=20, cex=0.5, col='blue')

var15 = file15[,'var']
var17 = file17[,'var']

theta15 = file15[,'theta']
theta17 = file17[,'theta']

plot(1:nrow(file15), (var15), pch=20, cex=0.5, col='grey')
points(1:nrow(file17), (var17), pch=20, cex=0.5, col='blue')

plot(density(log(var15)))
plot(density(theta15))
plot(density(theta17))

coplot(log(theta15[names(theta15)]) ~ log(theta17[names(theta15)]) | log(var15[names(theta15)]) * log(var17[names(theta15)]))
m = cbind(mFile15, mFile17)
head(m)
m = na.omit(m)
par(mfrow=c(1,1))
colnames(m) = paste(colnames(m), c(15, 15, 17, 17), sep='')
head(m)
coplot(logit(m[,'theta15']) ~ logit(m[,'theta17']) | log(m[,'var15']) * log(m[,'var17']))


x = as.numeric(rownames(m15[var.groups != '9' & m15[,'theta'] < 0.2,]))
seqcheck = refseq[[1]][x[10]:x[length(x)]]
writeXStringSet(DNAStringSet(seqcheck), 'Results/seq.fasta')



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
  getPosteriorPredict = function(theta, n=1){
    return(t(rmultinom(1000, n, theta)))
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
  return(list(theta=theta, var=var))
}



seq.1 = DNAString('AAAAAAAATTTTTTTCCCCCCCCGGGGGGGGG')
seq.1 = DNAString('AAAAAAAATTTTTTTTTTTTTTTTTTGC')
seq.1 = DNAString('AAAAAAAAAAATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTGC')
seq.1 = DNAString('AAAATGCCCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA')

seq = letterFrequency(seq.1, letters='ATGC', OR=0, as.prob=F)
getSequenceParameters(seq, 'A')

# temp.R

library(LearnBayes)
library(Biostrings)

##############################################################
seq.1 = DNAString('AAAAAAAATTTTTTTCCCCCCCCGGGGGGGGG')
seq.1 = DNAString('AAAAAAAATTTTTTTTTTTTTTTTTTGC')
seq.1 = DNAString('AAAAAAAAAAATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTGC')
seq.1 = DNAString('AAAATGCCCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA')
getAlpha = function(seq, prior=c(A=1/2, T=1/2, G=1/2, C=1/2)){
  #a = letterFrequency(seq, letters='ATGC', OR=0, as.prob=F)
  alpha = seq + prior #a + prior
  return(alpha)
}

getDirichletVariance = function(alpha){
  al0 = sum(alpha)
  denom = al0^2 * (al0+1)
  ret = sapply(alpha, function(x){
    return(x * (al0 - x) / denom)
  })
  return(ret)
}

seq = letterFrequency(seq.1, letters='ATGC', OR=0, as.prob=F)
getAlpha(seq)
sum(getDirichletVariance(getAlpha(seq)))
10 * log10(sum(getDirichletVariance(getAlpha(seq))))

10 * log10(sum(apply(getPosterior(getAlpha(seq), 1000), 2, var)))

a = getAlpha(seq)
p = getPosterior(a)
r = getPosteriorPredict(colMeans(p), 1000)
# get variance by adding variance of each binomial component of the posterior predictive data
(sum(apply(r, 2, var)))

# get posterior predictive values using multinomial random sample
getPosteriorPredict = function(theta, n=1){
  return(t(rmultinom(1000, n, theta)))
}

# get posterior theta from posterior dirichlet
getPosterior = function(alpha, n=1000){
  if(!require(LearnBayes)) stop('Package LearnBayes required')
  p = rdirichlet(n, alpha)
  colnames(p) = names(alpha)
  #m = colMeans(p)
  return(p)
}

getPosteriorPredict = function(theta, n=1){
  return(t(rmultinom(1000, n, theta)))
}

post = getPosterior(getAlpha(seq))
colMeans(post)

getBase.prob(getAlpha(seq), 'A')

##############################################################
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
  return(mean(all(p[,i] > p[,-i])))
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

sim.election=function()
{
  winner=rbinom(51,1,Obama.win.probs)
  sum(EV*winner)
}
sim.EV=replicate(1000,sim.election())

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



