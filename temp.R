
vMatch = matchPattern('NNN', refseq[[1]])






lt = lapply(lSignificant, function(x) {
  m = getTransitionMatrix(t(x[,c('A', 'T', 'G', 'C')]))
  m = round(m/sum(rowSums(m)), 2)
})
  



mBase = lSignificant[[1]][,c('A', 'T', 'G', 'C')]

mBase = temp[,c('A', 'T', 'G', 'C')]
mBase = t(mBase)

mTrans = matrix(0, 4, 4, dimnames = list(c('A', 'T', 'G', 'C'), c('A', 'T', 'G', 'C')))

temp = sort(mBase[,1], decreasing = T)
temp = mBase[,1]
i = names(which.max(temp))
mTrans[i,] = mTrans[i,] + temp



mDat = lMutation[[1]]
mDat = na.omit(mDat)
mDat = mDat[mDat[,'var.q'] != 3, ]
mDat = mDat[mDat[,'theta.q'] != 1, ]
head(mDat)
mDat = mDat[,c(2, 4)]
hist(mDat[,1])
hist(log(mDat[,1]))
mDat[,1] = log(mDat[,1])
head(mDat)
quantile(mDat[,1])
cut.pts = cut(mDat[,1], breaks = quantile(mDat[,1], 0:10/10), include.lowest = T)
table(cut.pts)
cut.pts = cut(mDat[,1], breaks = quantile(mDat[,1], 0:10/10), include.lowest = T, labels = c(1:10))
table(cut.pts)
mDat = cbind(mDat, cut.pts)
head(mDat)
plot(as.factor(mDat[,'cut.pts']), mDat[,'lambda.other'])
tapply(mDat[,'lambda.other'], mDat[,'cut.pts'], mean)
tapply(mDat[,'lambda.other'], mDat[,'cut.pts'], summary)

# at each cutoff model the gamma variable
plot.variance.cutoffs = function(lambda, beta=1, cut.pts){
  cut.pts = as.numeric(cut.pts)
  pos = max(cut.pts)
  for (i in 1:pos){
    # get the gamma variable
    ind = which(cut.pts == i)
    t = lambda[ind]
    r = range(t)
    s = seq(floor(r[1])-0.5, ceiling(r[2])+0.5, by=1)
    r[1] = floor(r[1])
    r[2] = ceiling(r[2])
    dg = dgamma(r[1]:r[2], mean(t), 1)
    df = table(round(t))
    df = df/sum(df)
    # which distribution can approximate the frequency of reactome terms
    hist(t, prob=T, sub=paste('Distribution of mutation rate at', i), breaks=s,
         xlab='Lambda', ylab='', ylim=c(0, max(dg, df)), main='Gamma Mutation Rate')
    # parameterized on the means
    lines(r[1]:r[2], dg, col='black', type='b')
    points(round(qgamma(0.95, mean(t), beta), 0), 0, pch=20, col='red')
    legend('topright', legend =c('Gamma'), fill = c('black'))
  }
}

plot.variance.cutoffs(mDat[,2], cut.pts = cut.pts)

mDat.var = mDat
mDat = lMutation[[1]]
mDat = na.omit(mDat)
mDat = mDat[mDat[,'var.q'] != 3, ]
mDat = mDat[mDat[,'theta.q'] != 1, ]
head(mDat)
mDat = mDat[,c(1, 4)]
hist(mDat[,1])
head(mDat)
quantile(mDat[,1])
cut.pts = cut(mDat[,1], breaks = quantile(mDat[,1], 0:10/10), include.lowest = T)
table(cut.pts)
cut.pts = cut(mDat[,1], breaks = quantile(mDat[,1], 0:10/10), include.lowest = T, labels = paste('t', 1:10))
table(cut.pts)
mDat = cbind(mDat, cut.pts)
head(mDat)
plot(as.factor(mDat[,'cut.pts']), mDat[,'lambda.other'])
tapply(mDat[,'lambda.other'], mDat[,'cut.pts'], mean)
tapply(mDat[,'lambda.other'], mDat[,'cut.pts'], summary)

mDat.both = cbind(mDat, mDat.var)
head(mDat.both)
colnames(mDat.both) = c('theta', 'lambda', 'cut.theta', 'var', 'lambda', 'cut.var')
mDat = mDat.both[,c(1, 2, 4, 3, 6)]
head(mDat)

coplot(mDat[,'cut.theta'] ~ mDat[,'lambda'] | as.factor(mDat[, 'cut.var']))


midpt = seq(0.05, 0.95, by = 0.1)
prior = c(1, 5.2, 8, 7.2, 4.6, 2.1, 0.7, 0.1, 0, 0)
prior = prior/sum(prior)
p = seq(0, 1, length=500)


plot.diagnostics.2 = function(mDat, ...){
  # remove the last quantile of the variance
  mDat = na.omit(mDat)
  mDat = mDat[mDat[,'var.q'] != 3, ]
  # remove first quantile of theta
  mDat = mDat[mDat[,'theta.q'] != 1, ]
  # plot theta
  plot(mDat[,'theta'], pch=20, cex=0.5, sub='Proportion of Reference', ylab='Theta', xlab='Sequence', ...)
  plot(logit(mDat[,'theta']), pch=20, cex=0.5, sub='Proportion of Reference', ylab='Logit Theta', xlab='Sequence', ...)
  # variance
  plot(mDat[,'var'], pch=20, cex=0.5, sub='Dirichlet Variance', ylab='Variance', xlab='Sequence', ...)
  plot(log(mDat[,'var']), pch=20, cex=0.5, sub='Dirichlet Variance', ylab='Log Variance', xlab='Sequence', ...)
  # plot the base rates
  plot(mDat[,'lambda.other'], logit(mDat[,'theta']), pch=20, cex=0.5, 
       sub='Mutation Rate ~ Gamma(lambda)', ylab='Logit Theta', xlab='Lambda', ...)
  plot(mDat[,'lambda.other'], mDat[,'theta'], pch=20, cex=0.5,
       sub='Mutation Rate ~ Gamma(lambda)', ylab='Theta', xlab='Lambda', ...)
  plot(mDat[,'lambda.base'], logit(mDat[,'theta']), pch=20, cex=0.5,
       sub='Reference Base Rate ~ Gamma(lambda)', ylab='Logit Theta', xlab='Lambda', ...)
  plot(mDat[,'lambda.base'], mDat[,'theta'], pch=20, cex=0.5,
       sub='Reference Base Rate ~ Gamma(lambda)', ylab='Theta', xlab='Lambda', ...)
  coplot(mDat[,'theta'] ~ mDat[,'lambda.other'] | as.factor(mDat[,'var.q']), columns = 3)
  ## plot the density and fit distribution, for mutation rate
  t = mDat[,'lambda.other']
  r = range(t)
  s = seq(floor(r[1])-0.5, ceiling(r[2])+0.5, by=1)
  r[1] = floor(r[1])
  r[2] = ceiling(r[2])
  #dn = dnbinom(r[1]:r[2], size = mean(t), mu = mean(t))
  #dp = dpois(r[1]:r[2], mean(t))
  dg = dgamma(r[1]:r[2], mean(t), 1)
  df = table(round(t))
  df = df/sum(df)
  # which distribution can approximate the frequency of reactome terms
  hist(t, prob=T, sub='Distribution of mutation rate', breaks=s,
       xlab='Lambda', ylab='', ylim=c(0, max(dg, df)), ...)
  # try negative binomial and poisson distributions
  # parameterized on the means
  lines(r[1]:r[2], dg, col='black', type='b')
  #lines(r[1]:r[2], dp, col='red', type='b')
  #lines(r[1]:r[2], dg, col='blue', type='b')
  points(round(qgamma(0.95, mean(t), 1), 0), 0, pch=20, col='red')
  legend('topright', legend =c('Gamma'), fill = c('black'))
  ## plot the lamda rates greater than cutoff
  ## NOTE: a possible error may need correcting, if there are no values over 0.95
  c = qgamma(0.95, mean(t), 1)
  i = which(mDat[,'lambda.other'] > c)
  # get limit of vector from the row names of mDat
  l = as.numeric(rownames(mDat)[nrow(mDat)])
  x = rep(0, times = l)
  rn = as.numeric(rownames(mDat)[i])
  x[rn] = mDat[i,'lambda.other']
  plot(x, pch=20, cex=0.5, sub='Significant Mutation Rate ~ Gamma(lambda)', ylab='Lambda', xlab='Sequence', 
       ylim=c(min(x[rn]-1), max(x[rn])),  ...)
}

getSignificantPositions = function(mDat){
  mDat = na.omit(mDat)
  mDat = mDat[mDat[,'var.q'] != 3, ]
  # remove first quantile of theta
  mDat = mDat[mDat[,'theta.q'] != 1, ]
  # get index of significant positions
  c = qgamma(0.95, mean(mDat[,'lambda.other']), 1)
  i = which(mDat[,'lambda.other'] > c)
  return(mDat[i,])
}






head(m15)
plot(m15[,'theta'])
plot(logit(m15[,'theta']))
plot(m15[,'var'])
plot(log(m15[,'var']))
plot(m15[,'lambda.base'])
plot(m15[,'lambda.other'])
plot(logit(m15[,'theta']))
plot(logit(m15[,'theta']))
plot(density(m15[,'lambda.base']))
plot(density(m15[,'lambda.other']))
plot(m15[,'lambda.other'], m15[,'theta'], pch=20, cex=0.5)
plot(m15[,'lambda.other'], logit(m15[,'theta']), pch=20, cex=0.5)
plot(m15[,'lambda.other'], m15[,'theta'], pch=20, cex=0.5)
plot(m15[,'lambda.base'], logit(m15[,'theta']), pch=20, cex=0.5)
plot(m15[,'lambda.base'], (m15[,'theta']), pch=20, cex=0.5)



m15 = f_getMutations(file.choose(), refseq)
m17 = f_getMutations(file.choose(), refseq)

m = cbind(m15, m17)
head(m)
m = na.omit(m)
par(mfrow=c(1,1))
colnames(m) = paste(colnames(m), c(15, 15, 17, 17), sep='')
head(m)
coplot(logit(m[,'theta15']) ~ logit(m[,'theta17']) | log(m[,'var15']) * log(m[,'var17']))

dfBoth = as.data.frame(m)

quantile(dfBoth$var15, 0:10/10)
var15.g = cut(dfBoth$var15, breaks = quantile(dfBoth$var15, 0:10/10), labels = 0:9, include.lowest = T)

quantile(dfBoth$var17, 0:10/10)
var17.g = cut(dfBoth$var17, breaks = quantile(dfBoth$var17, 0:10/10), labels = 0:9, include.lowest = T)

dfBoth$var15.g = var15.g
dfBoth$var17.g = var17.g
xtabs(~ var15.g + var17.g, data=dfBoth)
coplot(logit(dfBoth$theta15) ~ logit(dfBoth$theta17) | dfBoth$var15.g * dfBoth$var17.g)

# remove the first quantiles from theta


m15 = na.omit(mFile15)


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

plot(1:nrow(file15), (var15), pch=20, cex=0.5, col='grey')
points(1:nrow(file17), (var17), pch=20, cex=0.5, col='blue')

coplot(log(theta15[names(theta15)]) ~ log(theta17[names(theta15)]) | log(var15[names(theta15)]) * log(var17[names(theta15)]))


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
seq = c('A' = 1200, 'T' = 0, 'G' = 0, 'C' = 0)
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

# # get posterior predictive values using multinomial random sample
# getPosteriorPredict = function(theta, n=1){
#   return(t(rmultinom(1000, n, theta)))
# }

# get posterior theta from posterior dirichlet
getPosterior = function(alpha, n=1000){
  if(!require(LearnBayes)) stop('Package LearnBayes required')
  p = rdirichlet(n, alpha)
  colnames(p) = names(alpha)
  #m = colMeans(p)
  return(p)
}

getPosteriorPredict = function(theta, n=1){
  ret = t(apply(theta, 1, function(x) rmultinom(1, n, x)))
  colnames(ret) = colnames(theta)
  return(ret)
}

getPosteriorGamma = function(alpha, base, n=1000, rate=100){
  alpha.scale = (alpha/sum(alpha)) * rate
  i = which(names(alpha.scale) == base)
  alpha.new = c(alpha.scale[i], sum(alpha.scale[-i]))
  names(alpha.new) = c(base, 'Other')
  rg = sapply(seq_along(alpha.new), function(x) {
    return(rgamma(n, alpha.new[x], 1))
  })
  colnames(rg) = names(alpha.new)
  return(rg)
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

getBase.rate = function(alpha, base, scale=100){
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



