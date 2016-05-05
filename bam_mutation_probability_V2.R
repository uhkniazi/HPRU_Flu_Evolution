# Name: bam_mutation_probability_V2.R
# Auth: u.niazi@imperial.ac.uk
# Date: 04/05/2016
# Desc: gets a list of bam files and calculates a variance  and mutation probability for each file


##### source header files
source('header.R')
library(GenomicAlignments)

refseq = readDNAStringSet('Data_external/Reference_seq/Eng 195 concatenated by Daniel.fasta')
# load the sequence annotation object
load('Objects/lAnnotation.rds')


########## functions used in the script

plot.diagnostics = function(mDat, ...){
  if (!require(scatterplot3d)) stop('scatterplot3d library required')
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
  plot(mDat[,'lambda.base'], mDat[,'lambda.other'], pch=20, cex=0.5,
       sub='Ref rate vs Mutation rate', ylab='Lambda Mutation', xlab='Lambda Base', ...)
  # plot the variance, theta and mutation rate
  scatterplot3d(x=mDat[,'lambda.other'], z=logit(mDat[,'theta']), y=mDat[,'var'], pch=16, cex.symbols=0.5, scale.y=1, highlight.3d=T,
                angle=60, sub='Mutation Rate vs Theta and Variance', zlab='Logit Theta', xlab='Lambda', ylab='Variance', ... )
  ## plot the density and fit distribution, for mutation rate
  t = mDat[,'lambda.other']
  r = range(t)
  s = seq(max(0.5, floor(r[1]))-0.5, ceiling(r[2])+0.5, by=1)
  # calculate the mid points for histogram/discrete distribution
  h = hist(t, breaks=s, plot=F)
  dg = dgamma(h$mids, mean(t), 1)
  # which distribution can approximate the frequency
  hist(t, prob=T, sub='Distribution of mutation rate', breaks=s,
       xlab='Lambda', ylab='', ylim=c(0, max(dg, h$density)), ...)
  # parameterized on the means
  lines(h$mids, dg, col='black', type='b')
  points(qgamma(0.95, mean(t), 1), 0, pch=20, col='red')
  legend('topright', legend =c('Gamma'), fill = c('black'))
  ## plot the lamda rates greater than cutoff
  ## NOTE: a possible error may need correcting, if there are no values over 0.95
  c = qgamma(0.95, mean(t), 1)
  i = which(mDat[,'lambda.other'] > c)
  # get limit of vector from the row names of mDat, only if any significant positions
  if (length(i) > 0) {
    l = as.numeric(rownames(mDat)[nrow(mDat)])
    x = rep(0, times = l)
    rn = as.numeric(rownames(mDat)[i])
    x[rn] = mDat[i,'lambda.other']
    plot(x, pch=20, cex=0.5, sub='Significant Mutation Rate ~ Gamma(lambda)', ylab='Lambda', xlab='Sequence', 
         ylim=c(min(x[rn]-1), max(x[rn])),  ...)
  } # end if
  ## plot all the rates at different cutoffs of the theta
  mDat = mDat[,c('theta', 'lambda.other')]
  cut.pts = cut(mDat[,1], breaks = quantile(mDat[,1], 0:10/10), include.lowest = T, labels = c(1:10))
  # at each cutoff model the gamma variable
  # define function
  plot.lambda.cutoffs = function(lambda, beta=1, cut.pts){
    cut.pts = as.numeric(cut.pts)
    pos = max(cut.pts)
    for (i in pos:1){
      # get the gamma variable
      ind = which(cut.pts >= i)
      t = lambda[ind]
      r = range(t)
      s = seq(max(0.5, floor(r[1]))-0.5, ceiling(r[2])+0.5, by=1)
      # calculate the mid points for histogram/discrete distribution
      h = hist(t, breaks=s, plot=F)
      dg = dgamma(h$mids, mean(t), 1)
      # which distribution can approximate the frequency
      hist(t, prob=T, sub=paste('Distribution of mutation rate at', i), breaks=s,
           xlab='Lambda', ylab='', ylim=c(0, max(dg, h$density)), main='Gamma Mutation Rate')
      # parameterized on the means
      lines(h$mids, dg, col='black', type='b')
      points(qgamma(0.95, mean(t), 1), 0, pch=20, col='red')
      legend('topright', legend =c('Gamma'), fill = c('black'))
    }
  }
  plot.lambda.cutoffs(mDat[,'lambda.other'], cut.pts = cut.pts)
}

getSignificantPositions = function(mDat, p.cut=0.01, add=F, ...){
  # remove positions with no reads i.e. NA
  mDat = na.omit(mDat)
  # remove outlier with very high variance
  mDat = mDat[mDat[,'var.q'] != 3, ]
  # remove outlier, the first part of theta
  mDat = mDat[mDat[,'theta.q'] != 1, ]
  # split the sequence length into groups
  ivSequence = 1:nrow(mDat)
  names(ivSequence) = rownames(mDat)
  ## create breaks in the sequence i.e. bins
  bins = cut(ivSequence, breaks = 30, include.lowest = T)
  # split into regions
  lBins = split(ivSequence, bins)
  # calculate the p.values for each bin
  p.val = lapply(lBins, function(x){
    pgamma(mDat[x,'lambda.other'], shape = mean(mDat[x,'lambda.other']), 1, lower.tail = F)
  })
  # calculate adjusted p.value for each bin
  p.adj = lapply(p.val, p.adjust, 'BH')
  p.val = unlist(p.val)
  p.adj = unlist(p.adj)
  mDat = cbind(mDat, p.val, p.adj)
  i = which(mDat[,'p.val'] < p.cut)
  l = as.numeric(rownames(mDat)[nrow(mDat)])
  x = rep(0, times = l)
  # plot the frequencey and variance as well
  fr = rep(0, times = l)
  sig2 = rep(0, times = l)
  # use row names to index the positions as they map to sequence position
  rn = as.numeric(rownames(mDat)[i])
  # i is the index of the dataframe mapping to corresponding position rn in sequence
  x[rn] = mDat[i,'lambda.other']
  fr[as.numeric(rownames(mDat))] = round((1-mDat[,'theta'])*100, 2)
  sig2[as.numeric(rownames(mDat))] = mDat[,'var']
  # check if no significant genes
  if (!add && length(i) > 0) {
    plot(x, pch=20, cex=0.5, sub='Significant Mutation Rate ~ Gamma(lambda)', ylab='Lambda', xlab='Sequence', 
         ylim=c(min(x[rn]-1), max(x[rn])),  ...)
    # plot the theta - frequency
    ## choose colors
    col.fr = rep(1, length=length(fr))
    # set red colour for the significant positions
    col.fr[rn] = 2
    plot(fr, pch=20, cex=0.5, sub='Mutation Frequency', ylab='Frequency', xlab='Sequence', col=col.fr, ...)#, 
         #ylim=c(0, max(c(fr[rn], 1))),  ...)
    plot(sig2, pch=20, cex=0.5, sub='Variance', ylab='Variance', xlab='Sequence', col=col.fr, ...)#, 
  } else if (add && length(i) > 0) {
    points(x, pch=20, cex=0.5, ...)
    points(fr, pch=20, cex=0.5, ...)
  }
  return(mDat[i,])
}


f_getMutations = function(csBamfile, oDSRef){
  ## loading bam files
  #   oGAbam = readGAlignments(csBamfile)
  #   # reduce to get sequence length
  #   strand(oGAbam) = '*'
  #   oGRbam = reduce(as(oGAbam, 'GRanges'))
  #   rm(oGAbam); gc()
  # optionally bin the genome into smaller chunks if genome is too long
  # process chunks at a time
  # oGRbam.bin = f_split_GRanges_by_length(oGRbam, 1000)
  oPile = pileup(csBamfile, #scanBamParam = ScanBamParam(which = oGRbam), 
                 pileupParam = PileupParam(distinguish_strands = F, max_depth = 1000, min_base_quality = 30, 
                                           min_mapq = 0))
  # get the sequence
  seq = f_getSeq(oPile)
  # adjust the reference sequence size to the aligned pileup data
  #s = as.numeric(rownames(seq)[1])
  #oDSRef.sub = DNAStringSet(oDSRef[[1]], s)
  # get sequence parameters
  param = sapply(rownames(seq), function(x){
    getSequenceParameters(seq[x,], as.character(oDSRef[[1]][as.numeric(x)]))#, iNormalizingRate = iScale)
  })
  colnames(param) = rownames(seq)
  rownames(param) = c('theta', 'var', 'lambda.base', 'lambda.other', 'A', 'T', 'G', 'C')
  param = t(param)
  # create factors for quantiles
  param = na.omit(param)
  var.q = cut(param[,'var'],breaks = quantile(param[,'var'], c(0, 0.05, 0.95, 1)), include.lowest = T, 
                  labels = c('q.low', 'q,m', 'q.high'))
  theta.q = cut(param[,'theta'],breaks = quantile(param[,'theta'], c(0, 0.01, 0.95, 1)), include.lowest = T, 
                    labels = c('q.low', 'q.m', 'q.high'))
  param = cbind(param, var.q, theta.q)
  #return(t(param))
  ## reformat the data matrix 
  mFile = matrix(NA, nrow=width(oDSRef), ncol=ncol(param), dimnames=list(1:width(oDSRef), colnames(param)))
  m = match(rownames(param), rownames(mFile))
  mFile[m,] = param
  return(mFile)
}

csBamfiles = list.files('Data_external/Sam/', '*.bam$', full.names = T)

## calculate data for each bam file
lMutation = vector('list', length=length(csBamfiles))
lMutation = lapply(csBamfiles, f_getMutations, refseq)
csvSamples = gsub('Data_external/Sam//(\\w+)\\.bam', '\\1', csBamfiles)
names(lMutation) = csvSamples

temp = lapply(names(lMutation), function(x) plot.diagnostics(lMutation[[x]], main=x))

## get the positions in the sequence with 
## significant rates
lSignificant = lapply(names(lMutation), function(x) getSignificantPositions(lMutation[[x]], main=x, p.cut = 0.05))
names(lSignificant) = csvSamples
sapply(lSignificant, dim)
dfMutants = data.frame(Signif.Positions= do.call(rbind, endoapply(lSignificant, dim))[,-2])

# get only those samples with significant positions
dfMutants.sig = data.frame(Signif.Positions= dfMutants[dfMutants$Signif.Positions != 0,])
i = which(dfMutants$Signif.Positions != 0)
rownames(dfMutants.sig) = rownames(dfMutants)[i] 

## all the significant samples in one matrix together
mAllMutants.sig = matrix(0, nrow=width(refseq), ncol=nrow(dfMutants.sig), 
                         dimnames=list(1:width(refseq), rownames(dfMutants.sig)))
nm = colnames(mAllMutants.sig)
for(i in 1:ncol(mAllMutants.sig)){
  m = match(rownames(lSignificant[[nm[i]]]), rownames(mAllMutants.sig))
  mAllMutants.sig[m,i] = lSignificant[[nm[i]]][,'lambda.other']
}

# assign protein ids to positions
cvProteins = lAnnotation$ranges$ID
fProteins = rep(NA, length=nrow(mAllMutants.sig))
s = start(ranges(lAnnotation$ranges))
e = end(ranges(lAnnotation$ranges))
for (i in seq_along(s)){
  fProteins[s[i]:e[i]] = cvProteins[i]
}

dfAllMutants.sig = data.frame(mAllMutants.sig, protein=fProteins)
# remove positions with NNN
dfAllMutants.sig = na.omit(dfAllMutants.sig)
fProteins = factor(dfAllMutants.sig$protein, levels = cvProteins)
dfAllMutants.sig$protein = fProteins

## calculate theta for each protein in sample
f_getTheta = function(x, prot){
  # convert the positions with mutation and no mutation to a 
  # binomial variable
  x = ifelse(x > 0, 1, 0)
  ret = tapply(x, prot, function(y){
    # calculate binomial parameters
    n = length(y)
    return(round(sum(y)/n, 3))
  })
  return(unlist(ret))
}

## use the parent i.e. g14 to calculate a prior distribution
# lPriors = tapply(dfAllMutants.sig$G14_q10_sort, dfAllMutants.sig$protein, function(x) {
#   # convert the positions with mutation and no mutation to a 
#   # binomial variable
#   x = ifelse(x > 0, 1, 0)
#   # calculate binomial parameters
#   n = length(x)
#   s = sum(x)
#   f = n-s
#   # get the alpha and beta parameters for the beta prior
#   return(c(alpha=s, beta=f))
# })

## use g14, 13, 17 and 19 to calculate prior
g13 = f_getTheta(dfAllMutants.sig$G13_q10_sort, dfAllMutants.sig$protein)
g14 = f_getTheta(dfAllMutants.sig$G14_q10_sort, dfAllMutants.sig$protein)
g17 = f_getTheta(dfAllMutants.sig$G17_q10_sort, dfAllMutants.sig$protein)
g19 = f_getTheta(dfAllMutants.sig$G19_q10_sort, dfAllMutants.sig$protein)

mPrior = cbind(g13, g14, g17, g19)

mPriors = apply(mPrior, 1, function(x) unlist(getalphabeta(mean(x), var(x))))

#### make plots for all the experiments 
dfExp = data.frame(S100=dfAllMutants.sig$G11_q10_sort, S10=dfAllMutants.sig$G12_q10_sort, S0=dfAllMutants.sig$G13_q10_sort,
                    protein=dfAllMutants.sig$protein)
dfExp1 = dfAllMutants.sig
rownames(dfExp1) = rownames(dfAllMutants.sig)

mExp1 = NULL
for(i in 1:9){
  # read vector for protein i
  m = dfExp1[,i]
  names(m) = rownames(dfExp1)
  p = dfExp1[,'protein']
  # convert to binomial 1 0 variable
  m = ifelse(m > 0, 1, 0)
  # calculate the posterior theta for each protein
  ivPost = rep(NA, length=nlevels(p))
  names(ivPost) = levels(p)
  lev = levels(p)
  for (inner in 1:nlevels(p)){
    # get subvector for the current protein 
    m.sub = m[p == lev[inner]]
    # calculate posterior 
    trials = length(m.sub)
    suc = sum(m.sub)
    fail = trials - suc
    # get prior alpha beta
    #pr = lPriors[[lev[inner]]]
    pr = mPriors[,lev[inner]]
    po = rbeta(1000, suc+pr['alpha'], fail+pr['beta'])
    #po = rbeta(1000, suc+1, fail+1)
    ivPost[lev[inner]] = mean(po)
  }
  mExp1 = cbind(mExp1, ivPost)
  #   f = which(m == 0)
  #   m = m[-f]
  #   p = p[-f]
  #   # calculate lengths, i.e number of mutations
  #   mut.no = tapply(m, p, length)
  #   # normalize by dividing by length of the protein
  #   lp = as.numeric(table(dfExp1[,'protein']))
  #   mut.nor = mut.no/lp
  #   mExp1 = cbind(mExp1, mut.nor)
}
#colnames(mExp1) = colnames(dfExp1)[1:9]
cn = c('S100', 'S10', 'S0', 'W', 'P15', 'P16', 'P0', 'P18', 'P0')
colnames(mExp1) = cn
c = rainbow(nrow(mExp1))
barplot(mExp1, col=c, beside=T, las=2)
legend('topright', legend = rownames(mExp1), fill=c)

mExp1 = t(mExp1)
c = rainbow(nrow(mExp1))
barplot(mExp1, col=c, beside=T, las=2)
legend('topleft', legend = rownames(mExp1), fill=c)

# create a factor for samples to order the data
fSamples = c('100', '10', '0', '0', '5', '5', '0', '5', '0')
rownames(mExp1) = fSamples
fSamples = factor(fSamples, levels = c('0', '5', '10', '100'))

mExp1 = t(mExp1)

# plot the single passage experiment
mExp1.single = mExp1[,order(fSamples)]
mExp1.single = mExp1.single[,colnames(mExp1.single) %in% c('0', '10', '100')]
mExp1.single = t(mExp1.single)

# convert to a rate
mExp1.single = round(mExp1.single, 3)
c = rainbow(ncol(mExp1.single))
matplot(mExp1.single, type='b', pch=20, lty=1, xaxt='n', xlab='Dose Micromole', ylab='Theta', 
        main='Relationship between Dose and Theta for Each Gene', col=c, lwd=2)
axis(1, at = 1:nrow(mExp1.single), labels = rownames(mExp1.single))
legend('topleft', legend = colnames(mExp1.single), fill=c)

# plot the multiple passage
mExp1.mul = mExp1[,order(fSamples)]
mExp1.mul = mExp1.mul[,colnames(mExp1.mul) %in% c('0', '5')]
mExp1.mul = t(mExp1.mul)

# convert to a rate
mExp1.mul = round(mExp1.mul, 3)
c = rainbow(ncol(mExp1.mul))
matplot(mExp1.mul, type='b', pch=20, lty=1, xaxt='n', xlab='Dose Micromole', ylab='Theta', 
        main='Relation between Multiple passages and Theta for Each Gene', col=c, lwd=2)
axis(1, at = 1:nrow(mExp1.mul), labels = rownames(mExp1.mul))
legend('topleft', legend = colnames(mExp1.mul), fill=c)

### plot together
mExp1.all = mExp1[,order(fSamples)]
mExp1.all = t(mExp1.all)
# remove one of the odd samples
mExp1.all = mExp1.all[-7,]

# convert to a rate
mExp1.all = round(mExp1.all, 3)
c = rainbow(ncol(mExp1.all))
matplot(mExp1.all, type='b', pch=20, lty=1, xaxt='n', xlab='Dose Micromole', ylab='Theta', 
        main='Relationship between Dose and Theta for Each Gene', col=c, lwd=2)
axis(1, at = 1:nrow(mExp1.all), labels = rownames(mExp1.all))
legend('topleft', legend = colnames(mExp1.all), fill=c)



## mean for the treated samples
temp = mExp1[,colnames(mExp1) %in% c('100', '10', '5')]
colMeans(temp)
mean(colMeans(temp))


## bar plot to make figure from paper
np = mExp1.single[,'NP']
np = np[-6]

np0 = mean(np[1:4])
np10 = np[5]
mBar = cbind(np0, np10)
colnames(mBar) = c('0', '10')
barplot(mBar, ylim=c(0, 0.07), main='Theta for NP gene at 0 and 10 MMole', sub='Doses', ylab='Theta')


######### calculate transition matrices
lt = lapply(lSignificant, function(x) {
  m = getTransitionMatrix(t(x[,c('A', 'T', 'G', 'C')]))
  m = round(m/rowSums(m), 2)
})

# get G to A
g2a = lapply(lSignificant, function(x) {
  m = getTransitionMatrix(t(x[,c('A', 'T', 'G', 'C')]))
  m = round(m/rowSums(m), 2)
  m['G', 'A']
})
g2a = unlist(g2a)
names(g2a) = fSamples
#g2a = g2a[order(fSamples)]
g2a = tapply(g2a, fSamples, mean)
barplot(g2a, main='G to A transition probability', xlab='Doses')

c2t = lapply(lSignificant, function(x) {
  m = getTransitionMatrix(t(x[,c('A', 'T', 'G', 'C')]))
  m = round(m/rowSums(m), 2)
  m['C', 'T']
})
c2t = unlist(c2t)
names(c2t) = fSamples
#c2t = c2t[order(fSamples)]
c2t = tapply(c2t, fSamples, mean)
barplot(c2t, main='C to T transition probability', xlab='Doses')



############################################### older stuff needs cleanup

## experiment 1
dfExp1 = data.frame(S100=dfAllMutants.sig$G11_q10_sort, S10=dfAllMutants.sig$G12_q10_sort, S0=dfAllMutants.sig$G13_q10_sort,
                    protein=dfAllMutants.sig$protein)
rownames(dfExp1) = rownames(dfAllMutants.sig)

mExp1 = NULL
for(i in 1:3){
  # read vector for protein i
  m = dfExp1[,i]
  names(m) = rownames(dfExp1)
  p = dfExp1[,'protein']
  # convert to binomial 1 0 variable
  m = ifelse(m > 0, 1, 0)
  # calculate the posterior theta for each protein
  ivPost = rep(NA, length=nlevels(p))
  names(ivPost) = levels(p)
  lev = levels(p)
  for (inner in 1:nlevels(p)){
    # get subvector for the current protein 
    m.sub = m[p == lev[inner]]
    # calculate posterior 
    trials = length(m.sub)
    suc = sum(m.sub)
    fail = trials - suc
    # get prior alpha beta
    #pr = lPriors[[lev[inner]]]
    pr = mPriors[,lev[inner]]
    po = rbeta(1000, suc+pr['alpha'], fail+pr['beta'])
    #po = rbeta(1000, suc+1, fail+1)
    ivPost[lev[inner]] = mean(po)
  }
  mExp1 = cbind(mExp1, ivPost)
#   f = which(m == 0)
#   m = m[-f]
#   p = p[-f]
#   # calculate lengths, i.e number of mutations
#   mut.no = tapply(m, p, length)
#   # normalize by dividing by length of the protein
#   lp = as.numeric(table(dfExp1[,'protein']))
#   mut.nor = mut.no/lp
#   mExp1 = cbind(mExp1, mut.nor)
}
colnames(mExp1) = colnames(dfExp1)[1:3]
c = rainbow(nrow(mExp1))
barplot(mExp1, col=c, beside=T)
legend('topright', legend = rownames(mExp1), fill=c)

## experiment 2
dfExp2 = data.frame(G15=dfAllMutants.sig$G15_q10_sort, G16=dfAllMutants.sig$G16_q10_sort, G18=dfAllMutants.sig$G18_q10_sort,
                    protein=dfAllMutants.sig$protein)
rownames(dfExp2) = rownames(dfAllMutants.sig)

mExp2 = NULL
for(i in 1:3){
  m = dfExp2[,i]
  names(m) = rownames(dfExp2)
  p = dfExp2[,'protein']
  # remove 0s
  f = which(m == 0)
  m = m[-f]
  p = p[-f]
  # calculate lengths, i.e number of mutations
  mut.no = tapply(m, p, length)
  # normalize by dividing by length of the protein
  lp = as.numeric(table(dfExp2[,'protein']))
  mut.nor = mut.no/lp
  mExp2 = cbind(mExp2, mut.nor)
}
colnames(mExp2) = colnames(dfExp2)[1:3]
c = rainbow(nrow(mExp2))
barplot(mExp2, col=c, beside=T)
legend('topright', legend = rownames(mExp2), fill=c)

## plot both experiments together
mExp.join = cbind(mExp1, mExp2)
c = rainbow(nrow(mExp.join))
barplot(mExp.join, col=c, beside=T)
legend('topleft', legend = rownames(mExp.join), fill=c)


mAllMutants = matrix(NA, nrow=width(refseq), ncol=length(csvSamples), dimnames=list(1:width(refseq), csvSamples))

for(i in 1:ncol(mAllMutants)){
  # remove noisy areas
  mDat = na.omit(lMutation[[i]])
  mDat = mDat[mDat[,'var.q'] != 3, ]
  # remove first quantile of theta
  mDat = mDat[mDat[,'theta.q'] != 1, ]
  m = match(rownames(mDat), rownames(mAllMutants))
  mAllMutants[m,i] = mDat[,'lambda.other']
}

## various plots 
c = rainbow(ncol(mAllMutants.sig))
matplot(mAllMutants.sig, type='p', pch=20, cex=0.5, main='position of significant residues', col=c, ylim=c(1, max(mAllMutants.sig)))
plot.new()
legend('center', legend = colnames(mAllMutants.sig), fill=c)
colnames(mAllMutants) = gsub('(G\\d*).+', '\\1', colnames(mAllMutants))
boxplot(mAllMutants, main='All Mutation Rates', las=2, cex=0.5, ylim=c(0, 15))

# density of gamma rates
sapply(colnames(mAllMutants), function(n) {
  t = na.omit(mAllMutants[,n])
  r = range(t)
  s = seq(max(0.5, floor(r[1]))-0.5, ceiling(r[2])+0.5, by=1)
  # calculate the mid points for histogram/discrete distribution
  h = hist(t, breaks=s, plot=F)
  dg = dgamma(h$mids, mean(t), 1)
  # which distribution can approximate the frequency
  hist(t, prob=T, sub='Distribution of mutation rate', breaks=s,
       xlab='Lambda', ylab='', ylim=c(0, max(dg, h$density)), main=paste(n, 'Gamma Mutation Rate'))
  # parameterized on the means
  lines(h$mids, dg, col='black', type='b')
  points(qgamma(0.95, mean(t), 1), 0, pch=20, col='red')
  legend('topright', legend =c('Gamma'), fill = c('black'))
  })

fSamples = factor(c('Drug', 'Drug', 'Cont', 'Wild', 'Drug', 'Drug', 'Cont', 'Drug', 'Cont'), levels = c('Wild', 'Cont', 'Drug'))

mAllMutants.pool = na.omit(mAllMutants)

colnames(mAllMutants.pool) = fSamples

ivDrug = rowMeans(mAllMutants.pool[,fSamples == 'Drug'])
ivWild = (mAllMutants.pool[,'Wild'])
ivCont = rowMeans(mAllMutants.pool[,fSamples == 'Cont'])
cor(cbind(ivWild, ivDrug, ivCont))

mAllMutants.pool = cbind(ivDrug, ivWild, ivCont)

# density of gamma rates
sapply(colnames(mAllMutants.pool), function(n) {
  t = na.omit(mAllMutants.pool[,n])
  r = range(t)
  s = seq(max(0.5, floor(r[1]))-0.5, ceiling(r[2])+0.5, by=1)
  # calculate the mid points for histogram/discrete distribution
  h = hist(t, breaks=s, plot=F)
  dg = dgamma(h$mids, mean(t), 1)
  # which distribution can approximate the frequency
  hist(t, prob=T, sub='Distribution of mutation rate', breaks=s,
       xlab='Lambda', ylab='', ylim=c(0, max(dg, h$density)), main=paste(n, 'Gamma Mutation Rate'))
  # parameterized on the means
  lines(h$mids, dg, col='black', type='b')
  points(qgamma(0.95, mean(t), 1), 0, pch=20, col='red')
  legend('topright', legend =c('Gamma'), fill = c('black'))
})


fSamples = factor(c('Exp1.D', 'Exp1.D', 'Exp1.C', 
                    'Wild', 'Exp2.D', 'Exp2.D', 'Exp2.C', 'Exp2.D', 'Exp2.C'))

mAllMutants.pool = na.omit(mAllMutants)

colnames(mAllMutants.pool) = fSamples

boxplot(mAllMutants.pool)

ivExp1.D = rowMeans(mAllMutants.pool[,fSamples == 'Exp1.D'])
ivExp1.C = mAllMutants.pool[,'Exp1.C']
ivWild = (mAllMutants.pool[,'Wild'])
ivExp2.D = rowMeans(mAllMutants.pool[,fSamples == 'Exp2.D'])
ivExp2.C = rowMeans(mAllMutants.pool[,fSamples == 'Exp2.C'])

mAllMutants.pool = cbind(ivExp1.D, ivExp1.C, ivWild, ivExp2.D, ivExp2.C)
boxplot(mAllMutants.pool)
round(cor(mAllMutants.pool), 2)
# density of gamma rates
sapply(colnames(mAllMutants.pool), function(n) {
  t = na.omit(mAllMutants.pool[,n])
  r = range(t)
  s = seq(max(0.5, floor(r[1]))-0.5, ceiling(r[2])+0.5, by=1)
  # calculate the mid points for histogram/discrete distribution
  h = hist(t, breaks=s, plot=F)
  dg = dgamma(h$mids, mean(t), 1)
  # which distribution can approximate the frequency
  hist(t, prob=T, sub='Distribution of mutation rate', breaks=s,
       xlab='Lambda', ylab='', ylim=c(0, max(dg, h$density)), main=paste(n, 'Gamma Mutation Rate'))
  # parameterized on the means
  lines(h$mids, dg, col='black', type='b')
  points(qgamma(0.95, mean(t), 1), 0, pch=20, col='red')
  legend('topright', legend =c('Gamma'), fill = c('black'))
})
