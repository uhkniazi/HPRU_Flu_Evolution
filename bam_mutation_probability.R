# Name: bam_mutation_probability.R
# Auth: u.niazi@imperial.ac.uk
# Date: 23/02/2016
# Desc: gets a list of bam files and calculates a variance  and mutation probability for each file


##### source header files
source('header.R')
library(GenomicAlignments)

refseq = readDNAStringSet('Data_external/Reference_seq/Eng 195 concatenated by Daniel.fasta')

########## functions used in the script

plot.diagnostics = function(mDat, ...){
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
  mDat = na.omit(mDat)
  mDat = mDat[mDat[,'var.q'] != 3, ]
  # remove first quantile of theta
  mDat = mDat[mDat[,'theta.q'] != 1, ]
  # get index of significant positions
#   c = qgamma(0.95, mean(mDat[,'lambda.other']), 1)
#   i = which(mDat[,'lambda.other'] > c)
  p.val = pgamma(mDat[,'lambda.other'], shape = mean(mDat[,'lambda.other']), 1, lower.tail = F)
  p.adj = p.adjust(p.val, 'BH')
  mDat = cbind(mDat, p.val, p.adj)
  i = which(mDat[,'p.val'] < p.cut)
  l = as.numeric(rownames(mDat)[nrow(mDat)])
  x = rep(0, times = l)
  rn = as.numeric(rownames(mDat)[i])
  x[rn] = mDat[i,'lambda.other']
  # check if no significant genes
  if (!add && length(i) > 0) {
    plot(x, pch=20, cex=0.5, sub='Significant Mutation Rate ~ Gamma(lambda)', ylab='Lambda', xlab='Sequence', 
         ylim=c(min(x[rn]-1), max(x[rn])),  ...)
  } else if (add && length(i) > 0) {
    points(x, pch=20, cex=0.5, ...)
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
                                           min_mapq = 30))
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
                  labels = c('q5.low', 'q90', 'q5.high'))
  theta.q = cut(param[,'theta'],breaks = quantile(param[,'theta'], c(0, 0.05, 0.95, 1)), include.lowest = T, 
                    labels = c('q5.low', 'q90', 'q5.high'))
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

lapply(names(lMutation), function(x) plot.diagnostics(lMutation[[x]], main=x))
lSignificant = lapply(names(lMutation), function(x) getSignificantPositions(lMutation[[x]], main=x, p.cut = 0.05))
names(lSignificant) = csvSamples
dfMutants = data.frame(Signif.Positions= do.call(rbind, endoapply(lSignificant, dim))[,-2])

## all the samples in one matrix together
mAllMutants.sig = matrix(0, nrow=width(refseq), ncol=length(csvSamples), dimnames=list(1:width(refseq), csvSamples))

for(i in 1:ncol(mAllMutants.sig)){
  m = match(rownames(lSignificant[[i]]), rownames(mAllMutants.sig))
  mAllMutants.sig[m,i] = lSignificant[[i]][,'lambda.other']
}

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
