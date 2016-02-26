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
  ## plot the density and fit distribution, for mutation rate
  t = mDat[,'lambda.other']
  r = range(t)
  s = seq(floor(r[1])-0.5, ceiling(r[2])+0.5, by=1)
  r[1] = floor(r[1])
  r[2] = ceiling(r[2])
  # which distribution can approximate the frequency of reactome terms
  hist(t, prob=T, sub='Distribution of mutation rate', breaks=s,
       xlab='Lambda', ylab='', ...)
  # try negative binomial and poisson distributions
  # parameterized on the means
  dn = dnbinom(r[1]:r[2], size = mean(t), mu = mean(t))
  dp = dpois(r[1]:r[2], mean(t))
  lines(r[1]:r[2], dn, col='black', type='b')
  lines(r[1]:r[2], dp, col='red', type='b')
  legend('topright', legend =c('nbinom', 'poi'), fill = c('black', 'red'))
}


f_getMutations = function(csBamfile, oDSRef, iScale=1000){
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
                 pileupParam = PileupParam(distinguish_strands = F, max_depth = 1000))
  # get the sequence
  seq = f_getSeq(oPile)
  # adjust the reference sequence size to the aligned pileup data
  #s = as.numeric(rownames(seq)[1])
  #oDSRef.sub = DNAStringSet(oDSRef[[1]], s)
  # get sequence parameters
  param = sapply(rownames(seq), function(x){
    getSequenceParameters(seq[x,], as.character(oDSRef[[1]][as.numeric(x)]), iNormalizingRate = iScale)
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



