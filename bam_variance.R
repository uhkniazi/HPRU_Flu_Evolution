# Name: bam_variance.R
# Auth: u.niazi@imperial.ac.uk
# Date: 27/01/2016
# Desc: gets a list of bam files and calculates a variance for each file


source('header.R')
#source('~/Dropbox/Home/Data/R/My_Libraries/NGS_functions.R')

csBamfiles = list.files('Data_external/Sam/', '*.bam$', full.names = T)

## calculate variance for each bam file
lVariance = vector('list', length=length(csBamfiles))


f_getVariance = function(csBamfile){
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
  # get variance
  v = apply(seq, 1, function(x) getSequenceVariance(x, size = 1000))
  return(v)
}

lVariance = sapply(csBamfiles, f_getVariance)

sapply(lVariance, summary)
dir.create('Results')
pdf('Results/variance.pdf')
par(mfrow=c(3,1))
for (i in seq_along(lVariance)){
  plot(log(lVariance[[i]]), pch=20, main=names(lVariance)[i], xlab='sequence position', xaxt='n', ylab='log variance', 
       col='darkgrey', cex=0.5)
  pos.1 = range(as.numeric(names(lVariance[[i]])))[1]
  pos.2 = range(as.numeric(names(lVariance[[i]])))[2]
  axis(1, at = round(seq(pos.1, pos.2, length.out = 10)))
  abline(h = mean(log(lVariance[[i]])), lwd=2, col=2)
  abline(h = mean(log(lVariance[[i]]))+2*sd(log(lVariance[[i]])), lwd=2, col=2, lty=2)
  #lines(lowess(log(lVariance[[i]])), col=3, lwd=2)
  #abline(h = mean(log(lVariance[[i]]))-2*sd(log(lVariance[[i]])), lwd=2, col=2, lty=2)
}
plot.new()
par(mfrow=c(1,1))
m = sapply(lVariance, sum)
plot(m, xaxt='n', main='Total Variance in each sample', xlab='samples', ylab='sum of variance')
n = paste('G', 11:19)
axis(1, 1:length(m), labels = n)

l = sapply(lVariance, log)
names(l) = n
boxplot(l, main='log variance')
dev.off(dev.cur())


