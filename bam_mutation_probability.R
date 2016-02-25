# Name: bam_mutation_probability.R
# Auth: u.niazi@imperial.ac.uk
# Date: 23/02/2016
# Desc: gets a list of bam files and calculates a variance  and mutation probability for each file


##### source header files
source('header.R')
library(GenomicAlignments)

refseq = readDNAStringSet('Data_external/Reference_seq/Eng 195 concatenated by Daniel.fasta')


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
                 pileupParam = PileupParam(distinguish_strands = F, max_depth = 1000))
  # get the sequence
  seq = f_getSeq(oPile)
  # adjust the reference sequence size to the aligned pileup data
  #s = as.numeric(rownames(seq)[1])
  #oDSRef.sub = DNAStringSet(oDSRef[[1]], s)
  # get sequence parameters
  param = sapply(rownames(seq), function(x){
    getSequenceParameters(seq[x,], as.character(oDSRef[[1]][as.numeric(x)]))
  })
  colnames(param) = rownames(seq)
  rownames(param) = c('theta', 'var', 'rate.base', 'rate.other')
  param = t(param)
  # create factors for quantiles
  param = na.omit(param)
  var.q = cut(param[,'var'],breaks = quantile(param[,'var'], c(0, 0.05, 0.95, 1)), include.lowest = T, 
                  labels = c('q5.low', 'q90', 'q5.high'))
  theta.q = cut(param[,'theta'],breaks = quantile(param[,'theta'], c(0, 0.05, 0.95, 1)), include.lowest = T, 
                    labels = c('q5.low', 'q90', 'q5.high'))
#   rate.base.q = cut(param[,'rate.base'],breaks = quantile(param[,'rate.base'], c(0, 0.05, 0.95, 1)), include.lowest = T, 
#                 labels = c('q5.low', 'q90', 'q5.high'))
#   rate.other.q = cut(param[,'rate.other'],breaks = quantile(param[,'rate.other'], c(0, 0.05, 0.95, 1)), include.lowest = T, 
#                 labels = c('q5.low', 'q90', 'q5.high'))
#   
  param = cbind(param, var.q, theta.q)
  #return(t(param))
  ## reformat the data matrix 
  mFile = matrix(NA, nrow=width(oDSRef), ncol=ncol(param), dimnames=list(1:width(oDSRef), colnames(param)))
  m = match(rownames(param), rownames(mFile))
  mFile[m,] = param
  return(mFile)
}



