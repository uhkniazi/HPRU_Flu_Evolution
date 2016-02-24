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
  s = as.numeric(rownames(seq)[1])
  oDSRef.sub = DNAStringSet(oDSRef[[1]], s, width = nrow(seq))
  # get sequence parameters
  param = sapply(seq_along(1:nrow(seq)), function(x){
    getSequenceParameters(seq[x,], as.character(oDSRef.sub[[1]][x]))
  })
  colnames(param) = rownames(seq)
  rownames(param) = c('theta', 'var')
  return(param)
}


