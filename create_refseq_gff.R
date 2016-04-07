# Name: create_refseq_gff.R
# Auth: u.niazi@imperial.ac.uk
# Date: 07/04/2016
# Desc: create a gff file from the refseq


##### source header files
source('header.R')
library(GenomicAlignments)
library(rtracklayer)
library(annotate)

refseq = readDNAStringSet('Data_external/Reference_seq/Eng 195 concatenated by Daniel.fasta')

# find positions for NNN in the sequence
vMatch = matchPattern('NNN', refseq[[1]])
# convert views object to ranges
oIRmatch = ranges(vMatch)
# create a ranges object for the full sequence length
oIRfull = IRanges(1, width=width(refseq))
# get positions for subsequences
oIRseq = setdiff(oIRfull, oIRmatch)
# extract the sequences for these from refseq
oDNAseqs = DNAStringSet(refseq[[1]], start(oIRseq), end(oIRseq))

# lBlast = vector('list', length=length(oDNAseqs))
# for (i in 1:length(lBlast)){
#   lBlast[[i]] = blastSequences(oDNAseqs[[i]], hitListSize = 2, as = 'data.frame')
# }

writeXStringSet(oDNAseqs, 'Temp/seqs.fasta')
# after naming the sequences manually import the sequences again
oDNAseqs.import = readDNAStringSet('Temp/seqs.fasta')

oGRseqs = GRanges(seqnames='Eng_195_concatenated_by_Daniel A new nucleotide sequence entered manually', oIRseq)
oGRseqs$ID = names(oDNAseqs.import)
oGRseqs$seqs = as.character(oDNAseqs.import)
# export the gff
export(oGRseqs, 'Objects/virus.gff', format = 'gff3')
lAnnotation = list(seqs = oDNAseqs.import, ranges=oGRseqs)
save(lAnnotation, file='Objects/lAnnotation.rds')





