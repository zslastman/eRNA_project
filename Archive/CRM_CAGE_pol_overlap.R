setwd('/g/furlong/Harnett/TSS_CAGE_myfolder/')
#we can just pull out two strand specific vectors using this easily enough
bam_file<-'data/solexa/bam/chr2R.all.cage.68h.trimmed.quality.sort.05_15_2013.sort.bam'

#read in the reads with scanbam
reads <- scanBam(
  file = bam_file,
  param = ScanBamParam(
    what = c('pos', 'rname', 'strand', 'qwidth'),
    flag = scanBamFlag(isUnmappedQuery = FALSE)
  )
)[[1]]

#convert them into a GRanges object
reads.gr<-GRanges(reads[['rname']],IRanges(reads[['pos']],width=reads[['qwidth']]),strand=reads[['strand']])
cov1<-coverage(reads.gr)








#export them as two

#CRM_CAGE_Pol