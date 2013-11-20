# @ author Dermot Harnett, EMBL Heidelberg
# @date 16/5/2013
# @title Export coverage into wig format
########################################


# @Original author Tibor Pakozdi, EMBL Heidelberg
# @date Sat Apr 20 20:30:44 2013
# @title Export coverage into wig format
########################################

#usage: Rscript WigExport bamfile 

# add 'track name="..." type=wiggle_0' to top of wig file
# then gzip to compress

# set working directory based on operating system
setwd(dir='/g/furlong/Harnett')
# package list
library(IRanges)
library(multicore)
library(rtracklayer)
library(Rsamtools)

# dm3 chr limits
#save(sl,file='/g/furlong/Harnett/data_general/sl.melanogaster.object.R')
load('/g/furlong/Harnett/data_general/sl.melanogaster.object.R')

#for human
#load('database/rda_storage/sl_hsapiens_hg19.rda')

#Tibors exportWig function, slightly modified
exportWig <- function(bam_file = '', bin_width = 50L, shift_reads = 0L,outdir='./',col=c(200,0,0)) {
  
  cat(bam_file, '\n')
  
  # construct binning ranges
  bin_ranges <- IRangesList(
    #sapply(names(sl)[!grepl('M|U', names(sl))], function(tmp_chr) {
    sapply(names(sl)[!grepl('gl|hap', names(sl))], function(tmp_chr) {
      IRanges(start = seq(1, sl[[tmp_chr]] - bin_width, bin_width), width = bin_width)
    })
  )
  
  # load reads into R object
  reads <- scanBam(
    file = bam_file,
    param = ScanBamParam(
      what = c('pos', 'rname', 'strand', 'qwidth'),
      flag = scanBamFlag(isUnmappedQuery = FALSE)
    )
  )[[1]]
  
  # read width
  rwidth <- as.integer(reads[['qwidth']][1])
  reads[['qwidth']] <- NULL
  
  # correct rname vector if necessary
  reads[['rname']] <- as.character(reads[['rname']])
  if (!grepl('chr', reads[['rname']][1])) reads[['rname']] <- as.factor(paste('chr', reads[['rname']], sep = ''))
  
  # loop over all reasonable chromosomes
  chr_list <- names(sl)[!grepl('U|M', names(sl))]
  #chr_list <- names(sl)[!grepl('gl|hap', names(sl))]
  
  # bin limits
  seqlengths_bin <- sapply(bin_ranges, function(i) tail(end(i), 1))
  
  
  
  covobj <- do.call('c', mclapply(chr_list, function(tmp_chr) {
    ifchr<-reads[['rname']] == tmp_chr
    ifpos<-reads[['strand']][ifchr]=='+'
    pos<-reads[['pos']][ifchr]
    if (shift_reads > 0L) {
      
      read_ranges <- IRanges(
        start = ifelse(
          ifpos,
          pos + shift_reads,
          pos - shift_reads
        ),
        width = rwidth
      )
    } else {
      read_ranges <- IRanges(start = pos, width = rwidth)#add in strand to the iranegs object
    }
    
    read_score.p <- Rle(countOverlaps(bin_ranges[[tmp_chr]],read_ranges[ifpos], type = 'any'))
    read_score.n <- Rle(countOverlaps(bin_ranges[[tmp_chr]],read_ranges[!ifpos], type = 'any'))
    
    GRanges(
      seqnames = Rle(tmp_chr),
      seqlengths = seqlengths_bin,
      ranges = bin_ranges[[tmp_chr]],
      pscore = as.numeric(read_score.p),
      nscore = as.numeric(read_score.n),
      strand = '*'
    )
  }))
  
  
  
  
  covobj <- do.call('c', mclapply(chr_list, function(tmp_chr) {
    ifchr<-reads[['rname']] == tmp_chr
    ifpos<-reads[['strand']][ifchr]=='+'
    pos<-reads[['pos']][ifchr]
    if (shift_reads > 0L) {
      
      read_ranges <- IRanges(
        start = ifelse(
          ifpos,
          pos + shift_reads,
          pos - shift_reads
        ),
        width = rwidth
      )
    } else {
      read_ranges <- IRanges(start = pos, width = rwidth)#add in strand to the iranegs object
    }
    
    read_score.p <- Rle(countOverlaps(bin_ranges[[tmp_chr]],read_ranges[ifpos], type = 'any'))
    read_score.n <- Rle(countOverlaps(bin_ranges[[tmp_chr]],read_ranges[!ifpos], type = 'any'))
    
    GRanges(
      seqnames = Rle(tmp_chr),
      seqlengths = seqlengths_bin,
      ranges = bin_ranges[[tmp_chr]],
      pscore = as.numeric(read_score.p),
      nscore = as.numeric(read_score.n),
      strand = '*'
    )
  }))
  
  # normalize score
  #covobj@elementMetadata@listData$score <- score(covobj) / (1.35e+08 / sum(score(covobj)))
  covobj@elementMetadata@listData$score <- score(covobj) / (3.1e+09 / sum(score(covobj)))
  
  
  #define outputfile
  outfile=paste(sep='',outdir, strsplit(tail(strsplit(bam_file, '/')[[1]], 1), '[.]')[[1]][1], '.wig')
  #now add a color to the wig file
  write(file=outfile,x=paste0('track type=wiggle_0 color=',col[1],',',col[2],',',col[3]))
  # export wig
  export(covobj, outfile,append=T)
  
  
}

args <- commandArgs(trailingOnly=T)[1:5]
outfolder<-'/g/furlong/Harnett/TSS_CAGE_myfolder/data/solexa/wig/'

if(is.na( args[1])){
  message('usage - Rscript WigExport.R bamfile binsize outfolder readshift color.as.word')
  stop('no first argument')
}else{bam<-args[1]}

if(is.na(args[2])){
  message('no second argument,defaulting to bin width 10')
  bin<-10L
}else{bin<-as.numeric(args[2])}

if(is.na(args[3])){
  message(paste0('no third argument,defaulting to outfolder ',outfolder))
  
}else{outfolder<-as.numeric(args[3])}

if(is.na(args[4])){
  message('no fourth argument,defaulting to read shift 0')
  shift<-0
}else{shift<-as.numeric(args[4])}

if(is.na(args[5])){
  message('no fifth argument,defaulting to color black')
  wigcol<-'black'
}else{wigcol<-args[5]}


exportWig(col=col2rgb(wigcol),
          outdir=outfolder,
          shift_reads=0L,
          bin_width=bin,
          bam_file=bam)