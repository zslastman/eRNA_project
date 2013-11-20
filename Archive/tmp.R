
# bam2coverage ------------------------------------------------------------
#Read in a bam file and return the coverage as and RleList object - no normalization
bam2coveragePE<-function(bam,doshift=F,doresize=F,stranded=F,fragment.length=NA,use.spp=F,is.FirstMate){
  
  cat(paste0('attempting to read bam file  ',bam,'  \n'))
  #scan the bam file
  reads <- scanBam(
    file = bam,
    param = ScanBamParam(
      what = c('pos', 'rname', 'strand', 'qwidth','isize'),
      flag = scanBamFlag(isUnmappedQuery = FALSE,isFirstMate=T)
    )
  )[[1]]
  
  length(reads[[1]])
  
  #some isizes negative?
  
  
  message('..bam read..')
  #fix the names - necessary for some bams,by adding 'chr'
  if(!length(reads[['rname']][grepl(reads[['rname']],pattern='^2R$')])==0 ){
    reads[['rname']] <- Rle(paste0('chr',reads[['rname']]))
  }
  
  
  #and fixing the name of the mitochondrial genome
  if(!length(reads[['rname']][grepl(reads[['rname']],pattern='itochond')])==0 ){
    reads[['rname']][grepl(reads[['rname']],pattern='itochond')]<-  'chrM'
  }
  #make sure all the seqnames are okay
  stopifnot(unique(reads[['rname']]) %in% names(seqinfo(Dmelanogaster)))
  stopifnot(length(unique(reads[['qwidth']]))==1)
  #don't keep chrM or chrU
  reads.keep<-as.logical(reads[['rname']] %in% chrs.keep)
  #reads.keep<-as.logical(reads[['rname']] %in% 'chr2R')
  
  
  
  #now as a GRange object
  reads<-GRanges(reads[['rname']][reads.keep],
                 IRanges(reads[['pos']][reads.keep],
                         width=reads[['qwidth']][reads.keep]),
                 strand=Rle(as.character(reads[['strand']][reads.keep])),
                 seqinfo=si
  )
  taglength<-width(reads)[[1]]
  ifpos<-strand(reads)=='+'
  #coverageplot(peaks.pos$chr10[wpeaks[1]], peaks.neg$chr10[wpeaks[1]])
  
  
  #this calculates the cross correlation betwen strands to estimate the fragment length
  #fragment.length<-estimate.mean.fraglen(g[seqnames(g)=='chr2R'],method='correlation')
  
  #we'll just us ccf, and only on chr2L to save time (won't differ for other chrs)
  if((doshift || doresize)&(is.na(fragment.length))){
    message(paste0('...calculating fragment length....'))
    maxfraglength<-300
    fragment.length<-which.max(
      ccf(plot=F,
          as.vector(coverage(reads[strand(reads)=='+' & seqnames(reads) == 'chr2R'])[['chr2R']])[1000000:5000000],
          as.vector(coverage(reads[strand(reads)=='-' & seqnames(reads) == 'chr2R'])[['chr2R']])[1000000:5000000],
          lag.max=maxfraglength)$acf
    )
    fragment.length<-abs(fragment.length)
    message(paste0('.fragment length ',fragment.length))
    
    stopifnot(fragment.length>50)
    
  }
  
  
  #shift reads
  if(doshift==T){
    #move the negative ones towards the start, stopping at one
    reads[!ifpos]<-shift(reads[!ifpos],-pmin(round(fragment.length*0.5),start(reads[!ifpos])-1))
    reads[ifpos]<-shift(reads[ifpos],pmin(round(fragment.length*0.5),seqlengths(reads[ifpos])[as.character(seqnames(reads[ifpos]))]-end(reads[ifpos])))
  }
  #do resize
  if(doshift==F & doresize==T){
    message(paste0('..reads length resized... '))
    reads<-resize(reads,fragment.length)
  } 
  
  #and return the coverage
  if(stranded==F){return(coverage(reads))
  }else{
    return(#or the stranded coverage
      list(pos=coverage(reads[ifpos,]),
           neg=coverage(reads[!ifpos,])
      )
    )
  }  
}
