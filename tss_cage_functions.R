
# Functions for use in the eRNA/CAGE project ------------------------------
# Author: Dermot Harnett --------------------------------------------------

#source("http://bioconductor.org/biocLite.R")
#biocLite(pkgs='chipseq',lib.loc='~/Harnett/R')
#library(chipseq,lib.loc='~/Harnett/R/')
#library(chipseq)
.libPaths( c( "~/Harnett/R", .libPaths()) )
library(GenomicRanges)
library(multicore)
library(rtracklayer)
library(Rsamtools)
library(BSgenome.Dmelanogaster.UCSC.dm3)
#library(devtools)
#install_github("pryr")
#library(pryr)#useful just for the f() function
library(VGAM,lib.loc='~/Harnett/R')
library(igraph)
library(ggplot2)
library(reshape2)
library(LSD)
tps=c('tp24h','tp68h','tp1012h')


#library(snow)
#library(spp)

file.unprocessed.cage.tag.rles='data/objects/unprocessed.cage.tag.object.R'
file.cg.pl<-'data/objects/cg.pl.object.R'
file.cage.tag.rles<-'data/objects/cage.rle.object.R'
file.alltags<-'data/objects/alltags.object.R'
file.alltags.gr<-'data/objects/alltags.gr.object.R'
file.cage.tag.rles.rpgc<-'data/objects/cage.rpgc.object.R'
file.accession.df<-'data/objects/accession.df.object.R'
file.alltags.rpgc<-'data/objects/alltags.rpgc.object.R'
file.dnase.rles<-'data/objects/dnase.rles.object.R'
file.crmgrs<-'data/objects/crmgrs.object.R'
file.lambdas<-'data/objects/lambdas.object.R'
file.transcripts<-'data/objects/transcripts.object.R'
file.tss<-'data/objects/tss.object.R'
file.gene.model.data<-'data/objects/gene.model.data.object.R'
file.synonyms<-'data/objects/synonyms.object.R'
file.crm.list.binary<-'data/TF8008.txt'
file.cad3<-'data/cad3.object.R'
file.dnase.peaks<-'data/objects/dnase.peaks.R'

chrs.keep<-seqlevels(Dmelanogaster)[!grepl('U|M', seqlevels(Dmelanogaster))]
si<-seqinfo(Dmelanogaster)[chrs.keep]
genome.size<-sum(seqlengths(Dmelanogaster)[chrs.keep])
bigchrs<-chrs.keep[1:6]
#cluster <- makeCluster(8);

#These abbreviations are mainlmey for convenience int he interpreter, should be commented out
#for final draft
h<-head
n<-names
l<-length
cn<-colnames
rn<-rownames

pdftmp<-function(expr){
  pdf('tmp.pdf')
    expr
    dev.off()
  }

#function which, given a vector of integers and a single large vector,
#splits the large vector into chunks of the sizes specified. Can also do rles
chunkvect <- function(x, chunks, as.rle=F) {
  stopifnot(sum(chunks)==length(x))
  chunks=cumsum(chunks)
  chunkstarts=c(1,(chunks[-length(chunks)]+1))
  if(as.rle==T){func=Rle
  }else{func=identity}
  out=lapply(seq_along(chunks),function(n){
    func(x[chunkstarts[n]:chunks[n]])
  })
  names(out)=names(chunks)
  return(out)
}


qnormvect<-function(x){ x = qnorm((rank(x,ties.method='average')-0.5)/length(x));x}


##function to do power law normalization
fit.and.plnorm<-function(srles,lowerlim=200,refsize=10^8){
  ##function to do power law normalization
  fits<-sapply(names(srles),function(acc){

    sitecounts=as.vector(c(
      unlist(unname(srles[[acc]][['pos']][srles[[acc]][['pos']]!=0])),
      unlist(unname(srles[[acc]][['neg']][srles[[acc]][['neg']]!=0]))
    ))
    fit=power.law.fit(as.vector(sitecounts),xmin=lowerlim)#not that power.law.fit excludes the lower count numbers
    o=getOffset(fit$alpha,sum(sitecounts))#calculate the offset (second parameter)
    c(alpha=fit$alpha,offset=o)#now tack it onto the fit list and return
  })
  mean.alpha=mean(fits['alpha',])
  r.offset<-getOffset(mean.alpha,refsize)
  #4a now finally normalize all of our cage libraries to a common power law reference
  normed<-mapply(SIMPLIFY=F,names(srles),fits['alpha',],fits['offset',],FUN=function(acc,alpha,offset){

     lapply(srles[[acc]],function(x){ 

      #this is the normalization bit...
      beta = alpha/mean.alpha
      lambda = (r.offset/offset)^ (1/mean.alpha)
      x = x ^ beta
      x = x * lambda#do it like this to keep names
      return( x)

    })

  })
  return(normed)
}

#just ot make sure I'm going to try the normalization the way they have it written (wrong i think)


##function to do power law normalization
fit.and.plnorm.wrong<-function(srles,lowerlim=200,refsize=10^8){
  # stopifnot(exists(getOffset))
  require(VGAM,lib.loc='~/Harnett/R')

  ##function to do power law normalization
  fits<-sapply(names(srles),function(acc){

    sitecounts=as.vector(c(
      unlist(unname(srles[[acc]][['pos']][srles[[acc]][['pos']]!=0])),
      unlist(unname(srles[[acc]][['neg']][srles[[acc]][['neg']]!=0]))
    ))
    fit=power.law.fit(as.vector(sitecounts),xmin=lowerlim)#not that power.law.fit excludes the lower count numbers
    o=getOffset(fit$alpha,sum(sitecounts))#calculate the offset (second parameter)
    c(alpha=fit$alpha,offset=o)#now tack it onto the fit list and return
  })
  mean.alpha=mean(fits['alpha',])
  r.offset<-getOffset(mean.alpha,refsize)
  #4a now finally normalize all of our cage libraries to a common power law reference
  normed<-mapply(SIMPLIFY=F,names(srles),fits['alpha',],fits['offset',],FUN=function(acc,alpha,offset){

     lapply(srles[[acc]],function(x){ 

      #this is the normalization bit...
      beta=alpha/mean.alpha
      lambda = (r.offset/offset)^ (mean.alpha)
      x = x ^ beta
      x = x * lambda#do it like this to keep names
      return( x)

    })

  })
  return(normed)
}

# fit.and.plnorm(cg[1:4])

#define our power law  normalization function
pl.norm<-function(x,x.offset,x.alpha,ref.offset=r.offset,ref.alpha=mean.alpha){
  #take in values x and normalize them to a power law distribution of offset and slope given.
  #again, see  http://genomebiology.com/content/10/7/R79 for the math

  beta=x.alpha/ref.alpha
  lambda = (ref.offset/x.offset)^ (1/ref.alpha)
  x = x ^ beta
  x = x * lambda#do it like this to keep names
  return( x)

}

#function to generate negative binomial samples with a given mean and variance
rnbinom2<-function ( n,m,v ){ 
  s= ( m^2 )/ ( v-m )
  rnbinom ( n=n,mu=m,size=s )
 }
dnbinom2<-function( x , m , v ){
  s= ( m^2 )/ ( v-m )
  dnbinom (x,mu=m,size=s)
}

#function to derive empiracal pvalue for the coorrelation between two vectors via shuffling
vector.pvals.shuffle <-function ( a,b,iter=iterations ){ 
  cor.real<-cor ( a,b,method='s' )
  cors.shuffle<-replicate ( iter,{ cor ( a,sample ( b,length ( b ) ),method='s' ) } )
  mean ( cors.shuffle>cor.real )
 }#this does give uniform p values with normal input
lognz<-function(x){x[x!=0]<-log10(x[x!=0]);x }



#function to determine offset of power law distribution given total and a slope
getOffset<-function(alpha,T=1000000){
  require(VGAM)
  #see http://genomebiology.com/content/10/7/R79
  #we know the total number of tags T in a power law distribution is r0 * riemannZeta(alpha)
  #where r0 is the offset, and alpha is the slope of our power law. So we can work out the offset given T and alpha with...
  r0=T/zeta(alpha)
  return(r0)
  #confirmed this works using numbers from paper (slope of 1.25 and T=1000000 gives r0 of 217,623)
}


#strips mcols from a gr object
nomcols<-function(gr){mcols(gr)<-NULL;gr}#function that strips mcols from a GRanges object


combinegrs<-function(grlist){
  #functint hat accepts a list of grange objects
  metadata<-sapply(grlist,mcols)
  grs<-sapply(grlist,nomcols)
  grs<-sapply(names(grs),function(n){grs[[n]]$reg<-n;grs[[n]]})
  grs<-sort(do.call('c',unname(grs)))#gr object with all the original
  cols<-unique(unlist(sapply(metadata,colnames)))
  for(col in cols){mcols(grs)[[col]] <- Rle(NA)}
  for(col in cols){
    for(n in names(grlist)){
      if(!is.null(mcols(grlist[[n]])[[col]])){
        if(is.list(mcols(grlist[[n]])[[col]])){mcols(grs)[[col]]<-NA}
        mcols(grs)[[col]][grs$reg==n] <- mcols(grlist[[n]])[[col]]
      }
    }
  }
  return(grs)
  #extract all of our column names in our metadata
}

# Carries out Views on a GRange object ------------------------------------
GRViews<-function(rle,gr){
  stopifnot(gr==sort(gr))
  v=Views(rle[seqlevels(gr)],as(gr,'RangesList'))[seqlevels(gr)]
  v=v[sapply(v,function(x){length(x)!=0})]
  return(v)
}


# #expands a GRanges function for using it on a binned Rle Object ---------




get.cage.matrix<-function(gr,clist){
  #gr must be sorted for this to work properly
  stopifnot(gr==sort(gr))
  stopifnot(seqlevels(gr)==names(clist[[1]]))
  stopifnot(sapply(clist,function(c){names(c)==names(clist[[1]])}))
  gr.r<-as(gr,'RangesList')
  tmp<-mclapply(mc.cores=10,names(clist),function(acc){
    unlist(viewSums(Views(clist[[acc]][chrs.keep],gr.r[chrs.keep])))
  })
  tmp<-simplify2array(tmp)
  colnames(tmp)<-names(clist)
  tmp
}


viewOnBinned<-function(rlelist,gr){
  #expand the ranges
  rlelist<-as(rlelist,'SimpleRleList')
  chr<-names(seqlengths(gr))[1]
  binsize=ceiling(( seqlengths(gr)[[chr]])/ length(rlelist[[chr]]))
  start(gr)<-round(start(gr)/binsize)
  end(gr)<-round(end(gr)/binsize)
  seqlengths(gr)<-round(seqlengths(gr)/binsize)
  #drop chrs not in the list
  listchrs<-names(rlelist)
  gr<-gr[ seqnames(gr) %in% listchrs ]
  seqlevels(gr)<-seqlevels(gr)[seqlevels(gr)%in%listchrs]
  return(Views(rlelist,as(gr,'RangesList')))
  
}

# bam2coverage ------------------------------------------------------------
#Read in a bam file and return the coverage as and RleList object - no normalization
bam2coverage<-function(bam,doshift=F,doresize=F,stranded=F,fragment.length=NA,use.spp=F,is.FirstMateRead=NA,is.ProperPair=NA,use.size=F,maxsize=NA){
  cat(paste0('attempting to read bam file  ',bam,'  \n'))
  #scan the bam file
  reads <- scanBam(
    file = bam,
    param = ScanBamParam(
      what = c('pos', 'rname', 'strand', 'qwidth','isize'),
      flag = scanBamFlag(isUnmappedQuery = FALSE,isFirstMateRead=is.FirstMateRead,isProperPair=is.ProperPair)
    )
  )[[1]]

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
  if(!all(unique(reads[['rname']]) %in% names(seqinfo(Dmelanogaster)))){
    message(paste0('these chrs arent right',unique(reads[['rname']])))
    stop()
  }
  stopifnot(length(unique(reads[['qwidth']]))==1)
  #don't keep chrM or chrU
  reads.keep<-as.logical(reads[['rname']] %in% chrs.keep)
  #reads.keep<-as.logical(reads[['rname']] %in% 'chr2R')

  # keep isize forlater
  isize=reads$isize[reads.keep]

  #now as a GRange object
  reads<-GRanges(reads[['rname']][reads.keep],
             IRanges(reads[['pos']][reads.keep],
                     width=reads[['qwidth']][reads.keep]),
             strand=Rle(as.character(reads[['strand']][reads.keep])),
             seqinfo=si
  )
  taglength<-width(reads)[[1]]
  ifpos<-strand(reads)=='+'

  if(!is.na(maxsize)){
    reads=reads[isize<maxsize]
    isize=isize[isize<maxsize]
  }
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

  if(use.size){
    doresize=T
    fragment.length=abs(isize)
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


# . -----------------------------------------------------------------------


#function producing coverage plot fo a given range and rleList object
rlecovplot<-function(chr='chr2R',start=10020000,width=20000,rlelist1=chrom.rles.rpgc.sub.merge[[1]],rlelist2=F,seqinf=si,g=NA){
  if('GRanges' %in% is(g)){
    seqinfo(g)<-si
    g<-as(g,'RangesList')
  }else{g<-as(GRanges(chr,IRanges(start,start+width),seqinfo=si),'RangesList')}
 
  
  v1=Views(rlelist1,g)[[chr]]
  if('SimpleRleList' %in% is(rlelist2)){
    v2=Views(rlelist2,g)[[chr]]
    coverageplot(v1,v2)
  }else{
    coverageplot(v1)
  }  
}




# Based off Tibor's heatmap one liner --------------------------------------

heatmap.simple<-function(x,col1='white',col2='black',brks='quantile',filt=NA){
  if(!is.na(filt)){z=x[x>filt]
  }else{z<-x}
  if(brks=='quantile'){brks= quantile(z, probs = seq(0, 1, 0.01))}
  else if(brks=='uniform'){brks= seq(min(z),max(z),length.out=100)
  }else(warn('brks should be quantile, uniform'))
  image(t(x[nrow(x):1, ]), 
        breaks =brks,
        col = colorRampPalette(c(col1, col2))(length(brks ) - 1L),
        xaxt = 'none',
        yaxt = 'none')
}



# . -----------------------------------------------------------------------


#function that splits a GR into width 1 intervals
width1gr<-function(gr){
  
  #for each unique length bigger than 1
  extras<-sapply(unique(width(gr[width(gr)!=1])),function(w){
    #select those intervals equal to or larger than this
    tosplit<-gr[width(gr)>=w]
    #create width 1 intervals shifted forward by that amount
    shift(resize(tosplit,width=1,fix='start'),w-1)
  })
  gr<-resize(gr,width=1,fix='start')#resize the originals
  sort(c(gr,do.call('c',extras)))#concatenate the extras and sort
}


# from bioconductor website -----------------------------------------------


rollmeanRle <- function (x, k)
   {
     n <- length(x)
     cumsum(c(Rle(sum(window(x, 1, k))), window(x, k + 1, n) - window(x, 1, n - k)))/k
     }






# Clustering function -----------------------------------------------------


#This function returns a grange object of clusters when handed a SimpleRleList
#of tags. It ignores tag height, accept that it eliminates clusters whose maximum
#is lower than the given cutoff. Needs work...
Hcluster.cage.rle<-function(tags,cutoff=5){
  tags.gr<-as(tags,'GRanges')
  tags.gr<-tags.gr[tags.gr$score!=0]
  #we can just ignore anything lower than a certain threshold
  
  blocks.gr<-resize(tags.gr,width=301,fix='center')
  blocks.gr<-reduce(blocks.gr)
  
  #   tmp<-blocks.gr[sample(1:length(blocks.gr),replace=F,size=1000)]
  #   v<-Views(tags>0,as(blocks.gr,'RangesList'))
  #   cl<-viewApply(v,FUN=function(x){
  #     
  #     sites<-which(x)
  #     
  #     if(min(sites)+300 >= max(sites)){#and if all are within 300 bp return a cluster
  #       return(list(min(sites),
  #                   max(sites)))
  #     }else{
  #       #else do hierarchical clustering
  #       clusts<-cutree(hclust(dist(sites)),h=300)
  #       list(unique(clusts),function(clust){
  #         list(min(sites[clusts==clust]),
  #              max(sites[clusts==clust]),
  #              as.character(seqnames(blocks.gr)[block]))
  #       })
  #       
  #     }
  #     
  #   })                      
  #   
  #   
  ov<-findOverlaps(blocks.gr,tags.gr)
  #ov<-ov[1:10000,]
  blocks.pos<-do.call(cbind,mclapply(mc.cores=20,unique(ov@queryHits),function(block){
    sites<-ov@subjectHits[ov@queryHits==block]
    sites<-start(tags.gr)[sites]
    
    if(length(sites)==1){return(list(sites,#if we've only got one tag return a cluster
                                     sites,
                                     as.character(seqnames(blocks.gr)[block])))}
    
    if(min(sites)+300 >= max(sites)){#and if all are within 300 bp return a cluster
      return(list(min(sites),
                  max(sites),
                  as.character(seqnames(blocks.gr)[block])))
      
    }else{
      #else do hierarchical clustering
      clusts<-cutree(hclust(dist(sites)),h=300)
      
      
      sapply(unique(clusts),function(clust){
        list(min(sites[clusts==clust]),
             max(sites[clusts==clust]),
             as.character(seqnames(blocks.gr)[block]))
      })
      
    }
  }))
  
  blocks.gr<-GRanges(unlist(blocks.pos[3,]),IRanges(unlist(blocks.pos[1,]),unlist(blocks.pos[2,])),seqinfo=si)
  
  #now we merge clusters which are less than 50bp apart
  blocks.gr<-resize(blocks.gr,width=width(blocks.gr)+50,fix='start')
  blocks.gr<-reduce(blocks.gr)
  blocks.gr<-resize(blocks.gr,width=width(blocks.gr)-50,fix='start')
  
  #now mark blocks which have only one tagsite and are below the cutoff
  blocks.gr$max<-unlist(viewMaxs(Views(tags,as(blocks.gr,'RangesList'))))
  # blocks.gr$isolated<-distanceToNearest(blocks.gr)$distance
  #blocks.gr$singleton <- width(blocks.gr)==1 & blocks.gr$max<cutoff
  blocks.gr$singleton <- blocks.gr$max<cutoff
  
  blocks.gr<-  blocks.gr[!blocks.gr$singleton ]                        
  
  return(blocks.gr)
}


getBestSingleWindow <- function (reg, w, chrs.keep, big=bigchrs, cage=cg) {
   #make sure th
  stopifnot(all(seqnames(reg)%in%big))
  stopifnot(reg==sort(reg))
  # sapply(simplify=F,c(10,100,500),function(w){
  #make sure all our regions are at least the window size
  reg<-resize(reg,width=pmax(width(reg),w),fix='center')
  cat('.')
  #Do views on them
  maxwinds<-mclapply(mc.cores=10,cage,winds,FUN=function(acc,winds){
    #now get the windowed views for each window
    v<-unlist(viewApply(Views(acc[['pos']][big],winds[big]),FUN=function(x){ (runsum(x,w)) }))
    v.n<-unlist(viewApply(Views(acc[['neg']][big],winds[big]),FUN=function(x){ (runsum(x,w)) }))
    #get the bigger of the two strands
    ispos<-(sapply(v,max)>sapply(v.n,max))
    v<-ifelse( ispos,v,v.n)
    sapply(v,max)
  })
  return(simplify2array(maxwinds))
}


  #make sure th
get.best.window.mat <- function (reg, w, chrs.keep, big=bigchrs, cage=cg) {
   #make sure th
  stopifnot(all(seqnames(reg)%in%big))
  stopifnot(reg==sort(reg))
  # sapply(simplify=F,c(10,100,500),function(w){
  #make sure all our regions are at least the window size
  reg<-resize(reg,width=pmax(width(reg),w),fix='center')
  cat('.')
    #convert them to a Rangelist for Views
  winds<-as(reg,'RangesList')[big]
  winds<-winds[big]
  cat('.')
  # starts<-unlist(start(winds))
  #chrs<-names(starts) 
  #for each line  
  maxwinds<-mclapply(mc.cores=10,cage,winds,FUN=function(acc,winds){
    #now get the windowed views for each window
    v<-unlist(viewApply(Views(acc[['pos']][big],winds[big]),FUN=function(x){ (runsum(x,w)) }))
    v.n<-unlist(viewApply(Views(acc[['neg']][big],winds[big]),FUN=function(x){ (runsum(x,w)) }))
    #get the bigger of the two strands
    ispos<-(sapply(v,max)>sapply(v.n,max))
    v<-ifelse( ispos,v,v.n)
    sapply(v,max)
  })
  return(simplify2array(maxwinds))
}



get.best.chr.window.mat <- function (reg, w, chrs.keep, big=bigchrs, cage=cg) {
  
  # sapply(simplify=F,c(10,100,500),function(w){
  #make sure all our regions are at least the window size
  reg<-resize(reg,width=pmax(width(reg),w),fix='center')
  cat('.')
  #make sure th
  stopifnot(all(seqnames(reg)%in%big))
  stopifnot(reg==sort(reg))
  
  #convert them to a Rangelist for Views
  winds<-as(reg,'RangesList')[big]
  winds<-winds[big]
  cat('.')
  maxwinds<-mclapply(mc.cores=10,cage,winds,FUN=function(acc,winds){
    ispos<-unlist(viewApply(Views(acc[big],winds[big]),FUN=function(x){ max(runsum(x,w)) }))
  })
  return(simplify2array(maxwinds))
}


fill_density_plot<-function(dens,cutoff,colname='pink'){
  x1 <- which(dens$x >= log10(cutoff))[1]
  x2 <- length(dens$x)
  y1<-dens$y[x1]
  y2<-dens$y[x2]
  with(dens, polygon(x=c(x[c(x1,x1:x2,x2)]), y= c(0, y[x1:x2], 0), col=colname))
}


###Actually screw this just use DESeqs functions
deseqNormalize<-function(m,return.just.sizefactors=F){
#take in a matrix and perform DEseqs normalization on it. The columns should be big,
#if there's a list of matrices require they're all the same size, concatenate them,
#then spit them again at the end
  islist=F
  if(is(m,'list')){
    islist=T
    nrows<-sapply(m,nrow)
    splitfacts<-unlist(mapply(names(m),nrows,FUN=function(name,num){rep(name,num)}))
    m<-do.call(rbind,m)
  }



  require(DESeq)
  oldrownames<-rownames(m)
  rownames(m)<-1:nrow(m)#no duplicate names

  conditions=colnames(m)#
  cds=newCountDataSet(m,conditions)#create the DESeq object
  
  cds=estimateSizeFactors(cds)#returns object with size factors

  if(return.just.sizefactors){return(sizeFactors(cds))}

  m=counts(cds,normalized=T)#returns the normalized matrix

  rownames(m) =oldrownames#put the old rownames back in

  #join the list again if necessary
  if(islist){
    nrows=cumsum(nrows)
    starts=c(0,nrows)[1:length(nrows)]+1
    m=sapply(seq_along(nrows),function(n){  m[starts[n]:nrows[n],]})
    names(m)=names(nrows)
  }
  #return normalized matrix or list thereof
  return(m)

}




# Define functions for plotting each set ----------------------------------



ROCfunc <- function (pos,neg,main,sub) {
  scores<-stack(list(pos=pos,neg=neg))
  sumpred<-prediction(scores$values,scores$ind)
  perf<-performance(sumpred,'tpr','fpr')
  plot(perf,colorize=T,lwd=6)
  auc=slot(performance(sumpred,'auc'),'y.values')[[1]]
  aucstring=paste0('AUC = ',round(auc,4))
  text(x=0.90,y=0,labels=aucstring)
  
  posstring=paste0('Positive Set - ', length(pos))
  negstring=paste0('Negative Set - ', length(neg))
  text(x=0.85,y=0.2,labels=posstring)
  text(x=0.85,y=0.1,labels=negstring)
  
  title(main,sub)
  abline(coef=c(0,1),lty=3)
  grid(lwd=2)
}

PreRecfunc <- function (pos,neg,main,sub) {
  scores<-stack(list(pos=pos,neg=neg))
  sumpred<-prediction(scores$values,scores$ind)
  perf=performance(sumpred,'prec','rec')
  plot(perf,colorize=T,ylim=c(0.5,1))
  
  posstring=paste0('Positive Set - ', length(pos))
  negstring=paste0('Negative Set - ', length(neg))
  text(x=0.2,y=0.4,labels=posstring)
  text(x=0.2,y=0.35,labels=negstring)
  
  title(main,sub)
  cutoff<-max(perf@x.values[[1]][perf@y.values[[1]]>prec.cutoff],na.rm=T)
  
  grid(lwd=2)
  abline(h=prec.cutoff,lty=3)
  abline(v=cutoff,lty=3)
  
  perf<-performance(sumpred,'prec')
  cutoff<-min(perf@x.values[[1]][perf@y.values[[1]]>prec.cutoff],na.rm=T)
  cutstring=paste0('Cutoff at precision ',prec.cutoff,' - ',round(cutoff,4),' Normalized Tags')
  text(x=0.4,y=0.55,labels=cutstring)
  return(cutoff)
}
# list(pos=c(rep(T,30),rep(F,9)),neg=c(rep(F,100),rep(T,80),rep(F,7)))->above.cut.vects

Barplot.from.matrices<-function(posmat,negmat, logy=T,yname='Normalized Cage Signal 95% cl', xname='CRM',above.cut.vects=NA, barsize=2,cut=1,  bar=T,reorder=F  ){
  
  rownames(posmat)<-(1+nrow(negmat)):(nrow(posmat)+nrow(negmat))
  rownames(negmat)<-1:nrow(negmat)
  posmelt<-melt(posmat)
  negmelt<-melt(negmat)
  bothmelt<-rbind(posmelt,negmelt)
  bothmelt$set<-c(rep('Positive',nrow(posmelt)),rep('Negative',nrow(negmelt)))
  if(length(above.cut.vects)!=1){

    stopifnot(is(above.cut.vects,'list'))
    stopifnot(is(above.cut.vects,'vector'))

    bothmelt$abovecut<-c(above.cut.vects[[1]],above.cut.vects[[2]])

  }



  if(reorder==T){
    #NOW reorder the factor levels
    bothmelt$Var1 <-factor(bothmelt$Var1, levels = as.character(unique(bothmelt$Var1)))
    bothmelt$Var2 <-factor(bothmelt$Var2, levels = as.character(unique(bothmelt$Var2)))
    bothmelt$set <-factor(bothmelt$set, levels = as.character(unique(bothmelt$set)))
    
    means=by(bothmelt,bothmelt$Var1,FUN=function(df){mean(df$value)})
    names(means)<-unique(bothmelt$Var1)
    bothmelt<-bothmelt[order(means[as.character(bothmelt$Var1)]) ,]
    bothmelt<-bothmelt[order(bothmelt$set) ,]
    #NOW reorder the factor levels
    bothmelt$Var1 <-factor(bothmelt$Var1, levels = as.character(unique(bothmelt$Var1)))
    bothmelt$Var2 <-factor(bothmelt$Var2, levels = as.character(unique(bothmelt$Var2)))
  }
  
  bothmelt$logval <-bothmelt$value
  if(logy ){
    # m<-min(bothmelt$logval[ bothmelt$logval != 0 ])
    # bothmelt$logval <- bothmelt$logval / m 
    bothmelt$logval[ bothmelt$logval==0]<-1

  }

  #note that we take the logarithm of the logvals with ggplot later

  yscale=ifelse(logy,scale_y_log10(name=yname,breaks=c(1,10,100,1000,10000,100000)),scale_y_continuous(name=yname))

  if(length(above.cut.vects)!=1){
    return(
      #line graph of our se
      ggplot(data=bothmelt,aes(x=Var1,y=logval,color=set))+stat_summary(fun.data="mean_cl_boot",geom="errorbar",size=barsize)
      +
        yscale+
        #scale_y_continuous(name=yname)+
        scale_x_discrete(name=xname,labels='')+
        ggtitle('95% cl for mean of  normalized Signal')
    )
  }else{
    return(
       #line graph of our se
      ggplot(data=bothmelt,aes(x=Var1,y=logval,color=set))+stat_summary(fun.data="mean_cl_boot",geom="errorbar",size=barsize)+
        scale_y_log10(name=yname,breaks=c(1,10,100,1000,10000,100000))+
        #scale_y_continuous(name=yname)+
        scale_x_discrete(name=xname,labels='')+
        ggtitle('95% cl for mean of normalized Signal')
    )
  }  
}


# Barplot.from.matrices(posmat[1:10,],negmat[1:10,],bar=T)
Chrom.boxplots<-function(posmat,negmat, lowlim=F,yname='Mean Inp-Sub Chromatin Signals for crms above/below cutoff ', xname='Set',tit=''){
  
  rownames(posmat)<-(1+nrow(negmat)):(nrow(posmat)+nrow(negmat))
  rownames(negmat)<-1:nrow(negmat)
  posmelt<-melt(posmat)
  negmelt<-melt(negmat)
  bothmelt<-rbind(posmelt,negmelt)
  bothmelt$set<-c(rep('Above',nrow(posmelt)),rep('Below',nrow(negmelt)))
  bothmelt$nset<-c(rep(nrow(posmat),nrow(posmelt)),rep(nrow(negmat),nrow(negmelt)))

  #we need to do something iwth our zero/negative values to do a log scale....
  if(lowlim==F){
    lowlim=min(bothmelt$value[bothmelt$value>0])
  }
  bothmelt$value[bothmelt$value<=0] <- lowlim
  
  num.df <-  data.frame(set=c('Above','Below'),number=c(nrow(posmat),nrow(negmat)))
  num.df <- do.call(rbind, lapply(unique(bothmelt$Var2),function(v) cbind(num.df,v) ))
  num.df$Var2<-num.df$v
  num.df$yval<-min(bothmelt$value)


  #line graph of our se
  ggplot(data=bothmelt,aes(x=set,facet=Var2,y=value,color=Var2))+geom_boxplot()+
    scale_y_log10(name=yname)+    
    facet_wrap(~Var2)+
    #scale_y_continuous(name=yname)+
    scale_x_discrete(name=xname,labels=c('Above','Below'))+     
    geom_text(data=num.df ,aes(x=set,y=yval,label=number))+
    ggtitle(tit)+
    theme(strip.text=element_text(size=10))
}



dotplot.from.matrices<-function(posmat,negmat,yname='Normalized Cage Signal', xname='CRM'){
  
  rownames(posmat)<-(1+nrow(negmat)):(nrow(posmat)+nrow(negmat))
  rownames(negmat)<-1:nrow(negmat)
  posmelt<-melt(posmat)
  negmelt<-melt(negmat)
  bothmelt<-rbind(posmelt,negmelt)  
  bothmelt$CRM<-bothmelt$Var1
  bothmelt$set<-c(rep('Positive',nrow(posmelt)),rep('Negative',nrow(negmelt)))
  
      #line graph  our se
      ggplot(data=bothmelt,aes(x=CRM,y=value,color=set))+geom_point()+
      scale_y_log10(name=yname,breaks=c(1,10,100,1000,10000,100000))
       
}


simple.barplot.matrices<-function(posmat,negmat,yname='Normalized Cage Signal', xname='CRM',main.tit='tit'){
  
  rownames(posmat)<-(1+nrow(negmat)):(nrow(posmat)+nrow(negmat))
  rownames(negmat)<-1:nrow(negmat)
  posmelt<-melt(posmat)
  negmelt<-melt(negmat)
  bothmelt<-rbind(posmelt,negmelt)
  bothmelt$set<-c(rep('Positive',nrow(posmelt)),rep('Negative',nrow(negmelt)))
  bothmelt$CRM<-bothmelt$Var1
  bothmelt[['Normalized_Signal']]<-bothmelt$value

  
      #line graph  our se
      ggplot(data=bothmelt,aes(x=CRM,y=Normalized_Signal,color=set))+
      stat_summary(fun.data="mean_cl_boot",geom="errorbar")+
      ggtitle(main.tit)

      # scale_y_log10(name=yname,breaks=c(1,10,100,1000,10000,100000))
       
}


number.not.zero<-function(mat){apply(mat,1,function(x)sum(x!=0))}

#for now we'll just score the matrices by adding up the columns (libraries)
scorefunc<-function(x){rowMeans(x)}
#setn=names(possetlist)[6]




# dev.off()
# rlelist2=chrom.rles$PolII_6.8_R1
# rlecovplot(rlelist1=input.rles$Input_6.8_R1,rlelist2=chrom.rles$PolII_6.8_R1,g=resize(crmgrs[1],20000,fix='center'))