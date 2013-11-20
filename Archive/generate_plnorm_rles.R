###This is a script to generate RLE objects out of the bam files for our CAGE data

# @ author Dermot Harnett, EMBL Heidelberg
# @date 25/9/2013
# @title Export coverage into wig format
########################################
setwd(dir='/g/furlong/Harnett/TSS_CAGE_myfolder/')
source('src/tss_cage_functions.R')
library(VGAM,lib.loc='~/Harnett/R')
library(igraph,lib.loc='~/Harnett/R')
load(file.unprocessed.cage.tag.rles)

T=100000000#library size to normalize TO.

  #function to get slope of fitted alphas 
  fits<-sapply(names(cage.tag.rles)[],function(acc){
    sitecounts=sort(c(unlist(unname(cage.tag.rles[[acc]][['pos']])),unlist(unname(cage.tag.rles[[acc]][['neg']]))))
    sitecounts=sitecounts[sitecounts!=0]
    fit=power.law.fit(as.vector(sitecounts),xmin=200)#not that power.law.fit excludes the lower count numbers
    total.tags=sum(sitecounts)
    o=getOffset(fit$alpha,total.tags)#calculate the offset (second parameter)
    c(fit,offset=o,tagcount=total.tags)#now tack it onto the fit list and return
    })

  #average these
  mean.alpha=mean(unlist(fits['alpha',]))
  r.offset<-getOffset(mean.alpha,T)

  #now finally normalize each of our cage libraries
  cg.pl<-mapply(SIMPLIFY=F,names(cage.tag.rles),fits['alpha',],fits['offset',],FUN=function(acc,alpha,offset){
     lapply(cage.tag.rles[[acc]],function(srle){ 
      pl.norm(srle,x.alpha=alpha,x.offset=offset[[1]])
     })
  })

  #and save it
  save(cg.pl,file=file.cg.pl)
  
  message('CAGE Rle objects calculated, now summing reseqs')


  reseqs<-names(cg.pl)[grepl(names(cg.pl),pattern='reseq')]
  
  for(reseq in reseqs){
    orig<-gsub(reseq,pattern='_reseq',replacement='')
    cg.pl[[orig]]$pos<-cg.pl[[orig]]$pos+cg.pl[[reseq]]$pos
    cg.pl[[orig]]$neg<-cg.pl[[orig]]$pos+cg.pl[[reseq]]$neg
  }
  
  #now get rid of the reseqs  
  cg.pl<-cg.pl[!names(cg.pl) %in% reseqs]
  

  message('CAGE Rle objects calculated, now summing replicates')


  replicates<-names(cg.pl)[grepl(names(cg.pl),pattern='(+\\d)_?r2')]
  
  for(replicate in replicates){
    orig<-gsub(reseq,pattern='(+\\d)_?r2',replacement='\\1')
    cg.pl[[orig]]$pos<-cg.pl[[orig]]$pos+cg.pl[[replicate]]$pos
    cg.pl[[orig]]$neg<-cg.pl[[orig]]$pos+cg.pl[[replicate]]$neg
  }
  
  #now get rid of the replicates  
  cg.pl<-cg.pl[!names(cg.pl) %in% replicates]

  #dump chrs we don't want
  cg.pl<-rapply(cg.pl,how='replace',function(srle){srle[chrs.keep]})
  
  #sort the cage by the dataframex
  load(file.accession.df)
  cg.pl<-cg.pl[as.character(accession.df$accession)]

  message('summing cage tags')

  # Create summed CAGE Objects ----------------------------------------------
  alltags.pl<-cg.pl[[1]]
  for(n in 1:length(cg.pl[-1])){
    alltags.pl[['pos']]<-cg.pl[[n+1]][['pos']]+alltags.pl[['pos']]
    alltags.pl[['neg']]<-cg.pl[[n+1]][['neg']]+alltags.pl[['neg']]
    alltags.pl[['both']]<-cg.pl[[n+1]][['both']]+alltags.pl[['both']]
    
  }

  save(cg.pl,file=file.cg.pl)
  save(alltags.pl,file=file.alltags.pl)
  
   export(alltags.pl[['pos']],'data/solexa/wig/allcage.pl.pos.bedGraph')
   export(alltags.pl[['neg']],'data/solexa/wig/allcage.pl.neg.bedGraph')
   export(alltags.pl[['both']],'data/solexa/wig/allcage.pl.both.bedGraph')


for(acc in names(cg.pl)){
  export(cg.pl[[acc]][['pos']],paste0('data/solexa/wig/pos_indiv/',acc,'.plnorm.pos.bw'))
  export(cg.pl[[acc]][['neg']],paste0('data/solexa/wig/neg_indiv/',acc,'.plnorm.neg.bw'))
}

#now, get rid of the naems for online viewing

for(acc in 1:length(cg.pl)){
  cacc=paste0('c',acc)
  cat(cacc)
  export(cg.pl[[acc]][['pos']],paste0('data/solexa/wig/pos_indiv/',cacc,'.plnorm.pos.bw'))
  export(cg.pl[[acc]][['neg']],paste0('data/solexa/wig/neg_indiv/',cacc,'.plnorm.neg.bw'))
}
