setwd('~/Harnett/TSS_CAGE_myfolder/')
rm(list=ls())

source('src/tss_cage_functions.R')
library(reshape2)
library(ggplot2)
load(file.cage.tag.rles)
load(file.accession.df)
regs=list(
  pos=import(con='analysis/positive_8008.bed',asRangedData=F,seqinfo=si),
  neg=import(con='analysis/negative_8008.bed',asRangedData=F,seqinfo=si),
  int=import(con='analysis/random.intergenic.bed',asRangedData=F,seqinfo=si),
  cad3pos=import(con='analysis/cad3.pos.bed',asRangedData=F,seqinfo=si),
  cad3neg=import(con='analysis/cad3.neg.bed',asRangedData=F,seqinfo=si)
)
#define window sizes
window.sizes<-c('w25'=25,'w50'=50,'w100'=100,'w200'=200,'w300'=300,'w500'=500)
window.sizes<-window.sizes[c(1,6)]
#exclude the reseq libraries.
cg<-cage.tag.rles[!grepl('reseq',names(cage.tag.rles))]

#rm(cage.tag.rles)

window.sum.distributions<-sapply(simplify=F,regs,function(reg){
  sapply(simplify=F,window.sizes,function(w){
    # sapply(simplify=F,c(10,100,500),function(w){
    #make sure all our regions are at least the window size
    reg<-resize(reg,width=pmax(width(reg),w),fix='center')
    #make sure th
    #convert them to a Rangelist for Views
    winds<-as(reg,'RangesList')[chrs.keep]
    winds<-winds[bigchrs]
    starts<-unlist(start(winds))
    chrs<-names(starts)
    regnums<-1:length(as(winds,'GRanges'))
    
    #for each line  
    maxwinds<-mclapply(mc.cores=10,cg,function(acc){
      #now get the windowed views for each window
      v<-unlist(viewApply(Views(acc[['pos']][bigchrs],winds[bigchrs]),FUN=function(x){ (runsum(x,w)) }))
      v.n<-unlist(viewApply(Views(acc[['neg']][bigchrs],winds[bigchrs]),FUN=function(x){ (runsum(x,w)) }))
      #get the bigger of the two strands
      ispos<-(sapply(v,max)>sapply(v.n,max))
      v<-ifelse( ispos,v,v.n)
      #we can't just take the best window as that will give us only the leftmost.
      #return the middlevalue where the window can take multiple positions
      v.m<-sapply(v,function(x){
        m<-max(x)#find the best position
        
        if(m==0){return(c(NA))} #return NA if all are zero
        z<-which(x==m)#get positions with the best score
        return(
          ceiling((z[1]+rev(z)[1])/2))#return the middle such position
      })
      v.m<-starts+v.m-1##make this the chromosomal coordinate
      #      GRanges(chrs,IRanges(v.m,width=w),score=sapply(v,max))
      data.frame(pos=v.m,score=sapply(v,max),regnum=regnums) #output a data frame
    })
    
    
    cat('.')
    names(maxwinds)<-names(cg)
    maxwinds
  })
})


#melt into one huge dataframe
m<-melt(window.sum.distributions)

#m2<-melt(window.sum.distributions,level=2)
colnames(m)<-c('variable','value','acc','windowsize','region')
unique(m$variable)
unique(m$acc)
unique(m$windowsize)
unique(m$region)


z<-m[ m$variable=='score',]
z$pos<-m[ m$variable=='pos',]$value
z$regnum<-m[ m$variable=='regnum',]$value

unique(z$regnum)

#we can actually do the above by melting with id.vars as c('regnum') and then casting with ... ~ variable
#which is pointless now that I"ve done it anyway....


#library(ggplot2)
#caste(data=m,formula= ... ~ variable)
stopifnot(nrow(m[ m$variable=='pos',])==nrow(m[ m$variable=='score',]))
stopifnot(nrow(m[ m$variable=='pos',])==nrow(z))

#add a column of normalized values
v<-accession.df$intergenic.lib.size
names(v)<-accession.df$accession
z$val.libnorm<-z$value/v[ z$acc ]
#windowsize should be a numeric variable
z$windowsize<-as.numeric(substr(z$windowsize,2,100))

save(z,file='windowtable.robject.R')



