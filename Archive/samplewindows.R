#This script will produce roc curves using resamplings of the original pooled data set
#to see how many reads are necessary to discriminate our positive from our negative set.
setwd('~/Harnett/TSS_CAGE_myfolder/')
rm(list=ls())

source('src/tss_cage_functions.R')
library(reshape2)
library(ggplot2)
library(ROCR,lib.loc='~/Harnett/R')
load(file.alltags)
load(file.alltags.rpgc)

regs=list(
  pos=import(con='analysis/positive_8008.bed',asRangedData=F,seqinfo=si),
  neg=import(con='analysis/negative_8008.bed',asRangedData=F,seqinfo=si),
  #int=import(con='analysis/random.intergenic.bed',asRangedData=F,seqinfo=si),
  cad3pos=import(con='analysis/cad3.pos.bed',asRangedData=F,seqinfo=si),
  cad3neg=import(con='analysis/cad3.neg.bed',asRangedData=F,seqinfo=si)
)

#concatenate the regions with names for identification
regs<-sapply(names(regs),function(regn){
  reg<-regs[[regn]]
  mcols(reg)<-NULL
  reg$regnum<-1:length(reg);
  reg$region<-regn
  reg})
regs<-do.call(c,unname(regs))
regs<-sort(regs)
#define window sizes
window.sizes<-c('w25'=25,'w50'=50,'w100'=100,'w200'=200,'w300'=300,'w500'=500)

#function that takes an RleList and returns a new one with the counts resampled
#according to binomial probabilities 
sampleRleList<-function(srle,p){
  as(
    sapply(names(srle),function(chr){      
      chr<-srle[[chr]]
      sites<-chr>0
      subsample<-rbinom(n=length(as.vector(chr[sites])),prob=p,size=as.vector(chr[sites]))
      chr[sites]<-subsample
      chr
    }),'SimpleRleList')
}


samplesizes<-1/c(1000,100,20,10,5,2,1)
window.sizes=window.sizes
replist<-list(rep1='rep1',rep2='rep2',rep3='rep3')
#regs=regs[1:10]
#define an alltags object by sampling the old one
sample.reg.windowsize<-
sapply(simplify=F,replist,function(sample.replicate){  
  tmp<-mclapply(mc.cores=10,samplesizes,function(sampsize){#loop works
    alltags.s<-sapply(simplify=F,names(alltags),function(strand){
      sampleRleList(alltags[[strand]],sampsize)
    })#this works
    cat('.')
    sapply(simplify=F,window.sizes,function(w){
      #window the whole sampled Rle
      pos.windowed<-runsum(alltags.s[['pos']][bigchrs],w)
      neg.windowed<-runsum(alltags.s[['neg']][bigchrs],w)
      #make sure all our regions are at least the window size
      regs<-resize(regs,width=pmax(width(regs),w),fix='center')
      regs<-resize(regs,width=width(regs)-w+1,fix='start')#consider only windows inside the region
      #convert them to a Rangelist for Views
      winds<-as(regs,'RangesList')[chrs.keep]
      winds<-winds[bigchrs]#only the big chrs
      #get the views on the windowed rle
      v<-Views(pos.windowed,winds[bigchrs])[bigchrs]
      v.n<-Views(neg.windowed,winds[bigchrs])[bigchrs]
      #get the maximum 
      vm<-unlist(viewMaxs(v))
      v.nm<-unlist(viewMaxs(v.n))
      pmax(vm,v.nm)
    })
  })
  names(tmp)<-as.character(samplesizes)
  tmp
})
save(sample.reg.windowsize,file='data/objects/sample.reg.windowsize.object.R')
load(file='data/objects/sample.reg.windowsize.object.R')


rep='rep1'
w='w50'
s<-as.character(samplesizes)[3]



#now we should just be able to produces several roc plots for each sample size
for(rep in unique(names(sample.reg.windowsize))){
  for(w in unique(names(sample.reg.windowsize[[1]][[1]]))){
  #mfrow to a certan number
  par(mfcol=c(4,2))  
  for(sn in 1:length(samplesizes)){
    s<-as.character(samplesizes)[sn]
    allscores<-sample.reg.windowsize[[1]][[s]][[w]]
    scores<-stack(list(  pos=allscores[ regs$region=='pos' ],neg=allscores[ regs$region=='neg' ]))
    sumpred<-prediction(scores$values,scores$ind)
    perf<-performance(sumpred,'tpr','fpr')
    plot(perf,colorize=T,add= F)
    aucstring=paste0('AUC = ',round(performance(sumpred,'auc')@y.values[[1]],4))
    text(x=0.90,y=0,labels=aucstring)
    title(paste0('ROC curve for ',s,' of Summed Unnormalized Tags,8008 set'),sub=paste0('windowsize ',w))
    abline(coef=c(0,1),lty=3)
  }
}
}




    
# let's also produce density plots ----------------------------------------

#now we should just be able to produces several roc plots for each sample size
for(rep in unique(names(sample.reg.windowsize))){
  for(w in unique(names(sample.reg.windowsize[[1]][[1]]))){
    #mfrow to a certan number
    scores<-sapply(simplify=F, 1:length(samplesizes),function(sn){
      s<-as.character(samplesizes)[sn]
      allscores<-sample.reg.windowsize[[1]][[s]][[w]]
      scores<-stack(list(  pos=allscores[ regs$region=='cad3pos' ],neg=allscores[ regs$region=='cad3neg' ]))
      qplot(x=scores$value,color=scores$ind ,  geom='density')+scale_x_log10()
      
      
    })
    names(scores)<-samplesizes
    m<-melt(scores,id.vars=c('ind'))
    colnames(m)[4]<-'samplesize'
    m$samplesize<-factor(m$samplesize)
    head(m)
    qplot(x=m$value,color=ind ,facets=  ~ samplesize,  geom='density')+scale_x_log10()
    
  }
}
       
        




    sampleRleList<-function(srle,p){
  as(
    sapply(names(srle),function(chr){      
      chr<-srle[[chr]]
      sites<-chr>0
      subsample<-rbinom(n=length(as.vector(chr[sites])),prob=p,size=as.vector(chr[sites]))
      chr[sites]<-subsample
      chr
    }),'SimpleRleList')
}


    
#melt into one huge dataframe
# m<-melt(sample.reg.windowsize)
# 
# #m2<-melt(window.sum.distributions,level=2)
# h(m)
# colnames(m)<-c('variable','value','windowsize','region','samplesize')
# m$samplesize<-samplesizes[m$samplesize]
# z<-m[ m$variable=='score',]
# z$pos<-m[ m$variable=='pos',]$value
# z$regnum<-m[ m$variable=='regnum',]$value
# h(z)
# 
# 
# unique(z$regnum)

#now that we have z we can go through by replicate and 

