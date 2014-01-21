###This is a script to generate RLE objects out of the bam files for our CAGE data

# @ author Dermot Harnett, EMBL Heidelberg
# @date 18/10/2013
# @title Variance scatterplots
########################################
setwd(dir='/g/furlong/Harnett/TSS_CAGE_myfolder/')
source('src/tss_cage_functions.R')
library(VGAM)
load('data/objects/accession.df.full.object.R')
load('data/objects/tagseq.df.object.R')
load('data/objects/transcripts.object.R')
#cage and tagseq data
load('data/objects/all.cage.unprocessed.object.R')
cg=cage.tag.rles
rm(cage.tag.rles)
load('data/objects/all.tagseq.unprocessed.R')
ts=ts
rm(ts)
#tss annotation
load(file.tss)
#misc parameters
iterations=1000#number of iterations when doing empirical pvalues
#names of timepoints
tps=c('tp24h','tp68h','tp1012h')
accs=names(cg)#names of libraries
taccs=names(ts)#names of libraries
tss.gr=resize(sort(tss.gr),width=500,fix='center')

######2 prepare lists of matching accessions for each timeoint, so we compare only like accessions
accs=names(cg)
taccs=names(ts)
tmat=do.call(rbind,sapply(seq_along(accs),function(i){#for each cage library
  i=match(accs[i],accession.df$accession)
  #find the matching taglibrarys
  matches=which( tagseq.df$line==accession.df$line[i] & tagseq.df$timepoint==accession.df$timepoint[i] )
  if(length(matches)==0){return(NULL)}
  matrix(c(rep(accs[i],length(matches)),taccs[matches]),ncol=2)
}))#get a matrix a column of cage indices and a column of tagseq indicies

#use only these accessions - note duplications of names so we compare all the compbinations
#names of libraries
accs=tmat[,1]
taccs=tmat[,2]
accs=list(
  tp24 = accs[grepl('4h',accs)],
  tp68= accs[grepl('8h',accs)],
  tp1012= accs[grepl('12h',accs)]
  )
names(accs) = tps
#tagseq names
taccs =list(
  tp24 = taccs[grepl('4h',tmat[,1])],
  tp68= taccs[grepl('8h',tmat[,1])],
  tp1012= taccs[grepl('12h',tmat[,1])]
  )
names(taccs) = tps
# length(taccs[[2]])
# length(accs[[2]])
# taccs[[1]]
# accs[[1]]

######3 Load up the peaks for the 3' data as well.
tagseq.peaks=list(
  tp24=read.delim('/g/furlong/Harnett/TSS_CAGE_myfolder/data/Tagseq_peaks/peaks.polya.2h.default.gff.features.gff',header=F,comment.char='#'),
  tp68=read.delim('/g/furlong/Harnett/TSS_CAGE_myfolder/data/Tagseq_peaks/peaks.polya.6h.default.gff.features.gff',header=F,comment.char='#'),
  tp1012=read.delim('/g/furlong/Harnett/TSS_CAGE_myfolder/data/Tagseq_peaks/peaks.polya.10h.default.gff.features.gff',header=F,comment.char='#')
)
tagseq.peaks=sapply(tagseq.peaks,function(f){
  GRanges(paste0('chr',f$V1),IRanges(f$V4,f$V5),strand=f$V7,seqinfo=si)
})
names(tagseq.peaks) = tps
#create 500bp windows around the end of trnascripts for finding our tagseq peaks
ends.gr=resize(transcripts.gr,width=500,fix='end')
ends.gr=shift(ends.gr,250)



load(file.tss)

###### 4 Now get our matrices of counts
#Now go through each timepoint
alltssmat=list()
endmat=list()
for (tp in tps){
  #relevant rles to use for this tp
  tp.accs=accs[[tp]]
  #now get a matrix of strand specific values
  alltssmat[[tp]]=simplify2array(mclapply(mc.cores=10,unique(tp.accs),function(acc){
    stopifnot()
    v=rep(0,length(tss.gr))
    v[as.vector(strand(tss.gr)=='-')]=unlist(viewSums(GRViews(cg[[acc]]$neg,tss.gr)))
    v[as.vector(strand(tss.gr)=='+')]=unlist(viewSums(GRViews(cg[[acc]]$pos,tss.gr)))
    v
  }))
  colnames(alltssmat[[tp]])=unique(tp.accs)
}
dim(alltssmat[[1]])
#end
for (tp in tps){
   #select tss with only 1 peak at each tp
  uends = sort(tagseq.peaks[[tp]])
  tp.taccs=taccs[[tp]]
  #
  endmat[[tp]]=simplify2array(mclapply(mc.cores=10,unique(tp.taccs),function(acc){
    v=rep(0,length(uends))
    v[as.vector(strand(uends)=='-')]=unlist(viewSums(GRViews(ts[[acc]]$neg,uends[strand(uends)=='-'])))
    v[as.vector(strand(uends)=='+')]=unlist(viewSums(GRViews(ts[[acc]]$pos,uends[strand(uends)=='+'])))
    v
  }))
  colnames(endmat[[tp]])=unique(tp.taccs)
  #normalize our tagseq with just the library size for now
  gcov=with(tagseq.df,genome.coverage[accession=taccs])
  endmat[[tp]]=sweep(endmat[[tp]],MARGIN=2,STAT=gcov,FUN='*')
  #fix the names
  rownames(endmat[[tp]])=utss$TrID
}




#check if individual TSS have normal, negative binomial, or log(normally) distributed data
dummynormaldata=rnorm(100,1000,sd=1000)
dummyexpdata=rnorm(100,1,1)*
dummynegbinomdata=rnorm


#one way of doing this would be to fit distributions to our data, then do qqplots...
#likelihood function for normally distributed data
library(MASS)



normlik = function(d,mu,sig){-sum(log(dnorm(d,mu,sig)))}
nbinomlik = function(d,mu,sig){-sum(log(dnbinom2(d,mu,sig)))}
multnoiselik = function(d,mu,sig){  -sum(log(   dnorm(log(d),mu,sig)      )   )}

mle(minuslogl, start = formals(minuslogl), method = "BFGS",
    fixed = list(), nobs, ...)

#test individual row
ind=order(rowSums(alltssmat[[2]]),decreasing=T)[1]
test=alltssmat[[2]][ind , ]
# hist(test)
logfit=fitdistr(densfun='lognormal',x=test)
normfit=fitdistr(densfun='normal',x=test)
nbinomfit=fitdistr(densfun='negative binomial',x=test)
logfit$loglik-nbinomfit$loglik

#test test 
logliks=sapply(order(rowSums(alltssmat[[2]]),decreasing=T)[1:1000],     function(ind){
  test=alltssmat[[2]][ind , ]
    if(any(test==0)){return(NA)}
  logfit=fitdistr(densfun='lognormal',x=test)
  normfit=fitdistr(densfun='normal',x=test)
  logfit$loglik-nbinomfit$loglik
})
#test test but with our tagseq data
logliks.tagseq=sapply(order(rowSums(endmat[[1]]),decreasing=T)[1:1000],  function(ind){
  test=endmat[[1]][ind , ]
    if(any(test==0)){return(NA)}
  logfit=fitdistr(densfun='lognormal',x=test)
  normfit=fitdistr(densfun='normal',x=test)
  logfit$loglik-nbinomfit$loglik
})
logliks.tagseq

library(reshape2)
library(ggplot2)
pdf('Bfactorboxplots.pdf')
qplot(data= melt(list(logliks,logliks.tagseq)) ,y=value,x=L1,color=L1,geom='boxplot')
dev.off()
hist(test)


#Are our resequencing variants within poisson distance of eachother?

#define pairs of resequencing replicates

#for each of these pairs pull out the tss values

#plot the the theoretical MAs and the real MAs

