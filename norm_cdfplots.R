setwd('~/Harnett/TSS_CAGE_myfolder')

#Can't install CAGEr
#so use VGAM  and do the normalization ourselves
library(VGAM,lib.loc='~/Harnett/R')
library(igraph,lib.loc='~/Harnett/R')
library(ggplot2)
source('src/tss_cage_functions.R')
load('data/objects/all.cage.unprocessed.object.R')
load('data/objects/cg.libnorm.object.R')

u.cg<-cage.tag.rles
load('data/objects/cg.pl.object.R')
load(file.alltags)
load(file.tss)
accs=names(u.cg)

#create a folder for the results of this script
dir.create('analysis/Norm_comparison')

#get library sizes
sizefactor.df<-data.frame(acc=names(u.cg))
sizefactor.df$library.size<-sapply(u.cg[],function(lib)sum(as.numeric(sum(lib$neg))+as.numeric(sum(lib$pos))))



#get locations of top 1000 most highly tagged 500bp windows
rs=runsum(alltags$both,500)#running sum over windowsize500
topsites=slice(rs,lower=1000)
topsites<-as(topsites,'GRanges')
topsites
topval.lim=sort(unlist(runValue(rs)),dec=T)[1000]
topssites=as(rs>topval.lim,'GRanges')

#get cage for our tss and calculate the size factors using the top, middle and bottom quantiles.
tss.gr=tss.gr[!duplicated(tss.gr)]
cagesum<-unlist(viewSums(Views(alltags$both,as(tss.gr,'RangesList'))[unique(seqnames(tss.gr))]))
cagerank=order(cagesum,decreasing=T)#so rank 1 is the highest tss

#split tss into quantiles
topquants<-tss.gr[cagerank>quantile(cagerank,0.75)]
middlequants<-tss.gr[cagerank>quantile(cagerank,0.25) & cagerank<quantile(cagerank,0.75)]
bottomquants<-tss.gr[cagerank<quantile(cagerank,0.25)]

sizefactors.from.granges<-function(rle,gr){
	s<-sapply(rle,function(acc)unlist(viewSums(GRViews(acc,gr))))
  	sizefactors<-deseqNormalize(s,retur=T)#use DESeq method to get size factors
}

# sizefactor.df$topsites<-sizefactors.from.granges(sapply(cage.tag.rles,'[[',3),topssites)
sizefactor.df$topquant<-sizefactors.from.granges(sapply(u.cg,'[[',3),sort(topquants))
sizefactor.df$middlequants<-sizefactors.from.granges(sapply(u.cg,'[[',3),sort(middlequants))
sizefactor.df$bottomquants<-sizefactors.from.granges(sapply(u.cg,'[[',3),sort(bottomquants))

#Now plot the size factors
jpeg(h=600,w=960,'analysis/Norm_comparison/top_middle_sizefactors.jpeg')
qplot(x=sizefactor.df$topquant,sizefactor.df$middlequants,cex.main=8,main='DESeq sizefactors calculated using top vs. middle quantile CAGE signal TSS',cex.main=3,cex.lab=3)
dev.off()
jpeg(h=600,w=960,'analysis/Norm_comparison/bottom_middle_sizefactors.jpeg')
qplot(x=sizefactor.df$bottomquants,sizefactor.df$middlequants,cex.main=8,main='DESeq sizefactors calculated using bottom vs. middle quantile CAGE signal TSS',cex.lab=3)
dev.off()
jpeg(h=600,w=960,'analysis/Norm_comparison/libsize_middle_sizefactors.jpeg')
qplot(x=sizefactor.df$library.size,sizefactor.df$middlequants,cex.main=8,main='DESeq sizefactors calculated using library size vs. middle quantile CAGE signal TSS',cex.lab=3)
dev.off()


saccs=accs[c(sample(grep('68h',x=accs),10),sample(grep('24h',x=accs),10),sample(grep('1012h',x=accs),10))]
tpvect=c(rep('6-8hours',10),rep('2-4hours',10),rep('10-12hours',10))

return.rcdf.df<-function(srlpn,acc){# sitecounts<-unlist(viewSums(Views.gr(cage.tag.rles[[acc]][[3]],bin_ranges)))
	d=sort(c(unlist(unname(srlpn[['pos']])),unlist(unname(srlpn[['neg']]))))
	d=d[d!=0]
 	t=table(d)#the table, not counting zeroes
  	t=cumsum(t)
  	tot=t[length(t)]
  	t[]=c(0,t[-length(t)])
  	freq=tot-t
  	level=names(freq)
	cat('.')
	data.frame(freq,level=as.numeric(level),acc)
}

cg.de=mapply(SIMPLIFY=F,u.cg,sizefactor.df$middlequants,FUN=function(acc,sf){sapply(acc,function(x)x*sf)})

#RCDF graphs for non normalized and normalized data
#now let's do a graph with some of the library's and our reference
cumsum.df.u <-do.call(rbind,mapply(SIMPLIFY=F,u.cg[],saccs[],FUN=return.rcdf.df))
#Now for power law normalized data
cumsum.df.pl <-do.call(rbind,mapply(SIMPLIFY=F,cg.pl[],saccs[],FUN=return.rcdf.df))
#Now for the library size normalized data normalized data
cumsum.df.libnorm <-do.call(rbind,mapply(SIMPLIFY=F,cg.libnorm[],saccs[],FUN=return.rcdf.df))
#
cumsum.df.denorm <-do.call(rbind,mapply(SIMPLIFY=F,cg.de[],saccs[],FUN=return.rcdf.df))

dir.create('analysis/norm_cdfplots')

pdf('analysis/norm_cdfplots/cdfunp.pdf')
qplot(data=cumsum.df.u,x=level,y=freq,group=acc,log='xy',geom='line')
dev.off()

pdf('analysis/norm_cdfplots/cdfpl.pdf')
qplot(data=cumsum.df.pl<,x=level,y=freq,group=acc,log='xy',geom='line')
dev.off()

pdf('analysis/norm_cdfplots/cdflibnorm.pdf')
qplot(data=cumsum.df.libnorm<,x=level,y=freq,group=acc,log='xy',geom='line')
dev.off()

pdf('analysis/norm_cdfplots/cdfdeseq.pdf')
qplot(data=cumsum.df.denorm<,x=level,y=freq,group=acc,log='xy',geom='line')
dev.off()

object.size(cumsum.df.u/1000000)
