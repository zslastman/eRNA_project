setwd('~/Harnett/TSS_CAGE_myfolder')

#Can't install CAGEr
#so use VGAM  and do the normalization ourselves
library(VGAM,lib.loc='~/Harnett/R')
library(igraph,lib.loc='~/Harnett/R')
library(ggplot2)
source('src/tss_cage_functions.R')
load('data/objects/all.cage.unprocessed.object.R')
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




#RCDF graphs for non normalized and normalized data
#now let's do a graph with some of the library's and our reference
cumsum.df.u<-do.call(rbind,lapply(accs,function(acc){
	#get the summed reads in each genomic bin for this library
	# sitecounts<-unlist(viewSums(Views.gr(cage.tag.rles[[acc]][[3]],bin_ranges)))
	d=sort(c(unlist(unname(u.cg[[acc]][['pos']])),unlist(unname(u.cg[[acc]][['neg']]))))
	d=d[d!=0]
 	t=table(d)#the table, not counting zeroes
  	t=cumsum(t)
  	tot=t[length(t)]
  	t[]=c(0,t[-length(t)])
  	freq=tot-t
	cat('.')
	data.frame(freq,level=as.numeric(level),acc)
}))
#Now for power law normalized data
cumsum.df.pl<-do.call(rbind,lapply(accs,function(acc){
	#get the summed reads in each genomic bin for this library
	# sitecounts<-unlist(viewSums(Views.gr(cage.tag.rles[[acc]][[3]],bin_ranges)))
	d=sort(c(unlist(unname(cg.pl[[acc]][['pos']])),unlist(unname(cg.pl[[acc]][['neg']]))))
	d=d[d!=0]
 	t=table(d)#the table, not counting zeroes
  	t=cumsum(t)
  	tot=t[length(t)]
  	t[]=c(0,t[-length(t)])
  	freq=tot-t
	cat('.')
	data.frame(freq,level=as.numeric(level),acc)
}))
#Now for the DESeq normalized data





rle2rcdf<-function(d){
 #takes in a list of two srles and outputs the reverse cdf vector
  d=do.call('c',as.list(c(unname(d[[1]]),unname(d[[2]]))))
  
}

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

#function to get slope of fitted alphas 
fits<-sapply(names(u.cg)[],function(acc){
	sitecounts=sort(c(unlist(unname(u.cg[[acc]][['pos']])),unlist(unname(u.cg[[acc]][['neg']]))))
	sitecounts=sitecounts[sitecounts!=0]
	fit=power.law.fit(as.vector(sitecounts),xmin=200)#not that power.law.fit excludes the lower count numbers
	total.tags=sum(sitecounts)
	o=getOffset(fit$alpha,total.tags)#calculate the offset (second parameter)
	c(fit,offset=o,tagcount=total.tags)#now tack it onto the fit list and return
	# plfit1=power.law.fit(sitecounts)#not that power.law.fit excludes the lower count numbers
	# cat('.')
	# plfit1$alpha
})

#average these
mean.alpha=mean(unlist(fits['alpha',]))
r.offset<-getOffset(mean.alpha,T=1000000)

offsets=as.numeric(unlist(fits['offset',]))

#now let's plot the alphas vs. library size
qplot(x=offsets,y=unlist(fits['alpha',]),geom='point',color='red',
	main="Power law slope vs Library Size - Cage Libraries 6-8h",
	xlab='Offset',
	ylab='Power Law Slope')


jpeg(w=1920,h=1200,paste0('analysis/Norm_comparison','rev_loglog_cumulative_dens_prenorm','.jpeg'))
#Do a reverse cumulative plot, and add in a line for the reference distribution
qplot(data=cumsum.df,color=acc,y=freq,x=as.numeric(level),log='xy',geom='line',add=F,xlab='Tag count t',ylab='Num TSS >= t')	
dev.off()
#+ geom_abline(intercept=1000000,slope=-mean.alpha,)

reverse.cdf.plot<-function(x,add=F){
	freq<-cumsum(table(x))
	freq<-max(freq)-freq
	level<-as.numeric(names(freq))
	if(add){
	points(x=level,y=freq,type='l')
	}else{
	plot(type='l',main='reverse cdf plot',xlab='magnitude',ylab='reverse cumulative density',x=level,y=freq,log='xy')
	}
}
#some dummy power law data
reverse.cdf.plot(rpareto(100000,217000,1.25),add=F)
#does our 

#define our normalization function
pl.norm<-function(x,x.offset,x.alpha,ref.offset=r.offset,ref.alpha=mean.alpha){
	#take in values x and normalize them to a power law distribution of offset and slope given.
	#again, see  http://genomebiology.com/content/10/7/R79 for the math

	beta=x.alpha/ref.alpha
	lambda = (ref.offset/x.offset)^ (1/ref.alpha)
	x = x ^ beta
	x = x * lambda#do it like this to keep names
	return(  x)

}

alphas=as.numeric(fits['alpha',])
#now finally normalize each of our cage libraries
cage.tag.rles.plnorm<-mapply(SIMPLIFY=F,names(u.cg),alphas,offsets,FUN=function(acc,alpha,offset){
	 lapply(u.cg[[acc]],function(srle){ 
	 	pl.norm(srle,x.alpha=alpha,x.offset=offset[[1]])
	 })
})

#and save it
save(cage.tag.rles.plnorm,file=file.cage.tag.rles.plnorm)



#now plot our cdfs again


#
#Now let's do the cdf plots with the power law normalized libraries
cumsum.df<-do.call(rbind,lapply(names(cage.tag.rles.plnorm),function(acc){
	#get the summed reads in each genomic bin for this library
	# sitecounts<-unlist(viewSums(Views.gr(cage.tag.rles[[acc]][[3]],bin_ranges)))
	sitecounts=sort(c(unlist(unname(cage.tag.rles.plnorm[[acc]][['pos']])),unlist(unname(cage.tag.rles.plnorm[[acc]][['neg']]))))
	freq<-cumsum(table(sitecounts))
	freq<-max(freq)-freq
	level<-names(freq)
	cat('.')
	data.frame(freq,level=as.numeric(level),acc)

}))

jpeg(w=1920,h=1200,paste0('analysis/Norm_comparison','rev_loglog_cumulative_dens_prenorm','.jpeg'))
qplot(data=cumsum.df,color=acc,y=freq,x=as.numeric(level),log='xy',geom='line',add=F) + geom_abline(intercept=1000000,slope=-mean.alpha)
dev.off()

load(file.cage.tag.rles.rpgc)


#Now let's do the cdf plots with the RPGC normalized libraries
cumsum.df<-do.call(rbind,lapply(names(cage.tag.rles.rpgc),function(acc){
	#get the summed reads in each genomic bin for this library
	# sitecounts<-unlist(viewSums(Views.gr(cage.tag.rles[[acc]][[3]],bin_ranges)))
	sitecounts=sort(c(unlist(unname(cage.tag.rles.rpgc[[acc]][['pos']])),unlist(unname(cage.tag.rles.rpgc[[acc]][['neg']]))))
	freq<-cumsum(table(sitecounts))
	freq<-max(freq)-freq
	level<-names(freq)
	cat('.')
	data.frame(freq,level=as.numeric(level),acc)

}))
qplot(data=cumsum.df,color=acc,y=freq,x=as.numeric(level),log='xy',geom='line',add=F) + geom_abline(intercept=1000000,slope=-mean.alpha)



#Now let's do the cdf plots with the DESEQ normalized libraries
cumsum.df<-do.call(rbind,mapply(SIMPLIFY=F,names(u.cg)[],sizefactor.df$middlequants,FUN=function(acc,sizefactor){
	#get the summed reads in each genomic bin for this library
	# sitecounts<-unlist(viewSums(Views.gr(cage.tag.rles[[acc]][[3]],bin_ranges)))
	sitecounts=sort(c(unlist(unname(u.cg[[acc]][['pos']]*sizefactor)),unlist(unname(u.cg[[acc]][['neg']]*sizefactor))))
	freq<-cumsum(table(sitecounts))
	freq<-max(freq)-freq
	level<-names(freq)
	cat('.')
	data.frame(freq,level=as.numeric(level),acc)

}))
cumsum.df=data.frame(acc=cumsum.df[,'acc'])
qplot(data=cumsum.df,color=acc,y=freq,x=as.numeric(level),log='xy',geom='line',add=F) + geom_abline(intercept=1000000,slope=-mean.alpha)







#Show the effect this has on diffferent magnitudes of tss
load(file.tss)
load(file.crmgrs)
load(file.cad3)
source('src/tss_cage_functions.R')
w=100

#make list of grs
active.tss.gr<-tss.gr[ tss.gr$active68]
active.tss.gr<- active.tss.gr[ seqnames(active.tss.gr)%in%bigchrs]
gr=combinegrs(list(cad3=cad3.gr,crm=crmgrs,tss=active.tss.gr))
gr=sort(gr)

#get the matrix with info for each line on each of our grs
gr$cagemat=get.best.window.mat(reg=gr,w=300,chrs.keep,bigchrs,cage=cg[accs])

#is ther ea differential contribution from different libraries causing the problem?
rs=rowSums(gr$cagemat)
ncagemat=gr$cagemat/rs

qplot(data=melt(ncagemat),x=Var2,y=value,geom='boxplot',log='y',color=Var2,xlab=F,
                        ylab='Mean % of total reads at all CRMs/TSS',main='Average Contribution of normalized libraries to TSS/CRM \n power law Normalization')


#Test the success of the normalization methods by the boxplot showing contribution to each tss
#And also the pvalue normalization
