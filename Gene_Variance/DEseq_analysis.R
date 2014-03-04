#A Load libraries,functions,data
library(reshape2)
library(ggplot2)
setwd ( '/g/furlong/Harnett/TSS_CAGE_myfolder/' )
source ( 'src/tss_cage_functions.R' )
#load cage data and tagseq data
load(file.tss)
load(file.synonyms)
load('data/objects/accession.df.object.R')
load('data/objects/cg.object.R')
load('data/objects/cg.pl.object.R')

#timepoints
outfolder='analysis/gene_variance_all'
tps=c('tp24h','tp68h','tp1012h')
dir.create(outfolder)

####F Load info on gene function
#First we load the FBgn : Go file
fb_go.df <- read.delim(header=F,sep=' ','data/go_gene_association.short.txt',stringsAsFactors=F)
fb_go.df <- fb_go.df[,c(1,3,5)]
head(fb_go.df)
colnames(fb_go.df) <- c('FBgn','symb','GO')
names(fbgn2symb)=fbgn2symb=unique(fb_go.df$FBgn)
fbgn2symb=fb_go.df$symb[match(fbgn2symb,fb_go.df$FBgn)]
#now load the table linking fbgns and fbtrs
fbgn_fbtr_fbpp.df <- read.delim(header=F, comment.char='#',stringsAsFactors=F,sep='\t','/g/furlong/Harnett/TSS_CAGE_myfolder/data/fbgn_fbtr_fbpp_fb_2013_06.tsv.gz' )
colnames(fbgn_fbtr_fbpp.df) <- c('FBgn','FBtr','FBpp')
#Also load table of FBgns and gene info
#add
tss.gr$unp.cage = sapply(names(cg),function(acc){
	unlist(viewSums(GRViews(cg[[acc]]$both,resize(tss.gr,width=500,fix='center'))))
})
#add
tss.gr$tsscage = sapply(names(cg.pl),function(acc){
	unlist(viewSums(GRViews(cg.pl[[acc]]$both,resize(tss.gr,width=500,fix='center'))))
})
save(tss.gr)

#let's do this with dummy data first
#vector of gene signals


	


libsizes=libsizes/mean(libsizes)
for(n in seq_along(libsizes)){
	simcountmat[,n]=rpois(n=gn,lambda=exp(simcountmat[,n])*libsizes[n])
}








dobs=lapply(1:4,function(n){#estimate dispersions amongst
	reptypes=samplegroups[1:(n)]
	de.simcount=newCountDataSet(simcountmat,reptable[,reptypes])
	de.simcount=estimateSizeFactors(de.simcount)
	tmp2=estimateDispersions(de.simcount,'pooled',fitType='local')
})

#Now make a graph plotting our variance/mean for the different 
pdf('analysis/gene_variance_all/dispersionplots_simdata.pdf')

	xg = 10^seq(-0.5, 5, length.out = 100)
	maxdisp=max(unlist(sapply(dobs,function(x){max(fData(x)[[1]])})))

	plot(xg,seq(1,maxdisp,length.out=100),type='n',log='xy',ylim=c(0.001,10),main='Dispersion vs mean simulated data')

	linecols=c('red','green','blue','orange')
	for(n in seq_along(dobs)){
		linecol=linecols[ n]
		dob=dobs[[n]]
	    lines(xg, fitInfo(dob, name = 'pooled')$dispFun(xg), col = linecol, 
	        lwd = 4)
	}
	legend(x='bottomright',legend=c('Temporal','Genetic','Biological','Technical'),fill=linecols)

dev.off()


#############Now with the actual data

countmat=tss.gr[!duplicated(tss.gr)]$unp.cage
row.names(countmat)=tss.gr[!duplicated(tss.gr)]$TrID
samplegroups=c('timepoint','line','collection','prep','seq')
reptable=accession.df[accession.df$tissue=='embryo',]
reptable=reptable[,c('acc',samplegroups)]
countmat=countmat[,reptable$acc]
dim(countmat)
dim(reptable)
any(is.na(countmat))
any(countmat==0)

n=4

dobs=lapply(1:4,function(n){#estimate dispersions amongst
	reptypes=samplegroups[1:(n)]
	de.count=newCountDataSet(countmat,reptable[,reptypes])
	de.count=estimateSizeFactors(de.count)
	tmp2=estimateDispersions(de.count,'pooled',fitType='parametric')
})


#Now make a graph plotting our variance/mean for the different 
pdf('analysis/gene_variance_all/dispersionplots_actualdata.pdf')

	xg = 10^seq(-0.5, 5, length.out = 100)
	maxdisp=max(unlist(sapply(dobs,function(x){max(fData(x)[[1]])})))

	plot(xg,seq(1,maxdisp,length.out=100),type='n',log='xy',ylim=c(0.001,10),main='Dispersion vs mean simulated data')

	linecols=c('red','green','blue','orange')
	for(n in seq_along(dobs)){
		linecol=linecols[ n]
		dob=dobs[[n]]
	    lines(xg, fitInfo(dob, name = 'pooled')$dispFun(xg), col = linecol, 
	        lwd = 4)
	}
	legend(x='bottomright',legend=c('Temporal','Genetic','Biological','Technical'),fill=linecols)

dev.off()



#############Now with the actual data, and our normalization

countmat=tss.gr[!duplicated(tss.gr)]$tsscage
row.names(countmat)=tss.gr[!duplicated(tss.gr)]$TrID
samplegroups=c('timepoint','line','collection','prep','seq')
reptable=accession.df[accession.df$tissue=='embryo',]
reptable=reptable[,c('acc',samplegroups)]
countmat=countmat[,reptable$acc]
dim(countmat)
dim(reptable)
any(is.na(countmat))
any(countmat==0)
countmat=floor(countmat)

n=4
dobs=lapply(1:4,function(n){#estimate dispersions amongst
	reptypes=samplegroups[1:(n)]
	de.count=newCountDataSet(countmat,reptable[,reptypes])
#	de.count=estimateSizeFactors(de.count)
	tmp2=estimateDispersions(de.count,'pooled',fitType='parametric')
})


#Now make a graph plotting our variance/mean for the different 
pdf('analysis/gene_variance_all/dispersionplots_actualdata.pl.pdf')

	xg = 10^seq(-0.5, 5, length.out = 100)
	maxdisp=max(unlist(sapply(dobs,function(x){max(fData(x)[[1]])})))

	plot(xg,seq(1,maxdisp,length.out=100),type='n',log='xy',ylim=c(0.001,10),main='Dispersion vs mean actual data PLnorm')

	linecols=c('red','green','blue','orange')
	for(n in seq_along(dobs)){
		linecol=linecols[ n]
		dob=dobs[[n]]
	    lines(xg, fitInfo(dob, name = 'pooled')$dispFun(xg), col = linecol, 
	        lwd = 4)
	}
	legend(x='bottomright',legend=c('Temporal','Genetic','Biological','Technical'),fill=linecols)

dev.off()


max.finite<-function(x){max(x[is.finite(x)])}
tps=unique(reptable$timepoint)[c(2,1,3)]
#now for the different timepoints
tpdobs=lapply(tps,function(tp){#estimate dispersions amongst
	
	cmat= countmat[,reptable$timepoint==tp]
	de.count=newCountDataSet( cmat,reptable[reptable$timepoint==tp,] )

	#check what happens if we use only sample with a single collection
	tab=table(reptable$collection)
	collects2use=names(tab)[tab==1]
	cmat= countmat[,reptable$timepoint==tp & reptable$collection%in%collects2use]
	de.count=newCountDataSet( cmat,reptable[reptable$timepoint==tp & reptable$collection%in%collects2use,] )

	de.count=estimateSizeFactors(de.count)
	tmp2=estimateDispersions(de.count,'blind',fitType='parametric')
})

#Now make a graph plotting our variance/mean for the different 
pdf('analysis/gene_variance_all/dispersionplots_byTP_actualdata.pl.local.pdf')
	xg = 10^seq(-0.5, 6, length.out = 100)
	maxdisp=max(unlist(sapply(tpdobs,function(x){max.finite(fData(x)[[1]])})))
	plot(xg,seq(1,maxdisp,length.out=100),type='n',log='xy',ylim=c(0.1,1),
		main='Dispersion vs mean actual data PLnorm - 3 timepoints')

	linecols=c('red','green','blue')
	for(n in seq_along(tpdobs)){
		linecol=linecols[n]
		tpdob=tpdobs[[n]]
		px = rowMeans(counts(tpdob, normalized = TRUE))
		sel = (px > 0)
		px = px[sel]
	    py = fData(tpdob)[[1]][sel]	
    	ymin = 10^floor(log10(min(py[py > 0], na.rm = TRUE)) - 0.1)
    	
    	 points(px,pmax(py,ymin),col=linecol)
		
	    lines(xg, fitInfo(tpdob)$dispFun(xg), col = linecol, lwd = 1)
	}
	legend(x='bottomright',legend=tps,fill=linecols)

dev.off()

