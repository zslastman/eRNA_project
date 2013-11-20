setwd ( '/g/furlong/Harnett/TSS_CAGE_myfolder/' )
source ( 'src/tss_cage_functions.R' )

# load ( 'data/objects/crmgrs.object.R' )
# load ( 'data/objects/cad3.gr.object.R' )



######1 load up rles and metadata on them

#annotation and metadata
load(file.transcripts)
load('data/objects/accession.df.full.object.R')
load('data/objects/tagseq.df.object.R')
#cage and tagseq data
load('data/objects/all.cage.unprocessed.object.R')
cg=cage.tag.rles
rm(cage.tag.rles)
load('data/objects/all.tagseq.unprocessed.R')
ts=tagseq.rles
rm(tagseq.rles)
iterations=1000#number of iterations when doing empirical pvalues
#names of timepoints
tps=c('tp24h','tp68h','tp1012h')


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


###### 5 Now get our matrices of counts
#Now go through each timepoint
tssmat=list()
endmat=list()
for (tp in tps){
	#select tss with only 1 peak at each tp
	fov = findOverlaps(ends.gr,tagseq.peaks[[tp]])
	fov = fov[!duplicated(fov@queryHits) & !duplicated(fov@subjectHits),]
	#get the tss
	utss=resize(transcripts.gr[ fov@queryHits ],width=1,fix='start')
	utss=sort(resize(utss,width=500,fix='center'))
	#and en of transcripts from the peak data
	uends = sort(tagseq.peaks[[tp]][ fov@subjectHits ])
	#relevant rles to use for this tp
	tp.accs=accs[[tp]]
	tp.taccs=taccs[[tp]]
	#now get a matrix of strand specific values
	tssmat[[tp]]=simplify2array(mclapply(mc.cores=10,unique(tp.accs),function(acc){
		v=rep(0,length(utss))
		v[as.vector(strand(utss)=='-')]=unlist(viewSums(GRViews(cg[[acc]]$neg,utss[strand(utss)=='-'])))
		v[as.vector(strand(utss)=='+')]=unlist(viewSums(GRViews(cg[[acc]]$pos,utss[strand(utss)=='+'])))
		v
	}))
	colnames(tssmat[[tp]])=unique(tp.accs)
	#and for the ends
	endmat[[tp]]=simplify2array(mclapply(mc.cores=10,unique(tp.taccs),function(acc){
		v=rep(0,length(uends))
		v[as.vector(strand(uends)=='-')]=unlist(viewSums(GRViews(ts[[acc]]$neg,uends[strand(uends)=='-'])))
		v[as.vector(strand(uends)=='+')]=unlist(viewSums(GRViews(ts[[acc]]$pos,uends[strand(uends)=='+'])))
		v
	}))
	colnames(endmat[[tp]])=unique(tp.taccs)
	#fix the names
	rownames(tssmat[[tp]])=rownames(endmat[[tp]])=utss$TrID
	list()
}
dim(tssmat[[1]])
dim(endmat[[2]])
#We now have matrices with the rows representing genes at a given timepoint 



###### 6 Now get correlations of matching libraries

#get the correlation for each tss
gene.cors.tp=sapply(simplify=F,tps,function(tp){
	tp.accs=accs[[tp]]
	tp.taccs=taccs[[tp]]
		sapply(1:nrow(tssmat[[tp]]),function(i){
		#get the correlation over the comparable libraries
		#vector.pvals.shuffle(tssmat[i,accs],endmat[i,taccs])
		suppressWarnings(cor( tssmat[[tp]][i,tp.accs] ,endmat[[tp]][i,tp.taccs] ,method='s'))
	})
})

#and the pval as computed by shuffling - very 
gene.shuffle.pvals=sapply(simplify=F,tps,function(tp){
	tp.accs=accs[[tp]]
	tp.taccs=taccs[[tp]]
		mapply(a= tssmat[[tp]][1:3,tp.accs] , b = endmat[[tp]][1:3,tp.taccs],vector.pvals.shuffle)
		#get the correlation over the comparable libraries
		#vector.pvals.shuffle(tssmat[i,accs],endmat[i,taccs])
})

save('tmp.image.R')
#output the mean of all correlations

library(ggplot2)
library(LSD)
# dir.create('analysis/5primer_3prime_covariance')

for(tp in tps){
	#plot our correlations as a fuction of mean cage
	pdf(paste0('analysis/5primer_3prime_covariance/cage_tagseq_correlation_',tp,'.pdf'))
	heatscatter(x=log2(rowSums(tssmat[[tp]])),y=log2(rowSums(endmat[[tp]])),xlab=paste0('Log2 total cage signal vs. total tagseq signal ',tp),ylab='Log2 cage/tagseq spearman correlation',ncol=100,log='')
	dev.off()
	#plot our correlations as a fuction of mean cage
	pdf(paste0('analysis/5primer_3prime_covariance/correlation_ranked_by_cage_scatter_',tp,'.pdf'))
	heatscatter(x=log2(rowSums(tssmat[[tp]])),y=gene.cors.tp[[tp]][rank(rowSums(tssmat[[tp]]))],xlab='Log2 total cage signal',ylab='cage/tagseq spearman correlation',log='x',ylim=c(-1,1),cor=F)
	dev.off()
	#plot our correlations as a fuction of mean tagseq
	pdf(paste0('analysis/5primer_3prime_covariance/correlation_ranked_by_tagseq_scatter',tp,'.pdf'))
	heatscatter(x=log2(rowSums(endmat[[tp]])),y=gene.cors.tp[[tp]][rank(rowSums(endmat[[tp]]))],xlab='log2 total tagseq signal',ylab='cage/tagseq spearman correlation',log='x',ylim=c(-1,1),cor=F)
	dev.off()
}

#maybe show that different simulated data sets return different normalized methods as best via this method?
# gene.cors=gene.cors[1:length(gene.cors)]


save('tmp.image.R')

