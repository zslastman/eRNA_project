setwd ( '/g/furlong/Harnett/TSS_CAGE_myfolder/' )
source ( 'src/tss_cage_functions.R' )
rootfolder='/g/furlong/Harnett/TSS_CAGE_myfolder/'
outfolder='/g/furlong/Harnett/TSS_CAGE_myfolder/analysis/norm_comparison_vs_tagseq.R'
dir.create(outfolder,showWarnings=F)
	library(corrperm)

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
ts=ts
rm(ts)
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
filt.transcripts.gr = transcripts.gr[countOverlaps(transcripts.gr,transcripts.gr)==1]
ends.gr=resize(filt.transcripts.gr,width=500,fix='end')
ends.gr=shift(ends.gr,250)


#we'll perform steps 5 and onwards once for each normalization method
#b)sqrt(lib size)
#c)power law
#c.2)power law with the original equation
#d)qnorm((rank(reads) - 0.5)/length(reads), ties.method='average')
#e)qunatile normalization but counting sites with reads in ANY library



#Loop over all of these datasets
cormeanlist=list()
gene.shuffle.pvals=list()
gene.cors.tp=list()
setwd ( '/g/furlong/Harnett/TSS_CAGE_myfolder/' )
for(tagseq.data in list(
	'/g/furlong/Harnett/TSS_CAGE_myfolder/data/objects/ts.object.R',
	'/g/furlong/Harnett/TSS_CAGE_myfolder/data/objects/ts.pl.object.R'
	)){
	for(cage.data in list(
		paste0(rootfolder,'data/objects/cg.libnorm.object.R'),
		#paste0(rootfolder,'data/objects/cg.sqrt.object.R'),
		paste0(rootfolder,'data/objects/cg.pl.map.object.R')
		#paste0(rootfolder,'data/objects/cg.qn.object.R'),
		#paste0(rootfolder,'data/objects/cg.qn.allsites.object.R'),
		#paste0(rootfolder,'data/objects/all.cage.unprocessed.object.R')
	)){  #create a subfolder for each one
		setn= paste0(gsub('.*/(.*?).object.R$',x=cage.data,rep='\\1'),'.',gsub('.*/(.*?).object.R$',x=tagseq.data,rep='\\1'))
		#load the cage dataset,awkward but I forgot to save with RDS objects
		message('loading cage dataset')
		message(setn)

		cgenv<-new.env()
		load(file=cage.data,env=cgenv)
		cg=cgenv[[ls(cgenv)[1]]]
		rm(cgenv)

		tsenv<-new.env()
		load(file=tagseq.data,env=tsenv)
		ts=tsenv[[ls(tsenv)[1]]]
		rm(tsenv)

		message(length(cg))
		message(length(ts))
		#
		message('creating folder')
		setwd(outfolder)
		outsubfolder=paste0(outfolder,'/',setn)
		dir.create(outsubfolder,showWarnings=F)
		setwd(outsubfolder)
		#check our data:
		#list with all accs in ii
		stopifnot(is.list(cg))
		# stopifnot(all(names(cg)%in%unlist(accs)))
		stopifnot(all(sapply(cg,function(g){cat('.');'SimpleRleList' %in% is(g$pos)})))
		message('data checks out')
		# which(!sapply(cg,function(g){cat('.');'SimpleRleList' %in% is(g$both)}))
		# boths=sapply(names(cg),function(acc){cg[[acc]]$pos+cg[[acc]]$neg})
		###### 5 Now get our matrices of counts
		#Now go through each timepoint and just 
		tssmat=list()
		endmat=list()
		for (tp in tps){
			#select tss with only 1 peak at each tp
			fov = findOverlaps(ends.gr,tagseq.peaks[[tp]])
			fov = fov[!duplicated(fov@queryHits) & !duplicated(fov@subjectHits),]
			#get the tss
			utss=resize(filt.transcripts.gr[ fov@queryHits ],width=1,fix='start')
			utss=sort(resize(utss,width=500,fix='center'))
			#and en of transcripts from the peak data
			uends = sort(tagseq.peaks[[tp]][ fov@subjectHits ])
			stopifnot(length(uends)==length(utss))
			#relevant rles to use for this tp
			tp.accs=accs[[tp]]
			tp.taccs=taccs[[tp]]
			tp.accs = tp.accs[ tp.accs %in% names(cg) & tp.taccs %in% names(ts) ]
			tp.taccs = tp.taccs[ tp.accs %in% names(cg) & tp.taccs %in% names(ts) ]
			stopifnot(length(tp.accs)==length(tp.taccs))
			#now get a matrix of strand specific values
			libs=unique(tp.accs)
			tssmat[[tp]] = simplify2array(mclapply(mc.cores=10,libs,function(acc){
				v=rep(0,length(utss))
				v[as.vector(strand(utss)=='-')]=unlist(viewSums(GRViews(cg[[acc]]$neg,utss[strand(utss)=='-'])))
				v[as.vector(strand(utss)=='+')]=unlist(viewSums(GRViews(cg[[acc]]$pos,utss[strand(utss)=='+'])))
				stopifnot(length(v)==length(utss))
				v
			}))
			dim(tssmat[[tp]])
			colnames(tssmat[[tp]])=libs
			 #and for the ends
			libs=unique(tp.taccs)
			endmat[[tp]]=simplify2array(mclapply(mc.cores=10,libs,function(acc){
				v=rep(0,length(uends))
				v[as.vector(strand(uends)=='-')]=unlist(viewSums(GRViews(ts[[acc]]$neg,uends[strand(uends)=='-'])))
				v[as.vector(strand(uends)=='+')]=unlist(viewSums(GRViews(ts[[acc]]$pos,uends[strand(uends)=='+'])))
				v
			}))
			dim(tssmat[[tp]])
			dim(endmat[[tp]])
			colnames(endmat[[tp]])=libs
			stopifnot(nrow(tssmat[[tp]])==nrow(endmat[[tp]]))
			#now put in rownames
			rownames(tssmat[[tp]])=rownames(endmat[[tp]])=utss$TrID
			tokeep	= rowSums(tssmat[[tp]])!=0 & rowSums(endmat[[tp]])!=0 
			tssmat[[tp]]=tssmat[[tp]][ tokeep,]#get rid of our zeros
			endmat[[tp]]=endmat[[tp]][ tokeep,]#get rid of our zeros
				#normalize our tagseq with just the library size for now
			gcov=with(tagseq.df,genome.coverage[match(unique(tp.taccs),accession)])
			endmat[[tp]]=sweep(endmat[[tp]],MARGIN=2,STAT=gcov,FUN='/')
			stopifnot(nrow(tssmat[[tp]])==nrow(endmat[[tp]]))
				#fix the names
		}
	
		message('data matrices computed')
		sapply(tssmat,dim)
		sapply(endmat,dim)

		#quantile normalization across tss if we are using the unaltered data
		if(grepl('unprocessed',cage.data)){
			tssmat=sapply(names(tssmat),function(tp){
				apply(tssmat[[tp]],MARGIN=2,FUN=qnormvect)
			})
		}
		#We now have matrices with the rows representing genes at a given timepoint 
		###### 6 Now get correlations of matching libraries
		#get the correlation for each tss
		message('computing gene correlations')
		gene.cors.tp[[setn]]=sapply(simplify=F,tps,function(tp){
			tp.accs=accs[[tp]]
			tp.taccs=taccs[[tp]]
			tp.accs = tp.accs[ tp.accs %in% names(cg) & tp.taccs %in% names(ts) ]
			tp.taccs = tp.taccs[ tp.accs %in% names(cg) & tp.taccs %in% names(ts) ]
			tmat=tssmat[[tp]][,tp.accs]
			emat=endmat[[tp]][,tp.taccs]
			totest=which(rowSums(tssmat[[tp]])!=0 & rowSums(endmat[[tp]])!=0)
			stopifnot(c('integer','numeric')%in%is(totest))
			z=sapply(totest,function(i){
				a =	tmat[i,]
				b =	emat[i,]
				#get the correlation over the comparable libraries
				# cor.obj=cp.test(matrix(a),matrix(b))
				# c('pval'=cor.obj$S.p.value,'cor'=cor.obj$correlation)
				cor(a,b,meth='s')
			})
		})
	}
}

	# ####7 Now calculate teh pvalues for teh correlations via permutation
	# #and the pval as computed by shuffling - very 
	# gene.shuffle.pvals[[setn]]=sapply(simplify=F,tps,function(tp){
	# 	tp.accs=accs[[tp]]
	# 	tp.taccs=taccs[[tp]]
	# 	totest=which(rowSums(tssmat[[tp]])>0 & rowSums(endmat[[tp]])>0)
	# 	totest=sample(totest,min(1000,length(totest)))
	# 	sapply(totest,function(i){
	# 		a =	tssmat[[tp]][i,tp.accs]
	# 		b =	endmat[[tp]][i,tp.taccs]
	# 		if(all(a==0)|all(b==0)){return(NA)}
	# 		cp.test(matrix(a),matrix(b))$S.p.value
	#   	})
	# })
	#print summary stats to files
	cat(file='output.txt','mean correlation across all timepoints')
	cat(file='output.txt',mean(unlist(gene.cors.tp[[setn]]),na.rm=T))
	# cat(file='output.txt','mean pvals across all timepoints')
	# cat(file='output.txt',mean(unlist(gene.shuffle.pvals[[setn]]),na.rm=T))
	message('results done for')
	message(setn)
	# save('tmp.image.R')
	#output the mean of all correlations
	library(ggplot2)
	library(LSD)
	for(tp in tps){
		rst=rowSums(tssmat[[tp]])
		red=rowSums(endmat[[tp]])
		#cors
		cors=gene.cors.tp[[setn]][[tp]]
		#pvals
		# pvals=gene.shuffle.pvals[[setn]][[tp]]
		# pvals=-log10(pvals)
		# pvals[pvals==Inf]=max(pvals[pvals!=Inf])
		# names(pvals)=gsub(x=names(pvals),pattern='(FBtr\\d+).*',rep='\\1')
		# #plot our correlations as a fuction of mean cage
		pdf(paste0('tmpcage_tagseq_correlation_',tp,'.pdf'))
		heatscatter(x=rank(rst),y=rank(red),xlab=paste0('Rank total cage signal',tp),ylab='Rank total tagseq signal',ncol=100,log='')
		dev.off()
		#plot our correlations as a fuction of mean cage
		pdf(paste0('correlation_ranked_by_cage_scatter_',tp,'.pdf'))
		heatscatter(x=rank(rst)[names(cors)],y=cors,xlab='Rank in total cage signal',ylab='cage/tagseq spearman correlation',log='',ylim=c(-1,1),cor=F)
		dev.off()
		#plot our correlations as a fuction of mean tagseq
		pdf(paste0('correlation_ranked_by_tagseq_scatter',tp,'.pdf'))
		heatscatter(x=rank(red)[names(cors)],y=cors,xlab='Rank total tagseq signal',ylab='cage/tagseq spearman correlation',log='',ylim=c(-1,1),cor=F)
		dev.off()
		# #plot our pvals as a fuction of mean cage
		# pdf(paste0('correlation_ranked_by_cage_scatter_',tp,'.pdf'))
		# heatscatter(x=rank(rst)[names(pvals)],y= pvals,xlab='Rank in total cage signal',ylab='cage/tagseq correlation -log10(pval)',log='',cor=F)
		# dev.off()
		# #plot our pvals as a fuction of mean tagseq
		# pdf(paste0('correlation_ranked_by_tagseq_scatter',tp,'.pdf'))
		# heatscatter(x=rank(red)[names(pvals)],y= pvals,xlab='Rank total tagseq signal',ylab='cage/tagseq correlation -log10(pval)',log='',ylim=c(-1,1),cor=F)
		# dev.off()
	}
}}

save(gene.cors.tp,file='data/objects/cage_tag_cors.object.R')
save(gene.shuffle.pvals,file='data/objects/cage_tag_pvals.object.R')


#Now do some eCDF plots and plots between pairs of correlations
setwd(outfolder)
#Ecdf plots - just do this for each dataset

# for(setn in names(gene.cors.tp)){
# 	corvect=do.call('c',gene.cors.tp[[setn]])
# 	pdf(paste0('correllation_ecdf',setn,'.pdf'))
# 	plot(ecdf(corvect))
# 	dev.off()
# }
library(reshape2)
meltcors=melt(gene.cors.tp)
colnames(meltcors)=c('Correlation','Timepoint','Dataset')#name columns
meltcors$Dataset[meltcors$Dataset=="all.cage.unprocessed.object"]<-'TSS_quantnorm'
p=ggplot(meltcors, aes(x=Correlation, colour = Dataset)) + stat_ecdf() 
pdf(paste0('correllation_ecdf','allsets','.pdf'))
print(p)
dev.off()

#Now the scatterplot combinations
combmat=combn(names(gene.cors.tp),2)
sapply(1:ncol(combmat),function(i){
	x=do.call('c',gene.cors.tp[[combmat[1,i]]])
	y=do.call('c',gene.cors.tp[[combmat[2,i]]])
	x=x[names(x)%in%names(y)]
	y=y[names(y)%in%names(x)]
	stopifnot(length(x)==length(y))
	jpeg( paste0('corscatter',paste0(combmat[,i],collapse=''),'.jpeg',collapse='') )
	p=qplot(data=data.frame(),x,y,xlab=combmat[1,i],ylab=combmat[2,i])
	print(p)
	dev.off()
})

#Now compare library size norm and pl norm corrlation vs. signal on same graph with splines
# pdf(paste0('pl_vs_libnorm_corsplines',tp,'.pdf'))
pdf('/g/furlong/Harnett/TSS_CAGE_myfolder/analysis/norm_comparison_vs_tagseq.R/pl_vs_libnorm_corvsig.pdf')
sparpar=0.3
rst=rowSums(tssmat[[tp]])#rank in cage signal
red=rowSums(endmat[[tp]])
tp=tps[2]#timepoints
setn=names(gene.cors.tp)[1]#name of set - library normalization
cors=gene.cors.tp[[setn]][[tp]]
cors=cors[names(cors)%in%names(rst)]
rst=rank(rst)[names(cors)]#just the ones in both sets
#plot points
plot(cex=0.3,col='lightblue',x=rst,y=cors,,xlab='Rank in total cage signal',ylab='cage/tagseq spearman correlation',log='',ylim=c(-1,1))
#fit spline
lines(smooth.spline(spar=sparpar,rst,y=cors),lty=1,lwd=4,col='blue')
#now select data for power law normalizaed stuff
setn=names(gene.cors.tp)[2]
cors=gene.cors.tp[[setn]][[tp]]
cors=cors[names(cors)%in%names(rst)]
#plot points
points(cex=0.3,x=rst,y=cors,col='pink')
#fit spline
lines(smooth.spline(spar=sparpar,rst,y=cors),lty=1,lwd=4,col='red')
legend(x='bottomleft',fill=c('blue','red'),legend=c('library_normalization','powerlaw_normalization'))
dev.off()


#Now compare library size norm and pl norm corrlation vs. signal on same graph with splines
# pdf(paste0('pl_vs_libnorm_corsplines',tp,'.pdf'))
pdf('/g/furlong/Harnett/TSS_CAGE_myfolder/analysis/norm_comparison_vs_tagseq.R/pl_vs_libnorm_pltagseq.pdf')
sparpar=0.3
rst=rowSums(tssmat[[tp]])#rank in cage signal
red=rowSums(endmat[[tp]])
tp=tps[2]#timepoints
setn='cg.libnorm.ts.pl'#name of set - library normalization
cors=gene.cors.tp[[setn]][[tp]]
cors=cors[names(cors)%in%names(rst)]
rst=rank(rst)[names(cors)]#just the ones in both sets
#plot points
plot(cex=0.3,col='lightblue',x=rst,y=cors,main='powerlaw vs. library size normalization - power law normalized Tagseq data',xlab='Rank in total cage signal',ylab='cage/tagseq spearman correlation',log='',ylim=c(-1,1))
#fit spline
lines(smooth.spline(spar=sparpar,rst,y=cors),lty=1,lwd=4,col='blue')
#now select data for power law normalizaed stuff
setn='cg.pl.map.ts.pl'
cors=gene.cors.tp[[setn]][[tp]]
cors=cors[names(cors)%in%names(rst)]
#plot points
points(cex=0.3,x=rst,y=cors,col='pink')
#fit spline
lines(smooth.spline(spar=sparpar,rst,y=cors),lty=1,lwd=4,col='red')
legend(x='bottomleft',fill=c('blue','red'),legend=c('library_normalization','powerlaw_normalization'))
dev.off()


