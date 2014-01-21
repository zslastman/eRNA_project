setwd ( '/g/furlong/Harnett/TSS_CAGE_myfolder/' )
source ( 'src/tss_cage_functions.R' )
rootfolder='/g/furlong/Harnett/TSS_CAGE_myfolder/'
outfolder='/g/furlong/Harnett/TSS_CAGE_myfolder/analysis/cage_or_tagseq_vs_RNAseq'
dir.create(outfolder,showWarnings=F)

#load two cage datasets
# load(paste0(rootfolder,'data/objects/cg.libnorm.object.R'))
load(paste0(rootfolder,'data/objects/cg.pl.map.object.R'))
#load RNAseq
load(paste0(rootfolder,'data/objects/rna.seq.object.R'))
#and tagseq
load('data/objects/all.tagseq.unprocessed.R')
cg=cg.pl.map


load(file.tss)
load(file.transcripts)
# tss.gr=tss.gr[1:100,]

tss.gr=sort(tss.gr[!duplicated(tss.gr)])#remo

tss.gr$tsscage = sapply(names(cg),function(acc){
	unlist(viewSums(GRViews(cg[[acc]]$pos+cg[[acc]]$neg,resize(tss.gr,width=500,fix='center'))))
})
#and tagseq
#and load rna data
tss.gr$tsstagseq = cbind(sappy(tps,function(tp){sapply(names(ts),function(acc){


		unlist(viewSums(GRViews(ts[[acc]]$both,transcripts.gr,width=500,fix='center'))))
	})
})
#collapse the 
#and load rna data
tss.gr$c.rnaseq = sapply(names(c.rna.seq),function(acc){
	unlist(viewSums(GRViews(c.rna.seq[[acc]]$pos+c.rna.seq[[acc]]$neg,resize(tss.gr,width=500,fix='center'))))
})
tss.gr$lab.rnaseq = sapply(names(rna.seq),function(acc){
	unlist(viewSums(GRViews(rna.seq[[acc]],rna.seq[[acc]],resize(tss.gr,width=500,fix='center'))))
})

#normalize?
library(LSD)


unlist(viewSums(GRViews(ts[[acc]]$both,sort(resize(resize(transcripts.gr,width=1,fix='end'),width=500,fix='center')))))


#Correlation plots over the three timepoints
cage.tps=paste0('tp',gsub(colnames(tss.gr$tsscage),pat='.*_(\\d\\d?)(\\d\\d?h).*',rep='\\1\\2'))
tagseq.tps=paste0('tp',gsub(colnames(tss.gr$tsstagseq),pat='.*_(\\d\\d?)_(\\d\\d?h)_.*',rep='\\1\\2'))
stopifnot(cage.tps[1:320]%in%tps)
stopifnot(tagseq.tps%in%tps)
cage.tpsum=sapply(tps,function(tp){rowSums(tss.gr$tsscage[,cage.tps==tp])})
tagseq.tpsum=sapply(tps,function(tp){rowSums(tss.gr$tsstagseq[,tagseq.tps==tp])})
rna.tpsum=tss.gr$c.rnaseq[,c(2,4,1)]
colnames(rna.tpsum)=tps
pdf(paste0(outfolder,'/scatterplotsbytp.pdf'))
for(tp in tps){

	#Rank rank plots of summed cagev vs tagseq
	heatscatter(x=rank(cage.tpsum[,tp]),rank(tagseq.tpsum[,tp]),main=paste0('Rank Correlation of summed Signal 	',tp),cor=F,xlab='Rank Cage',ylab='Rank Tagseq')
	title(sub=paste0('spearman cor = ',cor(meth='s',cage.tpsum[,tp],tagseq.tpsum[,tp]) ))

	#Rank rank pots of summed cage vs rnaseq
	heatscatter(x=rank(rna.tpsum[,tp]),rank(cage.tpsum[,tp]),    main=paste0('Rank Correlation of Summed Signal ',tp),cor=F,xlab='Rank RNAseq',ylab='Rank TSS Cage')
	title(sub=paste0('spearman cor = ',cor(meth='s',rna.tpsum[,tp],cage.tpsum[,tp]) ))

	#Rank rank pots of summed tagseq vs rnaseq
	heatscatter(x=rank(rna.tpsum[,tp]),rank(tagseq.tpsum[,tp]),main=paste0('Rank Correlation of Summed Signal ',tp),cor=F,xlab='Rank Tagseq',ylab='Rank TSS Cage')
	title(sub=paste0('spearman cor = ',cor(meth='s',rna.tpsum[,tp],tagseq.tpsum[,tp]) ))

}
dev.off()

tss.gr$tsstagseq.o=tss.gr$tsstagseq
#normalize the tagseq
load('data/objects/tagseq.df.object.R')
tss.gr$tsstagseq=sweep(tss.gr$tsstagseq,2,tagseq.df$genome.coverage,FUN='*' )

#make matrices - 3 columns, one for each timepoint, with rows as TSS.
ctcors=sapply((1:nrow(cage.tpsum)),function(i){cor(meth='s',cage.tpsum[i,],tagseq.tpsum[i,])})
crcors=sapply((1:nrow(cage.tpsum)),function(i){cor(meth='s',cage.tpsum[i,],rna.tpsum[i,])})
trcors=sapply((1:nrow(cage.tpsum)),function(i){cor(meth='s',tagseq.tpsum[i,],rna.tpsum[i,])})

ma <- function(x,n=5){filter(x,rep(1/n,n), sides=2)}

#plot these correlations as a function of the mean signal
#Cage and Tagseq ranked by Cage
pdf('Cage_v_Tagseq_v_RNA_tpcor.pdf')
i=order(rowSums(cage.tpsum))
mtit='TSS\' Correlation between Cage and Tagseq over three timepoints'
qplot(geom='smooth',x=rank(rowSums(cage.tpsum)[i]),y=ma(ctcors[i],100),ylim=c(-1,1),ylab='Spearman\'s rho',xlab='Rank in Cage Expression',main=mtit)+geom_point(size=0.3)
#dev.off()
#Cage and RNAseq ranked by Cage
#pdf('tmp.pdf')
i=order(rowSums(cage.tpsum))
mtit='TSS\' Correlation between Cage and RNAseq over three timepoints'
qplot(geom='smooth',x=rank(rowSums(cage.tpsum)[i]),y=ma(crcors[i],100),ylim=c(-1,1),ylab='Spearman\'s rho',xlab='Rank in Cage Expression',main=mtit)+geom_point(size=0.3)
#dev.off()
#now but ranked by other
#pdf('tmp.pdf')
i=order(rowSums(tagseq.tpsum))
mtit='TSS\' Correlation between Cage and Tagseq over three timepoints'
qplot(geom='smooth',x=rank(rowSums(rna.tpsum)[i]),y=ma(ctcors[i],100),ylim=c(-1,1),ylab='Spearman\'s rho',xlab='Rank in Tagseq Expression',main=mtit)+geom_point(size=0.3)
#dev.off()
#
#pdf('tmp.pdf')
i=order(rowSums(rna.tpsum))
mtit='TSS\' Correlation between Cage and RNAseq over three timepoints'
qplot(geom='smooth',x=rank(rowSums(rna.tpsum)[i]),y=ma(crcors[i],100),ylim=c(-1,1),ylab='Spearman\'s rho',xlab='Rank in RNAseq Expression',main=mtit)+geom_point(size=0.3)
#dev.off()
#plot these correlations as a function of the mean signal
#Tagseq and Cage ranked by Tagseq
#pdf('tmp.pdf')
i=order(rowSums(tagseq.tpsum))
mtit='TSS\' Correlation between Tagseq and RNAseq over three timepoints'
qplot(geom='smooth',x=rank(rowSums(tagseq.tpsum)[i]),y=ma(trcors[i],100),ylim=c(-1,1),ylab='Spearman\'s rho',xlab='Rank in Tagseq Expression',main=mtit)+geom_point(size=0.3)
#dev.off()
#Tagseq and RNAseq ranked by Tagseq
#pdf('tmp.pdf')
mtit='TSS\' Correlation between Tagseq and RNAseq over three timepoints'
i=order(rowSums(rna.tpsum))
qplot(geom='smooth',x=rank(rowSums(rna.tpsum)[i]),y=ma(trcors[i],100),ylim=c(-1,1),ylab='Spearman\'s rho',xlab='Rank in RNAseq Expression',main=mtit)+geom_point(size=0.3)
dev.off()
