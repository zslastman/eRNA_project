#script to create plots showing agreement of cage, tagseq and RNAseq
setwd ( '/g/furlong/Harnett/TSS_CAGE_myfolder/' )
source ( 'src/tss_cage_functions.R' )
library(LSD)
dir.create(outfolder,showWarnings=F)

rootfolder='/g/furlong/Harnett/TSS_CAGE_myfolder/'
outfolder='/g/furlong/Harnett/TSS_CAGE_myfolder/analysis/cage_or_tagseq_vs_RNAseq'


##for now just use lib normalization for the ts
genes.gr$ts.pl=sweep(genes.gr$ts,2,tagseq.df$genome.coverage,FUN='*' )


#Correlation plots over the three timepoints
#get the relevant timepoints for each dataset
cage.tps=paste0('tp',gsub(colnames(genes.gr$cg.pl),pat='.*_(\\d\\d?)(\\d\\d?h).*',rep='\\1\\2'))
tagseq.tps=paste0('tp',gsub(colnames(genes.gr$ts.pl),pat='.*_(\\d\\d?)_(\\d\\d?h)_.*',rep='\\1\\2'))
stopifnot(cage.tps[1:320]%in%tps)
stopifnot(tagseq.tps%in%tps)
#create matrices with columns for each timepoint, rows for each gene
#and note that for now we can only compare 2 of our timepoints
cage.tpsum=sapply(tps[1:2],function(tp){rowSums(genes.gr$cg.pl[,cage.tps==tp])})
tagseq.tpsum=sapply(tps[1:2],function(tp){rowSums(genes.gr$ts.pl[,tagseq.tps==tp])})
rna.tpsum=genes.gr$c.rna.seq[,c(1,3)]
colnames(rna.tpsum)=tps[1:2]
pdf(paste0(outfolder,'/scatterplotsbytp_nomapfilt.pdf'))
for(tp in tps){

	#Rank rank plots of summed cagev vs tagseq
	heatscatter(x=rank(cage.tpsum[,tp]),rank(tagseq.tpsum[,tp]),main=paste0('Rank Correlation of summed Signal 	',tp),cor=F,xlab='Rank Cage',ylab='Rank Tagseq')
	title(sub=paste0('spearman cor = ',cor(meth='s',cage.tpsum[,tp],tagseq.tpsum[,tp]) ))
	#Rank rank pots of summed cage vs rnaseq
	heatscatter(x=rank(rna.tpsum[,tp]),rank(cage.tpsum[,tp]), main=paste0('Rank Correlation of Summed Signal ',tp),cor=F,xlab='Rank RNAseq',ylab='Rank TSS Cage')
	title(sub=paste0('spearman cor = ',cor(meth='s',rna.tpsum[,tp],cage.tpsum[,tp]) ))
	#Rank rank pots of summed tagseq vs rnaseq
	heatscatter(x=rank(rna.tpsum[,tp]),rank(tagseq.tpsum[,tp]),main=paste0('Rank Correlation of Summed Signal ',tp),cor=F,xlab='Rank Tagseq',ylab='Rank TSS Cage')
	title(sub=paste0('spearman cor = ',cor(meth='s',rna.tpsum[,tp],tagseq.tpsum[,tp]) ))

}
dev.off()

# #normalize the tagseq
# load('data/objects/tagseq.df.object.R')
# tss.gr$tsstagseq=sweep(tss.gr$tsstagseq,2,tagseq.df$genome.coverage,FUN='*' )



#make matrices - 3 columns, one for each timepoint, with rows as TSS.
ctcors=sapply((1:nrow(cage.tpsum)),function(i){cor(meth='s',cage.tpsum[i,],tagseq.tpsum[i,])})
crcors=sapply((1:nrow(cage.tpsum)),function(i){cor(meth='s',cage.tpsum[i,],rna.tpsum[i,])})
trcors=sapply((1:nrow(cage.tpsum)),function(i){cor(meth='s',tagseq.tpsum[i,],rna.tpsum[i,])})

ma <- function(x,n=5){filter(x,rep(1/n,n), sides=2)}

#plot these correlations as a function of the mean signal
#Cage and Tagseq ranked by Cage
pdf('Cage_v_Tagseq_v_RNA_tpcor_nomapfilt.pdf')
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




#Correlation plots over the three timepoints, but this time with the mappability filtering
#get the relevant timepoints for each dataset
cage.tps=paste0('tp',gsub(colnames(genes.gr$cg.mapfilt.pl),pat='.*_(\\d\\d?)(\\d\\d?h).*',rep='\\1\\2'))
stopifnot(cage.tps[1:320]%in%tps)
cage.tpsum=sapply(tps[1:2],function(tp){rowSums(genes.gr$cg.mapfilt.pl[,cage.tps==tp])})
pdf(paste0(outfolder,'/scatterplotsbytp_mapfilt.pdf'))
for(tp in tps[1:2]){

	#Rank rank plots of summed cagev vs tagseq
	heatscatter(x=rank(cage.tpsum[,tp]),rank(tagseq.tpsum[,tp]),main=paste0('Rank Correlation of summed Signal 	',tp),cor=F,xlab='Rank Cage-mapfiltered',ylab='Rank Tagseq')
	title(sub=paste0('spearman cor = ',cor(meth='s',cage.tpsum[,tp],tagseq.tpsum[,tp]) ))

	#Rank rank pots of summed cage vs rnaseq
	heatscatter(x=rank(rna.tpsum[,tp]),rank(cage.tpsum[,tp]), main=paste0('Rank Correlation of Summed Signal ',tp),cor=F,xlab='Rank RNAseq',ylab='Rank Gene Cage-mapfiltered')
	title(sub=paste0('spearman cor = ',cor(meth='s',rna.tpsum[,tp],cage.tpsum[,tp]) ))

	#Rank rank pots of summed tagseq vs rnaseq
	heatscatter(x=rank(rna.tpsum[,tp]),rank(tagseq.tpsum[,tp]),main=paste0('Rank Correlation of Summed Signal ',tp),cor=F,xlab='Rank Tagseq',ylab='Rank Gene Cage-maptilered')
	title(sub=paste0('spearman cor = ',cor(meth='s',rna.tpsum[,tp],tagseq.tpsum[,tp]) ))

}
dev.off()




#make matrices - 3 columns, one for each timepoint, with rows as TSS.
ctcors=sapply((1:nrow(cage.tpsum)),function(i){cor(meth='s',cage.tpsum[i,],tagseq.tpsum[i,])})
crcors=sapply((1:nrow(cage.tpsum)),function(i){cor(meth='s',cage.tpsum[i,],rna.tpsum[i,])})
trcors=sapply((1:nrow(cage.tpsum)),function(i){cor(meth='s',tagseq.tpsum[i,],rna.tpsum[i,])})

ma <- function(x,n=5){filter(x,rep(1/n,n), sides=2)}

#plot these correlations as a function of the mean signal
#Cage and Tagseq ranked by Cage
pdf('Cage_v_Tagseq_v_RNA_tpcor_mapfilt.pdf')
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





