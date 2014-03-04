setwd ( '/g/furlong/Harnett/TSS_CAGE_myfolder/' )
source ( 'src/tss_cage_functions.R' )
#annotation and metadata
load('data/objects/accession.df.full.object.R')
load('data/objects/tagseq.df.object.R')
load('data/objects/expr.mats.object.R')
#cage and tagseq data
tps = c('24h','68h','1012h')


###### 1 Now get correlations of matching libraries
load(file.accession.df)
load('data/objects/tagseq.df.object.R')
#get the right expression data and split by timepoint
cagemat = expr.mats$genes.gr$cg.pl
cagemat = sapply(tps,function(tp) cagemat[, accession.df$time == tp & accession.df$tissue=='embryo'])
tagmat = expr.mats$genes.gr$ts.pl
tagmat = sapply(tps,function(tp)tagmat[, tagseq.df$time ==tp & tagseq.df$tissue=='embryo'])#get the correlation for each tss


######2 prepare lists of matching accessions for each timeoint, so we compare only like accessions
accs=colnames(expr.mats$genes.gr$cg.pl)
taccs=colnames(expr.mats$genes.gr$ts.pl)

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


gene.cors.tp=sapply(simplify=F,tps,function(tp){
	tp.accs  = accs[[tp]]
	tp.taccs = taccs[[tp]]
		sapply(1:nrow(cagemat[[tp]]),function(i){
		#get the correlation over the comparable libraries
		#vector.pvals.shuffle(cagemat[i,accs],tagmat[i,taccs])
		suppressWarnings(cor( cagemat[[tp]][i,tp.accs] ,tagmat[[tp]][i,tp.taccs] ,method='s'))
	})
})

#and the pval as computed by shuffling - very 
gene.shuffle.pvals=sapply(simplify=F,tps,function(tp){
	tp.accs=accs[[tp]]
	tp.taccs=taccs[[tp]]
		mapply(a= cagemat[[tp]][1:3,tp.accs] , b = tagmat[[tp]][1:3,tp.taccs],vector.pvals.shuffle)
		#get the correlation over the comparable libraries
		#vector.pvals.shuffle(cagemat[i,accs],tagmat[i,taccs])
})

#output the mean of all correlations

library(ggplot2)
library(LSD)
# dir.create('analysis/5primer_3prime_covariance')

for(tp in tps){
	#plot our correlations as a fuction of mean cage
	pdf(paste0('analysis/5primer_3prime_covariance/cage_tagseq_correlation_',tp,'.pdf'))
	heatscatter(x=log2(rowSums(cagemat[[tp]])),y=log2(rowSums(tagmat[[tp]])),xlab=paste0('Log2 total cage signal vs. total tagseq signal ',tp),ylab='Log2 cage/tagseq spearman correlation',ncol=100,log='')
	dev.off()
	#plot our correlations as a fuction of mean cage
	pdf(paste0('analysis/5primer_3prime_covariance/correlation_ranked_by_cage_scatter_',tp,'.pdf'))
	heatscatter(x=log2(rowSums(cagemat[[tp]])),y= gene.cors.tp[[tp]],xlab='Log2 total cage signal',ylab='cage/tagseq spearman correlation',log='',ylim=c(-1,1),cor=F)
	dev.off()
	#plot our correlations as a fuction of mean tagseq
	pdf(paste0('analysis/5primer_3prime_covariance/correlation_ranked_by_tagseq_scatter',tp,'.pdf'))
	heatscatter(x=log2(rowSums(tagmat[[tp]])),y=gene.cors.tp[[tp]],xlab='log2 total tagseq signal',ylab='cage/tagseq spearman correlation',log='',ylim=c(-1,1),cor=F)
	dev.off()
}
for(tp in tps){
	#plot our correlations as a fuction of mean cage
	pdf(paste0('analysis/5primer_3prime_covariance/correlation_ranked_by_rankcage_scatter_',tp,'.pdf'))
	heatscatter(x=rank(rowSums(cagemat[[tp]])),y=gene.cors.tp[[tp]],xlab='rank total cage signal',ylab='cage/tagseq spearman correlation',log='',ylim=c(-1,1),cor=F)
	dev.off()
	#plot our correlations as a fuction of mean tagseq
	pdf(paste0('analysis/5primer_3prime_covariance/correlation_ranked_by_ranktagseq_scatter',tp,'.pdf'))
	heatscatter(x=rank(rowSums(tagmat[[tp]])),y=gene.cors.tp[[tp]],xlab='rank total tagseq signal',ylab='cage/tagseq spearman correlation',log='',ylim=c(-1,1),cor=F)
	dev.off()
}






#maybe show that different simulated data sets return different normalized methods as best via this method?
# gene.cors=gene.cors[1:length(gene.cors)]

#let's inspect the guys with weirdly high correlations - rank 800 to 1000 for cage at 68hrs will do this
highcorinds=order(rowSums(cagemat[[tp]]))[400:800]
highcorinds=which(gene.cors.tp[[tp]]>0.6 & rank(rowSums(cagemat[[tp]])) %in% 500:800  )
mean(gene.cors.tp[[tp]][highcorinds],na.rm=T)
#
j=1

doscatter<-function(i){print(qplot( cagemat[[tp]][i,tp.accs] ,tagmat[[tp]][i,tp.taccs],main='5prime vs 3prime scatterplot'))}
doscatter<-function(i){heatscatter( cagemat[[tp]][i,tp.accs] ,tagmat[[tp]][i,tp.taccs],main='5prime vs 3prime scatterplot')}

tp.accs  = accs[[tp]]
tp.taccs = taccs[[tp]]
i=highcorinds[j:(j+10)]
#
# ggsave(file='tmp.pdf',)
pdf('tmp.pdf')
for(n in i){doscatter(n)}
dev.off()
j=j+10


#Now do our correlation plot but also plot the nubmer of zeros
zeronum=apply(cagemat[[tp]],1,function(x){sum(x==0)})
summedsig=rank(rowSums(cagemat[[tp]]))
pdf('tmp.pdf')
heatscatter(x=summedsig,y=gene.cors.tp[[tp]],xlab='rank total cage signal',ylab='cage/tagseq spearman correlation',log='',ylim=c(-1,1))
par(new=T)
plot(x=summedsig,y=zeronum,xlab='',ylab='Number of zeroes',log='')
par(new=F)
heatscatter(x=zeronum,y= gene.cors.tp[[tp]],xlab='Number of zeroes',ylab='correlation',log='')
dev.off()


#Now do our correlation plot but also plot the nubmer of  ones

nonzeromin<-function(x){min(x[x!=0])}
singletagnorm=sapply(cg.pl.map,function(n){nonzeromin(n[[1]][[1]])})

pdf('tmp.pdf')
plot(x=rank(rowSums(cagemat[[tp]])),y=gene.cors.tp[[tp]],xlab='rank total cage signal',ylab='cage/tagseq spearman correlation',log='',ylim=c(-1,1))
par(new=T)
plot(x=rank(rowSums(cagemat[[tp]])),y= apply(cagemat[[tp]],1,function(x){sum(x%in%singletagnorm)}),xlab='',ylab='Number of zeroes',log='')
dev.off()

#run experiment to see if scalings can generate artificial correlations.
#re





