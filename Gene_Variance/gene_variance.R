#A Load libraries,functions,data
library(reshape2)
library(ggplot2)
setwd ( '/g/furlong/Harnett/TSS_CAGE_myfolder/' )
source ( 'src/tss_cage_functions.R' )
#load cage data and tagseq data
load ( '/g/furlong/Harnett/TSS_CAGE_myfolder/data/objects/unprocessed.cage.tag.object.R' )
load(file.tss)
load(file.synonyms)
cg=cage.tag.rles
rm(cage.tag.rles)
#timepoints
tps=c('tp24h','tp68h','tp1012h')

#B get count matrix data for some selected genes
message(' get count matrix data for some selected genes')

gene.tss=sapply(c('twi','bin','bap','mef2','tin','eve','beat-Ia'),function(symb){
	gene.tss=tss.gr[tss.gr$Gene==synonyms$gene_id[which(synonyms$syn==symb)]]
	gene.tss$name=paste(symb,letters[1:length(gene.tss)])
	gene.tss
})
gene.tss=do.call(c,unname(gene.tss))
gene.tss=sort(gene.tss)
mat=sapply(names(cg),function(acc){
	unlist(viewSums(GRViews(cg[[acc]]$both,resize(gene.tss,width=500,fix='center'))))
})
rownames(mat)=gene.tss$name
dummy data
mat=matrix(nrow=10,rnorm(rpois(1000,100)))
rownames(mat)<-letters[1:10]
Now given matrix we produce our p]ot

#C create plots for our selection of genes
dir.create('analysis/gene_variance/')
pdf('analysis/gene_variance/Selected_gene_boxplot.pdf')
qplot(data=melt(mat),color=gsub('^(\\w+).*','\\1',X1),x=X1,y=value,geom='boxplot',main='Non Normalized Expression for Selected Genes',ylab='Tags',log='y')+scale_color_discrete(name='TSS')
qplot(data=melt(mat),color=gsub('^(\\w+).*','\\1',X1),x=X1,y=value,geom='boxplot',main='Non Normalized Expression for Selected Genes',ylab='Tags',log='')+scale_color_discrete(name='TSS')
#calculate coefficient of variation
coef.var=apply(mat,1,function(row)(sd(row)/mean(row)))
#plot it
qplot(data=data.frame(),fill=rownames(mat),x=rownames(mat),y=coef.var, geom='bar',main='coefficient of variation',ylab='SD/Mean - non normalized tags')
#with and without log
qplot(data=melt(mat),color=gsub('^(\\w+).*','\\1',X1),x=X1,y=value,geom='boxplot',main='Non Normalized Expression for Selected Genes',ylab='Tags',log='y')+scale_color_discrete(name='TSS')
qplot(data=melt(mat),color=gsub('^(\\w+).*','\\1',X1),x=X1,y=value,geom='boxplot',main='Non Normalized Expression for Selected Genes',ylab='Tags',log='')+scale_color_discrete(name='TSS')
dev.off()

#D demonstrating variance of genes to Mat
#histograms
ggsave('/g/furlong/Harnett/TSS_CAGE_myfolder/analysis/gene_variance/Expression_Distributions.pdf',
ggplot(data=melt(mat),aes(x=value,facet=X1,fill=X1))+geom_histogram()+facet_wrap(~X1)
)
#with logs
ggsave('/g/furlong/Harnett/TSS_CAGE_myfolder/analysis/gene_variance/Log_Expression_Distributions.pdf',
ggplot(data=melt(mat),aes(x=value,facet=X1,fill=X1))+geom_histogram()+facet_wrap(~X1)+scale_x_log10()
)
#qqplots
ggsave('/g/furlong/Harnett/TSS_CAGE_myfolder/analysis/gene_variance/qqplots.pdf',
ggplot(data=melt(mat),aes(sample=value,facet=X1,fill=X1))+stat_qq()+facet_wrap(~X1)
)
#qqplots
ggsave('/g/furlong/Harnett/TSS_CAGE_myfolder/analysis/gene_variance/log2_qqplots.pdf',
ggplot(data=melt(mat),aes(sample=log2(value),facet=X1,fill=X1))+stat_qq()+facet_wrap(~X1)
)


####E     collect info for all transcripts
message('collect info for all transcripts')
#for now just use library normalized counts
cgenv<-new.env()
cage.data='/g/furlong/Harnett/TSS_CAGE_myfolder/data/objects/cg.libnorm.object.R'
load(file=cage.data,env=cgenv)
cg=cgenv[[ls(cgenv)[1]]]

quantnormovertss=F

tss.gr=sort(tss.gr[!duplicated(tss.gr)])#remo
alltssmat=sapply(names(cg),function(acc){
	unlist(viewSums(GRViews(cg[[acc]]$pos+cg[[acc]]$neg,resize(tss.gr,width=500,fix='center'))))
})
#qunatile normalize over TSS?
if(quantnormovertss){
		alltssmat=sapply(names(alltssmat),function(tp){
			apply(alltssmat[[tp]],MARGIN=2,FUN=qnormvect)
		})
	}
rownames(alltssmat)=tss.gr$TrID
coef.var=apply(alltssmat,1,function(row)(sd(row)/mean(row)))
gene.var=apply(alltssmat,1,function(row){sd(row)})

####F Load info on gene function
#load the annotation as a list of FBgns
tftab.df	= read.delim('data/tf_fbgns_biomart.txt')
FBgnsA		= as.character(tftab.df[,1])
FBgnsB 		= with(synonyms,syn[grepl('FBgn',syn,fixed=T)])
FBgnsB 		= FBgnsB[!FBgnsB %in% FBgnsA]
# #use our FBgn or other labels too look up tss.
# set.A.tss=with(synonyms,tss.gr[match(gene_id[match(FBgnsA,syn)],tss.gr$Gene)]
# set.B.tss=with(synonyms,tss.gr[match(gene_id[match(FBgnsB,syn)],tss.gr$Gene)]
# #Now get TrIDs from our genes
TrIDsA=sapply(FBgnsA,function(symb){
	tss.gr$TrID[which(tss.gr$Gene==synonyms$gene_id[which(synonyms$syn==symb)])]
})
TrIDsB=sapply(FBgnsB,function(symb){
	tss.gr$TrID[tss.gr$Gene==synonyms$gene_id[which(synonyms$syn==symb)]]
})

#select our 

#plot variation as a function of mean
pdf('analysis/gene_variance/variation_allgenes.pdf')

qplot(x=rowMeans(alltssmat),y=gene.var,ylab='Standard Deviation',xlab='Mean Expression',geom='point',log='x',main='Standard Deviation vs. Mean expression - All genes')
qplot(x=rowMeans(alltssmat),y=coef.var,ylab='Coefficient of Variation - SD/MEAN',xlab='Mean Expression',geom='point',log='x',main='Coefficient of Variation vs. Mean expression - All genes')

dev.off()







#load cage and 3' data, normalized





#get matrix of counts


#Compare sets using boxplots 
