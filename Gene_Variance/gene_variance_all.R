#A Load libraries,functions,data
library(reshape2)
library(ggplot2)
setwd ( '/g/furlong/Harnett/TSS_CAGE_myfolder/' )
source ( 'src/tss_cage_functions.R' )
#load cage data and tagseq data
load(file.tss)
load(file.synonyms)
cg=cage.tag.rles
rm(cage.tag.rles)
#timepoints
outfolder='analysis/gene_variance_all'
tps=c('tp24h','tp68h','tp1012h')
dir.create(outfolder)

####E     collect info for all transcripts
message('collect info for all transcripts')
for now just use   counts
cgenv<-new.env()
# cage.data='/g/furlong/Harnett/TSS_CAGE_myfolder/data/objects/cg.pl.object.R'
cage.data='/g/furlong/Harnett/TSS_CAGE_myfolder/data/objects/all.cage.unprocessed.object.R'
load(file=cage.data,env=cgenv)
cg=cgenv[[ls(cgenv)[1]]]


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
tss.gr$unp.cage=sapply(names(cage.tag.rles),function(acc){
	unlist(viewSums(GRViews(cage.tag.rles[[acc]]$both,resize(tss.gr,width=500,fix='center'))))
})

#We want a function that takes in a go term, and outputs FBtrs
golist=list(
	dna_binding='GO:0005575',#take go term
	regulatory_DNA_binding='GO:0000976',
translational_regulator = 'GO:0045182',
transcriptional_regulator = 'GO:0030528',
chromatin_binding='GO:0003682',
Structural_Molecule_Activity='GO:0005198',
Sequence_Specific_DNA_Binding='GO:0043565', 
regulatory_region_DNA_binding='GO:0000975',
olfactory_receptor_activity ='GO:0004984',
Gprotein_coupled_olfactory_receptor_activity='GO:0038022',
olfactory_receptor_binding='GO:0031849',
Biological_Adhesion='GO:0022610',
Metabolic_Process='GO:0008152',
Reproduction='GO:0000003',
Behavior_defective='GO:0007610')
go_tr_search<-function(goterms,fbgo=fb_go.df,fbtr=fbgn_fbtr_fbpp.df){
	fbgns = with(fbgo,FBgn[which(GO%in%goterms)])
	# stopifnot(length(fbgns)>0)
	fbtrs = with( fbtr , FBtr[ which(as.character(FBgn)%in%as.character(fbgns)) ]   )
	fbtrs
}
#now ge the TrIDs for each Go
gotrlist=sapply(golist,go_tr_search)
sapply(gotrlist,length)[]

#back up the matrix

alltssmat=tss.gr$tsscage

#very simple approach - we sort our TFs by mean, split into facets by GO term,
#and then produce boxplots for each one 
fbgn2symb=unique(fb_go.df$symb)#vector for converting fbgns to symbols
names(fbgn2symb)=unique(fb_go.df$FBgn)

go=names(gotrlist)[c]
tp=tps[1]

goodgos=names(gotrlist)[-c(1,4,8,10,11,12,15)]
for(go in goodgos){
	for(tp in tps){
		ingotrs=rownames(alltssmat)%in%gotrlist[[go]]
		mat=alltssmat[ingotrs,paste0('tp',accession.df$timepoint)==tp][,]


		means=rowMeans(mat)
		mat=mat[means>quantile(means,0),]
		dim(mat)
		#filter genes by 
		#take only the most expressed Transcript for each gene
		mdat=melt(mat)
		sorted=unique(mdat$Var1[order(means[as.character(mdat$Var1)])])
		mdat$Var1<-factor(mdat$Var1,levels=sorted,ordered=T)
		#Now replace the FBtr with gene symbols
		FBgns = fbgn_fbtr_fbpp.df$FBgn[match(mdat$Var1,fbgn_fbtr_fbpp.df$FBtr)]
		mdat$symb =fbgn2symb[FBgns]
		levelfbgns= fbgn_fbtr_fbpp.df$FBgn[match(sorted,fbgn_fbtr_fbpp.df$FBtr)]
		mdat$symb = factor( mdat$symb ,levels = unique(fbgn2symb[levelfbgns] ) )
		mdat$logval = log10(mdat$value)
		mdat$logval[mdat$logval==-Inf]=0
		dim(mdat)


		#only ones with reasonable high mean
		#now plot the boxplots
		ggsave(paste0(outfolder,'/expr_boxplots_log10_',tp,'_',go,'.pdf'),width=21,
			qplot(data=mdat[],geom='boxplot',x=symb,color=symb,y=logval,ylab='Log10 expression',main=paste0('Expression_Boxplots_ Go term:\n',go,'\ntimepoint ',tp))+
			scale_x_discrete(labels='')+
			guides(col=guide_legend(ncol=3,size=0.2))
			)
}
}
# 



#Let's write a function that calculates dispersion using a set of input vectors
log_disp<-function(alltssmat){
	logalltssmat=log(alltssmat)
	logalltssmat[logalltssmat==-Inf]<-0
	lognvar=apply(logalltssmat,1,var)
	tagcounts=rowSums(alltssmat)
	above.cutoff=tagcounts>100
	dispersion=lognvar - 1/tagcounts
}

DESeq_dis<-function(alltssmat){

}




disps=list()
for(tp in c('24h','68h','1012h')){
	allmat=tss.gr$tsscage
	#filter for timepoint
	alltssmat= allmat[,grepl(x=colnames(allmat),pattern=tp)]
	#what we want ot plot is not the variance but rather the variance in the count of normalized log counts
	logalltssmat=log10(alltssmat)
	logalltssmat[logalltssmat==-Inf]<-0
	#calculate varaince is log normal plots
	lognvar=apply(logalltssmat,1,var)

	#plot this as a function of mean
	pdf(paste0('analysis/gene_variance_all/var_mean',tp,'.pdf'))
	heatscatter(rowSums(logalltssmat),lognvar,xlab='log normalized count - mean',ylab='log normalized count - variance')
	dev.off()
	#use VDS as well
	#now subtract the poisson variance to get the dispersion
	tagcounts=rowSums(alltssmat)
	lcounts=log(tagcounts)
	above.cutoff=tagcounts>100
	dispersion=lognvar - 1/tagcounts
	disps[[tp]]=dispersion
	#now plot this as a function of mean
	pdf(paste0('analysis/gene_variance_all/var_minus_oneovermean',tp,'.pdf'))
	heatscatter(rowSums(logalltssmat),dispersion,xlab='log normalized count - mean',ylab='log normalized count - dispersion')
	dev.off()
}




#do boxplots with the dispersion
dispersion.list=sapply(gotrlist,function(trs){
	dispersion[names(dispersion) %in% trs & above.cutoff]
})
pdf('analysis/gene_variance_all/GO_disp_boxplots.pdf')
qplot(data=melt(dispersion.list),color=L1,x=L1,y=value,xlab='geneset',ylab='dispersion',geom='boxplot')
dev.off()

#Now also do boxplot with genes split by 'high,low and medium' since charles suggested it
lev=quantile(x=lcounts,probs=c(0.33,0.66,1))
names(lev)=c('low','medium','high')
expr_level=rep(0,length(lcounts))
for(lev.n in names(lev)){
	expr_level[lcounts<=lev[lev.n] & expr_level==0]<-lev.n
}
names(expr_level)=names(lcounts)
mdat=melt(dispersion.list)
mdat$FBtr=as.character(melt(sapply(dispersion.list,names))$value)

pdf('analysis/gene_variance_all/GO_disp_boxplots_split_by_expr_level.pdf')
qplot(data=mdat,color=L1,x=expr_level[mdat$FBtr],y=value,xlab='geneset',ylab='dispersion',geom='boxplot')
dev.off()


disp.mean.loess = loess(dispersion[above.cutoff] ~ lcounts[above.cutoff],span=2)
pdf('analysis/gene_variance_all/disp.mean.loess.pdf')
plot(lcounts[above.cutoff],dispersion[above.cutoff])
lo.pred=predict(disp.mean.loess,se=T)
i <- order(lcounts[above.cutoff])
lines(lcounts[above.cutoff][i],lo.pred$fit[i],pch=2, col="red")
dev.off()

#let's see how our variance look if we just calculate it for a single 




#select our 

#plot variation as a function of mean
pdf('analysis/gene_variance/variation_allgenes.pdf')

qplot(x=rowMeans(alltssmat),y=gene.var,ylab='Standard Deviation',xlab='Mean Expression',geom='point',log='x',main='Standard Deviation vs. Mean expression - All genes')
qplot(x=rowMeans(alltssmat),y=coef.var,ylab='Coefficient of Variation - SD/MEAN',xlab='Mean Expression',geom='point',log='x',main='Coefficient of Variation vs. Mean expression - All genes')

dev.off()

library(reshape2)
library(limma)
library(gllm,lib='~/Harnett/R')


#fit a model with no structure in the rows or columns
smallmat=matrix(exp(rnorm(100,3,0.5)),ncol=20,nrow=5)
meltmat=melt(smallmat)
smallmat
fit=gllm(y=meltmat$value,1:n



#Do some MA plots....

maPlotVects<-function(y,z,tit,xl='A',yl='M'){
	# z->>z
	# y->>y
	nonz= y!=0&z!=0
	y=y[nonz]
	z=z[nonz]
	m=log2(y)-log2(z)
	a=0.5*(log2(y)+log2(z))
	nonz=is.finite(m)
	plotMA(list(A=a,M=m),xlim=c(0,max(a)),ylim=c(-max(abs(m)),max(abs(m))),
		xlab=xl,ylab=yl,main=tit)
	big=a>quantile(unique(a),0.3)
	# mod=rlm(m[big]~a[big])
	# abline(mod,col='red',lty='dashed')
	# title(sub=paste0('slope = ',as.character(mod$coefficients[2])))

}
l1=order(accession.df$library.size[1:280],decreasing=T)[3]
l2=which.min(accession.df$library.size[1:280])
pdf('analysis/gene_variance_all/MA_unp.pdf')
maPlotVects(tss.gr$unp.cage[,l1],tss.gr$unp.cage[,l2],tit='Unprocessed Cage');dev.off()
pdf('analysis/gene_variance_all/MA_pl.pdf')
maPlotVects(tss.gr$tsscage[,l1],tss.gr$tsscage[,l2],tit='Normalized Cage')
dev.off()


	row(meltmat),~ rep(mean(meltmat$value),nrow(meltmat)) + meltmat$Var1 + meltmat$Var2)


for(tp in names(disps)){
dispersion=disps[[tp]]
#print expressed genes for dispersion rank analysis
expr = rowMeans(tss.gr$tsscage)>10
o=order(dispersion[expr],decreasing=T)
write.table(x=tss.gr$TrID[expr][o],file=paste0(tp,'disp_FBtrs.txt'),quote=F,row.names=F)
}


smallmat=alltssmat[1:3,c(1:2,141:142,301:302)]
smallmat[1,]=rnorm(6,300,30)
smallmat[2,]=rnorm(6,1000,120)
smallmat[3,]=rnorm(6,3000,300)
smallmat[3,]=rnorm(6,3000,300)
smallmat[3,]=rnorm(6,3000,300)
smallmat[3,]=rnorm(6,3000,300)
rownames(smallmat)
colnames(smallmat)
smallmat
smallmat[,5:6]=smallmat[,5:6]*2^rnorm(3*2,3,2)
timepoints = gsub(x=colnames(smallmat),pattern='.*(\\d+h)',rep='\\1')
names(timepoints)=colnames(smallmat)
meltmat=melt(smallmat)
meltmat$tp=timepoints[as.character(meltmat$Var2)]

fit=gllm(y=meltmat$value,1:nrow(meltmat),~ rep(mean(meltmat$value),nrow(meltmat)) + meltmat$Var1 + meltmat$Var2 + meltmat$Var1:meltmat$Var2)

var(fit$coefficients[ grepl(x= names(fit$coefficients) ,pattern='1012h') ])
var(fit$coefficients[ grepl(x= names(fit$coefficients) ,pattern='24h') ])
var(fit$coefficients[ grepl(x= names(fit$coefficients) ,pattern='68h') ])

#fit a model with no structure in the rows or columns
smallmat=matrix(exp(rnorm(10^2,3,0.5)),ncol=10,nrow=10)
meltmat=melt(smallmat)
fit=gllm(y=meltmat$value,1:nrow(meltmat),~ rep(mean(meltmat$value),nrow(meltmat)) + meltmat$Var1 + meltmat$Var2)



limm



#Simplest kind of modelling - what they did in the hourglass paper
#dummy data, a cell is described in the simplest case by it's 





