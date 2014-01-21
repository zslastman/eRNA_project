#A Load libraries,functions,data
setwd ( '/g/furlong/Harnett/TSS_CAGE_myfolder/' )
source ( 'src/tss_cage_functions.R' )
#load the annotation data with scores
load('data/objects/scored_regions.object.R')
#Also ollie's output files
library(rhdf5)
dir.create('analysis/gene_variance')
summaryHDF5 <- '/g/furlong/Harnett/TSS_CAGE_myfolder/data/olliver_processed_summary.hdf5'
dataHDF5 <- '/g/furlong/project/24_TSSCAGE/analysis/FULLDATA_HDF5/dataWithPC1-5qqnorm_PowerLaw.hdf5'
lsHDF5 <- h5ls(summaryHDF5,recursive=F)
lsdataHDF5 <- h5ls(dataHDF5,recursive=T)
head(lsHDF5)

var_decomp_tps = paste0('tp',h5read(dataHDF5,name='phenotype/col_header/dev_stage')[1:3])

#now get the variance data as per ollies email
Venv = h5read(summaryHDF5,name='Venv')
Vcis = h5read(summaryHDF5,name='Vcis')
Vtrans = h5read(summaryHDF5,name='Vtrans')
Vnoise = h5read(summaryHDF5,name='Vnoise')
Cnoise = h5read(summaryHDF5,name='Cnoise')
Ccis = h5read(summaryHDF5,name='Ccis')
Ctrans = h5read(summaryHDF5,name='Ctrans')


#Get the various processed variability components (no individual stages here....)
V = matrix(0,nrow=nrow(Vcis),ncol=5)
V[,1] = Venv
V[,2] = Vcis[,2]
V[,3] = Vcis[,1]-Vcis[,2]
V[,4] = Vtrans[,2]
V[,5] = Vtrans[,1]-Vtrans[,2]

genes.gr$geneDensity = getGeneDensity(genes.gr,w=200000)
#should sum to 1 within numeric tolerance
expect_true(all(identical(Vnoise[i,] +Vcis[i,]+Vtrans[i,]+Venv[] %% 1 < 0.001) ) )
expect_true(all(equals(Vnoise[i,] +Vcis[i,]+Vtrans[i,]+Venv[] ,1) ) )
#Now let's also get a hold of the covariance matrices
Cnoise = h5read(summaryHDF5,name='Cnoise')
Ccis = h5read(summaryHDF5,name='Ccis')
Ctrans = h5read(summaryHDF5,name='Ctrans')
#these are correlation matrices, it seems. Diagonals are therefor the variances 
i=sample(1:13000,size=1);Cnoise[,,i]
#it seems the ordering of the time stages is consistent but not chronological
#It should be, 68hours, 24 hours, 1012hours
colheader = h5read(dataHDF5,name='phenotype/col_header')

h5genes=h5read(summaryHDF5,name='geneID')
h5genes=simplify2array(strsplit(x=h5genes,split='_'))
h5genes=GRanges(h5genes[1,],IRanges(as.numeric(h5genes[2,]),as.numeric(h5genes[3,])),seqinfo=si)
match(h5genes,genes.gr)
#same geneset (YES)

###As a very rough and ready idea - we could just multiply the variances in the diagonal of each
#variance measure by the magnitude as measured in the power law data
#matrix with 
dim(Ccis)
dim(Ctrans)
dim(Cnoise)

cage.tpsum=sapply(tps[1:3],function(tp){rowSums(genes.gr$cg.pl[,cage.tps==tp])})

mat=matrix(1:9,ncol=3)
delta=  abs(row(mat)-col(mat))!=0#logical matrix for extracting diagonals
#Now go through the covariance matrices and subtract to get our stage specific variation
variance.bystage = sapply(simplify=F,list(Ccis,Ctrans,Cnoise),function(x){
	m=t(sapply(1:(dim(x)[3]),function(i){
	mat = x[,,i]
	d = diag(mat)
	d - max(c(d[delta],0))
	}))
	colnames(m)<-var_decomp_tps
	m
})
names(variance.bystage)=c('Cis','Trans','Noise')

str(variance.bystage)

variance.bystage=variance.bystage[,c(2,1,3)]

#now re-arrange to correct time stages

#now multiply those variances by the correct 

mdat = melt(variance.bystage)
unique(mdat$Var2)
unique(mdat$Var1)
head(mdat)
colnames(mdat)
#Now plot 


pdf('analysis/gene_variance/decomp_variance_boxplots.pdf')
qplot(data=mdat,color=factor(Var2),geom='boxplot',facet='L1',
	main='Divergence by Timepoint',x=factor(Var2),y=value,notch=T)+
	scale_x_discrete(name='Timepoint')+
	scale_y_log10(name='Divergence')
	#scale_y_continuous(name='Divergence',limits=c(0,0.1))
#now do teh means

	# pdf('analysis/gene_variance/decomp_variance_boxplots.pdf')
ggplot(data=mdat,aes(color=factor(Var2),x=factor(L1),y=value))+
ggtitle('Divergence by Timepoint')+
	stat_summary(fun.data='mean_cl_boot',geom='errorbar',position='jitter')+
	scale_x_discrete(name='Timepoint')+
	scale_y_log10(name='Divergence')
dev.off()

#checking ont he variance data
#the data may be normalized such that genes have equal variance however. We need to check
dim(Vnoise)
dim(Venv)
dim(Vcis)
dim(Vtrans)
i=1
Vnoise[i,i] +Vcis[i,i]+Vtrans[i,i]+Venv[i]
tmp
Vnoise[i,] +Vcis[i,]+Vtrans[i,]+Venv[]
Vnoise[i,] +Vcis[i,]+Vtrans[i,]+Venv[] == 1
Vnoise[i,] +Vcis[i,]+Vtrans[i,]+Venv[] -> tmp
tmp[13034]==1
is(tmp[13034])
class(tmp[13034])
tmp[13034] == as.numeric(1)
str(tmp[13034])
identical(tmp[13034], 1)
as.integer(tmp[13034]) == as.numeric(1)
tmp[13034] %% 1 



#Now we can do some plotting!
p = qplot(main = '' )
ggsave(filename='analysis/gene_variance/Variance_vs_Gene_Density',plot = p)

#cis variability vs gene density
#for stage specific and non stage specific
#trans variability vs gene density








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
