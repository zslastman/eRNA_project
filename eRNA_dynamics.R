#script examining changes in eRNA expression, and whether they dcoorelate with changes in enhancer activity
#1 load functions,datasets
setwd('/g/furlong/Harnett/TSS_CAGE_myfolder/')
source('src/tss_cage_functions.R')
source('src/generate_chromdata_cage.R')
load('data/objects/cg.pl.map.object.R')
load('data/objects/accession.df.full.object.R')
load(file.crm8008.gr)
load(file.cad3)
load(file.tss)

# for(i in seq_along(cg.pl.map)){cg.pl.map[[i]]$both = cg.pl.map[[i]]$pos  +cg.pl.map[[i]]$neg  }
# #exlcude our mesodermal libraries
 accession.df=accession.df[accession.df$tissue!='meso',]
 tps=accesion.df$timepoint

# cg=simplify2array(cg.pl.map)
# views=mclapply(mc.cores=10,cg['both',],function(srle){
# 	unlist(viewApply(GRViews(srle[[1]],gr),function(x)x))	#list of Rles
# })

# cad3.gr$cadcage = getBestWindowMat(cad3.gr,500,cg.pl.map)

#for now just use the cagecountmatlist thing from oldscript. Ugly but fine..
load('cagecountmatlist.object.R')
cad3.gr$cadcage = do.call(cbind,cagecountmatlist[c(2,5,8)])
accs = accession.df$acc[accession.df$acc%in% colnames(cad3.gr$cadcage)]
cad3.gr$cadcage = cad3.gr$cadcage[,accs]
tps=accession.df$timepoint[match(accs,accession.df$acc)]

#First we can do plots of fold change in total eRNA expression at the two timepoints
#first calculate mean for each cad3.gr
cad3.gr$tp24mean = rowMeans(cad3.gr$cadcage[,tps=='24h'])
cad3.gr$tp68mean = rowMeans(cad3.gr$cadcage[,tps=='68h'])
cad3.gr$tp1012mean = rowMeans(cad3.gr$cadcage[,tps=='1012h'])
cad3.gr$tp24tp68.lr = log(cad3.gr$tp68mean) - log(cad3.gr$tp24mean)
cad3.gr$tp68tp1012.lr = log(cad3.gr$tp1012mean) - log(cad3.gr$tp68mean)
cad3.gr$tp24tp68.la = 0.5*(log(cad3.gr$tp68mean) + log(cad3.gr$tp24mean))
cad3.gr$tp68tp1012.la = 0.5*(log(cad3.gr$tp1012mean) + log(cad3.gr$tp68mean))

#function to create a single factor from combinations
factorFromLogicals<-function(logl){
	#take in list of logical factors
	# logl=list(A=c(1,0,0,1,1,0),B=c(0,1,0,1,1,0),C=c(0,0,1,0,1,0))
	# logl=lapply(logl,as.logical)
	stopifnot(!is.null(names(logl)))#all named
	stopifnot(sapply(logl,is.logical))#logical
	R=length(logl)
	stopifnot(R<5)#logical
	perms=permutations(2,R,0:1,repeats=T)
	perms=lapply(1:nrow(perms),function(i){as.logical(perms[i,])})
	factvect=rep('',length(logl[[1]]))
	#get each binary combinatin of indices
#use these indices to create a list of names
	for(perm in perms){
		permname=paste(collapse='.',names(logl)[perm])
		log.vect=Reduce('+',logl[perm])==sum(perm)
		factvect[log.vect]=permname

	}
	factvect
#use these indices to combine logical vectors and insert name into vector
#now turn the whole thing into a factor, with levels ordered as the original names
}

activity_pattern = factorFromLogicals(list('Active.at.24.hrs'=cad3.gr$active24,'Active.at.68hrs'=cad3.gr$active68))
activity_pattern = factor(activity_pattern)
activity_pattern

#calculate fold change from 2-4, 4-8hrs, plot it as a function of mean

#now split our crms by their activity pattern and show fold change

library(ggplot2)
dir.create('analysis/eRNA_dynamics')

pdf('analysis/eRNA_dynamics/maplot2468.pdf')
print(
	qplot( colour= activity_pattern,x = cad3.gr$tp24tp68.la ,y=cad3.gr$tp24tp68.lr,
	main='MA plot for Cad eRNA - summed normalized Cage\n 2-4hrs vs 6-8hrs',geom='point',ylim=c(-5,5),ylab='M',xlab='A')+
	scale_color_manual(values=c('black','blue','red','green'))
)
dev.off()


pdf('analysis/eRNA_dynamics/densities2468.pdf',w=14,h=7)
print(
	qplot( fill= activity_pattern,alpha=0.3,x = cad3.gr$tp24tp68.lr,
	main='Density plot Cad eRNA - fold change in summed normalized Cage\n 2-4hrs vs 6-8hrs',
	geom='density',xlab='Fold change from 2-4 to 6-8 hrs')+
	scale_color_manual(values=c('black','blue','red','green'))
)
dev.off()

activity_pattern = factorFromLogicals(list('Active.at.68.hrs'=cad3.gr$active68,'Active.at.1012hrs'=cad3.gr$active1012))
activity_pattern = factor(activity_pattern)
activity_pattern

pdf('analysis/eRNA_dynamics/maplot681012.pdf')
print(
	qplot( colour= activity_pattern,x = cad3.gr$tp68tp1012.la ,y=cad3.gr$tp68tp1012.lr,
	main='MA plot for Cad eRNA - summed normalized Cage\n 6-8hrs vs 10-12hrs',geom='point',ylim=c(-5,5),ylab='M',xlab='A')+
	scale_color_manual(values=c('black','blue','red','green'))
)
dev.off()

pdf('analysis/eRNA_dynamics/densities681012.pdf')
print(
	qplot( colour= activity_pattern,x = cad3.gr$tp68tp1012.lr,
	main='Density plot Cad eRNA - fold change in summed normalized Cage\n 6-8hrs vs 10-12hrs',
	geom='density')+
	scale_color_manual(values=c('black','blue','red','green'))
)
dev.off()

#Now let's address what charles and ignacio said - given that a gene has pollII at one stage, but not K27ac, 
#does it tend to then have K27ac at the next stage?

#If A comes before B, then we see that A at 1 predicts B at 2 independantly of B at 1.
#on the other hand if B comes before A we don't see that
#

#What predicts 

h5createFile
