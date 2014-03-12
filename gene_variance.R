#A Load libraries,functions,data
setwd ( '/g/furlong/Harnett/TSS_CAGE_myfolder/' )
source ( 'src/tss_cage_functions.R' )
#load the annotation data with scores
load('data/objects/scored_regions.object.R')
load('data/objects/gene.annotation.object.R')
load('data/objects/accession.df.object.R')
load('data/objects/cg.mapfilt.pl.object.R')
library(rhdf5)

 #Also ollie's output files
dir.create('analysis/gene_variance')
summaryHDF5 <- '/g/furlong/Harnett/TSS_CAGE_myfolder/data/olliver_processed_summary.hdf5'
dataHDF5 <- '/g/furlong/project/24_TSSCAGE/analysis/FULLDATA_HDF5/dataWithPC1-5qqnorm_PowerLaw.hdf5'
lsHDF5 <- h5ls(summaryHDF5,recursive=F)
lsdataHDF5 <- h5ls(dataHDF5,recursive=F)
lsdataHDF5
head(lsHDF5)

h5ls(summaryHDF5,recursive=F,name='gene')

#pears
#1)libraries,constants,loading data
#2)Loading data from HDF5 file
#3)processing data from HDF5 file
#4)verifying data - correct gene annotations, etc etc, check if covar matrices normalized
#5)attribute variance scores to genes
#6)Plot variance distributions for classes of genes etc etc.



################################################################################
##### Loading Data from HDF5#########################################################################
################################################################################

#now get the variance data as per ollies email
Venv = h5read(summaryHDF5,name='Venv')
Vcis = h5read(summaryHDF5,name='Vcis')
Vtrans = h5read(summaryHDF5,name='Vtrans')
Vnoise = h5read(summaryHDF5,name='Vnoise')
Cnoise = h5read(summaryHDF5,name='Cnoise')
Ccis = h5read(summaryHDF5,name='Ccis')
Ctrans = h5read(summaryHDF5,name='Ctrans')
var_decomp_tps = paste0('tp',h5read(dataHDF5,name='phenotype/col_header/dev_stage')[1:3],'h')

#Get the various processed variability components (no individual stages here....)
V = matrix(0,nrow=ncol(Vcis),ncol=9)
V[,1] = Venv#variance due to environment
V[,2] = Vcis[2,]#shared cis
V[,3] = Vcis[1,]-Vcis[2,]#non shared cis
V[,4] = Vtrans[2,]#shared trans
V[,5] = Vtrans[1,]-Vtrans[2,]#non shared trans
V[,6] = Vnoise[2,]#shared noise
V[,7] = Vnoise[1,]-Vnoise[2,]#non-shared noise
V[,8] = Vcis[1,]+Vtrans[1,]+Vnoise[1,]#non developmental variance
V[,9] = Vcis[1,]+Vtrans[1,] / V[,8]#proportion of non-developmental variance which is due to genetic factors
rownames(V)=h5read(summaryHDF5,name='geneID')


############now look at the variance
mat=matrix(1:9,ncol=3)
delta=  abs(row(mat)-col(mat))!=0#logical matrix for extracting diagonals
gnum=dim(Ccis)[3]
#Now go through the covariance matrices and subtract to get our stage specific variation
var.spec = sapply(simplify=F,list(Ccis,Ctrans,Cnoise),function(x){
		m=t(sapply(1:(dim(x)[3]),function(i){
			mat = x[,,i]
			shared = min(mat)
			d = diag(mat) - shared
	}))
	colnames(m)<-var_decomp_tps
	rownames(m)=h5genes.summary.id
	m
})
#and the shared
var.shared = sapply(simplify=F,list(Ccis,Ctrans,Cnoise),function(x){
			m=t(sapply(1:(dim(x)[3]),function(i){
			mat = x[,,i]
			shared = min(mat)
	}))
	names(m)=h5genes.summary.id
	m
})
#and the combined for various tps
var.all = sapply(simplify=F,list(Ccis,Ctrans,Cnoise),function(x){
	t(sapply(1:(dim(x)[3]),function(i){
		d=diag(x[,,i]) 
		names(d)=var_decomp_tps
		d
	}))
})

var.all




#Plot the actual genetic variability...
qplot(V[,9],geom='density')





#also get the gene IDs
colheader = h5read(dataHDF5,name='phenotype/col_header')
h5genes.id = unique(colheader$geneID)
#note the summary object has different IDs
h5genes.summary.id = h5read(summaryHDF5,name='geneID')
h5genes.gr = sort(convertID2GRange(h5genes.id))


#now calculate absolute mean and sd so we can un-normalize our 
pheno = h5read(dataHDF5,name='phenotype')
dimnames(pheno$matrix) = list(gene.tp=NULL,line=NULL)
pheno.m= melt(pheno$matrix)
gene.sds=aggregate(pheno.m$value ,list(geneID=pheno$col_header$geneID[pheno.m$gene.tp] ), sd )
pheno.means=aggregate(pheno.m$value ,list(geneID=pheno$col_header$geneID[pheno.m$gene.tp] ), mean )

agby=list(dev=pheno$col_header$dev_stage,gene=pheno$col_header$geneID)
gene_dev_table=aggregate(x=rowMeans(pheno$matrix),by=agby,FUN=mean)
#now calculate absolute variances over all timepoints, 
genetic.vars=Vcis[1,]+Vtrans[1,] * pheno.sds$x[match(rownames(V),pheno.sds$geneID)]
names(genetic.vars) = rownames(V)
genetic.coefvar=genetic.vars / pheno.means$x[match(rownames(V),pheno.means$geneID)]



#now calculate absolute coefficients of variation variances over all timepoints, 



#also get the gene IDs
colheader = h5read(dataHDF5,name='phenotype/col_header')
h5genes.id = unique(colheader$geneID)
#note the summary object has different IDs
h5genes.summary.id = h5read(summaryHDF5,name='geneID')
h5genes.gr = sort(convertID2GRange(h5genes.id))
#get the phenotype (cage) info
##Find the matching genes for each of these peaks, then get an ID
h5match = match(h5genes.gr,resize(tss.gr,width=500,fix='center'))
h5.IDs = tss.gr$TrID[h5match]
hasID = !is.na(h5match)
#now print the table that David's scripts need
var.GOtab.df = data.frame(ID=h5.IDs,var = V[,5] + V[,3])
write.table(var.GOtab.df[hasID,],file='analysis/cisTransAllSumTab.txt',col.names=F,row.names=F,quote=F)


#var.spec###
#test system for order IDs works

h5read(summaryHDF5,name='lambda_c')

i=2
Ccis[,,i]

################################################################################
##### Processing Data from HDF5#################################################
################################################################################
#get phenotype data from my own objects, but using the hdf5 locations
# h5genes.gr$cg.pl = getStrandedExprMat(h5genes.gr,cg.pl)

h5genes.gr$cg.pl = getStrandedExprMat(h5genes.gr,cg.pl) 
h5genes.gr$cg.pl = h5genes.gr$cg.pl[,accession.df$tissue=='embryo']#check colnames
#Sum by timepoint
cage.tpsum=sapply(unique(accession.df$time),function(tp){rowSums(h5genes.gr$cg.pl[,accession.df$time[accession.df$tissue=='embryo']==tp])})
rownames(cage.tpsum) = h5genes.gr$id

#now aggregate the phenotype data in the hdf5 object similiarily
#Column haeders to parse the phenotype matrix
agby=list(dev=pheno$col_header$dev_stage,gene=pheno$col_header$geneID)
gene_dev_table=aggregate(x=rowMeans(pheno$matrix),by=agby,FUN=mean)
aggregate(x=,by=agby,FUN=mean)

meanstagetable = dcast(gene_dev_table,gene ~ dev,value.var='x')
#
dim(meanstagetable$gene )
length(h5genes.gr$id)
#
head(meanstagetable$gene )
head(h5genes.gr$id)
#now match the ordering in the hdf5 file to the sorted ordering in the GRange
thmatch = match(meanstagetable$gene, h5genes.gr$id)
#and merge them into the same table
mydat = cage.tpsum[thmatch,]
meanstagetable = meanstagetable[!is.na(thmatch),]
expect_identical( rownames(mydat) ,as.character(meanstagetable$gene) )
meanstagetable = cbind(meanstagetable,mydat) 
head(meanstagetable )
head(mydat )
####Plots comparing our own data and the hdf5 data
#Worryingly little correlation here
pdf('analysis/gene_variance/mydat_hdf5_correlation.pdf')
heatpairs(as.matrix(meanstagetable[,-1]),method='s')
dev.off()

pdf('analysis/gene_variance/mydat_hdf5_correlation_rank.pdf')
heatpairs(apply(as.matrix(meanstagetable[,-1]),2,rank),method='s')
dev.off()



############now look at the variance
mat=matrix(1:9,ncol=3)
delta=  abs(row(mat)-col(mat))!=0#logical matrix for extracting diagonals
gnum=dim(Ccis)[3]
#Now go through the covariance matrices and subtract to get our stage specific variation
var.spec = sapply(simplify=F,list(Ccis,Ctrans,Cnoise),function(x){
		m=t(sapply(1:(dim(x)[3]),function(i){
			mat = x[,,i]
			shared = min(mat)
			d = diag(mat) - shared
	}))
	colnames(m)<-var_decomp_tps
	rownames(m)=h5genes.summary.id
	m
})
var.shared = sapply(simplify=F,list(Ccis,Ctrans,Cnoise),function(x){
			m=t(sapply(1:(dim(x)[3]),function(i){
			mat = x[,,i]
			shared = min(mat)
	}))
	names(m)=h5genes.summary.id
	m
})
var.all = sapply(simplify=F,list(Ccis,Ctrans,Cnoise),function(x){
	t(sapply(1:(dim(x)[3]),function(i){
		d=diag(x[,,i]) 
		names(d)=var_decomp_tps
		d
	}))
})


names(var.all)=c('Cis','Trans','Noise')
names(var.spec)=c('Cis','Trans','Noise')
names(var.shared)=c('Cis','Trans','Noise')

var.spec.m = melt(var.spec)
var.shared.m = melt(var.shared)
var.all.m = melt(var.all)

expect_true(all(var.spec[[1]]>=0))

####set cutoffs on variance data - it seems to be bimodal
cutoffs=sapply(simplify=F,list(spec=var.spec,all=var.all),function(var.df){
	cutoffs=sapply(c('Cis','Trans'),function(ct){
		sapply(var_decomp_tps,function(tp){
			x=lognz(var.df[[ct]][,tp])
			#we need to iterate the cutoff fucntion, it fails sometimes
			for(i in 1:5){ try({cutoff = getBimodalCutoff( x,knum=2,proba=0.5)})} 
			cutoff
		})
	})
})



#create a logical column describing whether a gene is in the high or low variance group
var.spec.m$abovecut <- F
for(ct in c('Cis','Trans')){
	setinds = var.spec.m$L1==ct
	time=as.character(var.spec$Var2)
	var.spec.m$abovecut[setinds] = log10(var.spec.m$value[setinds]) > cutoffs$spec[time,ct]
}

#Now set the below cutoff values to zero
var.spec$value[var.spec$L1 %in% c('Cis','Trans') & var.spec$abovecut == F] <- 0  

#What are the proportions of zeroes in each catagory?
pdf('analysis/gene_variance/vardist_bar_violinplots.pdf')
qplot(geom='violin',data=var.spec.m[var.spec.m$L1!='Noise',],color=L1,x=Var2,y=value,log='y',ylab="variance",xlab='Timepoint')
qplot(geom='violin',data=var.shared.m[,],x=L1,y=value,log='y',ylab="variance")
ggplot(var.spec[var.spec$abovecut==T,], aes(Var2,fill=Var2)) + geom_bar() +
  facet_wrap(~ L1)
dev.off()


#############Produce Plots comparing the variance to the Gene Density
head(var.spec.m)
colnames(var.spec.m) <- c('gene_id','timepoint','value','CisOrTrans','abovecut')
h5genes.gr$genedensity50 = getGeneDensity(h5genes.gr,genes=genes.gr,windowsize=50000)
h5genes.gr$prop_gen_var = V[match(h5genes.gr$id,rownames(V)),9]

#var.spec.m$genedensity50 = h5genes.gr$genedensity50[ match(var.spec.m$gene_id,h5genes.gr$id) ]
#
pdf('analysis/gene_variance/density_vs_variance_heatscatter.pdf')
heatscatter(x=h5genes.gr$genedensity50,y=h5genes.gr$prop_gen_var,main='Gene Density vs proportional genetic variance',xlab='gene density 50kb',ylab='ratio of genetic variance to noise')
dev.off()

for(tp in var_decomp_tps){
	varmeasure = #Cis/Trans/Both,spec/shared/all
	maintit=paste0('Gene Density vs', varname , tp )
	heatscatter(maintit,x=gene_density,y=varmeasure) 

}
dev.off()
var.spec


# qplot(geom='density',log='x',x=rowSums(var.spec$Cis),main='Distribution of Cis Variance Scores')
	# qplot(geom='density',log='x',x=rowSums(var.spec$Trans),main='Distribution of Trans Variance Scores')
	# qplot(geom='density',log='x',x=rowSums(var.spec$Noise),main='Distribution of Noise Variance Scores')







d = structure(list(category = structure(c(2L, 2L, 1L, 1L, 1L, 1L, 
1L, 1L, 1L, 3L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 3L, 3L, 3L, 3L, 
3L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L), .Label = c("4X4", 
"HATCH", "SEDAN"), class = "factor"), L1 = structure(c(1L, 
1L, 1L, 1L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 
3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 
3L, 3L, 3L, 3L, 3L, 3L, 3L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 1L, 1L, 
1L), .Label = c("FIXED", "FREE", "MOBILE"), class = "factor"), 
    value = c(6440.32, 6287.22, 9324, 7532, 7287.63, 6827.27, 
    6880.48, 7795.15, 7042.51, 2708.41, 1373.69, 6742.87, 7692.65, 
    7692.65, 8116.56, 7692.65, 7692.65, 7692.65, 7962.65, 8116.56, 
    5691.12, 2434, 8343, 7727.73, 7692.65, 7721.15, 1944.38, 
    6044.23, 8633.65, 7692.65, 7692.65, 8151.65, 7692.65, 7692.65, 
    2708.41, 3271.45, 3333.82, 1257.48, 6223.13, 7692.65, 6955.46, 
    7115.46, 7115.46, 7115.46, 7115.46, 6955.46, 7615.46, 2621.21, 
    2621.21, 445.61)), .Names = c("category", "L1", "value"
), class = "data.frame", row.names = c(NA, -50L))

library(ggplot2)

p = ggplot(d, aes(x=d$value, fill=d$category)) + geom_density(alpha=.3)
p + facet_grid(d$L1 ~ .)


plot(1:3)
qplot(data=d, x=value, fill=category,geom='density',alpha=.3) +facet_grid(L1 ~ .)



pdf('analysis/gene_variance/gene_density_vs_var_plots.jpeg')

ggplot(data=var.spec,aes(x=genedensity50,y=log10(value))) + 
	facet_grid(CisOrTrans ~ .) +
 	stat_density2d(aes(fill=..level..), geom="polygon",h=0.5) +
  scale_fill_gradient(low="blue", high="white")
dev.off()
for(ct in unique(var.spec$L1)){
	d=var.spec[var.spec$L1==ct,]
	heatscatter(main=paste0(ct, ' vs Gene Density'),x=d$genedensity50,y=log10(d$value))

}

#


head(var.spec[var.spec$L1=='Trans',])
varstagetable = var.spec[var.spec$L1=='Trans',]
varstagetable = dcast(varstagetable, Var1 ~ Var2)
colnames(varstagetable) <- paste('Variance',colnames(varstagetable))
var.mean.table = as.matrix(cbind(varstagetable[,-1],meanstagetable[sum2datmatch,-c(1,5,6,7)]))
heatpairs(main = 'Magnitude Vs. Trans Variance, decomp data',method='s',var.mean.table)
#








colnames(var.spec) = c('Gene','Timepoint','Value','VarianceType')
head(var.spec)
var.spec$Timepoint = factor(var.spec$Timepoint,levels=tps)
#Now plot 
pdf('analysis/gene_variance/decomp_variance_boxplots.pdf')
qplot(data=var.spec,color=factor(Timepoint),geom='boxplot',
	main='Divergence by Timepoint',x=factor(Timepoint),y=Value,notch=T)+
	scale_x_discrete(name='Timepoint')+
	scale_y_log10(name='Divergence')+
	facet_wrap(~VarianceType)
dev.off()
	#scale_y_continuous(name='Divergence',limits=c(0,0.1))
#now do teh means

#Now the mean variances at each time split by timepoint
pdf('analysis/gene_variance/decomp_variance_boxplots.pdf')
print(
	ggplot(data=var.spec,aes(color=factor(Timepoint),x=factor(Timepoint),y=Value))+
		ggtitle('Divergence by Timepoint')+
		stat_summary(fun.data='mean_cl_boot',geom='errorbar')+
		scale_x_discrete(name='Timepoint')+
		scale_y_log10(name='Divergence')+
		facet_wrap(~VarianceType)
	)
dev.off()






genes.gr$geneDensity = getGeneDensity(genes.gr,w=200000)
genes.gr$geneDensity = getGeneDensity(genes.gr,w=200000)

