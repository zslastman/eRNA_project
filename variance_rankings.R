#A Load libraries,functions,data
setwd ( '/g/furlong/Harnett/TSS_CAGE_myfolder/' )
source ( 'src/tss_cage_functions.R' )
#load the annotation data with scores
load('data/objects/scored_regions.object.R')
load('data/objects/gene.annotation.object.R')
load('data/objects/accession.df.object.R')
#load('data/objects/cg.mapfilt.pl.object.R')
library(rhdf5)

#Ollie's output files
outfolder = 'analysis/variance_rankings' 
dir.create(outfolder)
summaryHDF5 <- '/g/furlong/Harnett/TSS_CAGE_myfolder/data/olliver_processed_summary.hdf5'
dataHDF5 <- '/g/furlong/project/24_TSSCAGE/analysis/FULLDATA_HDF5/dataWithPC1-5qqnorm_PowerLaw.hdf5'
lsHDF5 <- h5ls(summaryHDF5,recursive=F)
lsdataHDF5 <- h5ls(dataHDF5,recursive=F)


####################################
#and aggregate the data from the HDF data file
####################################
var_decomp_tps = paste0('tp',h5read(dataHDF5,name='phenotype/col_header/dev_stage')[1:3],'h')
colheader = h5read(dataHDF5,name='phenotype/col_header')
h5.data.genes.id = unique(colheader$geneID)#original order - aggregate files in this order
h5datgenes.gr = sort(convertID2GRange(h5.data.genes.id))#sorted GRanges object
#now calculate absolute mean and sd so we can un-normalize our variance data
pheno = h5read(dataHDF5,name='phenotype')#get phenotype data
pheno$col_header$dev_stage=paste0('tp',pheno$col_header$dev_stage,'h')#convert to my format (letter in fron tavoids problems)
dimnames(pheno$matrix ) = list(gene.tp=NULL,line=NULL)
pheno.m= melt(pheno$matrix)
agbygene=list(geneID=pheno$col_header$geneID[pheno.m$gene.tp] )
agbygenedev=list(geneID=pheno$col_header$geneID[pheno.m$gene.tp],dev= pheno$col_header$dev_stage[pheno.m$gene.tp] )
gene.sds=aggregate(pheno.m$value ,agbygene, sd )
gene.means=aggregate(pheno.m$value ,agbygene, mean )
gene.dev.sds=aggregate(pheno.m$value ,agbygenedev, sd )
gene.dev.means=aggregate(pheno.m$value ,agbygenedev, mean )
#and also get the means for each timestage
meanstagetable = dcast(gene.dev.means,geneID ~ dev,value.var='x')
sdstagetable = dcast(gene.dev.sds,geneID ~ dev,value.var='x')#for rescaling to within timestage sd
inc.dev.sds = gene.sds$x[sumord]#for unscaling the cross timestage sd
inc.dev.means = gene.means$x[sumord]#for unscaling the cross timestage sd
#now match the peaks to actual genes.
dn=distanceToNearest(h5sumgenes.gr,genes.gr)
h5sumgenes.gr$Gene=genes.gr[dn$subjectHits]$id
h5sumgenes.gr$Gene[dn$distance>500] <- NA
#indexes for re-ordering the data h5 genes to the UNSORTED summary h5 order - 
sumord = match(rownames(V),gene.sds$geneID)#order genes origianlly in 
sumord = sumord[!is.na(sumord)]


####################################
#Now get variance decomp data from the summary file
####################################
#read summary object
#note the summary object has different IDs
h5genes.summary.id = h5read(summaryHDF5,name='geneID')#original order
h5sumgenes.gr = sort(convertID2GRange(h5genes.summary.id))#sorted GRanges object
Venv = h5read(summaryHDF5,name='Venv')
Vcis = h5read(summaryHDF5,name='Vcis')
Vtrans = h5read(summaryHDF5,name='Vtrans')
Vnoise = h5read(summaryHDF5,name='Vnoise')
Cnoise = h5read(summaryHDF5,name='Cnoise')
Ccis = h5read(summaryHDF5,name='Ccis')
Ctrans = h5read(summaryHDF5,name='Ctrans')
var_decomp_tps = paste0('tp',h5read(dataHDF5,name='phenotype/col_header/dev_stage')[1:3],'h')
#Variance across all timestages...
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
rownames(V)=h5genes.summary.id
############now look at the variance
mat=matrix(1:9,ncol=3)
delta=  abs(row(mat)-col(mat))!=0#logical matrix for extracting diagonals
gnum=dim(Ccis)[3]#number of genes in summary object
#Now go through the covariance matrices and subtract to get our stage specific variation
var.spec = setNames(nm=c('Cis','Trans','Noise'),
    sapply(simplify=F,list(Ccis,Ctrans,Cnoise),function(x){
      m=t(sapply(1:(dim(x)[3]),function(i){
        mat = x[,,i]
        shared = min(mat)
        d = diag(mat) - shared
    }))
    colnames(m)<-var_decomp_tps
    rownames(m)=h5genes.summary.id
    m
  })
)
#and the shared
var.shared = setNames(nm=c('Cis','Trans','Noise'),
    m=sapply(simplify=F,list(Ccis,Ctrans,Cnoise),function(x){
        m=t(sapply(1:(dim(x)[3]),function(i){
        mat = x[,,i]
        shared = min(mat)
    }))
    names(m)=h5genes.summary.id
    m
  })
)

#and the combined for various tps
var.all = setNames(nm=c('Cis','Trans','Noise'),
    sapply(simplify=F,list(Ccis,Ctrans,Cnoise),function(x){
    m=t(sapply(1:(dim(x)[3]),function(i){
      d=diag(x[,,i]) 
      names(d)=var_decomp_tps
      d
    }))
    rownames(m)=h5genes.summary.id
    m
  })
)
#now, we need to scale the variances by multiplying them by the sd for all stages, and then dividing by the
#sd for that stage
var.all.rescaled = var.all
for(s in 1:length(var.all) ){
  for(tp in colnames(var.all[[s]]))  {
    var.all.rescaled[[s]][,tp] = (var.all[[s]][,tp] * gene.sds[sumord,]$x ) / sdstagetable[sumord,tp]
  }
}
####The variance data is bimodal - set cutoffs on it to define two groups
cutoffs=sapply(simplify=F,list(spec=var.spec,all=var.all),function(var.df){
  cutoffs=sapply(c('Cis','Trans'),function(ct){
    sapply(var_decomp_tps,function(tp){
      x=lognz(var.df[[ct]][,tp])
      #we need to iterate the cutoff fucntion, it fails sometimes
      for(i in 1:10){ try({cutoff = getBimodalCutoff( x,knum=2,proba=0.5)})} 
      cutoff
    })
  })
})

var.spec.m = melt(var.spec)
var.shared.m = melt(var.shared)
var.all.m = melt(var.all)
head(var.spec.m)
expect_true(all(var.spec[[1]]>=0))

#create a logical column describing whether a gene is in the high or low variance group
var.spec.m$abovecut <- F
for(ct in c('Cis','Trans')){
  setinds = var.spec.m$L1==ct
  time=as.character(var.spec.m$Var2)
  var.spec.m$abovecut[setinds] = log10(var.spec.m$value[setinds]) > cutoffs$spec[time,ct]
}
#create a logical column describing whether a gene is in the high or low variance group
var.all.m$abovecut <- F
for(ct in c('Cis','Trans')){
  setinds = var.all.m$L1==ct
  time=as.character(var.all.m$Var2)
  var.all.m$abovecut[setinds] = log10(var.all.m$value[setinds]) > cutoffs$all[time,ct]
}








###############################################################################
#Now fit a linear model to the data using Limma
##################################################
#Do this with only a selection of the most expressed genes.
#Genes which have reasonably high expression at all timestages
#The top 100 of these
rs=rowSums(meanstagetable[,2:4])
rv=apply(meanstagetable[,2:4],1,sd)
inds=which(rs>quantile(rs,0.2) & rv < quantile(rv,0.2))#consistently high peaks
mypeaks = h5datgenes.gr[]
#or if we're just using all of them
mypeaks = h5sumgenes.gr
my.df = accession.df[accession.df$tissue=='embryo' & accession.df$time=='68h',]

load('data/objects/cg.mapfilt.pl.object.R')
peakcounts.mapfilt=getStrandedExprMat(gr=mypeaks,cg.mapfilt.pl[my.df$accession])
colnames(peakcounts.mapfilt)=my.df$accession
row.names(peakcounts.mapfilt)=mypeaks$id#Get the expression matrix for these
save(peakcounts.mapfilt,file='data/objects/peakcounts.mapfilt.object.R')

# load('data/objects/cg.pl.object.R')
# peakcounts.nomapfilt=getStrandedExprMat(gr=mypeaks,cg.pl[my.df$accession])
# colnames(peakcounts.nomapfilt)=my.df$accession
# row.names(peakcounts.nomapfilt)=mypeaks$id#Get the expression matrix for these
# save(peakcounts.nomapfilt,file='data/objects/peakcounts.object.R')

# peakcounts=getStrandedExprMat(gr=mypeaks,cg[my.df$accession])
# colnames(peakcounts)=my.df$accession
# row.names(peakcounts)=mypeaks$id#Get the expression matrix for these
# save(peakcounts,file='data/objects/peakcounts.object.R')

peakcounts = peakcounts.mapfilt
expect_is(peakcounts,'matrix')


cage.tps = accession.df$timepoint[ match(colnames(peakcounts),accession.df$accession)]
cage.tpsum=sapply(unique(cage.tps),function(tp){rowSums(peakcounts[,cage.tps==tp])})

#1a)Remove mesodermal only stuff,Remove the 'bad libraries' mentioned by ignacio
samplegroups=c('timepoint','line','collection','prep','seq')
reptable=my.df[!my.df$isbadsample & !my.df$tissue=='meso',c('accession',samplegroups)]
reptable <- data.frame(lapply(reptable, as.character), stringsAsFactors=T)
# #Include only lines with 2 in each timepoint
# hasmultiple.collections=sapply(unique(reptable$line),function(line){
#   length(unique(reptable[reptable$line==line,]$collection))>1
# })
# reptable=reptable[ hasmultiple.collections[as.character(reptable$line)] ,]

#and fit leaner model
require(limma)
require(edgeR)
peakcounts=peakcounts[,as.character(reptable$accession)]
peakcounts.ob <- DGEList(counts=peakcounts, genes=rownames(peakcounts))
peakcounts.ob = calcNormFactors(peakcounts.ob)
v <- voom(peakcounts.ob,model.matrix( ~ reptable$line ),plot=F)
fit <- lmFit(v,design)
vfit <- eBayes(fit,trend=T)

#Now scale the peakcounts object in the same way as our variance decomp pipeline does.
sds = apply(peakcounts,1,sd)
means = apply(peakcounts,1,mean)
scaled.peakcounts = peakcounts
for(i in 1:nrow(peakcounts)){
  scaled.peakcounts[i,] = (peakcounts[i,] - means[i]) / sds[i]
}
#and fit linear model
peakcounts.ob <- DGEList(counts=scaled.peakcounts, genes=rownames(scaled.peakcounts))
peakcounts.ob = calcNormFactors(peakcounts.ob)
v <- voom(peakcounts.ob,model.matrix( ~ reptable$line ),plot=F)
fit <- lmFit(v,design)
vfit.scaled <- eBayes(fit,trend=T)
#this design matrix is too big I think. 

#Fit a random effects model
require(nlme)
n=20
inds=order(rowSums(peakcounts),decreasing=T)[1:n]#best rows of peakcounts
peakcounts.filt = peakcounts[inds,]
peakcounts.filt.m = melt(peakcounts.filt)
colnames(peakcounts.filt.m) = c('Gene','accession','value')
peakcounts.filt.m$collection = factor(reptable$collection[match(peakcounts.filt.m$accession,reptable$accession)])
peakcounts.filt.m$line = factor(reptable$line[match(peakcounts.filt.m$accession,reptable$accession)])
peakcounts.filt.m$value = as.numeric(as.character( peakcounts.filt.m$value))

head(peakcounts.filt)
head(peakcounts.filt.m,n=100)
lmfit=lm(data=peakcounts.filt.m,value ~ Gene * line)
lmefit = lme(data=peakcounts.filt.m,fixed = value ~ Gene * line , random = ~ collection*Gene)


###############################################################################
#Plots showing correspondence between Limma and variance decomposition
##################################################
# Comparing the data based on the input in these objects to the data in mine
load('data/objects/peakcounts.object.R')
peakcounts = peakcounts.mapfilt
#summary data for timepoints
cage.tps = accession.df$timepoint[ match(colnames(peakcounts),accession.df$accession)]
cage.tpsum=sapply(unique(cage.tps),function(tp){rowSums(peakcounts[,cage.tps==tp])})
sumord2 = match(rownames(peakcounts,meanstagetable$geneID)#hdf5 ordering to peakcounts
h5cagemean68 = meanstagetable$tp68h[sumord2] 

#plots showing summary of input and output from limma
pdf(paste0(outfolder,'/','limma_input_output.pdf'))
#Summary of variance from Limma,
qplot(main='Limma output',limmavar,geom='density',log='',fill='blue',xlab='genetic var from limma fit')
qplot(main='Limma output',limmavar,geom='density',log='x',fill='blue',xlab='genetic var from limma fit')
# compared to mean
heatscatter(x= log10( apply(vfit$coef,1,var) ) ,y= log10( cage.tpsum[,] ) ,log='',main='Calculated Genetic Variance from Limma vs. Kinship Decomp 68h',xlab = 'Log10 Variance in Coeffs' , ylab = 'Log10 Mean')
# for the scaled data
heatscatter(x= log10( apply(vfit.scaled$coef,1,var) ) ,y= log10(cage.tpsum[,]) ,log='',main='Calculated Genetic Variance from Limma with scaling vs. Kinship Decomp 68h',xlab = 'Log10 Variance in Coeffs' , ylab = 'Log10 Mean')
# Now dividing our variance measure by the mean
heatscatter(x= log10( apply(vfit$coef,1,var) / cage.tpsum[,] ) ,y= log10(cage.tpsum[,] ) ,log='',main='Calculated Genetic Variance from Limma vs. Kinship Decomp 68h',xlab = 'Log10 Variance in Coeffs' , ylab = 'Log10 Mean')
# for the scaled data
heatscatter(x= log10( apply(vfit.scaled$coef,1,var) / cage.tpsum[,] ) ,y= log10( cage.tpsum[,] ) ,log='',main='Calculated Genetic Variance from Limma with scaling vs. Kinship Decomp 68h',xlab = 'Log10 Variance in Coeffs' , ylab = 'Log10 Mean')
#using the stduncaled in the object
______

dev.off()

###############################################################################
#plots showing summary of input and output from variance Decomp
###############################################################################
prop.genetic = (var.all$Cis+var.all$Trans) / var.all$Cis+var.all$Trans+var.all$Noise
abs.genetic = var.all$Cis+var.all$Trans
pdf( paste0(outfolder,'/','hd5_input_output.pdf') )
#First - I should be able to calculate the proportion of the variance which is due to dev stages myself...
#calculate variance amongst stages - approximately the variance between the means at each stage
MyStageVar = setNames(nm= meanstagetable[,1] ,apply(meanstagetable[,-1],1,var) )
#proportion vairance due to stages
MyStageVar = MyStageVar / (gene.sds^2)
#now get value in summary hdf5
H5StageVar=as.vector(Venv)
#and order the dataHDF5 stuff correctly
MyStageVar = MyStageVar[ sumord ]
#Compare values in scatterplot
heatscatter(x= log10( MyStageVar ) ,y= log10( H5StageVar ) ,method='s',log='' ,main='Comparison of my estimated variance due to stages\n and the H5 estimate',xlab = 'Log10 MyStageVar' , ylab = 'Log10 H5StageVar')
#done
dev.off()
#Now look at distribution of 

#different tps and different sources
#Basic plotting of variance data,shows bimodality
pdf(file='analysis/gene_variance/VarDensity_stagespec_Densityplot.pdf')
  
  for(tp in var_decomp_tps){
    transvar = var.spec$Trans[,tp]
    cisvar = var.spec$Cis[,tp]
    maintit = paste0(tp,' Distribution of Specific Cis/Trans Variance Scores')
    print(qplot(geom='point',alpha=0.1,log='xy',y=transvar ,x=cisvar,main=maintit) +
    geom_hline(yintercept=10^ cutoffs$spec[tp,'Trans'] )+
    geom_vline(xintercept=10^ cutoffs$spec[tp,'Cis']) )
  }
  
  for(tp in var_decomp_tps){
    transvar = var.all$Trans[,tp]
    cisvar = var.all$Cis[,tp]
    maintit = paste0(tp,' Distribution of Total Cis/Trans Variance Scores')
    print(qplot(geom='point',alpha=0.1,log='xy',y=transvar ,x=cisvar,main=maintit) +
    geom_hline(yintercept=10^cutoffs$all[tp,'Trans'])+
    geom_vline(xintercept=10^cutoffs$all[tp,'Cis']))
  }

  transvar = var.shared$Trans
  cisvar = var.shared$Cis
  maintit = paste0(' Distribution of shared Cis/Trans Variance Scores')
  print(qplot(geom='point',alpha=0.1,log='xy',y=transvar ,x=cisvar,main=maintit))
  qplot(geom='density',data=var.shared.m,log='x',x=as.numeric(value),color=factor(L1),main='Stage Specific Variance Sources') 
  qplot(geom='density',data=var.spec.m,log='x',x=as.numeric(value),color=factor(L1),main='Shared Specific Variance Sources')
dev.off()

#mean vs variance
#Proportions using cutoffs
#What are the proportions of zeroes in each catagory?
pdf('analysis/gene_variance/vardist_bar_violinplots.pdf')
qplot(geom='violin',data=var.spec.m[var.spec.m$L1!='Noise',],color=L1,x=Var2,y=value,log='y',ylab="variance",xlab='Timepoint')
qplot(geom='violin',data=var.shared.m[,],x=L1,y=value,log='y',ylab="variance")
ggplot(var.spec[var.spec$abovecut==T,], aes(Var2,fill=Var2)) + geom_bar() +
  facet_wrap(~ L1)
dev.off()

dev.off()

###############################################################################
#plots howing comparisons between variance decomp and limma
###############################################################################
pdf(paste0(outfolder,'/','hd5_input_output.pdf'))
#Plots showing - data comparisons - do they match up
heatscatter(x= log10( cage.tpsum[,] ) ,y= log10( meanstagetable$tp68h[sumord2] ) ,log='',main='My input to LImma vs. Input from data hdf5',xlab = 'Log10 Cage - mine' , ylab = 'Log10 Cage - hdf5')
#plots showing variances - do they match up



decompvar = abs.genetic[,var_decomp_tps[1]]
names(limmavar) = rownames(fit$coefficients)
limmavar = limmavar[ match( names(decompvar),names(limmavar)) ]

qplot(
  main = 'Proportional Genetic Variance According to linear Model Vs. Proportional Genetic variance\n According to kinship matrix Decomposition',
  x=(limmavar),xlab='Genetic Variance - Limma',
  y=(decompvar),ylab='Genetic Variance - Kinship Decomposition ',
  log='xy'
  geom='point'
)

dev.off()










tp = var_decomp_tps[[1]]
repped_coeffs = colnames(fitob$coeff)









