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
head(gene.dev.sds)
#and also get the means for each timestage
sumord = match(rownames(V),gene.sds$geneID)#order genes origianlly in 
sumord = sumord[!is.na(sumord)]
all(gene.sds$geneID[sumord]==rownames(V))
meanstagetable = dcast(gene.dev.means,geneID ~ dev,value.var='x')
sdstagetable = dcast(gene.dev.sds,geneID ~ dev,value.var='x')#for rescaling to within timestage sd
head(sdstagetable)
inc.dev.sds = gene.sds$x[sumord]#for unscaling the cross timestage sd
inc.dev.means = gene.means$x[sumord]#for unscaling the cross timestage sd
identical( gene.sds$geneID[sumord],rownames(V)[])#for unscaling the cross timestage sd
#now match the peaks to actual genes.
dn=distanceToNearest(h5sumgenes.gr,genes.gr)
h5sumgenes.gr$Gene=genes.gr[dn$subjectHits]$id
h5sumgenes.gr$Gene[dn$distance>500] <- NA




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
for(s in 1:length(var.all) ){
  for(tp in colnames(var.all[[s]]))  {
    var.all[[s]][,tp] = (var.all[[s]][,tp] * gene.sds[sumord,]$x ) / sdstagetable[sumord,tp]
  }
}



###############################################################################
#Now calculate the technical noise using the information from replicates and limma
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

load('data/objects/cg.pl.object.R')
peakcounts.nomapfilt=getStrandedExprMat(gr=mypeaks,cg.pl[my.df$accession])
colnames(peakcounts.nomapfilt)=my.df$accession
row.names(peakcounts.nomapfilt)=mypeaks$id#Get the expression matrix for these
save(peakcounts.nomapfilt,file='data/objects/peakcounts.object.R')

peakcounts = peakcounts.mapfilt

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
#Calculate the proportion of variation which is genetic for each timepoint
prop.genetic = (var.all$Cis+var.all$Trans) / var.all$Cis+var.all$Trans+var.all$Noise
abs.genetic = var.all$Cis+var.all$Trans
str(prop.genetic)
pdf(paste0(outfolder,'/','hd5_input_output.pdf'))
#summary of variance
#all tps
#different tps and different sources
#mean vs variance
#Proportions using cutoffs

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






###############################################################################
#Plots comparing individual features to our variances
###############################################################################
load('data/objects/scored.annotation.R')
load('data/objects/tfgrlist.object.R')
load('data/objects/motif.overlap.object.R')#get motif overlap
load('data/objects/motifs.gr.object.R')
load('data/objects/dnase.peaks.object.R')



# functions for boxplot labelling
give.n <- function(x){
  return(c(y = median(x)*1.05, label = length(x))) 
  # experiment with the multiplier to find the perfect position
}
mean.n <- function(x){
  return(c(y = median(x)*0.97, label = round(mean(x),2))) 
  # experiment with the multiplier to find the perfect position
}

genetic.vars = list(
  'Proportional' = (var.all$Cis[,tp] + var.all$Trans[,tp]) / var.all$Noise[,tp],
  'Cis' = var.all$Cis[,tp] ,
  'Trans' = var.all$Trans[,tp] ,
  'Genetic' = var.all$Cis[,tp] + var.all$Trans[,tp],
  'Non developmental' = var.all$Cis[,tp] + var.all$Trans[,tp] + var.all$Noise[,tp],
  'Limma' = apply(fit$coeff,1,sd)#limma only uses replicated ones anyway...
  # 'Limma_replicated_only' = apply(fit$coeff[,repped_coeffs],1,sd)
)
windowsizes = c('50kb'=50000,'25kb'=25000,'12.5kb'=12500)

gn = 'Limma'
w=windowsizes[[2]]

for(gn in names(genetic.vars)){
genetic.var = genetic.vars[[gn]]


pdf(paste0('analysis/variance_rankings/indiv_features_vs_Var_',gn,' w=',w,'_.pdf'))
# pdf('analysis/variance_rankings/density_vs_variance_heatscatter.pdf')
h5sumgenes.gr$genedensity50 = getGeneDensity(h5sumgenes.gr,genes=genes.gr,windowsize=w)
heatscatter(x=h5sumgenes.gr$genedensity50,y=log10(genetic.var),main=paste0('Gene Density vs ',gn,' windowsize ',w),xlab='gene density',ylab=paste0('Log10',gn))
# dev.off()
# Nearby TF peaks
h5sumgenes.gr$tfpeakdensity = rowSums(sapply( tfgrlist[['6-8']] , function(tfpeaks.gr){
  getGeneDensity(h5sumgenes.gr,genes=tfpeaks.gr,windowsize=w)
}))
# pdf('analysis/variance_rankings/tfpeakdensity_vs_variance_heatscatter.pdf')
heatscatter(x=h5sumgenes.gr$tfpeakdensity,y=log10(genetic.var),main=paste0('Local big5 TF peak density vs ',gn,' windowsize ',w),xlab='tf peak density',method='s',ylab=gn)
# dev.off()
# Nearby Dnase Peaks
# pdf('analysis/variance_rankings/tfpeakdensity_vs_variance_heatscatter.pdf')
# load('data/objects/dnase.peaks.object.R')
h5sumgenes.gr$dnasepeakdensity = getGeneDensity(h5sumgenes.gr,genes=dnase.peaks,windowsize=w)
heatscatter(x=h5sumgenes.gr$dnasepeakdensity,y=log10(genetic.var),main=paste0('Local Dnase Peak Density vs ',gn,' windowsize ',w),xlab='Dnase Peak density',method='s',ylab=gn)
# dev.off()
dev.off()

allmotifs.gr<-import('/g/furlong/Harnett/TSS_CAGE_myfolder/data/fimo_out_JASPAR_CORE_insects_02/fimo.gff',asRangedData=F)
allmotifs.gr = allmotifs.gr[countOverlaps(allmotifs.gr,dnase.peaks)>0]

# Promotor Type
gene.motif.overlap = matrix(0,ncol=ncol(motif.overlap$tss.gr),nrow=length(genes.gr))
colnames(gene.motif.overlap) = colnames(motif.overlap$tss.gr)
rownames(gene.motif.overlap) = genes.gr$id
# the below code sums all tss for the gene
for(i in length(tss.gr$Gene)){
  tss.ov = motif.overlap$tss.gr[i,]#get the motifs overlapping that gene
  gid = tss.gr$Gene[i]#get the corresponding gene
  gene.motif.overlap[gid,]=gene.motif.overlap[gid,] + tss.ov
}

pdf(paste0('analysis/variance_rankings/motifs_vs_Var_',gn,' w=',w,'_.pdf'))

#or just calculate overlap directly on the peaks....
motifs=unique(motifs.gr$name)
thresh=0.0001
h5sumgenes.gr$motif.overlap=sapply(motifs,function(motif){
      mots=motifs.gr[motifs.gr$name==motif & motifs.gr$pvalue < thresh]
      countOverlaps(resize(h5sumgenes.gr,width=500,fix='center'),mots)
})
for(motif in motifs){
h5sumgenes.gr$mot.overlap = h5sumgenes.gr$motif.overlap[ ,motif ] > 0 
# qplot(geom='errorbar',log='y', y = genetic.var, x = h5sumgenes.gr$TATA.overlap , fill= h5sumgenes.gr$TATA.overlap )
d=data.frame(motif=factor(h5sumgenes.gr$mot.overlap) , 'GeneticVariance' = genetic.var)
print(ggplot(d, aes(factor(motif),GeneticVariance,color=motif)) +
  stat_summary(fun.data='mean_cl_boot',geom='errorbar',size=4) +
  stat_summary(fun.data = give.n, geom = "text", fun.y = median) +
  ggtitle(paste0(motif,' presence vs ',gn)))
  scale_y_continuous(limits=c(0,2))
}

dev.off()



#could also skip tss with no expression, or just use the main one for each gene.
#now show boxplots for various sets.
expr.pattern = rep(NA,length(genes.gr))
expr.pattern[genes.gr$endoderm] = 'endoderm'
expr.pattern[genes.gr$mesoderm] = 'mesoderm'
expr.pattern[genes.gr$ectoderm] = 'ectoderm'
expr.pattern[genes.gr$just.ub] = 'ubiquitious'
h5sumgenes.gr$expr.pattern = expr.pattern[match(h5sumgenes.gr$Gene,genes.gr$id)]

pdf(paste0('analysis/variance_rankings/indivexpr_features_vs_Var_',gn,' w=',w,'_.pdf'))

# Expression catagory (from BDGP)
qplot(log = 'y', main = 'expression pattern vs. Genetic variability' ,
 x = h5sumgenes.gr$expr.pattern, fill = h5sumgenes.gr$expr.pattern, y = genetic.var, geom = 'Violin' )

expr.tab=c(table(h5sumgenes.gr$expr.pattern),'NA'=sum(is.na(h5sumgenes.gr$expr.pattern)))
expr.tab=data.frame(pattern=names(expr.tab),freq=expr.tab)

d=data.frame('Pattern'=factor(h5sumgenes.gr$expr.pattern) , 'GeneticVariance' = genetic.var)

print(
ggplot(d, aes(factor(Pattern),GeneticVariance)) +
  stat_summary(data=d,fun.data='mean_cl_boot',geom='errorbar',size=4, colour = c('yellow','red','orange','blue','black')) +
  stat_summary(fun.data = give.n, geom = "text", fun.y = median) 
)

dev.off()

}
#Finally, take the number of motifs (any )

#Now let's output files for Davids analysis suite.

for(gn in names(genetic.vars)){
  genetic.var = genetic.vars[[gn]]

  hasID = !is.na(h5sumgenes.gr$Gene)
  IDs = genes.gr[match(h5sumgenes.gr$Gene[hasID],genes.gr$id)]$acc
  var.GOtab.df = data.frame(gene_ID = IDs, score = genetic.var[hasID])
  write.table(var.GOtab.df[hasID,],file=paste0('analysis/GO_BDGP_analysis/go_var_table_',gn,'_w',w,'_.txt'),sep='\t',col.names=T,row.names=F,quote=F)
}






