
library(rhdf5)
summaryHDF5 <- '/g/furlong/Harnett/TSS_CAGE_myfolder/data/olliver_processed_summary.hdf5'
Ccis = h5read(summaryHDF5,name='Ccis')
############now look at the variance
examplematmat=matrix(1:9,ncol=3)
delta=  abs(row(mat)-col(mat))!=0#logical matrix for extracting diagonals - true when NOT diagonal
gnum=dim(Ccis)[3]
#Now go through the covariance matrices and subtract to get our stage specific variation
m=t(sapply((1:gnum),function(i){
      mat = Ccis[,,i]#for each gene
      offdiag=mat[delta]#get off diagonal elements
      shared = max(0,min(offdiag))#now the maximum of these and the 
      specific = diag(mat) - shared
}))
m=t(m)#so that the columns are timepoints
mean( m < 0 ) #shows that about 15% of the stage specific variances are negative ?!


var.shared = sapply(simplify=F,list(Ccis,Ctrans,Cnoise),function(x){
    m=t(sapply(1:(dim(x)[3]),function(i){
      shared = max(0,min(mat[delta]))
  }))
})
names(var.spec)=c('Cis','Trans','Noise')
var.spec.m = melt(var.spec)
expect_true(all(var.spec[[1]]>0))

####set cutoffs on variance data - it seems to be bimodal
cutoffs=sapply(c('Cis','Trans'),function(ct){
  sapply(var_decomp_tps,function(tp){
    x=log10(var.spec[[ct]][,tp])
    #we need to iterate the cutoff fucntion, it fails sometimes
    for(i in 1:5){ try({cutoff = getBimodalCutoff( x,knum=2,proba=0.5)})} 
    cutoff
  })
})

#Basic plotting of variance data,shows bimodality
pdf(file='analysis/gene_variance/VarDensity_Densityplot.pdf')
  
  for(tp in var_decomp_tps){
    transvar = var.spec$Trans[,tp]
    cisvar = var.spec$Cis[,tp]
    maintit = paste0(tp,' Distribution of Cis/Trans Variance Scores')
    print(qplot(geom='point',alpha=0.1,log='xy',y=transvar ,x=cisvar,main=maintit)+
    geom_hline(yintercept=10^cutoffs[tp,'Trans'])+
    geom_vline(xintercept=10^cutoffs[tp,'Cis']))
  }
  qplot(geom='density',data=var.spec,log='x',x=as.numeric(value),color=factor(L1),main='Variance Sources')
dev.off()

