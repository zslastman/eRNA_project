setwd ( '/g/furlong/Harnett/TSS_CAGE_myfolder/' )
source ( 'src/tss_cage_functions.R' )
library(corrperm)
load ( 'data/objects/cg.pl.object.R' )
load ( 'data/objects/crm8008.gr.object.R' )
load ( 'data/objects/tssgr.R' )
load ( 'data/objects/cad3.gr.object.R' )
load( file.crm8008.gr)
# load ( 'data/objects/alltags.rpgc.object.R' )

#function that defines what's too low to use and what isn't
accs=names ( cg.pl )
crm8008.gr=


# #make list of grs
# tss.gr<-tss.gr[!duplicated(tss.gr)]
# tss.gr<-tss.gr[seqnames(tss.gr)%in% bigchrs]
# tss.gr=sort(tss.gr)
# active.tss.gr<-tss.gr[  tss.gr$active68 ]
# active.tss.gr<- active.tss.gr[  seqnames ( active.tss.gr )%in%bigchrs ]
# gr=combinegrs ( list ( cad3=cad3.gr,crm=crm8008.gr,tss=active.tss.gr ) )



##Cad target interactions
trs=cad3.gr$transcripts
names(trs) <- as.character(1:length(trs))
trs=melt(trs)
trs$tssind=match(trs$value,tss.gr$TrID)
trs=trs[ !is.na(trs$tssind) , ]

cadtarget=list(as.numeric(trs$L1) ,trs$tssind)

#nearest gene interactions
near.match=distanceToNearest(crm8008.gr,tss.gr)
near.match=near.match[near.match$distance!=0,]
n.dists=near.match$distance
near.match=list(near.match$queryHits,near.match$subjectHits)#get nearest tss.
#different chr interactions
rand.far=list(rep(NA,length(crm8008.gr)),rep(NA,length(crm8008.gr)))
for (  chr in  unique ( seqnames ( crm8008.gr )  ) ){ 
  #list of chromosomes not on this one
  diffchr <-as.vector (  seqnames ( crm8008.gr ) != chr  )
  inds = which(seqnames(crm8008.gr)==chr)
  #randomly pick from tss not in this chromosomes
  targs=sample(tss.gr[seqnames(tss.gr)!=chr],length(inds))
  #now convert out targets to indices in tss.gr
  tinds=match(targs,tss.gr)
  rand.far[[1]][inds]=inds
  rand.far[[2]][inds]=tinds   
}
  
cols=colnames(cad3.gr$cad3cage)[colnames(cad3.gr$cad3cage)%in%colnames(tss.gr$tsscage)]
tss.gr$tsscage = tss.gr$tsscage[,cols]
cad3.gr$cad3cage = cad3.gr$cad3cage[,cols]
crm8008.gr$crmcage = crm8008.gr$crmcage[,cols]

iterations=1000
#now calculate the pvalues for tss-gene interactions
vector.pvals.shuffle <-function ( a,b,iter=1000 ){ 
  cor.real<-cor ( a,b,method='s' )
  cors.shuffle<-replicate ( iter,{ cor ( a,sample ( b,length ( b ) ),method='s' ) } )
  mean ( cors.shuffle>cor.real )
 }
#To get our correlations we input our two matrices and our 
get.cor<-function(cmat,tmat,indices){
  stopifnot(is.matrix(cmat))
  stopifnot(is.matrix(tmat))
  stopifnot(is.list(indices))
  mapply(indices[[1]],indices[[2]],FUN=function(a,b){
     a=cmat[a,]
     b=tmat[b,]
   # cor.obj=cp.test(nrep=1000,matrix(a),matrix(b))$S.p.value
   cor(a,b,meth='spearman')
  })
}
#To get our correlations we input our two matrices and our 
get.pvalmat<-function(cmat,tmat,indices){
  stopifnot(is.matrix(cmat))
  stopifnot(is.matrix(tmat))
  stopifnot(is.list(indices))
  mapply(indices[[1]],indices[[2]],FUN=function(a,b){
     a=cmat[a,]
     b=tmat[b,]
   # cor.obj=cp.test(nrep=1000,matrix(a),matrix(b))$S.p.value
   vector.pvals.shuffle(a,b)
  })
}
p.c = get.pvalmat( cad3.gr$cad3cage ,tss.gr$tsscage,cadtarget)
cat('.')
p.r = get.pvalmat( crm8008.gr$crmcage ,tss.gr$tsscage,rand.far)
cat('.')
p.n = get.pvalmat( crm8008.gr$crmcage ,tss.gr$tsscage,near.match)

save.image('tmp.image.cor.object.R')


dir.create('analysis/crm_tss_covariance')


pdf('analysis/crm_tss_covariance/covar_histograms.pdf')
#Visualize Correlations using histograms
par ( mfrow=c ( 2,2 ) )
#plot ( density ( unlist ( pvalmat ),na.rm=T ),xlim=c ( 0,1 ) )
hist ( unlist ( p.r ),xlim=c ( 0,1 ),breaks=50,main='p-values using random gene on other chr',xlab='pval' )
# title ( sub= paste0 ( mean ( unlist ( p.r )<0.05,na.rm=T ) )
hist ( unlist ( p.n ),xlim=c ( 0,1 ),breaks=50,main='p-values using nearest gene',xlab='pval' )
hist ( unlist ( p.c ),xlim=c ( 0,1 ),breaks=50,main='p-values using the cad3 target assignments',xlab='pval' )

dev.off()


#now condition on above-cutoff eRNA expression - i.e. a mean
p.c.filt=p.c[cadcutvect[cadtarget[[1]]]]
p.r.filt=p.r[cutvect[cadtarget[[1]]]]
p.n.filt=p.n[cutvect[cadtarget[[1]]]]


pdf('analysis/crm_tss_covariance/covar_filt_histograms.pdf')
#Visualize Correlations using histograms
par ( mfrow=c ( 2,2 ) )
#plot ( density ( unlist ( pvalmat ),na.rm=T ),xlim=c ( 0,1 ) )
hist ( unlist ( p.r.filt ),xlim=c ( 0,1 ),breaks=50,main='p-values using random gene on other chr',xlab='pval' )
# title ( sub= paste0 ( mean ( unlist ( p.r )<0.05,na.rm=T ) )
hist ( unlist ( p.n.filt ),xlim=c ( 0,1 ),breaks=50,main='p-values using nearest gene',xlab='pval' )
hist ( unlist ( p.c.filt ),xlim=c ( 0,1 ),breaks=50,main='p-values using the cad3 target assignments',xlab='pval' )

dev.off()

#now analyze the distributions using qqplots
pdf(paste0('analysis/crm_tss_covariance/' ,'QQplot_byMapFilter','.pdf'))
samplist=list(Nearest_Neighbor=p.n,Annotated_CAD_Target_gene=p.c,Random_Match=p.r,uniform_pvals=runif(min=0,max=1,n=10000))
samplist=sapply(simplify=F,samplist,sort)
sampquants=lapply(samplist,function(x)ppoints(length(x)))
samplist=stack(samplist)
samplist$unifquant=stack(sampquants)$values
samplist$nlpvals=-log10(samplist$values)
samplist$nlquant=-log10( samplist$unifquant )
# qplot(data=samplist,x=unifquant ,y=pvals,color=ind,geom='point')
qplot(data=samplist,x=nlquant ,y=nlpvals,color=ind,geom='point',xlab='-log10 uniform quantiles',ylab='-log10 Pvalues')+
scale_colour_discrete(name='Pair Type')
dev.off()

png(paste0(OUTDIR, 'QQplot_byMapFilter','.png'), height=5, width=5, res=400, units='in', type='Xlib')
i1=1:100000
plot(-log10((seq(1, length(allP)) - 0.5)/length(allP))[i1], -log10(sort(allP, decreasing=F))[i1], pch=3, cex=0.03, xlab='-log10 Uniform Quantiles', ylab='-log10 P-value Quantiles', xlim=c(0,-log10(1/length(allP))), ylim=c(0, max(-log10(allP))))
i2=seq(100000, length(allP), 500)
points(-log10((seq(1, length(allP)) - 0.5)/length(allP))[i2], -log10(sort(allP, decreasing=F))[i2], pch=3, cex=0.03)
abline(0,1, col='black', lty=2)



pdf('analysis/crm_tss_covariance/correlation_dist_scatterplot.pdf')
#Visualize Correlations using histograms
dists=distanceToNearest(crm8008.gr,tss.gr)$distance
dists=dists[dists!=0]
cors=get.cor(crm8008.gr$crmcage ,tss.gr$tsscage,near.match)
heatscatter(cor=F,x=log10(dists),y=cors,xlab='Log10 distance from gene to TSS',ylab='Spearmann correlation',main='Correlation of CRMs to the nearest gene vs. distance')
dev.off()

#




#now plot these pvalues as a function of the crms expression

inds=cadtarget
pvals=p.c
crminds=inds[[1]]]
crmexpr=rowMeans(crm8008.gr$crmcage)

pdf('analysis/crm_tss_covariance/cadtarget_pval_vs_crm_expr.pdf')
qplot(x=crmexpr,pvals= -log10(pvals),ylab='-ve log10 pvalue',xlab='crm exprr')
qplot(x=log10(crmexpr),pvals= -log10(pvals),ylab='-ve log10 pvalue',xlab='Log10 crm exprr')
dev.off()


inds=rand.far
pvals=p.c
crminds=inds[[1]]]
crmexpr=rowMeans(crm8008.gr$crmcage)
pdf('analysis/crm_tss_covariance/fartarget_pval_vs_crm_expr.pdf')
qplot(x=crmexpr,pvals= -log10(pvals),ylab='-ve log10 pvalue',xlab='crm exprr')
qplot(x=log10(crmexpr),pvals= -log10(pvals),ylab='-ve log10 pvalue',xlab='Log10 crm exprr')
dev.off()


inds=near.match
pvals=p.n
crminds=inds[[1]]]
crmexpr=rowMeans(crm8008.gr$crmcage)
pdf('analysis/crm_tss_covariance/neartarget_pval_vs_crm_expr.pdf')
qplot(x=crmexpr,pvals= -log10(pvals),ylab='-ve log10 pvalue',xlab='crm exprr')
qplot(x=log10(crmexpr),pvals= -log10(pvals),ylab='-ve log10 pvalue',xlab='Log10 crm exprr')
dev.off()




# #let's look at mean/variance over the cage mat
# cagemat<-gr$cagemat
# cagemat=cagemat[ 9000:10000, ]



# cagemat<-cagemat[ gr$region!='tss', ]
# rsums=rowSums ( cagemat )
# cagemat.sort<-cagemat[ order ( rsums ), ]
# vars<-coef.vars<-apply ( cagemat.sort,1,var )

# coef.vars<-apply ( cagemat.sort,1,function ( x ){ sd ( x )/mean ( x ) } )

# qplot ( x=rsums,y=coef.vars,log='xy',geom='2Ddensity' )

# #vars
# plot ( log10 (  ( vars ) ) )


# cagemat<-sapply ( 1:nrow ( cagemat ),function ( n ){ cagemat[ n, ]/rsums[ n ] } )
# cagemat<-t ( cagemat )
# cagemat<-cagemat[ !grepl ( NaN,cagemat[ ,1 ] ), ]
# cagemat<-cagemat[ !grepl ( Inf,cagemat[ ,1 ] ), ]
# cagemat
# plot ( colSums ( cagemat,na.rm=T ) )
# ggplot2::qplot ( x=colSums ( cagemat,na.rm=T ),geom='point' )
# #now let's look at how qvalue behave under reasonable assumptions...








# Why are pvalue so high? -------------------------------------------------

#is ther ea differential contribution from different libraries causing the problem?
ncagemat<-t ( sapply ( 1:nrow ( gr$cagemat ),function ( r ){ gr$cagemat[ r, ]/sum ( gr$cagemat[ r, ] ) } ) )

v<-accession.df$run
names ( v )<-accession.df$accession
d=melt ( ncagemat )
d$run=v[ d$Var2 ]

head ( d )
qplot ( data=d,x=Var2,y=value,geom='boxplot',log='y',color=v[ Var2 ],xlab=F,
                        ylab='Mean % of total reads at all CRMs/TSS',main='Average Contribution of normalized libraries to TSS/CRM' )
#same graph but sort grouped by lane
qplot ( data=d,x=Var2,y=value,geom='boxplot',log='y',color=run,xlab='',
      ylab='Mean % of total reads at all CRMs/TSS',main='Average Contribution of normalized libraries to TSS/CRM' )

#check if pca extraction helps
cagemat=gr$cagemat[ 10000:11000, ]
ncagemat<-t ( sapply ( 1:nrow ( cagemat ),function ( r ){ cagemat[ r, ]/sum ( cagemat[ r, ] ) } ) )

v<-accession.df$run
names ( v )<-accession.df$accession
d=melt ( ncagemat )
d$run=v[ d$Var2 ]

qplot ( data=d,x=Var2,y=value,geom='boxplot',log='y',color=run,xlab='',
      ylab='Mean % of total reads at all CRMs/TSS',main='Average Contribution of normalized libraries to TSS/CRM' )
cagemat.pca<-subtract.pcas ( cagemat,3 )
dim ( cagemat.pca )
ncagemat.pca<-t ( sapply ( 1:nrow ( cagemat.pca ),function ( r ){ cagemat.pca[ r, ]/sum ( cagemat.pca[ r, ] ) } ) )

v<-accession.df$run
names ( v )<-accession.df$accession
d=melt ( ncagemat.pca )
d$run=v[ d$Var2 ]

qplot ( data=d,x=Var2,y=value,geom='boxplot',log='y',color=run,xlab='',
      ylab='Mean % of total reads at all CRMs/TSS',main='Average Contribution of normalized libraries to TSS/CRM' )




rnbinom2<-function ( n,m,v ){ 
  s= ( m^2 )/ ( v-m )
  rnbinom ( n=n,mu=m,size=s )
 }


sim.pvals=replicate ( n=10000,{ 
  vector.pvals ( rnbinom2 ( 100,100,200 ),rnbinom2 ( 100,100,200 ) )
 } )
hist ( sim.pvals )


data ( hedenfalk )
pvals<-hedenfalk
pvals<-sim.pvals

hist ( pvals )
qobj<-qvalue ( pvals )
hist ( qobj$qvalues )
qsummary ( qobj )
qplot ( qobj )




#simulating our data


#each column has an associated library size
#each each row has an associated signal strength
#So our best bet might be to just 

libsizes=accession.df$library.size[ 1:10 ]
signals=1:10000



#subsampling to give equal sizes
sampleRleList<-function ( srle,p ){ 
  as ( 
    sapply ( names ( srle ),function ( chr ){       
      chr<-srle[[ chr ]]
      sites<-chr>0
      subsample<-rbinom ( n=length ( as.vector ( chr[ sites ] ) ),prob=p,size=as.vector ( chr[ sites ] ) )
      chr[ sites ]<-subsample
      chr
     } ),'SimpleRleList' )
 }

#first the boxplot showing unequal contributions for TSS.
#Now for crms

#repeat these but with the new normalization

#now subsample to see if this fixes the problem

#


