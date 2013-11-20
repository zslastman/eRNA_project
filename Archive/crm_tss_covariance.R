setwd ( '/g/furlong/Harnett/TSS_CAGE_myfolder/' )
source ( 'src/tss_cage_functions.R' )
load ( file.cage.tag.rles.plnorm )
load ( 'data/objects/crmgrs.object.R' )
load ( 'data/objects/tssgr.R' )
load ( 'data/objects/cad3.gr.object.R' )

# load ( 'data/objects/alltags.rpgc.object.R' )

#function that defines what's too low to use and what isn't
cutoff.func<-function ( m ){ rowSums ( m )/ncol(m)<0.16 }

window.size=100
libnum=length ( cage.tag.rles.plnorm )
cg=cage.tag.rles.plnorm
accs=names ( cg )

#function to derive the p value for a linera model between two vectors
vector.pvals <-function ( a,b,weights=NULL ){ 
  model=lm ( a~b )
  s=summary ( model )
  pval=pf ( q=s$fstatistic[ 1 ],df1=s$fstatistic[ 2 ],df2=s$fstatistic[ 3 ],lower.tail=F )
  return ( pval )
 }

#function to derive empiracal pvalue for the coorrelation between two vectors via shuffling
vector.pvals.shuffle <-function ( a,b,iter=1000 ){ 
  cor.real<-cor ( a,b,method='s' )
  cors.shuffle<-replicate ( iter,{ cor ( a,sample ( b,length ( b ) ),method='s' ) } )
  mean ( cors.shuffle>cor.real )
 }#this does give uniform p values with normal input

#make list of grs
active.tss.gr<-tss.gr[  tss.gr$active68 ]
active.tss.gr<- active.tss.gr[  seqnames ( active.tss.gr )%in%bigchrs ]
gr=combinegrs ( list ( cad3=cad3.gr,crm=crmgrs,tss=active.tss.gr ) )
gr=sort ( gr )

#get the matrix with info for each line on each of our grs
gr$cagemat=get.best.window.mat ( reg=nomcols ( gr ),w=window.size,chrs.keep,bigchrs,cage=cg[ accs ] )


#we should exclude anything below our cutoff
gr$above.cutoff = cutoff.func ( gr$cagemat ) 





# Target-Gene mapping -----------------------------------------------------

#we need a list mapping indexes to other indexes
#simplest version is just neighbors
#another simple version is just all against all
#first get the random set from other chrs
num=3
rand.far=sapply ( 1:length ( gr ),function ( i ){ 
  if (  gr$reg[ i ]=='tss' ){ return ( NULL ) }
  diffchr<-as.vector ( seqnames ( gr )!=seqnames ( gr )[ i ] )
  difftype<-as.vector (   ( gr$reg=='tss' ) !=  ( gr$reg[ i ]=='tss' )   )
  sample (  ( 1:length ( gr ) )[ diffchr & difftype ],num )#sample from stuff not on the same chr, and compare crms to genes and vica versa
 } )


rep = 3

is.tss  <-as.vector (  gr$reg == 'tss'  )

gr$near.match = NA
gr$rand.far.match = NA
gr$cad.target = NA

#for each chr make random pics from other chrs - of stated length
#now pply these to th
for (  chr in  unique ( seqnames ( gr )  ) ){ 
  #potential tss to choose f:rom
  diffchr <-as.vector (  seqnames ( gr ) != chr  )
  #how many crms on this chr
  compnum = sum ( !is.tss & !diffchr )
  #randomly pick tss on different chr
  randpicks <- sample ( seq_along ( gr )[  diffchr & is.tss  ],size = compnum * rep,rep=T )
  randpicks <- split ( randpicks , rep ( 1:compnum, each = rep ) )#list of vectors of length rep
  gr$rand.far.match[  !is.tss & !diffchr  ] <- randpicks
 }

#we want the nearest tss to each crm - we can do this with the regnum
gr$near[  gr$reg!='tss'  ] <- nearest ( gr[ gr$reg!='tss' ],gr[ gr$reg=='tss' ] )#find neighbor in tss
gr$near[  gr$reg!='tss'  ] <- which ( gr$reg=='tss' )[ gr$near[  gr$reg!='tss'  ]]#convert to index over whole gr


#also use the cad 
gr$cad.target <- NA
incad <- gr$reg=='cad3'

gr$cad.target[incad] =  rapply ( gr$transcripts[ incad ], how='replace' ,function ( tr ){ 
  m=which ( gr$TrID==tr ) 
  if(length(m) == 0){ return( NA ) 
  }else{m}
})

# targetlists=c ( rand.far,near,cad.target )
# targetlist=cad.target

pvalfunc=vector.pvals.shuffle

get.pvalmat<-function ( targetlist,cagemat,pvalfunc=vector.pvals.shuffle,abovecut=gr$above.cutoff ){ 
  if ( length ( abovecut )>1 ){ cutfilter=abovecut
   }else{ cutfilter=T }
  s<-which ( !is.na ( targetlist ) & sapply ( targetlist,is.vector ) & sapply ( targetlist,is.integer ) )
  t=targetlist
  #make sure our targetlist is correct format
  t<-mclapply ( mc.cores=20,s,function ( i ){ 
      sapply ( targetlist[[ i ]],function ( j ){ 
        a=cagemat[ i, ]
        b=cagemat[ j, ]
        if ( all ( a==0 )|all ( b==0 ) ){ return ( NA ) }
# 
#         
        p=tryCatch ( { 
          pvalfunc ( a,b )
                   },error=function ( e ){  stop ( paste0 ( i,'...',j )  ) }
           )
        cat ( '...j' )
        return ( p )
       } )
   } )
  
  t#output the pvalues
 }

#apply to our different target lists
p.c=get.pvalmat ( gr$cad.target,gr$cagemat )
p.n=get.pvalmat ( gr$near,gr$cagemat )
p.r=get.pvalmat ( gr$rand.far,gr$cagemat )

par ( mfrow=c ( 2,2 ) )
#plot ( density ( unlist ( pvalmat ),na.rm=T ),xlim=c ( 0,1 ) )
hist ( unlist ( p.r ),xlim=c ( 0,1 ),breaks=50,main='p-values using random gene on other chr',xlab='pval' )
# title ( sub= paste0 ( mean ( unlist ( p.r )<0.05,na.rm=T ) )
hist ( unlist ( p.n ),xlim=c ( 0,1 ),breaks=50,main='p-values using nearest gene',xlab='pval' )
hist ( unlist ( p.c ),xlim=c ( 0,1 ),breaks=50,main='p-values using the cad3 target assignments',xlab='pval' )








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


