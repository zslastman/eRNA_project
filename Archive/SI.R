#This script is for


setwd('/g/furlong/Harnett/TSS_CAGE_myfolder/')
source('src/load_annotations.R')
source('src/generate_Rle.R')
source('src/generate_chromdata_cage.R')
#source('src/normalization.R')
library(reshape2)
library(ggplot2)

signal.cutoff<-100
load(file.cage.tag.rles)
load(file.alltags)
accs<-names(cage.tag.rles)

#shape indices


#define the ranges we want to look at

#first make a shortened set to do the trial runs on
#grs=list(crm8008=crm8008.gr,CAD3=cad3.gr,TSS=tss.gr)
# cad3.gr=cad3.gr[1:20]
# accs<-names(cage.tag.rles)[1:3]
# tss.gr<-tss.gr[ tss.gr$TrID%in%unlist(cad3.gr$transcripts)]


grs=list(crm8008=crm8008.gr,CAD3=cad3.gr,TSS=tss.gr)


#get a giant grange
gr<-combinegrs(grs)
#add an 'NA' if it doesn't


gr<-resize(gr,500,fix='center')

gr<-keepSeqlevels(gr,chrs.keep)
polrlelist<-chrom.rles.rpgc.sub.merge[[1]]

#make our crms the same size
gr<-resize(gr,width=500,fix='center')

# collect statistics on chromatin and overall cage ------------------------

#get the view of PolII occupancy for each granges
polviews<-Views(polrlelist,as(gr,'RangesList'))[unique(seqnames(gr))]

#check again which ones are intergenic
intergenic <- distanceToNearest(gr,transcripts.gr)$distance>intergenic_dist

#now re-do the polymerase and other chromatin views Views
polviews<-Views(polrlelist,as(gr,'RangesList'))[unique(seqnames(gr))] 
polsums<-unlist( viewSums(polviews) )
polmaxs<-unlist(viewMaxs(polviews))

K4me1views<-Views(chrom.rles.rpgc.sub.merge$H3K4me1_6.8,as(gr,'RangesList'))[unique(seqnames(gr))] 
K4me1sums<-unlist( viewSums(K4me1views) )
K4me1maxs<-unlist(viewMaxs(K4me1views))

K4me3views<-Views(chrom.rles.rpgc.sub.merge$H3K4me3_6.8,as(gr,'RangesList'))[unique(seqnames(gr))] 
K4me3sums<-unlist( viewSums(K4me3views) )
K4me3maxs<-unlist(viewMaxs(K4me3views))

K4ratio=K4me3maxs/K4me1maxs

#get the max out of our views object and enter it into the granges object
gr@elementMetadata[['Pollmax.pos']]<- unlist( viewWhichMaxs(polviews) )
gr@elementMetadata[['Pollsum']]<- unlist( viewSums(polviews) )
#note that we DONT re do our polymerase maximum - we want want the original peak, proximal to the crm

#get the summed CAGE tag info
strandviews<-sapply(alltags,function(strand){Views(strand,as(gr,'RangesList'))[unique(seqnames(gr))]})
strandsums<-sapply(strandviews,function(strand){unlist(viewSums(strand))})
strandmaxs<-sapply(strandviews,function(strand){unlist(viewMaxs(strand))})
gr$plusisgreater<-apply(strandmaxs,1,function(x)x[1]>=x[2])
gr$allsum<-strandsums[,1]+strandsums[,2]
gr$domsum<-pmax(strandsums[,1],strandsums[,2])


#now we split our gr object again
grlist<-sapply(unique(gr$region),function(reg){gr[gr$region==reg]})

#get a list of matrices for each region with the positive, negative, summed, dominant and weaker strand for
#each one
matlist<-sapply(simplify=F,grlist,function(gr){
  mcols(gr)<-NULL
  #now get matrices with columns for each library
  grlibmats<-sapply(simplify=F,c('pos','neg'),function(strand){
    simplify2array(mclapply(mc.cores=10,accs,function(acc){
      unlist(viewSums(Views(cage.tag.rles[[acc]][[strand]],as(gr,'RangesList')))[ chrs.keep])
    }))
  })
  stopifnot(is.matrix(grlibmats[['pos']]))
  
  grlibmats$both <- grlibmats$pos+grlibmats$neg
  
  #defining dominant strand
  grlibmats$dominant<-grlibmats$neg
  grlibmats$weak<-grlibmats$pos
  grlibmats$dominant[gr$plusisgreater,]<-grlibmats$pos[gr$plusisgreater,]
  grlibmats$weak[gr$plusisgreater,]<-grlibmats$neg[gr$plusisgreater,]
  grlibmats$plusigreater<-gr$plusisgreater
  
  grlibmats
})


#function that takes in a grange object and a 




##1)his section assigns enhancers to their targets

#1a)The cad enhancers

cadtargettrs<-sapply(cad3.gr$transcripts,function(trl){
  sapply(trl,function(tr){which(tss.gr$TrID==tr)})
  })#target gene

#also get the nearest TSS for each CAD
dist<-distanceToNearest(cad3.gr,tss.gr)
cad3.gr$nearest.tss<-dist$subjectHits
cad3.gr$dist.to.nearest.tss<-dist$distance


#1b) also for the 8008
dist<-distanceToNearest(crm8008.gr,tss.gr)
cad.nearest.tss.gr<-tss.gr[dist$subjectHits]
cad.nearest.tss.gr$dist<-dist$distance

#looking at bimodality
sapply(1:length(grlist),function(n){ grlist$shape.indexes<-apply(matlist$dominant,1,shape.index)
                                     
strandratio<-rowSums(grlibmats$dominant)/rowSums(grlibmats$weak)
strandratio[is.infinite(strandratio)]<-0
strandratio[is.nan(strandratio)]<-0




#calculate shape index for CRMs
shape.index<-function(v){
  v<-v[v!=0]
  s<-sum(v)
  2+sum((v/s)*log2(v/s))
}


tmp<-sapply(1:length(grlist),function(n){ grlist$shape.indexes<<-apply(matlist$dominant,1,shape.index)

                                     
                                     
save.image('SI.image.R')






#add target info to the cads
grlist$CAD3$target.tss<-cadtargettrs
#and ids to the tss
grlist$TSS$FBtr<-tss.gr$TrID

# filter out weak cage signal ---------------------------------------------

#to start with just use the summed reads
#filter with these as late as possible
grlist$CAD3$abovecut<-grlist$CAD3$domsum>signal.cutoff
grlist$crm8008$abovecut<-grlist$crm8008$domsum>signal.cutoff

# grlist$CAD3<-grlist$CAD3[grlist$CAD3$domsum>signal.cutoff]
# grlist$crm8008<-grlist$crm8008[ grlist$crm8008$domsum>signal.cutoff]


#get nearest non target transcripts for cad
nontargets<-grlist$TSS[-unlist(grlist$CAD3$target.tss)]#exclude targetted trs
cad.nearest.nontarget<-nontargets[nearest(cad3.gr,nontargets)]
nomatch<-sum(is.na(cadtargettrs))

#number of matching tss per enhancer
pdf('analysis/cad3_assigned_gene_num.pdf')
hist(sapply(cadtargettrs,function(cadtargettrs){if(is.na(cadtargettrs)){0}else{length(cadtargettrs)}  }),main='distribution of tss matches for cad3 enhancers',
     xlab='number of matching tss',breaks=20,sub=paste0('note that ',nomatch,' of ',length(cad3.gr),' enhancers have no gene assigned'))
dev.off()

#Do shape CAD CRMs have similiar shape indices to their targets?

#go through each CRM and get the shape index of it's target
cad.target.si<-do.call(rbind,lapply(1:length(grlist$CAD3),function(n){
  if(length(grlist$CAD3$target.tss[[n]][[1]])==0){return(NULL)}#skip cases where there's no match
  if(!grlist$CAD3$abovecut[n]){return(NULL)}#skip cases where the signal is weak
  cbind(grlist$CAD3.strong$shape.indexes[n], grlist$TSS$shape.indexes[ grlist$CAD3$target.tss[[n]] ])
}))

#as data frame with two columns
cad.target.si<-as.data.frame(cad.target.si)
colnames(cad.target.si)<-c('CRM','Target')

qplot( cad.target.si$CRM ,cad.target.si$Target)

l<-lm( cad.target.si$CRM ~cad.target.si$Target)
summary(l)




# now for 8008 crms -------------------------------------------------------

#and the nearest ones to the 8008 crms
crm.nearest<-grlist$TSS[nearest(grlist$crm8008,grlist$TSS)]$regnum


#go through each CRM and get the shape index of it's target
crm.target.si<-do.call(rbind,lapply(1:length(crm.nearest),function(n){
  if(length(crm.nearest[n])==0){return(NULL)}#skip cases where there's no match
  if(!grlist$crm8008$abovecut[n]){return(NULL)}#skip cases where the signal is weak
  cbind(grlist$crm8008$shape.indexes[n], grlist$TSS$shape.indexes[crm.nearest[[n]] ])
}))

crm.target.si<-as.data.frame(crm.target.si)
colnames(crm.target.si)<-c('CRM','Target')

#plot this
qplot( crm.target.si$CRM ,crm.target.si$Target)
qplot( crm.target.si$CRM ,crm.target.si$Target)


l<-lm( crm.target.si$CRM ~crm.target.si$Target)

summary(l)



# Now see if we can detect correlations in eRNA/TSS signal ----------------

  


#go through each CRM 
cad.target.si<-do.call(rbind,lapply(1:length(grlist$CAD3),function(n){
  if(length(grlist$CAD3$target.tss[[n]][[1]])==0){return(NULL)}#skip cases where there's no match
  if(length(grlist$CAD3$[[n]][[1]])==0){return(NULL)}#skip cases where there's no match
  

  #skip cases where the crm is very weak
  #go through all targets
  #skipping cases where the target is weak
  #perform the regression, outputting crm number, tss number, and pvalue
  cbind(grlist$CAD3$shape.indexes[n], grlist$TSS$shape.indexes[ grlist$CAD3$target.tss[[n]] ])

}))


v1<-rnorm(10,10)
v2<-rnorm(10,10)

plot(v1,v2)
l<-lm(v1~v2)
s<-summary(l)
s$
  


cage.mat<-getBestWindowMat(crm8008.gr[1:10],100,cage=cage.tag.rles.rpgc)
cage.mat
jpeg('tmp.jpeg')
print(plot(cage.mat[6,],cage.mat[7,]))
a=cage.mat[6,]
b=cage.mat[7,]

dev.off()

vector.pvals <-function(a,b,weights=NULL){
  model=lm(a~b)
  s=summary(model)
  pval=pf(q=s$fstatistic[1],df1=s$fstatistic[2],df2=s$fstatistic[3],lower.tail=F)
  return(pval)
}

                                          
vector.pvals(a,b)