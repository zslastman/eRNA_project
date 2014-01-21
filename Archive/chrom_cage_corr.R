load(file.alltags)
load(file.alltags.rpgc)
source('src/generate_chromdata_cage.R')

regs=list(
  pos=import(con='analysis/positive_8008.bed',asRangedData=F,seqinfo=si),
  neg=import(con='analysis/negative_8008.bed',asRangedData=F,seqinfo=si),
  #int=import(con='analysis/random.intergenic.bed',asRangedData=F,seqinfo=si),
  cad3pos=import(con='analysis/cad3.pos.bed',asRangedData=F,seqinfo=si),
  cad3neg=import(con='analysis/cad3.neg.bed',asRangedData=F,seqinfo=si)
)



# install.packages('doBy',lib='~/Harnett/R')
 install.packages('multcomp',lib='~/Harnett/R')
library(doBy,lib.loc='~/Harnett/R')
source('')
#What level of tags actually predicts polII?
gr<-crm8008.gr
gr<-keepSeqlevels(gr,chrs.keep)
polrlelist<-chrom.rles.rpgc.sub.merge$PolII_6.8

names(polrlelist)
seqlevels(gr)
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

K4ratio=K4me3sums/K4me1sums
K4peak.l<-0< (countOverlaps(gr,modK4me3)+
                countOverlaps(gr,chrompeaks[['K4me3_6-8h']]))

#get the CAGE tag info
strandviews<-sapply(alltags,function(strand){Views(strand,as(gr,'RangesList'))[unique(seqnames(gr))]})
strandsums<-sapply(strandviews,function(strand){unlist(viewSums(strand))})
strandmaxs<-sapply(strandviews,function(strand){unlist(viewMaxs(strand))})
cagesums<-pmax(strandsums[,1],strandsums[,2])
cagemax<-pmax(strandmaxs[,1],strandmaxs[,2])

qplot(x=cagesums,y=K4me3sums)
qplot(x=cagesums,y=K4me1sums)
qplot(x=cagesums,y=K4ratio)

cagesums





