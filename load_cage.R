###This is a script to generate RLE objects out of the bam files for our CAGE data

# @ author Dermot Harnett, EMBL Heidelberg
# @date 16/5/2013
# @title Export coverage into wig format
########################################
#This script loads up the the cage data form teh bams, parses the filenames
#into sensible metadata. It then uses Jack's mappability data to create
#'mapfiltered' cage data, and then it power law normalizes both (using
#the unfiltered data to calculate parameters for the normalization
#it then creates merged data files which consider technical replicates together
#and finally it creates sumed files which have the toals for each timesstage

setwd(dir='/g/furlong/Harnett/TSS_CAGE_myfolder/')
source('src/tss_cage_functions.R')
library(VGAM,lib.loc='~/Harnett/R')
library(igraph,lib.loc='~/Harnett/R')
library(ggplot2)
library(testthat)
#variables
refsize=100000000


################################################################################
#1 read in the bam files as rles
cage_bam_files.68 = list.files('/g/furlong/Harnett/24_TSSCAGE/data/solexa/bam/',full.names=T,pattern='_68h.*trimmed.*quality.bwa.sort.bam$')
cage_bam_files.24 = list.files('/g/furlong/Harnett/24_TSSCAGE/data/solexa/bam/',full.names=T,pattern='24h.*trimmed.*quality.bwa.sort.bam$')
cage_bam_files.1012 = list.files('/g/furlong/Harnett/24_TSSCAGE/data/solexa/bam/',full.names=T,pattern='1012h.*trimmed.*quality.bwa.sort.bam$')
cage_bam_files.meso = list.files('/g/furlong/Harnett/24_TSSCAGE/data/solexa/bam/',full.names=T,pattern='splitmeso.*sort.*.bam$')
cage_bam_files = c(cage_bam_files.68,cage_bam_files.24,cage_bam_files.1012)


################################################################################
#2 perpare metadata dataframe from the filenames
accession.df<-data.frame(
  accession=gsub(pattern='.*(split[meso]?\\w\\w\\w?_\\w+_\\d\\d?\\d\\d?h.*?).trimmed.*',replacement='\\1',x=cage_bam_files,perl=T),
  sample=gsub(pattern='.*split(\\d\\d?\\w)_\\w+_\\d\\d?\\d\\d?h.*?.trimmed.*',replacement='\\1',x=cage_bam_files,perl=T),
  run=gsub(pattern='.*split(\\d\\d?)\\w_\\w+_\\d\\d?\\d\\d?h.*?.trimmed.*',replacement='\\1',x=cage_bam_files,perl=T),
  lane=gsub(pattern='.*split\\d\\d?(\\w)_\\w+_\\d\\d\\d?\\d?h.*?.trimmed.*',replacement='\\1',x=cage_bam_files,perl=T),
  reseq= grepl(pattern='reseq',x=cage_bam_files,perl=T),
  line=gsub(pattern='.*split\\w\\d?\\w_(\\d+).*?\\d\\d\\d?\\d?h.*?.trimmed.*',replacement='\\1',x=cage_bam_files,perl=T),
  timepoint=gsub(pattern='.*split\\w\\d?\\w_\\w+_(\\d\\d\\d?\\d?h).*?.trimmed.*',replacement='\\1',x=cage_bam_files,perl=T),
  tissue='embryo',
  stringsAsFactors=F
  )  
#add data for mesodermal cage libraries
accession.df=rbind(accession.df,data.frame(
  accession=gsub(pattern='.*splitmeso_(\\w\\d\\d?_r\\d).*trimmed.*',replacement='\\1',x=cage_bam_files.meso,perl=T),
  sample = NA,
  run=NA,
  lane=NA,
  reseq=F,
  line='BiTS_line',
  timepoint=gsub(pattern='.*meso(\\d\\dh).*',replacement='\\1',x=cage_bam_files.meso,perl=T),
  tissue='meso',
  stringsAsFactors=F
  ))
#convert all our line information to the RAL #
accession.df$RAL=!grepl('25\\d\\d\\d',x=accession.df$line)
linetable = read.table('data/line_name_table.txt',header=T)
bloom2ral = linetable[,2]
names(bloom2ral) = linetable[,1]
accession.df$line[ !accession.df$RAL ] = bloom2ral[accession.df$line[ !accession.df$RAL]]
#now insert the correct line numbers into our accessions
inds=!is.na(accession.df$sample)
accession.df$acc[inds] = paste0('split',accession.df$sample[inds],'_',accession.df$line[inds],'_',accession.df$timepoint[inds],sep='')
# #now mark out duplicates as reseqs
dup=duplicated(accession.df$acc)
accession.df[dup,]$acc = paste0(accession.df[dup,]$acc,'_reseq')
accession.df[dup,]$reseq = T

head(accession.df)

################################################################################
#3 Read in data on replicates
message('Reading in data on replicates')
#tp - line - col - prep - seq
accession.df$collection = 1:nrow(accession.df)
accession.df$prep = 1:nrow(accession.df)
accession.df$seq = 1:nrow(accession.df)
#get info on bad samples and replicates
replicatetable = read.table('data/replicate_pairs.txt',header=T,sep='\t')
badsamples=as.character(unlist(replicatetable[as.logical(replicatetable$IsBad),1:2]))
accession.df$isbadsample = accession.df$sample %in% badsamples
#go through each accession, and if it's a 'reseq', change it's collection and prep to the match
for(i in 1:nrow(accession.df)){
  if(accession.df[i,]$reseq){
      reseq<-accession.df$accession[i]
      orig<-gsub(reseq,pattern='_reseq',replacement='')
      o=which(accession.df$accession==orig)
      stopifnot(length(o)>0 & accession.df$line[i]==accession.df$line[o])
      accession.df$collection[i]=accession.df$collection[o]
      accession.df$prep[i]=accession.df$collection[o]
      # accession.df$
  }
}
#go through each technical replicate and change the collections to be the same
for(i in 1:nrow(replicatetable)){
  if(replicatetable[i,]$ReplicateType=='Technical'){
      r1<-replicatetable$Sample.ID.r1[i]
      r2<-replicatetable$Sample.ID.r2[i]
      i=which(accession.df$sample==r1)
      o=which(accession.df$sample==r2)
      stopifnot(length(o)>0 & accession.df$line[i]==accession.df$line[o])
      accession.df$collection[i]=accession.df$prep[o]
  }
}
save(accession.df,file='data/objects/accession.df.object.R')
#18E and 27H are bio reps
#13B and 27B are techreps
accession.df[accession.df$sample%in%c('18E','27H','13B','27B'),]


################################################################################
#4 Read in the bam files as RLE lists --------------------------------------
message('reading bam files....')
cg<-mclapply(mc.cores=5,mc.cleanup=T,c(cage_bam_files),function(x)bam2coverage(x,doshift=F,doresize=T,stranded=T,fragment.length=1))
# save(cg,file='data/objects/cg.object.R')
#now for the mesodermal
cg<-c(cg,mclapply(mc.cores=5,mc.cleanup=T,c(cage_bam_files.meso),function(x)bam2coverage(x,doshift=F,doresize=T,stranded=T,fragment.length=1)))
#and name our Rle object appropriately
names(cg) <-accession.df$accession 
#create a third strand for each library with the other two summed
for(acc in 1:length(cg)){cg[[acc]][['both']]<-cg[[acc]]$pos+cg[[acc]]$neg}
#calculate size of  library
accession.df$library.size <- sapply(names(cg),function(acc){sum(as.numeric(sum(cg[[acc]]$both)))})
accession.df$genome.coverage<-(accession.df$library.size/sum(seqlengths(si)))
#####Now incorporate Jack's mapability data,setting non mappabable sites to zero
cg.mapfilt <- cg
load('data/objects/allmap.object.R')
nallmap=!allmap
cg.mapfilt=mclapply(mc.cores=5,names(cg),function(acc){
  sapply(c('pos','neg'),function(s){
    cg[[acc]][[s]][nallmap]<-0
    cat('.')
    cg[[acc]][[s]]
  })
})
save(cg.mapfilt,file='data/objects/cg.mapfilt.object.R')
save(cg,file='data/objects/cg.object.R')
#load('data/objects/cg.object.R')
#load('data/objects/accession.df.full.object.R')

# Create boxplot
# Density plot
# Save plot as pdf
pdf(file="fileName.pdf")
plot(density(rnorm(100),na.rm=TRUE, data=dataName, legend=T, xlab="xLabel", ylab="yLabel", main="main label here")
dev.off()

# Save plot as pdf
pdf(file="fileName.pdf")
plotHere
dev.off()
# Save plot as pdf
pdf(file="fileName.pdf")
plotHere
dev.off()

# Create boxplot
boxplot(DV~IV, data=dataName, horizontal=FALSE, legend=T, xlab="xLabel", ylab="yLabel", ylim=c(0,100), main="main label here")

# Density plot
plot(density(DV,na.rm=TRUE, data=dataName, legend=T, xlab="xLabel", ylab="yLabel", main="main label here")

################################################################################
## 5  do the power law normalization 
fits<-sapply(accession.df$accession,function(acc){
  sitecounts=sort(c(unlist(unname(cg[[acc]][['pos']])),unlist(unname(cg[[acc]][['neg']]))))
  sitecounts=sitecounts[sitecounts!=0]
  fit=power.law.fit(as.vector(sitecounts),xmin=200)#not that power.law.fit excludes the lower count numbers
  total.tags=sum(sitecounts)
  o=getOffset(fit$alpha,total.tags)#calculate the offset (second parameter)
  c(fit,offset=o,tagcount=total.tags)#now tack it onto the fit list and return
})
accession.df$alpha<-as.numeric(fits['alpha',])
mean.alpha=mean(accession.df$alpha)
r.offset <-getOffset(mean.alpha,refsize)
accession.df$offset<-as.numeric(fits['offset',])
#4a now finally normalize all of our cage libraries to a common power law reference
cg.pl<-mapply(SIMPLIFY=F,names(cg),accession.df$alpha,accession.df$offset,FUN=function(acc,alpha,offset){
   lapply(cg[[acc]],function(srle){ 
    pl.norm(srle,x.alpha=alpha,x.offset=offset[[1]])
  })
})
for(acc in names(cg.pl)){cg.pl[[acc]][['both']]<-cg.pl[[acc]]$pos+cg.pl[[acc]]$neg}
names(cg.pl)<-accession.df$acc
save(cg.pl,file='data/objects/cg.pl.object.R')
#and also the pr-map filtered stuff
cg.mapfilt.pl<-mapply(SIMPLIFY=F,names(cg.mapfilt),accession.df$alpha,accession.df$offset,FUN=function(acc,alpha,offset){
   lapply(cg.mapfilt[[acc]],function(srle){ 
    pl.norm(srle,x.alpha=alpha,x.offset=offset[[1]])
  })
})
for(acc in names(cg.mapfilt.pl)){cg.mapfilt.pl[[acc]][['both']]<-cg.mapfilt.pl[[acc]]$pos+cg.mapfilt.pl[[acc]]$neg}
names(cg.mapfilt.pl)<-accession.df$acc
save(cg.mapfilt.pl,file='data/objects/cg.mapfilt.pl.object.R')

################################################################################

accession.df$accession %in% names(cg.pl)
accession.df$acc%in% names(cg.pl)

splitaccs=with(accession.df,split(acc,paste0(timepoint,tissue)))
splitaccs[[3]] %in% names(cg.pl)


#6 create a summed alltags object for each timepoint
message('summing cage tags')
#split our accs into a list of grouped libraries
splitaccs=with(accession.df,split(accession,paste0(timepoint,tissue)))
#now use this list to sum up the libraries
alltaglist=mclapply(mc.cores=4,splitaccs,function(accs){
  alltags=list(
      pos=Reduce('+',sapply(cg[[acc]],'[[','pos')),
      neg=Reduce('+',sapply(cg[[acc]],'[[','neg'))
    )
  alltags$both=alltags$pos+alltags$neg#and add a both rle for each one.
  alltags
})
#export bigwigs for these
for(set in names(alltaglist)){
    export(alltaglist[[set]]$pos,paste0('/g/furlong/Harnett/TSS_CAGE_myfolder/data/solexa/wig/allcage.',set,'.pos.bw'))
    export(alltaglist[[set]]$neg,paste0('/g/furlong/Harnett/TSS_CAGE_myfolder/data/solexa/wig/allcage.',set,'.neg.bw'))
}
#Also one which just has the sum of ALL sites
save(alltaglist,file='data/objects/alltaglist.object.R')
##And for the pl normalized 
#now use this list to sum up the libraries
alltaglist.pl=mclapply(mc.cores=4,splitaccs,function(accs){
  alltags=list(
      pos= Reduce('+',sapply(cg.pl[accs],'[[','pos')) ,
      neg=Reduce('+',sapply(cg.pl[accs],'[[','neg'))
    )
  alltags$both=alltags$pos+alltags$neg#and add a both rle for each one.
  alltags = lapply(alltags,'/',length(accs))
})
#export bigwigs for these
for(set in names(alltaglist.pl)){
    export( alltaglist.pl[[set]]$pos,paste0('/g/furlong/Harnett/TSS_CAGE_myfolder/data/solexa/wig/allcage.',set,'.pl.pos.bw'))
    export(alltaglist.pl[[set]]$neg,paste0('/g/furlong/Harnett/TSS_CAGE_myfolder/data/solexa/wig/allcage.',set,'.pl.neg.bw'))
    export(alltaglist.pl[[set]]$both,paste0('/g/furlong/Harnett/TSS_CAGE_myfolder/data/solexa/wig/allcage.',set,'.pl.both.bw'))

}
#Also one which just has the sum of ALL sites
save(alltaglist.pl,file='data/objects/alltaglist.pl.object.R')



################################################################################
##7now create rles with our replciates and reseqs merged
#now create lists, but 
cg.mapfilt.merge=cg.mapfilt
message('CAGE Rle objects calculated, now summing reseqs')
reseqs<-accession.df$acc[grepl(accession.df$acc,pattern='reseq')]
for(reseq in reseqs){
  orig<-gsub(reseq,pattern='_reseq',replacement='')
  cg.mapfilt.merge[[orig]]$pos<-cg.mapfilt.merge[[orig]]$pos+cg.mapfilt.merge[[reseq]]$pos
  cg.mapfilt.merge[[orig]]$neg<-cg.mapfilt.merge[[orig]]$pos+cg.mapfilt.merge[[reseq]]$neg
}
#now get rid of the reseqs  
cg.mapfilt.merge<-cg.mapfilt.merge[!names(cg.mapfilt.merge) %in% reseqs]
#now the replicates
message('... now summing replicates')
replicates<-accession.df$accession[grepl(accession.df$accession,pattern='(+\\d)_?r2')]  
for(replicate in replicates){
  orig<-gsub(reseq,pattern='(+\\d)_?r2',replacement='\\1')
  cg.mapfilt.merge[[orig]]$pos<-cg.mapfilt.merge[[orig]]$pos+cg.mapfilt.merge[[replicate]]$pos
  cg.mapfilt.merge[[orig]]$neg<-cg.mapfilt.merge[[orig]]$pos+cg.mapfilt.merge[[replicate]]$neg
}
#now get rid of the replicates  
cg.mapfilt.merge<-cg.mapfilt.merge[!names(cg.mapfilt.merge) %in% replicates]
save(cg.mapfilt.merge,file='data/objects/cg.mapfilt.merge.object.R')

################################################################################
#7b and for the normalized data
cg.mapfilt.merge.pl <- cg.mapfilt.pl
message('CAGE Rle objects calculated, now summing reseqs')
reseqs<-names(cg.mapfilt.merge.pl)[grepl(names(cg.mapfilt.merge.pl),pattern='reseq')]

for(reseq in reseqs){
  orig<-gsub(reseq,pattern='_reseq',replacement='')
  cg.mapfilt.merge.pl[[orig]]$pos<-cg.mapfilt.merge.pl[[orig]]$pos+cg.mapfilt.merge.pl[[reseq]]$pos
  cg.mapfilt.merge.pl[[orig]]$neg<-cg.mapfilt.merge.pl[[orig]]$pos+cg.mapfilt.merge.pl[[reseq]]$neg
}

#now get rid of the reseqs  
cg.mapfilt.merge.pl<-cg.mapfilt.merge.pl[!names(cg.mapfilt.merge.pl) %in% reseqs]

message('CAGE Rle objects calculated, now summing replicates')

replicates<-names(cg.mapfilt.merge.pl)[grepl(names(cg.mapfilt.merge.pl),pattern='(+\\d)_?r2')]  
for(replicate in replicates){
  orig<-gsub(reseq,pattern='(+\\d)_?r2',replacement='\\1')
  cg.mapfilt.merge.pl[[orig]]$pos<-cg.mapfilt.merge.pl[[orig]]$pos+cg.mapfilt.merge.pl[[replicate]]$pos
  cg.mapfilt.merge.pl[[orig]]$neg<-cg.mapfilt.merge.pl[[orig]]$pos+cg.mapfilt.merge.pl[[replicate]]$neg
}

  #now get rid of the replicates  
cg.mapfilt.merge.pl<-cg.mapfilt.merge.pl[!names(cg.mapfilt.merge.pl) %in% replicates]
save(cg.mapfilt.merge.pl,file='data/objects/cg.mapfilt.pl.merge.object.R')




